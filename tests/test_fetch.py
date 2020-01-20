import json
import pytest
from os import path
from Bio import SeqIO
from crispy_service import fetch, utils
import requests

@pytest.fixture
def req(monkeypatch):
    """A fake replacement for the requests library"""
    class FakeReqCtx(object):
        def __init__(self, status_code, text):
            self.text = text
            self.status_code = status_code
        def raise_for_status(self):
            if self.status_code >= 400:
                raise requests.HTTPError()
    class FakeRequests(object):
        def __init__(self):
            self.retvals = []
            self.urls = []
        def set_return(self, value):
            self.retvals.append(value)
        def get(self, url):
            self.urls.append(url)
            status_code, text = self.retvals.pop(0)
            return FakeReqCtx(status_code, text)
    req = FakeRequests()
    monkeypatch.setattr(requests, 'get', req.get)
    return req


def get_testfile(filename):
    return path.join(path.dirname(__file__), filename)


def test_grab_index(req):
    req.set_return((220, 'test'))
    ret = fetch.grab_index('example')
    assert ret == 'test'
    assert req.urls[0] == 'https://antismash.secondarymetabolites.org/upload/example/index.html'


def test_grab_index_error(req):
    req.set_return((440, ''))
    with pytest.raises(requests.HTTPError):
        fetch.grab_index('example')


def test_find_dl_link():
    text = 'test\n<li><a href="NC_003888.3.zip">Download all results</a></li><li><a href="NC_003888.3.geneclusters.xls">Download XLS overview file</a></li><li><a href="NC_003888.3.final.embl">Download EMBL summary file</a></li><li><a href="NC_003888.3.final.gbk">Download GenBank summary file</a></li><li><a href="biosynML.xml" target="_blank">Download BiosynML file</a></li><li><a href="metabolicModel/antiSMASH_model_with_template_sco.xml">Download metabolic model SBXML file</a></li></ul>\ntest'
    link = fetch.find_dl_link(text)
    assert link == 'NC_003888.3.final.gbk'


def test_find_dl_link_broken_link():
    text = 'test\n<li><a href=NC_003888.gbk> Download GenBank summary file</a></li>'
    link = fetch.find_dl_link(text)
    assert link == 'NC_003888.json'


def test_find_dl_link_prefers_json():
    text = 'test\n<li><a href="NC_003888.gbk">Download GenBank summary file</a></li>'
    link = fetch.find_dl_link(text)
    assert link == 'NC_003888.json'


def test_download_gbk(req):
    req.set_return((220, "I'm totally a genbank file"))
    ret = fetch.download_gbk('example', 'file.gbk')
    assert ret == "I'm totally a genbank file"
    assert req.urls[0] == 'https://antismash.secondarymetabolites.org/upload/example/file.gbk'


def test_download_gbk_error(req):
    req.set_return((440, ''))
    with pytest.raises(requests.HTTPError):
        fetch.download_gbk('example', 'file.gbk')


def test_genome_json():
    rec = SeqIO.read(get_testfile('nisin.final.gbk'), 'genbank')
    ret = fetch.genome_json(rec)
    bp_orfs = [f for f in rec.features if f.type == 'CDS']
    orfs = utils.jsonify_orfs(bp_orfs)
    expected = {
        'description': 'Lactococcus lactis subsp. lactis nisin biosynthetic gene cluster, complete sequence',
        'length': 15016,
        'id': 'HM219853.1',
        'organism': 'Lactococcus lactis subsp. lactis',
        'clusters': [
            {'start': 0, 'end': 15016, 'name': 'Region 1', 'type': 'lantipeptide',
             'description': 'Nisin_A_biosynthetic_gene_cluster (100% of genes show similarity)'},
        ],
        'orfs': orfs,
    }
    assert ret == expected


def test_convert_json():
    with open(get_testfile('nisin.json'), 'r') as handle:
        as_json = json.load(handle)

    as_orfs = [f for f in as_json['records'][0]['features'] if f['type'] == "CDS"]
    orfs = utils.convert_json_orfs(as_orfs)
    expected = {
        'description': 'Lactococcus lactis subsp. lactis nisin biosynthetic gene cluster, complete sequence',
        'length': 15016,
        'id': 'HM219853.1',
        'organism': 'Lactococcus lactis subsp. lactis',
        'clusters': [
            {'start': 0, 'end': 15016, 'name': 'Region 1', 'type': 'lanthipeptide',
             'description': 'Nisin A'},
        ],
        'orfs': orfs,
    }
    ret = fetch.convert_json(as_json)
    assert ret == expected


def test_fetch(tmpdir, req):
    req.set_return((220, '<a href="final.gbk">Download GenBank summary file</a>'))
    with open(get_testfile('nisin.final.gbk'), 'r') as fh:
        text = fh.read()
    rec = SeqIO.read(get_testfile('nisin.final.gbk'), 'genbank')
    req.set_return((220, text))

    json_rec, filename = fetch.fetch('test', str(tmpdir))
    #FIXME: change this once creating the json record is implemented
    assert json_rec == fetch.genome_json(rec)
    assert filename == 'input.gbk'
    assert req.urls[0] == 'https://antismash.secondarymetabolites.org/upload/test/index.html'
    assert req.urls[1] == 'https://antismash.secondarymetabolites.org/upload/test/final.gbk'


def test_fetch_create_dir(tmpdir, req):
    req.set_return((220, '<a href="final.gbk">Download GenBank summary file</a>'))
    with open(get_testfile('nisin.final.gbk'), 'r') as fh:
        text = fh.read()
    rec = SeqIO.read(get_testfile('nisin.final.gbk'), 'genbank')
    req.set_return((220, text))

    json_rec, filename = fetch.fetch('test', path.join(str(tmpdir), 'nonexistent'))
    #FIXME: change this once creating the json record is implemented
    assert json_rec == fetch.genome_json(rec)
    assert filename == 'input.gbk'
    assert req.urls[0] == 'https://antismash.secondarymetabolites.org/upload/test/index.html'
    assert req.urls[1] == 'https://antismash.secondarymetabolites.org/upload/test/final.gbk'
