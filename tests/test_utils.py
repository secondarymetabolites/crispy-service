from os import path
from Bio import SeqIO
import json

from crispy_service import utils


def get_testfile(filename):
    return path.join(path.dirname(__file__), filename)


def test_jsonify_orfs():
    rec = SeqIO.read(get_testfile('results.gbk'), 'genbank')
    orfs = [f for f in rec.features if f.type == 'CDS']
    assert len(orfs) == 11
    expected = [
        {'start': 827, 'end': 1001, 'strand': 1, 'id': 'nisA', 'gene': 'nisA', 'protein_id': 'ADJ56352.1'},
        {'start': 1108, 'end': 4090, 'strand': 1, 'id': 'nisB', 'gene': 'nisB', 'protein_id': 'ADJ56353.1'},
        {'start': 4100, 'end': 5903, 'strand': 1, 'id': 'nisT', 'gene': 'nisT', 'protein_id': 'ADJ56354.1'},
        {'start': 5895, 'end': 7140, 'strand': 1, 'id': 'nisC', 'gene': 'nisC', 'protein_id': 'ADJ56355.1'},
        {'start': 7136, 'end': 7874, 'strand': 1, 'id': 'nisI', 'gene': 'nisI', 'protein_id': 'ADJ56356.1'},
        {'start': 7875, 'end': 9924, 'strand': 1, 'id': 'nisP', 'gene': 'nisP', 'protein_id': 'ADJ56357.1'},
        {'start': 9992, 'end': 10679, 'strand': 1, 'id': 'nisR', 'gene': 'nisR', 'protein_id': 'ADJ56358.1'},
        {'start': 10671, 'end': 12015, 'strand': 1, 'id': 'nisK', 'gene': 'nisK', 'protein_id': 'ADJ56359.1'},
        {'start': 12113, 'end': 12791, 'strand': 1, 'id': 'nisF', 'gene': 'nisF', 'protein_id': 'ADJ56360.1'},
        {'start': 12792, 'end': 13521, 'strand': 1, 'id': 'nisE', 'gene': 'nisE', 'protein_id': 'ADJ56361.1'},
        {'start': 13507, 'end': 14152, 'strand': 1, 'id': 'nisG', 'gene': 'nisG', 'protein_id': 'ADJ56362.1'},
    ]
    ret = utils.jsonify_orfs(orfs)
    assert ret == expected


def test_convert_json_ofs():
    with open(get_testfile("nisin.json"), 'r') as handle:
        as_json = json.load(handle)

    orfs = [f for f in as_json['records'][0]['features'] if f['type'] == "CDS"]
    assert len(orfs) == 11
    expected = [
        {'start': 827, 'end': 1001, 'strand': 1, 'id': 'nisA', 'gene': 'nisA', 'protein_id': 'ADJ56352.1'},
        {'start': 1108, 'end': 4090, 'strand': 1, 'id': 'nisB', 'gene': 'nisB', 'protein_id': 'ADJ56353.1'},
        {'start': 4100, 'end': 5903, 'strand': 1, 'id': 'nisT', 'gene': 'nisT', 'protein_id': 'ADJ56354.1'},
        {'start': 5895, 'end': 7140, 'strand': 1, 'id': 'nisC', 'gene': 'nisC', 'protein_id': 'ADJ56355.1'},
        {'start': 7136, 'end': 7874, 'strand': 1, 'id': 'nisI', 'gene': 'nisI', 'protein_id': 'ADJ56356.1'},
        {'start': 7875, 'end': 9924, 'strand': 1, 'id': 'nisP', 'gene': 'nisP', 'protein_id': 'ADJ56357.1'},
        {'start': 9992, 'end': 10679, 'strand': 1, 'id': 'nisR', 'gene': 'nisR', 'protein_id': 'ADJ56358.1'},
        {'start': 10671, 'end': 12015, 'strand': 1, 'id': 'nisK', 'gene': 'nisK', 'protein_id': 'ADJ56359.1'},
        {'start': 12113, 'end': 12791, 'strand': 1, 'id': 'nisF', 'gene': 'nisF', 'protein_id': 'ADJ56360.1'},
        {'start': 12792, 'end': 13521, 'strand': 1, 'id': 'nisE', 'gene': 'nisE', 'protein_id': 'ADJ56361.1'},
        {'start': 13507, 'end': 14152, 'strand': 1, 'id': 'nisG', 'gene': 'nisG', 'protein_id': 'ADJ56362.1'},
    ]
    ret = utils.convert_json_orfs(orfs)
    assert ret == expected


def test_jsonify_orfs_error_paths():

    rec = SeqIO.read(get_testfile('results.gbk'), 'genbank')
    bp_orfs = [f for f in rec.features if f.type == 'CDS']
    assert len(bp_orfs) == 11

    cds = bp_orfs[0]
    cds.qualifiers['locus_tag'] = ['fake_tag']
    cds = bp_orfs[1]
    del cds.qualifiers['gene']
    cds.qualifiers['protein_id'] = ['fake_protein']
    cds = bp_orfs[2]
    cds.location = None
    cds = bp_orfs[-2]
    del cds.qualifiers['gene']
    del cds.qualifiers['protein_id']
    cds = bp_orfs[-1]
    del cds.qualifiers['gene']
    del cds.qualifiers['protein_id']

    expected = [
        {'start': 827, 'end': 1001, 'strand': 1, 'id': 'fake_tag', 'gene': 'nisA', 'protein_id': 'ADJ56352.1',
         'locus_tag': 'fake_tag'},
        {'start': 1108, 'end': 4090, 'strand': 1, 'id': 'fake_protein', 'protein_id': 'fake_protein'},
        {'start': 5895, 'end': 7140, 'strand': 1, 'id': 'nisC', 'gene': 'nisC', 'protein_id': 'ADJ56355.1'},
        {'start': 7136, 'end': 7874, 'strand': 1, 'id': 'nisI', 'gene': 'nisI', 'protein_id': 'ADJ56356.1'},
        {'start': 7875, 'end': 9924, 'strand': 1, 'id': 'nisP', 'gene': 'nisP', 'protein_id': 'ADJ56357.1'},
        {'start': 9992, 'end': 10679, 'strand': 1, 'id': 'nisR', 'gene': 'nisR', 'protein_id': 'ADJ56358.1'},
        {'start': 10671, 'end': 12015, 'strand': 1, 'id': 'nisK', 'gene': 'nisK', 'protein_id': 'ADJ56359.1'},
        {'start': 12113, 'end': 12791, 'strand': 1, 'id': 'nisF', 'gene': 'nisF', 'protein_id': 'ADJ56360.1'},
        {'start': 12792, 'end': 13521, 'strand': 1, 'id': 'ORF_000001'},
        {'start': 13507, 'end': 14152, 'strand': 1, 'id': 'ORF_000002'},
    ]
    ret = utils.jsonify_orfs(bp_orfs)
    assert ret == expected
