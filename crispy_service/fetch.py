import json
import os
import re
import requests
from os import path
from io import StringIO
from Bio import SeqIO
from . import utils

BASE_URL = 'https://antismash.secondarymetabolites.org/upload/{id}'

def grab_index(asID):
    """Get the index.html file for an antiSMASH ID"""
    index_url = '{}/index.html'.format(BASE_URL.format(id=asID))
    index = requests.get(index_url)
    index.raise_for_status()
    return index.text


def find_dl_link(text):
    """Find the correct download link in the index.html page text"""
    pattern = re.compile(r'href=["]?(?P<link>[a-zA-Z0-9_.-]+)["]?>[ ]?Download GenBank summary file')
    match = pattern.search(text)
    link = match.group('link')
    if link.endswith('final.gbk'):
        return link
    return link[:-3] + 'json'


def download_gbk(asID, link):
    """Download a GenBank file for a given antiSMASH id and filename"""
    download_url = '{}/{}'.format(BASE_URL.format(id=asID), link)
    req = requests.get(download_url)
    req.raise_for_status()
    return req.text


def genome_json(seq_rec):
    """generate a JSON datastructure from a SeqRecord"""
    # Support both as5 and as4 annotations
    bp_regions = [f for f in seq_rec.features if f.type == 'region' or f.type == 'cluster']
    bp_orfs = [f for f in seq_rec.features if f.type == 'CDS']

    json_rec = {}

    json_rec['length'] = len(seq_rec)
    json_rec['id'] = seq_rec.id
    json_rec['description'] = seq_rec.description
    if 'organism' in seq_rec.annotations:
        organism = seq_rec.annotations['organism']
    else:
        organism = "Unknown organism"
    json_rec['organism'] = organism

    json_regions = []
    region_counter = 1
    for bp_region in bp_regions:
        region = {}
        region['start'] = int(bp_region.location.start)
        region['end'] = int(bp_region.location.end)
        region['type'] = ", ".join(bp_region.qualifiers['product'])

        if 'knownclusterblast' in bp_region.qualifiers:
            description = bp_region.qualifiers['knownclusterblast'][0].split('\t')[-1]
        else:
            description = '-'
        region['description'] = description

        region['name'] = "Region {}".format(region_counter)
        region_counter += 1

        json_regions.append(region)
    json_rec['clusters'] = json_regions  # TODO: Rename this in the JSON as well

    json_rec['orfs'] = utils.jsonify_orfs(bp_orfs)

    return json_rec


def convert_json(as_json):
    """Generate JSON datastructures from an antiSMASH JSON file."""
    seq_rec = as_json['records'][0]
    as_regions = [f for f in seq_rec['features'] if f['type'] == 'region']
    as_orfs = [f for f in seq_rec['features'] if f['type'] == 'CDS']

    json_rec = {}

    json_rec['length'] = len(seq_rec['seq']['data'])
    json_rec['id'] = seq_rec['id']
    json_rec['description'] = seq_rec['description']
    if 'organism' in seq_rec['annotations']:
        organism = seq_rec['annotations']['organism']
    else:
        organism = "Unknown organism"
    json_rec['organism'] = organism

    json_regions = []
    region_counter = 1
    for as_region in as_regions:
        region = {}
        region['start'], region['end'], _ = utils.parse_serialised_location(as_region['location'])
        region['type'] = ", ".join(as_region['qualifiers']['product'])

        region['description'] = utils.parse_description_from_json(as_json, region_counter)

        region['name'] = "Region {}".format(region_counter)
        region_counter += 1

        json_regions.append(region)

    json_rec['clusters'] = json_regions  # TODO: Rename this in the JSON as well

    json_rec['orfs'] = utils.convert_json_orfs(as_orfs)

    return json_rec


def fetch(asID, dirname):
    """Fetch a file from antiSMASH, parse it, save it and return JSON datastructure"""
    index = grab_index(asID)
    link = find_dl_link(index)
    if link.endswith(".gbk"):
        gbk_text = download_gbk(asID, link)
        gbk_rec = list(SeqIO.parse(StringIO(gbk_text), 'genbank'))[0]
        json_rec = genome_json(gbk_rec)
    else:
        as_json = json.loads(download_gbk(asID, link))
        gbk_rec = utils.json_to_gbk(as_json)
        json_rec = convert_json(as_json)

    if not path.exists(dirname):
        os.mkdir(dirname)
    filename = 'input.gbk'
    SeqIO.write([gbk_rec], path.join(dirname, filename), 'genbank')

    return json_rec, filename
