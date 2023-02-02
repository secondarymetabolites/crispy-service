from Bio import SeqIO
from os import path

from crispy_service import utils
from crispy_service import annotations as ann


# TODO: Move to a shared location
def get_testfile(filename):
    return path.join(path.dirname(__file__), filename)


def test_json_annotations():
    rec = SeqIO.read(get_testfile('results.gbk'), 'genbank')
    bp_orfs = [f for f in rec.features if f.type == 'CDS']
    orfs = utils.jsonify_orfs(bp_orfs)
    scan_results = [
        {
            'start': 25, 'end': 48, 'strand': 1,
            'sequence': 'CGTTGACTTTCCATTATGCT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        {
            'start': 95, 'end': 118, 'strand': 1,
            'sequence': 'TAACAAATAATTGGCATATC', 'pam': 'GGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        {
            'start': 829, 'end': 852, 'strand': 1,
            'sequence': 'GAGTACAAAAGATTTTAACT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        {
            'start': 1165, 'end': 1188, 'strand': -1,
            'sequence': 'GTAAAACTCCGTTTATCGTT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        {
            'start': 86, 'end': 109, 'strand': -1,
            'sequence': 'CCAATTATTTGTTAAGAATA', 'pam': 'CGA',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 1,
        },
    ]
    expected = {
        'name': 'Cluster 1',
        'orfs': orfs,
        'grnas': {
            'CY00000001': {
                'id': 'CY00000001',
                'start': 25, 'end': 48, 'strand': 1,
                'orf': '-',
                'can_edit': False,
                'changed_aas': {},
                'sequence': 'CGTTGACTTTCCATTATGCT', 'pam': 'TGG',
                '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
            },
            'CY00000002': {
                'id': 'CY00000002',
                'start': 95, 'end': 118, 'strand': 1,
                'orf': '-',
                'can_edit': False,
                'changed_aas': {},
                'sequence': 'TAACAAATAATTGGCATATC', 'pam': 'GGG',
                '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
            },
            'CY00000003': {
                'id': 'CY00000003',
                'start': 829, 'end': 852, 'strand': 1,
                'orf': 'nisA',
                'can_edit': {'AtoG': True, 'CtoT': True},
                'changed_aas': {'AtoG': ['T3A', 'K4G'], 'CtoT': ['T3I']},
                'sequence': 'GAGTACAAAAGATTTTAACT', 'pam': 'TGG',
                '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
            },
            'CY00000004': {
                'id': 'CY00000004',
                'start': 1165, 'end': 1188, 'strand': -1,
                'orf': 'nisB',
                'can_edit': {'AtoG': True, 'CtoT': True},
                'changed_aas': {'AtoG': ['F26P'], 'CtoT': ['R24Q', 'S25N']},
                'sequence': 'GTAAAACTCCGTTTATCGTT', 'pam': 'TGG',
                '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
            },
            'CY00000005': {
                'id': 'CY00000005',
                'start': 86, 'end': 109, 'strand': -1,
                'orf': '-',
                'can_edit': False,
                'changed_aas': {},
                'sequence': 'CCAATTATTTGTTAAGAATA', 'pam': 'CGA',
                '0bpmm': 0, '1bpmm': 0, '2bpmm': 1,
            },
        },
    }

    ret = ann.json_annotations(scan_results, rec)
    assert ret == expected


def test_get_name_from_cluster():
    rec = SeqIO.read(get_testfile('results.gbk'), 'genbank')
    name = ann.get_name_from_cluster(rec)
    assert name == "Cluster 1"


def test_generate_orf_entry():
    rec = SeqIO.read(get_testfile('results.gbk'), 'genbank')
    bp_orfs = [f for f in rec.features if f.type == 'CDS']
    assert len(bp_orfs) == 11
    expected = utils.jsonify_orfs(bp_orfs)
    ret = ann.generate_orf_entry(rec)
    assert ret == expected


def test_generate_grna_entries():
    rec = SeqIO.read(get_testfile('results.gbk'), 'genbank')
    grnas = [
        {
            'start': 25, 'end': 48, 'strand': 1,
            'sequence': 'CGTTGACTTTCCATTATGCT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        {
            'start': 95, 'end': 118, 'strand': 1,
            'sequence': 'TAACAAATAATTGGCATATC', 'pam': 'GGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        {
            'start': 829, 'end': 852, 'strand': 1,
            'sequence': 'GAGTACAAAAGATTTTAACT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        {
            'start': 1165, 'end': 1188, 'strand': -1,
            'sequence': 'GTAAAACTCCGTTTATCGTT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        {
            'start': 86, 'end': 109, 'strand': -1,
            'sequence': 'CCAATTATTTGTTAAGAATA', 'pam': 'CGA',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 1,
        },
    ]
    expected = {
        'CY00000001': {
            'id': 'CY00000001',
            'start': 25, 'end': 48, 'strand': 1,
            'orf': '-',
            'can_edit': False,
            'changed_aas': {},
            'sequence': 'CGTTGACTTTCCATTATGCT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        'CY00000002': {
            'id': 'CY00000002',
            'start': 95, 'end': 118, 'strand': 1,
            'orf': '-',
            'can_edit': False,
            'changed_aas': {},
            'sequence': 'TAACAAATAATTGGCATATC', 'pam': 'GGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        'CY00000003': {
            'id': 'CY00000003',
            'start': 829, 'end': 852, 'strand': 1,
            'orf': 'nisA',
            'can_edit': {'AtoG': True, 'CtoT': True},
            'changed_aas': {'AtoG': ['T3A', 'K4G'], 'CtoT': ['T3I']},
            'sequence': 'GAGTACAAAAGATTTTAACT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        'CY00000004': {
            'id': 'CY00000004',
            'start': 1165, 'end': 1188, 'strand': -1,
            'orf': 'nisB',
            'can_edit': {'AtoG': True, 'CtoT': True},
            'changed_aas': {'AtoG': ['F26P'], 'CtoT': ['R24Q', 'S25N']},
            'sequence': 'GTAAAACTCCGTTTATCGTT', 'pam': 'TGG',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 0,
        },
        'CY00000005': {
            'id': 'CY00000005',
            'start': 86, 'end': 109, 'strand': -1,
            'orf': '-',
            'can_edit': False,
            'changed_aas': {},
            'sequence': 'CCAATTATTTGTTAAGAATA', 'pam': 'CGA',
            '0bpmm': 0, '1bpmm': 0, '2bpmm': 1,
        },
    }
    ret = ann.extend_grna_entries(grnas, rec)
    assert ret == expected
