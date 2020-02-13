from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    CompoundLocation,
)
from Bio.SeqRecord import SeqRecord
import logging

_UNIQUE_ID = 1
IDENT_PREFERENCE = ['locus_tag', 'gene', 'protein_id']


def jsonify_orfs(bp_orfs):
    """Convert BioPython ORF records to JSON-compatible data structure"""
    json_orfs = []
    for bp_orf in bp_orfs:
        if bp_orf.location is None:
            logging.warning("Skipping CDS with invalid location %s", bp_orf)
            continue
        orf = {}
        orf['start'] = int(bp_orf.location.start)
        orf['end'] = int(bp_orf.location.end)
        orf['strand'] = int(bp_orf.location.strand)
        orf['id'] = get_ident(bp_orf)
        for key in IDENT_PREFERENCE:
            if key in bp_orf.qualifiers:
                orf[key] = bp_orf.qualifiers[key][0]

        json_orfs.append(orf)

    return json_orfs


def convert_json_orfs(as_json_orfs):
    """Convert antiSMASH JSON ORF records into our format."""
    json_orfs = []
    for as_orf in as_json_orfs:
        orf = {}
        orf['start'], orf['end'], orf['strand'] = parse_serialised_location(as_orf['location'])
        orf['id'] = get_ident(as_orf)
        for key in IDENT_PREFERENCE:
            if key in as_orf['qualifiers']:
                orf[key] = as_orf['qualifiers'][key][0]

        json_orfs.append(orf)
    return json_orfs


def parse_serialised_location(location):
    """Parse the antiSMASH style serialised location."""
    # Possible locations look like
    # [23:42](+) -> 23, 42, 1
    # join{[<17:19](-), [23:>42](-)} -> 17, 42, -1
    # TODO: Replace this by the logic in secmet once we can use AGPL code

    def parse_strand(loc_part):
        strand_char = loc_part[-2]
        if strand_char == '+':
            return 1
        elif strand_char == '-':
            return -1
        else:
            return 0

    def force_number(raw_coord):
        if raw_coord[0] in ('<', '>'):
            raw_coord = raw_coord[1:]
        return int(raw_coord)

    def parse_coords(loc_part):
        raw_coords = loc_part.split('(')[0]
        raw_start, raw_end = raw_coords[1:-1].split(":")
        return force_number(raw_start), force_number(raw_end)

    loc_parts = []

    if location.startswith("join"):
        loc_parts = [part for part in location[5:-1].split(", ")]
    else:
        loc_parts.append(location)

    parsed_locs = []
    part = loc_parts[0]
    start, end = parse_coords(part)
    strand = parse_strand(part)

    for part in loc_parts[1:]:
        new_start, new_end = parse_coords(part)
        new_strand = parse_strand(part)
        assert new_strand == strand

        start = min(start, new_start)
        end = max(end, new_end)

    return start, end, strand


def parse_description_from_json(as_json, region_nr):
    """Describe a region by parsing the best knownclusterblast hit."""
    fallback = "-"
    seq_rec = as_json['records'][0]

    if 'antismash.modules.clusterblast' not in seq_rec['modules']:
        return fallback

    if 'knowncluster' not in seq_rec['modules']['antismash.modules.clusterblast']:
        return fallback

    for res in seq_rec['modules']['antismash.modules.clusterblast']['knowncluster']['results']:
        if res['total_hits'] < 1:
            continue
        if res['region_number'] == region_nr:
            return res['ranking'][0][0]['description']

    return fallback


def get_ident(bp_orf):
    """Get an identifier for an ORF, generating one if needed"""
    global _UNIQUE_ID
    ident = None
    if hasattr(bp_orf, 'qualifiers'):
        qualifiers = bp_orf.qualifiers
    else:
        qualifiers = bp_orf['qualifiers']
    for key in IDENT_PREFERENCE:
        if key in qualifiers:
            ident = qualifiers[key][0]
            break
    if ident is None:
        ident = "ORF_{:06}".format(_UNIQUE_ID)
        _UNIQUE_ID += 1

    return ident


def json_to_gbk(as_json):
    """Generate a GenBank record from an antiSMASH JSON file."""
    as_rec = as_json['records'][0]
    seq = Seq(as_rec['seq']['data'], generic_dna)
    features = []
    for json_feature in as_rec['features']:
        start, end, strand = parse_serialised_location(json_feature['location'])
        loc = FeatureLocation(start, end, strand)
        feature = SeqFeature(loc, type=json_feature['type'], qualifiers=json_feature['qualifiers'])
        features.append(feature)

    record = SeqRecord(seq, id=as_rec['id'], name=as_rec['name'], description=as_rec['description'],
                       features=features, annotations=as_rec['annotations'], dbxrefs=as_rec['dbxrefs'])
    return record


def locations_overlap(first, second) -> bool:
    """Returns True if the two provided locations overlap"""
    if isinstance(first, CompoundLocation):
        return any(locations_overlap(part, second) for part in first.parts)
    if isinstance(second, CompoundLocation):
        return any(locations_overlap(first, part) for part in second.parts)
    return first.end > second.start and second.end > first.start
