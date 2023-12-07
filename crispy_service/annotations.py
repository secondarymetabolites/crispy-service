"""CRISPy-web related annotation logic."""

from typing import Dict, List, Union

from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

from . import utils
from .best import (
    AtoG,
    BestEditWindow,
    CtoT,
)


CRISPR_BEST_MODES = [AtoG, CtoT]


def json_annotations(grnas, region, best_size=7, best_offset=13):
    """Return a JSON-compatible version of the target region"""
    json_region = {}

    json_region['name'] = get_name_from_cluster(region)
    json_region['orfs'] = generate_orf_entry(region)
    json_region['grnas'] = extend_grna_entries(grnas, region, best_size, best_offset)

    return json_region


def get_name_from_cluster(region: SeqRecord) -> str:
    """Get the region name from the first cluster's description"""
    name = ''
    bp_clusters = [f for f in region.features if f.type == 'cluster']
    if bp_clusters:
        for note in bp_clusters[0].qualifiers.get('note', []):
            if note.startswith('Cluster number: '):
                name = 'Cluster{}'.format(note.split(':')[-1])
    return name


def generate_orf_entry(region: SeqRecord) -> List[Dict[str, Union[str, int]]]:
    """Create a list of JSONified ORF records"""
    bp_orfs = [f for f in region.features if f.type == 'CDS']
    json_orfs = utils.jsonify_orfs(bp_orfs)
    return json_orfs


def extend_grna_entries(grnas, region, best_size=7, best_offset=13):
    """Extend each record in a list of JSONified gRNA records"""
    result = {}
    idx = 1
    for grna in grnas:
        _id = 'CY{:08d}'.format(idx)
        grna['id'] = _id
        idx += 1

        grna['orf'] = '-'
        grna['can_edit'] = {}
        grna['changed_aas'] = {}
        location = FeatureLocation(grna["start"], grna["end"], grna["strand"])
        for feature in region.features:
            if feature.type == 'CDS' and BestEditWindow.overlaps(feature, location, best_size, best_offset):
                grna['orf'] = utils.get_ident(feature)
                grna['can_edit'] = {}
                for mode in CRISPR_BEST_MODES:
                    if BestEditWindow.can_edit(region, location, mode, best_size, best_offset):
                        edit_window = BestEditWindow(region, feature, location, mode, best_size, best_offset)
                        grna['changed_aas'][mode.name()] = list(map(str, edit_window.get_mutations()))
                        grna['can_edit'][mode.name()] = bool(len(grna['changed_aas'][mode.name()]))
        result[_id] = grna
    return result
