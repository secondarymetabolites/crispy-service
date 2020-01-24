"""CRISPy-web related annotation logic."""
from . import utils
from .best import (
    AtoG,
    BestEditWindow,
    CtoT,
)
from Bio.SeqFeature import ExactPosition, FeatureLocation, SeqFeature


CRISPR_BEST_MODES = [AtoG, CtoT]


def annotate_grnas(feats, target_region):
    for feat_list in feats:
        for feat in feat_list:
            zero_bp_mismatches, one_bp_mismatches, two_bp_mismatches = feat[3][:3]

            # find_offtarg reports the self-match as well, so number of true off-targets is one less
            zero_bp_mismatches -= 1

            target_region.features.append(
                SeqFeature(FeatureLocation(ExactPosition(feat[0]), ExactPosition(feat[1])),
                           type="gRNA", strand=feat[2], qualifiers={
                    "0bpmm": [str(zero_bp_mismatches)],
                    "1bpmm": [str(one_bp_mismatches)],
                    "2bpmm": [str(two_bp_mismatches)],
                })
            )

    return target_region


def json_annotations(grnas, region, best_size=7, best_offset=13):
    """Return a JSON-compatible version of the target region"""
    json_region = {}

    json_region['name'] = get_name_from_cluster(region)
    json_region['orfs'] = generate_orf_entry(region)
    json_region['grnas'] = extend_grna_entries(grnas, region, best_size, best_offset)

    return json_region


def get_name_from_cluster(region):
    """Get the region name from the first cluster's description"""
    name = ''
    bp_clusters = [f for f in region.features if f.type == 'cluster']
    if bp_clusters:
        for note in bp_clusters[0].qualifiers.get('note', []):
            if note.startswith('Cluster number: '):
                name = 'Cluster{}'.format(note.split(':')[-1])
    return name


def generate_orf_entry(region):
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
        grna['can_edit'] = False
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
