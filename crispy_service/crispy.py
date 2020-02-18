from typing import Dict, List, Tuple, Union

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from nearmiss import Searcher


def build_comparison_text(records: List[SeqRecord], window_size: int) -> str:
    """ Builds a single string containing forward and reverse strands of all
        record sequences, separated by a large enough section of '$' characters
        to prevent the uniqueness window from matching any of the neighbouring
        sequences.

        Arguments:
            records: a dictionary of record id to SeqRecord instance
            window_size: the size of the unique window

        Returns:
            a string of the combined sequences
    """
    seqs = []
    separator = "$" * window_size
    for record in records:
        seqs.append(str(record.seq))
        seqs.append(str(record.seq.reverse_complement()))
    return separator.join(seqs)


def crispy_scan(haystack: List[SeqRecord], needle: SeqRecord, pam: str = "GG",
                unique_size: int = 13, full_size: int = 23, threads: int = -1,
                ) -> List[Tuple[int, List[int]]]:
    if unique_size < 1:
        raise ValueError("unique size cannot be below 1")
    if full_size < unique_size:
        raise ValueError("full size cannot be below unique size")

    def build_json_base(location, seq_section, result) -> Dict[str, Union[str, int]]:
        base = {
            'start': location.start,
            'end': location.end,
            'strand': location.strand,
            'sequence': str(seq_section[:-3]),
            'pam': str(seq_section[-3:]),
            'all_hits': result,  # new to JSON, for handy sorting
            '0bpmm': result[0] - 1,  # remove self-hit
        }
        # add remaining mismatch info
        for i, val in enumerate(result[1:]):
            base['{}bpmm'.format(i+1)] = val
        return base

    # set the size of the window to the unique size
    # and shift one back since in the previous system it skipped a leading N
    before_window = (-unique_size - 1, -1)

    final_result = []

    comparison_text = build_comparison_text(haystack, unique_size)

    idx = 0
    for strand in [1, -1]:
        if strand == -1:
            searcher = Searcher(str(needle.seq.reverse_complement()))
        else:
            searcher = Searcher(str(needle.seq))
        results = searcher.find_repeat_counts(target=pam, before_window=before_window,
                                              other_text=comparison_text, threads=threads)

        for pam_start, result in sorted(results.items(), key=lambda x: x[1]):
            # set the window location, accounting for strand
            if strand == -1:
                start = len(needle.seq) - pam_start - len(pam)
                end = start + full_size
            else:
                start = pam_start - full_size + len(pam)
                end = pam_start + len(pam)
            # skip anything for which the full window shown would be truncated
            if start < 0 or end >= len(needle.seq):
                continue

            location = FeatureLocation(start, end, strand)
            seq = location.extract(needle.seq)
            final_result.append(build_json_base(location, seq, result))
            idx += 1

    # order by lowest hits, then by start position
    final_result.sort(key=lambda x: (x["all_hits"], x["start"]))

    return final_result
