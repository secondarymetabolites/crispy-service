#!/usr/bin/env python3
"""Scan a sequence for CRISPR PAMs from the command line"""

import argparse
import logging
import multiprocessing

from Bio import SeqIO

from crispy_service.crispy import crispy_scan
from crispy_service.annotations import extend_grna_entries
from crispy_service import version


def main():
    parser = argparse.ArgumentParser(description='CRISPy scan CLI')
    parser.add_argument('-d', '--debug', dest='debug',
                        action='store_true', default=False,
                        help='print diagnostic messages while running')
    parser.add_argument('-t', '--threads', dest='threads',
                        type=int, default=multiprocessing.cpu_count(),
                        help='Number of threads to use. (default: %(default)s)')
    parser.add_argument('-V', '--version', action="version", version=version.__version__)
    parser.add_argument('--start',
                        type=int, default=0,
                        help="first base coordinate of region to analyse")
    parser.add_argument('--end',
                        type=int, default=-1,
                        help="last base coordinate of region to analyse")
    parser.add_argument("--pam", dest="pam",
                        default="GG",
                        help="PAM sequence to use. (default: %(default)s)")
    parser.add_argument("--unique-size", dest="unique_size",
                        type=int, default=13,
                        help="Spacer bps that need to be unique. (default: %(default)s)")
    parser.add_argument("--full-size", dest="full_size",
                        type=int, default=23,
                        help="Full size of the spacer in bp. (default: %(default)s)")
    parser.add_argument("--cbest-size", dest="best_size",
                        type=int, default=7,
                        help="Size of the cBEST window. (default: %(default)s)")
    parser.add_argument("--cbest-offset", dest="best_offset",
                        type=int, default=13,
                        help="Distance from PAM that the cBEST window starts at. (default: %(default)s)")
    parser.add_argument("--cbest-mode",
                        choices=("none", "CtoT", "AtoG"), default="none",
                        help="cBEST mode to run in. (default: %(default)s)")
    parser.add_argument("--stop-only",
                        action="store_true", default=False,
                        help="Show only cBEST mutations leading to STOP codons.")
    parser.add_argument('filename',
                        help="Filename to analyse")
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(format='%(levelname)s %(asctime)s: %(message)s', level=logging.INFO)

    # TODO: Once the rest of the code handles record lists, use the full list.
    record = list(SeqIO.parse(args.filename, 'genbank'))[0]
    if args.end == -1:
        args.end = len(record)
    target_region = record[args.start:args.end]

    results = crispy_scan(
        [record],
        target_region,
        args.pam,
        args.unique_size,
        args.full_size,
        threads=args.threads)

    grnas = list(extend_grna_entries(results, target_region, args.best_size).values())

    print("Start\tEnd\tStrand\tORF\tSequence\tPAM\tmutations\t0bp mismatches\t1bp mismatches\t2bp mismatches")
    for rna in grnas:

        if args.cbest_mode != "none":
            if not rna["can_edit"].get(args.cbest_mode):
                continue
            if args.stop_only:
                found = False
                for changed_aa in rna["changed_aas"][args.cbest_mode]:
                    if changed_aa[-1] == "*":
                        found = True

                if not found:
                    continue

        print(
            rna["id"],
            rna["start"],
            rna["end"],
            rna["strand"],
            rna["orf"],
            rna["sequence"],
            rna["pam"],
            ", ".join(rna["changed_aas"].get(args.cbest_mode, [])),
            rna["0bpmm"],
            rna["1bpmm"],
            rna["2bpmm"],
            sep="\t")


if __name__ == '__main__':
    main()
