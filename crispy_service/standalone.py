#!/usr/bin/env python3
"""Scan a sequence for CRISPR PAMs from the command line"""

import argparse
import logging
import json

from Bio import SeqIO

from crispy_service.crispy import crispy_scan
from crispy_service.annotations import json_annotations
from crispy_service import version


def main():
    parser = argparse.ArgumentParser(description='CRISPy scan CLI')
    parser.add_argument('-d', '--debug', dest='debug',
                        action='store_true', default=False,
                        help='print diagnostic messages while running')
    parser.add_argument('-t', '--threads', dest='threads',
                        type=int, default=2,
                        help='Number of threads to use.')
    parser.add_argument('-V', '--version', action="version", version=version.__version__)
    parser.add_argument('--start',
                        type=int, default=0,
                        help="first base coordinate of region to analyse")
    parser.add_argument('--end',
                        type=int, default=-1,
                        help="last base coordinate of region to analyse")
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

    results = crispy_scan([record], target_region, "GG", 13, 20, threads=args.threads)

    json_regions = json_annotations(results, target_region, 7, 13)

    print(json.dumps(json_regions))


if __name__ == '__main__':
    main()
