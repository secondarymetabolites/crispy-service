#!/usr/bin/env python3
"""Prepare sequences for CRISPy processing"""

import json
import redis
import argparse
import logging
import requests
import time
from os import path
from Bio import SeqIO
from crispy_models import Queue, Session
from crispy_service.fetch import (
    convert_json,
    fetch,
    genome_json,
)
from crispy_service.utils import json_to_gbk
from crispy_service import version


def main():  # noqa: C901  # no, flake8, it's not too complex
    parser = argparse.ArgumentParser(description='CRISPy sequence preparation service')
    parser.add_argument('-d', '--debug', dest='debug',
                        action='store_true', default=False,
                        help='print diagnostic messages while running')
    parser.add_argument('-q', '--queue', dest='queue',
                        default='redis://127.0.0.1:6379/0',
                        help='URI of the Queue server')
    parser.add_argument('-u', '--upload-dir', dest='upload_dir',
                        default='../uploads',
                        help='Path to directory where the session files are stored')
    parser.add_argument('-V', '--version', action="version", version=version.__version__)
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(format='%(levelname)s %(asctime)s: %(message)s', level=logging.INFO)

    UPLOAD_PATH = path.abspath(args.upload_dir)

    db = redis.Redis.from_url(args.queue, decode_responses=True)
    queue = Queue(db, 'prepare')
    while True:
        job_key = queue.next()
        if job_key is None:
            time.sleep(1)
            continue
        logging.info('processing {}'.format(job_key))
        job_id = int(job_key.split(':')[-1])
        job = Session(db, session_id=job_id)
        dirname = path.join(UPLOAD_PATH, str(job_id))
        if job.asid:
            logging.info('grabbing file for id {}'.format(job.asid))
            try:
                json_rec, filename = fetch(job.asid, dirname)
                job.filename = filename
            except requests.HTTPError:
                job.state = 'error'
                job.error = 'Downloading sequence file failed. Is the ID spelled correctly?'
                continue
            except Exception:
                job.state = 'error'
                job.error = 'Failed to load the sequence'
                continue
        else:
            logging.info('parsing uploaded file {}'.format(job.filename))
            full_path = path.join(dirname, job.filename)
            if job.filename.endswith('.gbk') or job.filename.endswith('.gb'):
                try:
                    rec = list(SeqIO.parse(full_path, 'genbank'))[0]
                except (IndexError, ValueError, AssertionError, AttributeError, UnicodeDecodeError):
                    job.state = 'error'
                    job.error = 'Parsing the uploaded file failed. Please check it is a valid GenBank file'
                    continue

                try:
                    json_rec = genome_json(rec)
                except (TypeError):
                    job.state = 'error'
                    job.error = 'Invalid input file. Please check the file contains a nucleotide sequence.'
                    continue
            elif job.filename.endswith('.json'):
                with open(full_path, 'r') as handle:
                    as_json = json.load(handle)
                json_rec = convert_json(as_json)
                gbk_record = json_to_gbk(as_json)
                job.filename = 'input.gbk'
                SeqIO.write([gbk_record], path.join(dirname, job.filename), 'genbank')
            else:
                job.state = 'error'
                job.error = 'Invalid input file'
                continue

        job.genome = json_rec
        job.state = 'loaded'

        logging.info('done with {}'.format(job_key))


if __name__ == '__main__':
    main()
