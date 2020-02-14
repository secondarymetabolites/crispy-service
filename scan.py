#!/usr/bin/env python3
"""Scan a sequence for CRISPR PAMs"""

import redis
import time
import argparse
import logging
from os import path
from Bio import SeqIO
from crispy_models import Queue, Session
from crispy_service.crispy import crispy_scan
from crispy_service.annotations import json_annotations
from crispy_service import version


def main():
    parser = argparse.ArgumentParser(description='CRISPy scan service')
    parser.add_argument('-d', '--debug', dest='debug',
                        action='store_true', default=False,
                        help='print diagnostic messages while running')
    parser.add_argument('-q', '--queue', dest='queue',
                        default='redis://127.0.0.1:6379/0',
                        help='URI of the Queue server')
    parser.add_argument('-u', '--upload-dir', dest='upload_dir',
                        default='../uploads',
                        help='Path to directory where the session files are stored')
    parser.add_argument('-t', '--threads', dest='threads',
                        type=int, default=2,
                        help='Number of threads to use.')
    parser.add_argument('-V', '--version', action="version", version=version.__version__)
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(format='%(levelname)s %(asctime)s: %(message)s', level=logging.INFO)

    UPLOAD_PATH = path.abspath(args.upload_dir)

    db = redis.Redis.from_url(args.queue, decode_responses=True)
    queue = Queue(db, 'scan')
    while True:
        job_key = queue.next()
        if job_key is None:
            time.sleep(5)
            continue

        logging.info('processing {}'.format(job_key))
        job_id = int(job_key.split(':')[-1])
        job = Session(db, session_id=job_id)
        dirname = path.join(UPLOAD_PATH, str(job_id))

        try:
            record = SeqIO.parse(path.join(dirname, job.filename), 'genbank')
            target_region = record[0][job.from_coord:job.to_coord]

            results = crispy_scan(record, target_region, args.threads, job.pam, job.uniq_size, job.full_size)

            json_region = json_annotations(record, results, job.best_size, job.best_offset)
            job.region = json_region
            job.state = 'done'
            logging.info('done with {}'.format(job_key))
        except ValueError as e:
            job.state = 'error'
            job.error = e.message
            logging.info('job {} failed: {!r}'.format(job_key, e.message))
        except Exception as err:
            job.state = 'error'
            job.error = str(err)
            raise


if __name__ == '__main__':
    main()
