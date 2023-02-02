# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
""" Clean up old CRISPy-web jobs """

from argparse import ArgumentParser
import logging
from datetime import datetime, timedelta
from pathlib import Path
import shutil

from crispy_models import Session
import redis


def main():
    """ command line wrangling """
    parser = ArgumentParser()
    parser.add_argument("-a", "--age", dest="age", type=int, default=7,
                        help="Maximum age of job in days (default: %(default)s).")
    parser.add_argument("-d", "--debug", dest="debug",
                        action="store_true", default=False,
                        help="print diagnostic messages while running")
    parser.add_argument("-q", "--queue", dest="queue",
                        default="redis://127.0.0.1:6379/0",
                        help="URI of the Queue server")
    parser.add_argument("-u", "--upload-dir", dest="upload_dir",
                        default="../uploads",
                        help="Path to directory where the session files are stored")
    args = parser.parse_args()

    setup_logging(args.debug)
    age = timedelta(days=args.age)
    run(args.queue, Path(args.upload_dir).resolve(), age)


def run(queue: str, upload_dir: Path, age: timedelta):
    """ Clean up old CRISPy-web jobs """
    db = redis.Redis.from_url(queue, encoding="utf-8", decode_responses=True)
    jobs = db.keys("crispy:session:*")
    now = datetime.now()
    for key in jobs:
        job_id = int(key.split(":")[-1])
        job = Session(db, session_id=job_id)
        if not job.last_changed_datetime < now - age:
            continue
        logging.debug("%s => %s, older than %s deleting", key, job.last_changed, age)
        delete_job_files(job_id, upload_dir)
        db.delete(key)


def delete_job_files(job_id: int, upload_dir: Path):
    """ Delete the job's files if they exist """

    jobdir = upload_dir / str(job_id)
    if not jobdir.exists():
        logging.debug("No directory %s, skipping", jobdir)
        return

    shutil.rmtree(jobdir)


def setup_logging(debug: bool):
    """ Set up logging """
    log_level = logging.WARNING
    if debug:
        log_level = logging.DEBUG

    log_format = "%(levelname)-8s %(asctime)s   %(message)s"
    date_format = "%d/%m %H:%M:%S"
    logging.basicConfig(format=log_format, datefmt=date_format)
    logger = logging.getLogger()
    logger.setLevel(log_level)


if __name__ == "__main__":
    main()
