# CRISPy-web backend service

This repository contains the backend service logic for [CRISPy-web](https://crispy.secondarymetabolites.org/).

The same logic can be used in a local program `crispy-standalone` also installable from this repository.

## Installation

Clone the sources from GitHub:

```
git clone https://github.com/secondarymetabolites/crispy-service
cd crispy-service
```

Ideally, create and load a python virtual environment. E.g.

```
python3 -m venv env
source env/bin/activate
```

(or use virtualenvmanager or conda or the like).

**Note:** Supported Python versions are 3.9, 3.10, and 3.11

Then, install crispy-standalone and all its dependencies:

```
pip install .
```

## Usage

```
crispy-standalone NC_003888.3.gbk --start 246525 --end 271084 --cbest-mode CtoT --stop-only > sco_cluster_3_cbest_stops.csv
```

See `crispy-standalone --help` for all possible options.

**Note:** `crispy-standalone` output is tab-separated, but using `.csv` as the file ending is recommended if you plan on opening the files in Excel or other spreadsheet tools.

## License

This code is released under the GNU Affero General Public license.
See [`LICENSE.txt`](LICENSE.txt) for details.
