#!/usr/bin/env bash
#
# Run
#
#   $ release/ci_release_script.sh version prevversion releasedir
#
# Example:
#
#   $ release/ci_release_script.sh 1.9rc1 1.8 release-1.9rc1

release/aptinstall.sh

python3 -m venv release/venv_main
source release/venv_main/bin/activate

pip install -U pip wheel
pip install -r release/requirements.txt
pip install -e .

python release/releasecheck.py $1 $2 release-$1
