#!/usr/bin/env bash
#
# Run
#
#   $ release/ci_release_script.sh version prevversion
#
# Example:
#
#   $ release/ci_release_script.sh 1.9rc1 1.8

doc/aptinstall.sh

pip install --upgrade setuptools pip wheel build
pip install -r doc/requirements.txt
pip install .

python release/releasecheck.py $1 $2 release-$1
