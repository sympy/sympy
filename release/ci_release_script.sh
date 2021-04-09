#!/usr/bin/env bash

release/aptinstall.sh

python3 -m venv release/venv_main
source release/venv_main/bin/activate

pip install -U pip wheel
pip install -r release/requirements.txt
pip install -e .

python release/releasecheck.py $1 release-$1
