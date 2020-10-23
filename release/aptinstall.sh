#!/usr/bin/env bash

sudo apt install\
	antlr4\
	libgfortran5\
	python3-venv\
	python3-pip\
	python3-gmpy2\
	#

python3 -m venv venv_main
. venv_main/bin/activate

pip install -U pip
pip install -r requirements.txt
