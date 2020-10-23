#!/usr/bin/env bash

set -o errexit

sudo apt install\
	antlr4\
	libgfortran5\
	python3-venv\
	python3-pip\
	python3-gmpy2\
	texlive-latex-recommended\
	texlive-fonts-recommended\
	texlive-fonts-extra\
	texlive-xetex\
	latexmk\
       	dvipng\
       	librsvg2-bin\
   	imagemagick\
	inkscape\
	libcanberra-gtk-module\
       	docbook2x\
	graphviz\
	#

python3 -m venv release/venv_main
source release/venv_main/bin/activate

pip install -U pip wheel
pip install -r release/requirements.txt
pip install -e .
