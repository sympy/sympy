#!/usr/bin/env bash

set -o errexit

sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update

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
	gfortran\
	#

sudo apt install\
	python3.6\
	python3.6-venv\
	python3.7\
	python3.7-venv\
	python3.9\
	python3.9-venv\
	# python3.8 is installed by default in 20.04

python3 -m venv release/venv_main
source release/venv_main/bin/activate

pip install -U pip wheel
pip install -r release/requirements.txt
pip install -e .
