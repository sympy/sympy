#!/usr/bin/env bash
#
# Install non-Python dependencies needed to build the docs on Ubuntu. A
# requirements.txt file is provided to install the Python dependencies.
#
# From a fresh install of Ubuntu:
#
#    $ doc/aptinstall.sh
#    $ pip install -r doc/requirements.txt
#
# Then you should be set to build the docs.

set -o errexit

sudo apt-get update
sudo apt-get install\
  texlive-latex-recommended\
  texlive-fonts-recommended\
  texlive-fonts-extra\
  texlive-xetex\
  latexmk\
  dvipng\
  librsvg2-bin\
  imagemagick\
  chromium-browser\
  libcanberra-gtk-module\
  docbook2x\
  graphviz\
  #

sudo dot -c
