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

browser_packages=(chromium-browser)
for browser in google-chrome chrome chromium chromium-browser; do
  if command -v "$browser" >/dev/null 2>&1 && "$browser" --version >/dev/null 2>&1; then
    browser_packages=()
    echo "Using existing browser: $browser"
    break
  fi
done

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
  "${browser_packages[@]}"\
  docbook2x\
  graphviz\
  #

sudo dot -c
