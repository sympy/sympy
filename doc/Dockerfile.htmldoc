###############################################################################
# This Dockerfile can be used to build an image where the HTML version of the
# docs for SymPy can be built.
#
# BASIC USAGE
# ===========
#
# If SYMPY_ROOT is the directory where the sympy repository lives, then you can
# build the image with
#
#   $ cd SYMPY_ROOT/doc
#   $ docker build -f Dockerfile.htmldoc -t sympy_htmldoc .
#
# Once the image is built, you can build the docs at any time (and from any
# directory) using
#
#   $ docker run --rm -v /absolute/path/to/SYMPY_ROOT:/sympy sympy_htmldoc
#
# (substitute the actual absolute filesystem path to SYMPY_ROOT).
#
# The documentation will be built in SYMPY_ROOT's doc/_build/html directory.
#
# LIVEHTML SERVER
# ===============
#
# Alternatively, you can use the image to run the "livehtml" server. For
# example, use
#
#   $ docker run --rm -it \
#       -v /absolute/path/to/SYMPY_ROOT:/sympy \
#       -p 8000:80 \
#       sympy_htmldoc live
#
# and then navigate your browser to localhost:8000. You can use a different
# port by changing the 8000 in the command.
#
# This will automatically detect changes in the doc sources, rebuild, and
# update the page in the browser.
#
# When finished, you can stop the server with ctrl-c in the terminal.
#
# If you want to instead run the server in detached mode, you can use
#
#   $ docker run --rm -d --name=sympy-livehtml \
#       -v /absolute/path/to/SYMPY_ROOT:/sympy \
#       -p 8000:80 \
#       sympy_htmldoc live
#
# and then
#
#   $ docker stop sympy-livehtml
#
# when you are finished.
#
###############################################################################

FROM python:3.8.12-slim-buster

RUN apt-get update
RUN apt-get install -y make librsvg2-bin imagemagick graphviz git

COPY requirements.txt /tmp
RUN pip install --upgrade pip
RUN pip install -r /tmp/requirements.txt

RUN echo '#!/bin/bash \n\
case $1 in \n\
  live) \n\
    cd /sympy/doc; make livehtml LIVEHOST=0.0.0.0 LIVEPORT=80 \n\
    ;; \n\
  *) \n\
    cd /sympy/doc; make html \n\
    ;; \n\
esac \n\
' > /usr/local/bin/makehtml.sh
RUN chmod +x /usr/local/bin/makehtml.sh

ENTRYPOINT ["makehtml.sh"]
