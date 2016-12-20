How to Build Documentation
==========================

To make the html documentation, install the prerequisites, e.g. on
Debian/Ubuntu (similarly for other distributions)::

    apt-get install python-sphinx texlive-latex-recommended dvipng librsvg2-bin imagemagick docbook2x

and do::

    make html

and to view it, do::

    epiphany _build/html/index.html
