How to Build Documentation
==========================

Debian/Ubuntu
-------------

To make the html documentation, install the prerequisites::

    apt-get install python-sphinx texlive-latex-recommended dvipng librsvg2-bin imagemagick docbook2x graphviz

and do::

    make html

and to view it, do::

    epiphany _build/html/index.html

Fedora
------

Fedora (and maybe other RPM based distributions), install the prerequisites::

    dnf install python3-sphinx librsvg2 ImageMagick docbook2X texlive-dvipng-bin texlive-scheme-medium librsvg2-tools

After that, run::

    make html

If you get **mpmath** error, install python3-mpmath package::

    dnf install python3-mpmath

And view it at::

    _build/html/index.html
