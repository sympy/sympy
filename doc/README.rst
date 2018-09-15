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

Windows 10
----------

Some dependencies like dvipng or docbook2x are not available for windows.

For Windows 10, however, the Windows Subsystem for Linux is available, and you
can install Ubuntu shell on top of it after following up the tutorial below.

    https://github.com/MicrosoftDocs/WSL/blob/live/WSL/install-win10.md

In your command prompt, simply run `ubuntu` to transfer to linux terminal,
and follow up the Debian/Ubuntu tutorial above to install all the dependencies,
and then you can run `make html` to build.
(Note that you would also have to install `make` via `apt-get install make`)

If you want to change directory in your prompt to your working folder of sympy
in windows file system, you can prepend `cd /mnt/` to your file path in windows,
and run in your shell to navigate to the folder.
(Also note that linux uses `/` instead of `\` for path)

This method would provide better compatability than cygwin or msys2,
and more convenience than a virtual machine, if you partially need linux
environment for your work.

However this method is only viable for Windows 10 64-bit users.
