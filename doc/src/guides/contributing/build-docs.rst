==========================
Build the Documentation
==========================


Start by installing the required dependencies for the documentation.

Required dependencies
^^^^^^^^^^^^^^^^^^^^^^

Debian/Ubuntu
~~~~~~~~~~~~~~~

For Debian/Ubuntu::

   apt-get install python3-sphinx texlive-latex-recommended dvipng librsvg2-bin imagemagick docbook2x graphviz
   python -m pip install sphinx-math-dollar sphinx-reredirects myst-parser linkify-it-py

If you get mpmath error, install python-mpmath package::

   apt-get install python-mpmath

If you get matplotlib error, install python-matplotlib package::

   apt-get install python-matplotlib

Fedora
~~~~~~~~

For Fedora (and maybe other RPM-based distributions), install the
prerequisites::

   dnf install python3-sphinx librsvg2 ImageMagick docbook2X texlive-dvipng-bin
   texlive-scheme-medium librsvg2-tools
   python -m pip install sphinx-math-dollar sphinx-reredirects myst-parser linkify-it-py

If you get mpmath error, install python3-mpmath package::

   dnf install python3-mpmath

If you get matplotlib error, install python3-matplotlib package::

   dnf install python3-matplotlib

Mac
~~~~

For Mac, first install homebrew: https://brew.sh/

Then install these packages with homebrew::

   brew install imagemagick graphviz docbook librsvg

Install these packages with either pip or conda::

   python -m pip install mpmath matplotlib sphinx sphinx-math-dollar sphinx-reredirects myst-parser linkify-it-py

Or::

   conda install -c conda-forge mpmath matplotlib sphinx sphinx-math-dollar sphinx-reredirects myst-parser linkify-it-py

Making your Sphinx build successful on the Windows system is tricky because
some dependencies like ``dvipng`` or ``docbook2x`` are not available.

Windows 10
~~~~~~~~~~~~

For Windows 10, however, the Windows Subsystem for Linux can be a possible
workaround solution, and you can install Ubuntu shell on your Windows system
after following the tutorial below:

https://github.com/MicrosoftDocs/WSL/blob/live/WSL/install-win10.md

In your command prompt, run ``ubuntu`` to transfer to Linux terminal, and
follow the Debian/Ubuntu tutorial above to install the dependencies, and then
you can run ``make html`` to build. (Note that you also have to install
``make`` via ``apt-get install make``.)

If you want to change the directory in your prompt to your working folder of
SymPy in the Windows file system, you can prepend ``cd /mnt/`` to your file
path in Windows, and run in your shell to navigate to the folder. (Also note
that Linux uses ``/`` instead of ``\`` for file paths.)

This method provides better compatibility than Cygwin or MSYS2 and more
convenience than a virtual machine if you partially need a Linux environment
for your workflow, however this method is only viable for Windows 10 64-bit
users.

Build the Docs
^^^^^^^^^^^^^^^

The documentation can be built by running the ``makefile`` in the ``doc``
subdirectory.

To start, in your preferred web browser, use the drop down menu and select
“open file” to navigate into the sympy/doc folder saved on your computer. In
the doc folder, select the _build folder, then the html folder, and in the html
folder, open the index.html file.

To build the HTML documentation, run::

   cd doc

   make html

This builds a local version of the documentation in ``doc/_build/html`` in your
web browser.

Open ``_build/html/index.html``.
