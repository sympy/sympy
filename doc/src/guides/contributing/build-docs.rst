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
   python -m pip install sphinx-math-dollar

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
   python -m pip install sphinx-math-dollar

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

   python -m pip install mpmath matplotlib sphinx sphinx-math-dollar sphinx-reredirects

Or::

   conda install -c conda-forge mpmath matplotlib sphinx sphinx-math-dollar sphinx-reredirects

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

PDF Documentation
^^^^^^^^^^^^^^^^^

.. note::

   It is not necessary for the majority of contributors to build the PDF
   documentation. The PDF documentation will be built automatically on GitHub
   Actions on pull requests. PDF documentation for each release is included on
   the `GitHub releases page <https://github.com/sympy/sympy/releases>`_.

   If the PDF documentation build fails on GitHub Actions, 99% of the time
   this is due to bad LaTeX math formatting. Double check that any math you
   have added is formatted correctly, and make sure you use \`\`double
   backticks\`\` for code (\`single backticks\` will render as math, not
   code). See the resources in the :ref:`style guide
   <style_guide_latex_recommendations>` for tips on formatting LaTeX math.

Building the PDF documentation requires a few extra dependencies. First you
will need to have a TeXLive installation that includes XeLaTeX and latexmk.
You will also need to have Chrome or Chromium installed, as it is used to
convert some SVG files for the PDF.

On Ubuntu, you can install these with::

    apt-get install chromium-browser texlive texlive-xetex texlive-fonts-recommended texlive-latex-extra latexmk lmodern

On Mac, you can use::

    brew install texlive

and also make sure the `Google Chrome browser
<https://www.google.com/chrome/>`_ is installed in ``/Applications``.

To build the pdf docs run::

    cd doc

    make latexpdf

The resulting PDF will be in::

    _build/latex/sympy-<version>.pdf

where ``<version>`` is the SymPy version (e.g., ``sympy-1.10.dev.pdf``).
