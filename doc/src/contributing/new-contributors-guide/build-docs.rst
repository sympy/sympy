.. _build-the-documentation:

==========================
Building the Documentation
==========================


Start by installing the required dependencies for the documentation.

Required dependencies
^^^^^^^^^^^^^^^^^^^^^^

You can either install the dependencies locally on your machine, or you can
build a Docker image containing them.

Docker
~~~~~~

If you have `Docker <https://docs.docker.com/engine/install/>`_, then instead of
following the OS-specific installation instructions below, you may choose to
build a Docker image::

   cd doc

   docker build -f Dockerfile.htmldoc -t sympy_htmldoc .

If you choose this option, you can now skip down to the "Build the Docs"
section below.

Debian/Ubuntu
~~~~~~~~~~~~~~~

For Debian/Ubuntu::

   apt-get install python3-sphinx texlive-latex-recommended dvipng librsvg2-bin imagemagick docbook2x graphviz

Install pip using::

   sudo apt install python3-pip

However, you can also create a virtual environment and use pip in it using::

   python3 -m venv /path/to/my/venv  # create the venv

Then activate it using::

   source /path/to/my/venv/bin/activate  # need to rerun this each time you open a new terminal

After installing pip through either of the two methods given above, run::

   python -m pip install -r doc/requirements.txt

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

   python -m pip install -r doc/requirements.txt

If you get mpmath error, install python3-mpmath package::

   dnf install python3-mpmath

If you get matplotlib error, install python3-matplotlib package::

   dnf install python3-matplotlib

Mac
~~~~

For Mac, first install homebrew: https://brew.sh/

Then install these packages with homebrew::

   brew install imagemagick graphviz docbook librsvg

Install the docs dependencies with either pip or conda::

   python -m pip install -r requirements.txt

Or::

   conda install -c conda-forge --file requirements.txt

Making your Sphinx build successful on the Windows system is tricky because
some dependencies like ``dvipng`` or ``docbook2x`` are not available.

Windows 10
~~~~~~~~~~~~

For Windows 10, however, the Windows Subsystem for Linux can be a possible
workaround solution, and you can install Ubuntu shell on your Windows system
after following the tutorial below:

https://learn.microsoft.com/en-us/windows/wsl/install

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

or

Follow `instruction <https://chocolatey.org/install>`_ to install Chocolatey

Install make and other dependencies::

   choco install make graphviz rsvg-convert imagemagick

Install python dependencies::

   pip install -r doc/requirements.txt

Build the Docs
^^^^^^^^^^^^^^^

Docker
~~~~~~

If you chose to build using Docker, and followed the instructions above to
build the ``sympy_htmldoc`` image, then you can build the docs with::

   docker run --rm -v /absolute/path/to/sympy:/sympy sympy_htmldoc

(Be sure to substitute the actual absolute filesystem path to sympy!) This
command can be run from any directory.

Local Installation
~~~~~~~~~~~~~~~~~~

If you chose to follow OS-specific instructions above and installed the
required dependencies locally, the documentation can be built by running the
``makefile`` in the ``doc`` subdirectory::

   cd doc

   make html

View the Docs
^^^^^^^^^^^^^

Once you have built the docs, the generated files will be found under
``doc/_build/html``. To view them in your preferred web browser, use the drop
down menu and select “open file”, navigate into the ``sympy/doc/_build/html``
folder, and open the ``index.html`` file.


Auto-Rebuild with the Live Server
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The instructions given above told you how to build the docs once, and load them
in the browser. After you make changes to the document sources, you'll have to
manually repeat the build step, and reload the pages in the browser.

There is an alternative approach that sets up a live server, which will monitor
the docs directory, automatically rebuild when changes are detected, and
automatically reload the page you are viewing in the browser.

If you want to use this option, the procedure again depends on whether you are
using Docker, or a local installation.

Docker
~~~~~~

To start the live server with Docker, you can use::

   docker run --rm -it \
        -v /absolute/path/to/sympy:/sympy \
        -p 8000:80 \
        sympy_htmldoc live

and then navigate your browser to ``localhost:8000``. You can use a different
port by changing the ``8000`` in the command. Again, be sure to substitute the
actual absolute filesystem path to sympy.

When finished, you can stop the server with ``ctrl-c`` in the terminal.

Alternatively, you may run the server in detached mode, using::

   docker run --rm -d --name=sympy-livehtml \
        -v /absolute/path/to/sympy:/sympy \
        -p 8000:80 \
        sympy_htmldoc live

and then stop it with::

   docker stop sympy-livehtml


Local Installation
~~~~~~~~~~~~~~~~~~

If you installed the build dependencies locally, then simply use::

   cd doc

   make livehtml

to start the server. Your web browser should then automatically open a new tab,
showing the index page of the SymPy docs.

When you are finished, you can use ``ctrl-c`` in the terminal to stop the
server.


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

    brew install --cask chromium

    brew tap homebrew/cask-fonts

    brew install font-dejavu

On Windows 10, you can use::

    choco install chromium strawberryperl miktex dejavufonts

If DejaVu fonts are not installed in ``C:\Windows\Fonts``, then open
``~\AppData\Local\Microsoft\Windows\Fonts``, select all DejaVu fonts,
right-click and click ``Install for all users``.

To build the pdf docs run::

    cd doc

    make pdf

The resulting PDF will be in::

    _build/latex/sympy-<version>.pdf

where ``<version>`` is the SymPy version (e.g., ``sympy-1.10.dev.pdf``).
