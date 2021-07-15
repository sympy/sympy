.. _devsetup:

================================
Development Environment Setup
================================

The first step to contributing to the code base is creating your development environment.

Git Setup
-----------

SymPy is available on `GitHub <https://github.com/micropython/micropython>`_ and uses
`Git <https://git-scm.com>`_ for source control. The workflow is such that
code is pulled and pushed to and from the main repository. Install the respective version
of Git for your operating system to start development.

.. note::
   Refer to the installation instructions in
   the `Git installation instructions <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.
   Learn about the basic git commands in this `Git Handbook <https://guides.github.com/introduction/git-handbook/>`_
   or any other sources on the internet.

Get the SymPy Code
-------------------

It is recommended practice to create a fork of the SymPy project for your development purposes. Create your own fork of the SymPy project (if you have not yet). Go to the SymPy GitHub repository:

.. code-block:: bash

   https://github.com/sympy/sympy

You will now have a fork at <https://github.com/<your-user-name>/sympy>.

Then, nn your machine browse to where you would like to store SymPy, and clone (download) the latest code from SymPy's original repository (about 77 MiB):

.. code-block:: bash

   $ git clone https://github.com/<your-user-name>/sympy

You must `configure the remote repositories <https://git-scm.com/book/en/v2/Git-Basics-Working-with-Remotes>`_ for collaboration with the upstream project:

.. code-block:: bash

   $ cd sympy
   $ git remote add upstream https://github.com/sympy/sympy

After the configuration, your setup should be similar to this:

.. code-block:: bash

   $ git remote -v
   origin	https://github.com/<your-user-name>/sympy (fetch)
   origin	https://github.com/<your-user-name>/sympy (push)
   upstream	https://github.com/sympy/sympy (fetch)
   upstream	https://github.com/sympy/sympy (push)

For further development, it is recommended
to create a development branch.

.. code-block:: bash

    $ git checkout -b dev-branch

The new branch can be of any name.

Virtual Environment Setup
---------------------------

You may want to take advantage of using virtual environments to isolate your development version of SymPy from any system wide installed versions, e.g. from ``apt-get install python-sympy``.

We recommend using ``conda`` to create a virtual environment:

.. code-block:: bash

    $ conda create -n sympy-dev python=3 mpmath flake8

You now have a environment that you can use for testing your development copy of SymPy. For example, clone your SymPy fork from Github:

.. code-block:: bash

    $ git clone git@github.com:<your-github-username>/sympy.git
    $ cd sympy

Now activate the environment:

.. code-block:: bash

    $ conda activate sympy-dev


Run the Tests
--------------

There are several ways of running Sympy tests but the easiest is to use the ``bin/test`` script, consult 'the wiki details on running tests <https://github.com/sympy/sympy/wiki/Running-tests>`_.

The script takes a number of options and arguments and then passes them to ``sympy.test(*paths, **kwargs)``. Run ``bin/test --help`` for all supported arguments.

Run all tests by using the command:

.. code-block:: bash

    $ bin/test

To run tests for a specific file, use:

.. code-block:: bash

    $ bin/test test_basic

Where ``test_basic`` is from file ``sympy/core/basic.py``.

To run tests for modules, use:

.. code-block:: bash

   $  bin/test /core /utilities

This will run tests for the ``core`` and ``utilities`` modules.

Similary, run quality tests with:

.. code-block:: bash

    $ bin/test code_quality

Build the Documentation
------------------------

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

Debugging
----------

To start sympy in debug mode set the SYMPY_DEBUG variable. For instance in a unix-like system you would do

    $ SYMPY_DEBUG=True bin/isympy

or in Windows

    > set SYMPY_DEBUG=True
    > python bin/isympy

Now just use for example the `limit()` function. You will get a nice printed tree, which is very useful for debugging.
