.. _installation:

Installation
------------

The SymPy CAS can be installed on virtually any computer with Python
2.6 or above. SymPy does require `mpmath`_ Python library to be
installed first.  The current recommended method of installation is
directly from the source files.  Alternatively, executables are
available for Windows, and some Linux distributions have SymPy
packages available.

SymPy officially supports Python 2.6, 2.7, 3.2, 3.3, 3.4, and PyPy.

Source
======

SymPy currently recommends that users install directly from the source files.
You will first have to download the source files via the archive. Download the
latest release (tar.gz) from the `downloads site`_ and open it with your
operating system's standard decompression utility.

After the download is complete, you should have a folder called "sympy". From
your favorite command line terminal, change directory into that folder and
execute the following::

    $ python setup.py install

Alternatively, if you don't want to install the package onto your computer, you
may run SymPy with the "isympy" console (which automatically imports SymPy
packages and defines common symbols) by executing within the "sympy" folder::

    $ ./bin/isympy

You may now run SymPy statements directly within the Python shell::

    >>> from __future__ import division
    >>> from sympy import *
    >>> x, y, z, t = symbols('x y z t')
    >>> k, m, n = symbols('k m n', integer=True)
    >>> f, g, h = symbols('f g h', cls=Function)
    >>> diff(x**2/2, x)
    x

Git
===

If you are a developer or like to get the latest updates as they come, be sure
to install from git. To download the repository, execute the following from the
command line::

    $ git clone git://github.com/sympy/sympy.git

Then, execute either the `setup.py` or the `bin/isympy` scripts as demonstrated
above.

To update to the latest version, go into your repository and execute::

    $ git pull origin master

If you want to install SymPy, but still want to use the git version, you can run
from your repository::

    $ setupegg.py develop

This will cause the installed version to always point to the version in the git
directory.

Anaconda
========

Although SymPy does not have any hard dependencies, many nice features are
only enabled when certain libraries are installed.  For example, without
Matplotlib, only simple text-based plotting is enabled.  With the IPython
notebook or qtconsole, you can get nicer `\LaTeX` printing by running
``init_printing()``.  An easy way to get all these libraries in addition to
SymPy is to install `Anaconda <http://continuum.io/downloads>`_, which is
a free Python distribution from Continuum Analytics that includes SymPy,
Matplotlib, IPython, NumPy, and many more useful packages for scientific
computing.

Other Methods
=============

An installation executable (.exe) is available for Windows users at the
`downloads site`_. In addition, various Linux distributions have SymPy
available as a package. Others are strongly encouraged to download from source
(details above).

Run SymPy
=========

After installation, it is best to verify that your freshly-installed SymPy
works. To do this, start up Python and import the SymPy libraries::

    $ python
    >>> from sympy import *

From here, execute some simple SymPy statements like the ones below::

    >>> x = Symbol('x')
    >>> limit(sin(x)/x, x, 0)
    1
    >>> integrate(1/x, x)
    log(x)

For a starter guide on using SymPy effectively, refer to the :ref:`tutorial`.

Questions
=========

If you have a question about installation or SymPy in general, feel free to
visit our chat on `Gitter`_. In addition, our `mailing list`_ is an excellent
source of community support.

If you think there's a bug or you would like to request a feature, please open
an `issue ticket`_.

.. _downloads site: https://github.com/sympy/sympy/releases
.. _Gitter: https://gitter.im/sympy/sympy
.. _issue ticket: https://github.com/sympy/sympy/issues
.. _mailing list: https://groups.google.com/forum/#!forum/sympy
.. _mpmath: http://mpmath.org/
