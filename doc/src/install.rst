.. _installation:

Installation
------------

The SymPy CAS can be installed on virtually any computer with Python 2.6 or
above. SymPy does require `mpmath`_ Python library to be installed first.  The
current recommended method of installation is through Anaconda, which includes
mpmath, as well as several other useful libraries.  Alternatively, executables
are available for Windows, and some Linux distributions have SymPy packages
available.

SymPy officially supports Python 2.6, 2.7, 3.2, 3.3, 3.4, and PyPy.

Anaconda
========

`Anaconda <http://continuum.io/downloads>`_ is a free Python a free Python
distribution from Continuum Analytics that includes SymPy, Matplotlib,
IPython, NumPy, and many more useful packages for scientific computing. This
is recommended because many nice features of SymPy are only enabled when
certain libraries are installed.  For example, without Matplotlib, only simple
text-based plotting is enabled.  With the IPython notebook or qtconsole, you
can get nicer `\LaTeX` printing by running ``init_printing()``.

If you already have Anaconda and want to update SymPy to the latest version,
use::

    conda update sympy

Git
===

If you wish to contribute to SymPy or like to get the latest updates as they
come, install SymPy from git. To download the repository, execute the
following from the command line::

    git clone git://github.com/sympy/sympy.git

To update to the latest version, go into your repository and execute::

    git pull origin master

If you want to install SymPy, but still want to use the git version, you can run
from your repository::

    setupegg.py develop

This will cause the installed version to always point to the version in the git
directory.

Other Methods
=============

An installation executable (.exe) is available for Windows users at the
`downloads site`_. In addition, various Linux distributions have SymPy
available as a package. You may also install SymPy from source or using pip.

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

Mpmath
======

Versions of SymPy prior to 0.7.7 included `mpmath`_, but it now depends on it as
an external dependency.  If you installed SymPy with Anaconda, it will already
include mpmath. Use::

  conda install mpmath

to ensure that it is installed.

If you do not wish to use Anaconda, you can use ``pip install mpmath``.

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
.. _mailing list: http://groups.google.com/group/sympy
.. _mpmath: http://mpmath.org/
