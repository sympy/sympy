
.. include:: ../definitions.def

Installing, configuring and running SymPy
=========================================

The easiest way to get SymPy is to type ``sudo pip install sympy``,
which will download the latest version of the package from PyPI and
install it.  If you want to get the source and install it manually,
visit `this <http://code.google.com/p/sympy>`_ page and download the
latest tarball from *Featured Downloads* section, or use the following
direct link:

.. parsed-literal::

    $ :input:`wget http://sympy.googlecode.com/files/sympy-0.7.2.tar.gz`
    $ :input:`tar -xz -C sympy --strip-components 1 -f sympy-0.7.2.tar.gz`

You will also find an installer for Windows there. An alternative way is
to clone SymPy's `git <http://www.git-scm.org>`_ repository from
`GitHub <http://github.com/sympy/sympy>`_:

.. parsed-literal::

    $ :input:`git clone git://github.com/sympy/sympy.git`

Note that this will give you development version of SymPy, so certain features
may not work properly, APIs are unstable and documentation may be outdated. In
return you get the cutting edge features before SymPy is officially released and
you can easily participate in SymPy's development process.

To use it, issue:

.. parsed-literal::

    $ :input:`cd sympy`
    $ :input:`python`
    Python 2.7.2+ (default, Oct  4 2011, 20:06:09)
    [GCC 4.6.1] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> from sympy import *
    >>> var('x')
    x
    >>> diff(sin(x), x)
    cos(x)

If you want SymPy to be available globally, you can install it using
``./setup.py install``. SymPy is available in major Linux distributions,
so you can install it also using package manager of your distribution
(for example in Ubuntu install ``python-sympy`` and in Gentoo install
``dev-python/sympy``). Note that distributions often contain very outdated
versions of packages, so make sure you don't install some ancient version
of SymPy on your system.

By default, SymPy doesn't have any dependencies besides Python. The following
versions of Python are supported: 2.5, 2.6, 2.7, 3.2 and 3.3. Version 2.4 was
originally supported but was dropped in SymPy 0.7.1. Python 3.0 and 3.1
aren't officially supported due to their limitations, but it's likely that
you will be able to run SymPy there and use most of its features. SymPy also
can be run on PyPy 1.6+. If you use Jython, IronPython or other less popular
interpreter, you may find out SymPy not working properly (Jython) or you
won't be able to run it at all (IronPython).

Certain modules may require additional dependencies to provide more
features or improve speed of computations. For example, :mod:`sympy.polys`
and :mod:`sympy.mpmath` can take advantage of gmpy library to radically
improve speed of otherwise pure Python library. If gmpy is not available,
those modules fall back automatically to pure Python implementation of
arithmetic routines. Other optional dependencies are IPython, Matplotlib,
NumPy, SciPy, Cython, Pyglet, LaTeX distribution and more.

SymPy in Python/IPython
-----------------------

Sessions in standard Python's interpreter and IPython look very similar,
just the banner and prompt look differently, for example:

.. parsed-literal::

    $ :input:`python`
    Python 2.7.2+ (default, Oct  4 2011, 20:06:09)
    [GCC 4.6.1] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import sympy
    >>> x = sympy.Symbol('x')
    >>> sympy.integrate(3*x**2)
    x**3
    >>> sympy.init_printing(use_unicode=True, no_global=True)
    >>> sympy.integrate(3*x**2)
     3
    x

Interactive SymPy (``isympy``)
------------------------------

For users' convenience, SymPy's distribution includes a simple script called
isympy. This script is installed by SymPy's installers and is also available in
``bin/isympy`` in source distribution. isympy uses either IPython (if available) or
standard Python's interpreter with readline support. On startup isympy sets
up the environment to make interaction with SymPy more pleasant. It enables
new division, imports everything from :mod:`sympy`, injects a few commonly
used symbols into the global namespace, and initializes the pretty printer.

Here is an example session with isympy:

.. parsed-literal::

    $ :input:`bin/isympy`
    IPython console for SymPy 0.7.2 (Python 2.7.1-64-bit) (ground types: gmpy)

    These commands were executed:
    >>> from __future__ import division
    >>> from sympy import *
    >>> x, y, z, t = symbols('x y z t')
    >>> k, m, n = symbols('k m n', integer=True)
    >>> f, g, h = symbols('f g h', cls=Function)

    Documentation can be found at http://www.sympy.org

    In [1]: integrate(3*x**2, x)
    Out[1]:
     3
    x

    In [2]: %quit
    Do you really want to exit ([y]/n)? y
    Exiting ...

Command-line arguments
~~~~~~~~~~~~~~~~~~~~~~

There are a variety of command-line options supported by isympy:

``-h``, ``--help``
    show help
``-c CONSOLE``, ``--console=CONSOLE``
    select type of interactive session: ``ipython``, ``python``. Default is
    ``ipython`` if IPython is installed, otherwise, ``python``.
``-p PRETTY``, ``--pretty=PRETTY``
    setup pretty printing: ``unicode``, ``ascii`` or ``no``. Default is ``unicode``
    if the terminal supports it, otherwise, ``ascii``.
``-t TYPES``, ``--types=TYPES``
    setup ground types: ``gmpy``, ``python`` or ``sympy``. Default is ``gmpy`` if
    it's installed, otherwise ``python``.
``-o ORDER``, ``--order=ORDER``
    setup ordering of terms: ``[rev-]lex``, ``[rev-]grlex``, ``[rev-]grevlex`` or
    ``old``. Default is ``lex``.
``-q``, ``--quiet``
    print only version information at startup
``-C``, ``--no-cache``
    disable caching

Environment variables
---------------------

``SYMPY_USE_CACHE``
    By default SymPy caches all computations. If this is undesirable, for
    example due to limited amount of memory, set this variable to ``no``
    to disable caching. Note that some operations will run much slower with
    the cache off. Setting this variable to ``no`` is equivalent to running
    isympy with ``-C`` option.
``SYMPY_GROUND_TYPES``
    SymPy is a pure Python library, however to improve the speed of computations
    it can take advantage of gmpy library to speedup coefficient arithmetics
    (also known as ground domain arithmetics). Ground types are set automatically,
    so if gmpy is not available, it simply won't be used. However, if gmpy is
    available but for some reason it is undesirable to use it, set this variable
    to ``python``, to disable usage of gmpy. Use ``-t`` or ``--type`` option to
    achieve the same in isympy.

Running the test suite
----------------------

To verify that SymPy works properly on your computer, you can run SymPy's
test suite. This is done either with ``bin/test`` command or :func:`test`
in an interactive session. For example, to test :mod:`sympy.core` issue:

.. parsed-literal::

    $ :input:`bin/test sympy/core`
    ============================= test process starts ==============================
    executable:   /usr/bin/python  (2.7.2-final-0)
    architecture: 64-bit
    ground types: gmpy
    random seed:  65858271

    sympy/core/tests/test_arit.py[48] ...f..........................................
    ..                                                                          [OK]
    sympy/core/tests/test_assumptions.py[28] ............................       [OK]
    sympy/core/tests/test_basic.py[10] ..........                               [OK]
    sympy/core/tests/test_cache.py[1] .                                         [OK]
    sympy/core/tests/test_complex.py[13] .............                          [OK]
    sympy/core/tests/test_containers.py[5] .....                                [OK]
    sympy/core/tests/test_count_ops.py[2] ..                                    [OK]
    sympy/core/tests/test_diff.py[6] ......                                     [OK]
    sympy/core/tests/test_equal.py[6] ......                                    [OK]
    sympy/core/tests/test_eval.py[8] .....f..                                   [OK]
    sympy/core/tests/test_eval_power.py[13] .............                       [OK]
    sympy/core/tests/test_evalf.py[24] ........................                 [OK]
    sympy/core/tests/test_expand.py[6] ......                                   [OK]
    sympy/core/tests/test_expr.py[59] ..............................................
    .............                                                               [OK]
    sympy/core/tests/test_exprtools.py[4] ....                                  [OK]
    sympy/core/tests/test_facts.py[11] ...........                              [OK]
    sympy/core/tests/test_functions.py[27] .....fff...................          [OK]
    sympy/core/tests/test_logic.py[11] ...........                              [OK]
    sympy/core/tests/test_match.py[26] ...f......................               [OK]
    sympy/core/tests/test_numbers.py[46] ...........................................
    ...                                                                         [OK]
    sympy/core/tests/test_operations.py[4] ....                                 [OK]
    sympy/core/tests/test_priority.py[5] .....                                  [OK]
    sympy/core/tests/test_relational.py[7] .......                              [OK]
    sympy/core/tests/test_sets.py[18] ..................                        [OK]
    sympy/core/tests/test_subs.py[30] ..............................            [OK]
    sympy/core/tests/test_symbol.py[9] ....X....                                [OK]
    sympy/core/tests/test_sympify.py[26] ...f......................             [OK]
    sympy/core/tests/test_truediv.py[3] ...                                     [OK]
    sympy/core/tests/test_var.py[5] .....                                       [OK]

    ________________________________ xpassed tests _________________________________
    sympy/core/tests/test_symbol.py:

    tests finished: 453 passed, 7 expected to fail, 1 expected to fail but passed,
    in 6.30 seconds

This tells us that all standard tests in :mod:`sympy.core`'s  pass (dots).
In case of failure, ``.`` would change to ``F`` and ``OK`` to ``FAIL``
(additionally all failures would be colored in red and listed at the end of
output from SymPy's test utility). Non-standard tests are those marked with
``f`` and ``X`` characters. The former means that the test was supposed to
fail and failed (XFAIL), whereas the later means that the test was supposed
to fail but passed (XPASS).

To run the whole test suite issue ``bin/test`` or :func:`test` without any
arguments. Running the whole test suite takes more than ten minutes on
Pentium-M 1.6 GHz and less than 5 minutes on Xeon 3.0 GHz (one core).

There is another test utility in SymPy, ``bin/doctest``, which verifies
examples in docstrings and documentation. If you are going to contribute
to SymPy, make sure you run both ``bin/test`` and ``bin/doctest`` before
submitting a `pull request <http://help.github.com/send-pull-requests>`_.

SymPy in web browsers
---------------------

SymPy is available in the following web applications:

* SymPy Live (http://live.sympy.org)
* SymPy Gamma (http://gamma.sympy.org)
* Sage Notebook (http://www.sagenb.org)

SymPy Live was developed specifically for SymPy. It is a simple web shell
that looks similar to isympy under standard Python's interpreter. SymPy
Live uses Google App Engine as computational backend.

IPython Notebook and QTConsole (IPython 0.12+) have native support for SymPy.
Start IPython with either ``ipython notebook`` or ``ipython qtconsole`` and
issue ``%load_ext sympy.interactive.ipythonprinting`` magic command. This
will give you nice pretty printing via LaTeX.
