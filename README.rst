SymPy
=====

A Python library for symbolic mathematics.

http://sympy.org/

See the AUTHORS file for the list of authors.

And many more people helped on the SymPy mailinglist, reported bugs, helped
organize SymPy's participation in the Google Summer of Code, the Google Highly
Open Participation Contest, wrote and blogged about SymPy...

License: New BSD License (see the LICENSE file for details)
covers all files in the sympy repository unless stated otherwise.

0. Download
-----------

::

    $ git clone git://github.com/sympy/sympy.git

For other options (tarballs, debs, etc.), see the web page of SymPy.

1. Documentation and usage
--------------------------

Everything is at:

http://docs.sympy.org/

You can generate everything at the above site in your local copy of SymPy by::

    $ cd doc
    $ make html
    $ epiphany _build/html/index.html # Or your preferred web browser

If you don't want to read that, here is a short usage:

From this directory, start python and::

    >>> from sympy import Symbol, cos
    >>> x = Symbol('x')
    >>> e = 1/cos(x)
    >>> print e.series(x, 0, 10)
    1 + (1/2)*x**2 + (5/24)*x**4 + (61/720)*x**6 + (277/8064)*x**8 + O(x**10)

SymPy also comes with a console that is a simple wrapper around the
classic python console (or ipython when available) that loads the
sympy namespace and executes some common commands for you.

To start it, issue::

    $ bin/isympy

from this directory if SymPy is not installed or simply::

    $ isympy

if SymPy is installed somewhere in your ``PATH``.

3. Installation
---------------

To install SymPy, simply run::

    $ python setup.py install

If you install it system-wide, you may need to prefix the previous command with ``sudo``::

    $ sudo python setup.py install

4. Tests
--------

To execute all tests, run::

    $./setup.py test

in the current directory.

For more fine-grained running of tests, use ``bin/test`` or respectively
``bin/doctest``. The master branch is automatically tested by Travis CI,
the results can be seen here:

.. image:: https://secure.travis-ci.org/sympy/sympy.png?branch=master
    :target: http://travis-ci.org/sympy/sympy

To test pull requests, use `sympy-bot <https://github.com/sympy/sympy-bot>`_.

5. Usage in Python 3
-------------------

SymPy also supports Python 3. Currently, this is implemented by maintaining a
Python 2 compatible codebase and running our own 2to3 script on it. Run it with::

    $ bin/use2to3

When ran, it will create a new directory, py3k-sympy, which holds a Python 3
compatible version of the code. SymPy can then be used normally with Python 3
from that directory (installation, interactive shell, etc.).

6. Clean
--------

To clean everything (thus getting the same tree as in the repository)::

    $./setup.py clean

7. Brief History
----------------

SymPy was started by Ondrej Certik in 2005, he wrote some code during the
summer, then he wrote some more code during the summer 2006. In February 2007,
Fabian Pedregosa joined the project and helped fixed many things, contributed
documentation and made it alive again. 5 students (Mateusz Paprocki, Brian
Jorgensen, Jason Gedge, Robert Schwarz and Chris Wu) improved SymPy incredibly
during the summer 2007 as part of the Google Summer of Code. Pearu Peterson
joined the development during the summer 2007 and he has made SymPy much more
competitive by rewriting the core from scratch, that has made it from 10x to
100x faster. Jurjen N.E. Bos has contributed pretty printing and other patches.
Fredrik Johansson has wrote mpmath and contributed a lot of patches. Since
then, a lot more people have joined the development and some people have also
left. You can see the full list in doc/src/aboutus.rst, or online at:

http://docs.sympy.org/dev/aboutus.html#sympy-development-team

For people that don't want to be listed there, see the git history.


8. Citation
-----------

To cite SymPy in publications use::

    SymPy Development Team (2013). SymPy: Python library for symbolic mathematics
    URL http://www.sympy.org.

A BibTeX entry for LaTeX users is::

    @Manual{,
    title = {SymPy: Python library for symbolic mathematics},
    author = {{SymPy Development Team}},
    year = {2013},
    url = {http://www.sympy.org},
    }

SymPy is BSD licensed, so you are free to use it whatever you like, be it
academic, commercial, creating forks or derivatives, as long as you copy the BSD
statement if you redistribute it (see the LICENSE file for details).
That said, although not required by the SymPy license, if it is convenient for
you, please cite SymPy when using it in your work and also consider
contributing all your changes back, so that we can incorporate it and all of us
will benefit in the end.
