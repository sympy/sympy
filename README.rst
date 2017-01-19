SymPy
=====

|pypi version| |pypi download| |Build status| |Gitter Badge| |Zenodo Badge|

.. |pypi version| image:: https://img.shields.io/pypi/v/sympy.svg
   :target: https://pypi.python.org/pypi/sympy
.. |pypi download| image:: https://img.shields.io/pypi/dm/sympy.svg
   :target: https://pypi.python.org/pypi/sympy
.. |Build status| image:: https://secure.travis-ci.org/sympy/sympy.svg?branch=master
   :target: http://travis-ci.org/sympy/sympy
.. |Gitter Badge| image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/sympy/sympy
   :target: https://gitter.im/sympy/sympy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
.. |Zenodo Badge| image:: https://zenodo.org/badge/18918/sympy/sympy.svg
   :target: https://zenodo.org/badge/latestdoi/18918/sympy/sympy

A Python library for symbolic mathematics.

http://sympy.org/

See the AUTHORS file for the list of authors.

And many more people helped on the SymPy mailing list, reported bugs, helped
organize SymPy's participation in the Google Summer of Code, the Google Highly
Open Participation Contest, Google Code-In, wrote and blogged about SymPy...

License: New BSD License (see the LICENSE file for details) covers all files
in the sympy repository unless stated otherwise.

Our mailing list is at
https://groups.google.com/forum/?fromgroups#!forum/sympy.

We have community chat at `Gitter <https://gitter.im/sympy/sympy>`_. Feel free
to ask us anything there. We have a very welcoming and helpful community.


Download
--------

Get the latest version of SymPy from
https://pypi.python.org/pypi/sympy/

To get the git version do

::

    $ git clone git://github.com/sympy/sympy.git

For other options (tarballs, debs, etc.), see
http://docs.sympy.org/dev/install.html.

Documentation and usage
-----------------------

Everything is at:

http://docs.sympy.org/

You can generate everything at the above site in your local copy of SymPy by::

    $ cd doc
    $ make html

Then the docs will be in `_build/html`. If you don't want to read that, here
is a short usage:

From this directory, start python and::

    >>> from sympy import Symbol, cos
    >>> x = Symbol('x')
    >>> e = 1/cos(x)
    >>> print e.series(x, 0, 10)
    1 + x**2/2 + 5*x**4/24 + 61*x**6/720 + 277*x**8/8064 + O(x**10)

SymPy also comes with a console that is a simple wrapper around the
classic python console (or IPython when available) that loads the
sympy namespace and executes some common commands for you.

To start it, issue::

    $ bin/isympy

from this directory if SymPy is not installed or simply::

    $ isympy

if SymPy is installed.

Installation
------------

SymPy has a hard dependency on the `mpmath <http://mpmath.org/>`
library (version >= 0.19).  You should install it first, please refer to
the mpmath installation guide:

https://github.com/fredrik-johansson/mpmath#1-download--installation

To install SymPy itself, then simply run::

    $ python setup.py install

If you install it system-wide, you may need to prefix the previous command with ``sudo``::

    $ sudo python setup.py install

See http://docs.sympy.org/dev/install.html for more information.

Contributing
------------

We welcome contributions from anyone, even if you are new to open
source. Please read our `introduction to contributing
<https://github.com/sympy/sympy/wiki/Introduction-to-contributing>`_. If you
are new and looking for some way to contribute a good place to start is to
look at the issues tagged `Easy to Fix
<https://github.com/sympy/sympy/issues?q=is%3Aopen+is%3Aissue+label%3A%22Easy+to+Fix%22>`_.

Please note that all participants of this project are expected to follow our
Code of Conduct. By participating in this project you agree to abide by its
terms. See `CODE_OF_CONDUCT.md <CODE_OF_CONDUCT.md>`_.

Tests
-----

To execute all tests, run::

    $./setup.py test

in the current directory.

For more fine-grained running of tests or doctest, use ``bin/test`` or
respectively ``bin/doctest``. The master branch is automatically tested by
Travis CI.

To test pull requests, use `sympy-bot <https://github.com/sympy/sympy-bot>`_.

Usage in Python 3
-----------------

SymPy also supports Python 3. If you want to install the latest version in
Python 3, get the Python 3 tarball from
https://pypi.python.org/pypi/sympy/

To install the SymPy for Python 3, simply run the above commands with a Python
3 interpreter.

Clean
-----

To clean everything (thus getting the same tree as in the repository)::

    $ ./setup.py clean

You can also clean things with git using::

    $ git clean -Xdf

which will clear everything ignored by ``.gitignore``, and::

    $ git clean -df

to clear all untracked files.  You can revert the most recent changes in git
with::

    $ git reset --hard

WARNING: The above commands will all clear changes you may have made, and you
will lose them forever. Be sure to check things with ``git status``, ``git
diff``, ``git clean -Xn`` and ``git clean -n`` before doing any of those.

Bugs
----

Our issue tracker is at https://github.com/sympy/sympy/issues.  Please report
any bugs that you find.  Or, even better, fork the repository on GitHub and
create a pull request.  We welcome all changes, big or small, and we will help
you make the pull request if you are new to git (just ask on our mailing list
or Gitter).

Brief History
-------------

SymPy was started by Ondřej Čertík in 2005, he wrote some code during the
summer, then he wrote some more code during the summer 2006. In February 2007,
Fabian Pedregosa joined the project and helped fixed many things, contributed
documentation and made it alive again. 5 students (Mateusz Paprocki, Brian
Jorgensen, Jason Gedge, Robert Schwarz and Chris Wu) improved SymPy incredibly
during the summer 2007 as part of the Google Summer of Code. Pearu Peterson
joined the development during the summer 2007 and he has made SymPy much more
competitive by rewriting the core from scratch, that has made it from 10x to
100x faster. Jurjen N.E. Bos has contributed pretty printing and other patches.
Fredrik Johansson has written mpmath and contributed a lot of patches.

SymPy has participated in every Google Summer of Code since 2007. You can see
https://github.com/sympy/sympy/wiki#google-summer-of-code for full details.
Each year has improved SymPy by bounds. Most of SymPy's development has come
from Google Summer of Code students.

In 2011, Ondřej Čertík stepped down as lead developer, with Aaron Meurer, who
also started as a Google Summer of Code student, taking his place. Ondřej
Čertík is still active in the community, but is too busy with work and family
to play a lead development role.

Since then, a lot more people have joined the development and some people have
also left. You can see the full list in doc/src/aboutus.rst, or online at:

http://docs.sympy.org/dev/aboutus.html#sympy-development-team

The git history goes back to 2007, when development moved from svn to hg.  To
see the history before that point, look at http://github.com/sympy/sympy-old.

You can use git to see the biggest developers.  The command::

     $ git shortlog -ns

will show each developer, sorted by commits to the project.  The command::

     $ git shortlog -ns --since="1 year"

will show the top developers from the last year.

Citation
--------

To cite SymPy in publications use

    Meurer A, Smith CP, Paprocki M, Čertík O, Kirpichev SB, Rocklin M, Kumar A,
    Ivanov S, Moore JK, Singh S, Rathnayake T, Vig S, Granger BE, Muller RP,
    Bonazzi F, Gupta H, Vats S, Johansson F, Pedregosa F, Curry MJ, Terrel AR,
    Roučka Š, Saboo A, Fernando I, Kulal S, Cimrman R, Scopatz A. (2017) SymPy:
    symbolic computing in Python. *PeerJ Computer Science* 3:e103
    https://doi.org/10.7717/peerj-cs.103

A BibTeX entry for LaTeX users is

.. code-block:: none

    @article{10.7717/peerj-cs.103,
     title = {SymPy: symbolic computing in Python},
     author = {Meurer, Aaron and Smith, Christopher P. and Paprocki, Mateusz and \v{C}ert\'{i}k, Ond\v{r}ej and Kirpichev, Sergey B. and Rocklin, Matthew and Kumar, AMiT and Ivanov, Sergiu and Moore, Jason K. and Singh, Sartaj and Rathnayake, Thilina and Vig, Sean and Granger, Brian E. and Muller, Richard P. and Bonazzi, Francesco and Gupta, Harsh and Vats, Shivam and Johansson, Fredrik and Pedregosa, Fabian and Curry, Matthew J. and Terrel, Andy R. and Rou\v{c}ka, \v{S}t\v{e}p\'{a}n and Saboo, Ashutosh and Fernando, Isuru and Kulal, Sumith and Cimrman, Robert and Scopatz, Anthony},
     year = 2017,
     month = jan,
     keywords = {Python, Computer algebra system, Symbolics},
     abstract = {
                SymPy is an open source computer algebra system written in pure Python. It is built with a focus on extensibility and ease of use, through both interactive and programmatic applications. These characteristics have led SymPy to become a popular symbolic library for the scientific Python ecosystem. This paper presents the architecture of SymPy, a description of its features, and a discussion of select submodules. The supplementary material provide additional examples and further outline details of the architecture and features of SymPy.
             },
     volume = 3,
     pages = {e103},
     journal = {PeerJ Computer Science},
     issn = {2376-5992},
     url = {https://doi.org/10.7717/peerj-cs.103},
     doi = {10.7717/peerj-cs.103}
    }

SymPy is BSD licensed, so you are free to use it whatever you like, be it
academic, commercial, creating forks or derivatives, as long as you copy the
BSD statement if you redistribute it (see the LICENSE file for details).  That
said, although not required by the SymPy license, if it is convenient for you,
please cite SymPy when using it in your work and also consider contributing
all your changes back, so that we can incorporate it and all of us will
benefit in the end.
