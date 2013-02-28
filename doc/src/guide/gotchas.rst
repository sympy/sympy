
.. include:: ../definitions.def

Gotchas and pitfalls
====================

SymPy is being written in and runs under `Python <http://www.python.org/>`_,
a general purpose programming language, so there are a few things that may
be quite different from what can be experienced in other symbolic mathematics
or computer algebra systems like Maple or Mathematica. These are some of the
gotchas and pitfalls that you may encounter when using SymPy.

``1/3`` is not a rational number
--------------------------------

Users of classical symbolic mathematics systems like Maple or Mathematica,
are accustomed to typing ``1/3`` and get the rational number one over three. In
SymPy this gives either ``0`` or a floating point number, depending on whether
we use old or new division. This is considered most disturbing difference
between SymPy and other mathematical systems.

First, this strange behavior comes from the fact that Python is a
general purpose programming language  and for a very long time it didn't
have support for rational numbers in the standard library. This changed
in Python 2.6, where the :class:`Fraction` class was introduced, but it would
be anyway unusual for Python to make ``/`` return a rational number.

To construct a rational number in SymPy, one can use :class:`Rational`
class::

    >>> r = Rational(1, 3)
    >>> r
    1/3

    >>> type(r)
    <class 'sympy.core.numbers.Rational'>

    >>> int(r)
    0
    >>> float(r)
    0.333333333333

    >>> r.evalf()
    0.333333333333333

There are also other ways::

    >>> Integer(1)/3
    1/3
    >>> S(1)/3
    1/3

``S`` is SymPy's registry of singletons. It implements the ``__call__`` method,
which is a shorthand for :func:`sympify`. Using ``S`` is the most concise
way to construct rational numbers. The last way is to pass a string with
``1/3`` to :func:`sympify`::

    >>> sympify("1/3")
    1/3
    >>> type(_)
    <class 'sympy.core.numbers.Rational'>

:func:`sympify` implements a :mod:`tokenize`--based preparser that puts
Python's numeric types in envelopes consisting of SymPy's numeric class
constructors.

You can also avoid this problem by not typing ``int/int`` when other
terms are involved. For example, write ``2*x/3`` instead of ``2/3*x``.
And you can type ``sqrt(x)`` instead of ``x**Rational(1, 2)``, as the
two are equivalent.

``^`` is not exponentiation operator
------------------------------------

SymPy uses the same default arithmetic operators as Python. Most of these,
like ``+``, ``-``, ``*`` and ``/``, are standard. There are, however, some
differences when comparing SymPy to standalone mathematical systems. One
of the differences is lack of implied multiplication, to which Mathematica
users may be accustomed::

    >>> var('x')

    >>> 2*x
    2*x

    >>> 2x
    Traceback (most recent call last):
    ...
    SyntaxError: invalid syntax

    >>> 2 x
    Traceback (most recent call last):
    ...
    SyntaxError: invalid syntax

More importantly, Python uses ``**`` to denote exponentiation, whereas
other mathematical systems use ``^`` operator. Notable exceptions to
this rule are Axiom and Maple, which allow both, though most users may
not be aware of this. For example in Mathematica, ``**`` operator is
used for non-commutative multiplication. So in Sympy the following
expression is perfectly valid::

    >>> (x + 1)**2
           2
    (x + 1)

    >>> type(_)
    <class 'sympy.core.power.Pow'>

but using ``^``::

    >>> (x + 1)^2
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand type(s) for ^: 'Add' and 'int'

gives use :exc:`TypeError`. For users' convenience, :func:`sympify` converts
``^`` to ``**`` by default in a string::

    >>> sympify("(x + 1)^2")
           2
    (x + 1)

    >>> type(_)
    <class 'sympy.core.power.Pow'>

People who what pure Python behavior of :func:`sympify` can disable this
automatic conversion by passing ``convert_xor=False`` to it.

``=`` is not comparison operator
--------------------------------

The equals sign (``=``) is the assignment operator in Python, not equality
operator. In other many mathematical systems, ``=`` is used for comparing
values and/or for constructing equalities, but with SymPy you have to use
``==`` for the former and ``Eq(x, y)`` for the later. Note that instances
of :class:`Eq` class, in boolean context, collapse to ``==``::

    >>> var('x,y')

    >>> x == y
    False

    >>> Eq(x, y)
    x = y
    >>> bool(_)
    False

Why you shouldn't write ``10**-1000``
-------------------------------------

Symbolic mathematics systems are expected to work with expressions of
arbitrary size, limited only by the size of available memory. Python
supports arbitrary precision integers by default, but allows only fixed
precision floats. Thus you can write::

    >>> 10**-10
    1e-10

but::

    >>> 10**-1000
    0.0

is not what we expect. To overcome this, we have to make the base an
instance of SymPy's floating point type::

    >>> Float(10.0)**-1000
    1.00000000000000e-1000

Note that we can't write simply ``Float(10)``, because SymPy automatically
converts this to an instance of :class:`Integer` class and thus::

    >>> type(Float(10)**-1000)
    <class 'sympy.core.numbers.Rational'>

Of course we could issue::

    >>> (Float(10)**-1000).evalf()
    1.00000000000000e-1000

but this it is neither readable, nor efficient.

You can also pass the entire number as a string to :class:`Float`. If you
do this, you must use the scientific notation syntax::

    >>> Float("1e-1000")
    1.00000000000000e-1000

Finally, we note that it is preferable to use exact (i.e., rational)
numbers when the values of the numbers are exactly known. Many parts of
SymPy work better when rational numbers are used instead of floating
point numbers. This is because rational numbers do not suffer from some
of the problems of floating point numbers, like rounding errors.

This is especially the case for exponents::

    >>> factor(x**2.0 - 1)
    x**2.0 - 1

    >>> factor(x**2 - 1)
    (x - 1)*(x + 1)

The first expression is not factored because the factorization only
holds for the exponent of `2` *exactly*. This problem can also come
up when using floating point coefficients::

    >>> solve([2*x + y**2, y - x], [x, y])
    [(-2, -2), (0, 0)]

    >>> solve([2.0*x + y**2, y - x], [x, y])
    Traceback (most recent call last):
    ...
    DomainError: can't compute a Groebner basis over RR

Here, the algorithm for solving systems of polynomial equations relies
on computing a |groebner| basis (see the :ref:`groebner-bases` section
below for more information on these). But the algorithm for computing
this currently does not support floating point coefficients, so
:func:`solve` fails in that case.

How to deal with limited recursion depth
----------------------------------------

Very often algorithms in symbolic mathematics and computer algebra are
highly recursive in nature. This can be a problem even for relatively
small inputs in SymPy, because Python interpreters set a limit on the
depth of recursion. Suppose we want to compute, manipulate and print the
following function composition:

.. math::

    \underbrace{(f \circ f \circ \ldots \circ f)}_{1000}(x)

Computing this isn't a problem::

    >>> f = Function('f')
    >>> x = Symbol('x')

    >>> u = x

    >>> for i in xrange(1000):
    ...     u = f(x)
    ...

    >>> type(u)
    f

However, if we try to get the number of all subexpressions of ``u`` that
contain ``f``, we get the following error::

    >>> len(u.find(f))
    Traceback (most recent call last):
    ...
    RuntimeError: maximum recursion depth exceeded while calling a Python object

The same happens when we try to print ``u``::

    >>> len([ c for c in str(u) if c == 'f' ])
    Traceback (most recent call last):
    ...
    RuntimeError: maximum recursion depth exceeded while calling a Python object

Python provides, at least partially, a solution to this problem by
allowing the user to relax the limit on recursion depth::

    >>> import sys
    >>> sys.setrecursionlimit(1050)

    >>> len(u.find(f))
    1000

To print ``u`` we have to relax the limit even more::

    >>> len([ c for c in str(u) if c == 'f' ])
    Traceback (most recent call last):
    ...
    RuntimeError: maximum recursion depth exceeded while calling a Python object

    >>> sys.setrecursionlimit(5500)

    >>> len([ c for c in str(u) if c == 'f' ])
    1000

This should be a warning about the fact that often it is possible to
perform computations with highly nested expressions, but it is not
possible to print those expressions without relaxing the recursion depth
limit. SymPy never uses ``sys.setrecursionlimit`` automatically, so
it's users responsibility to relax the limit whenever needed.

Unless you are using a highly nested expression like the one above, you
generally won't encounter this problem, as the default limit of 1000 is
generally high enough for the most common expressions.

Expression caching and its consequences
---------------------------------------

To improve speed of computations, SymPy by default caches all intermediate
subexpressions. The difference is easily visible when running tests:

.. parsed-literal::

    $ :input:`SYMPY_USE_CACHE=yes bin/test sympy/integrals/tests/test_risch.py`
    ============================= test process starts ==============================
    executable:   /usr/bin/python2.6  (2.6.6-final-0)
    architecture: 64-bit
    ground types: gmpy

    sympy/integrals/tests/test_risch.py[20] .....ffff...........                [OK]

    ======= tests finished: 16 passed, 4 expected to fail, in 28.18 seconds ========

    $ :input:`SYMPY_USE_CACHE=no bin/test sympy/integrals/tests/test_risch.py`
    ============================= test process starts ==============================
    executable:   /usr/bin/python2.6  (2.6.6-final-0)
    architecture: 64-bit
    ground types: gmpy

    sympy/integrals/tests/test_risch.py[20] .....ffff...........                [OK]

    ======= tests finished: 16 passed, 4 expected to fail, in 64.82 seconds ========

(note the time needed to run the tests at the end of the each test run)
and in interactive sessions:

.. parsed-literal::

    $ :input:`bin/isympy -q`
    IPython console for SymPy 0.7.1 (Python 2.7.1-64-bit) (ground types: gmpy)

    In [1]: f = (x-tan(x)) / tan(x)**2 + tan(x)

    In [2]: %time integrate(f, x);
    CPU times: user 0.46 s, sys: 0.00 s, total: 0.46 s
    Wall time: 0.49 s

    In [4]: %time integrate(f, x);
    CPU times: user 0.24 s, sys: 0.00 s, total: 0.24 s
    Wall time: 0.25 s

    $ :input:`bin/isympy -q -C`
    IPython console for SymPy 0.7.1 (Python 2.7.1-64-bit) (ground types: gmpy, cache: off)

    In [1]: f = (x-tan(x)) / tan(x)**2 + tan(x)

    In [2]: %time integrate(f, x);
    CPU times: user 1.82 s, sys: 0.00 s, total: 1.82 s
    Wall time: 1.84 s

    In [4]: %time integrate(f, x);
    CPU times: user 1.82 s, sys: 0.00 s, total: 1.82 s
    Wall time: 1.83 s

(``-C`` is equivalent to setting ``SYMPY_USE_CACHE="no"``).

The main consequence of caching is that SymPy can use a lot of resources
in certain situations. One can use :func:`clear_cache` to reduce memory
consumption:

.. sourcecode:: ipython

    In [6]: from sympy.core.cache import clear_cache

    In [7]: clear_cache()

    In [8]: %time integrate(f, x);
    CPU times: user 0.46 s, sys: 0.00 s, total: 0.46 s
    Wall time: 0.47 s

As caching influences computation times, any benchmarking must be performed
with cache off. Otherwise those measurements will be either inaccurate or
completely wrong (measuring how fast SymPy can retrieve data from cache,
rather than actual computing times):

.. sourcecode:: ipython

    $ bin/isympy -q
    IPython console for SymPy 0.7.1 (Python 2.7.1-64-bit) (ground types: gmpy)

    In [1]: %timeit sin(2*pi);
    10000 loops, best of 3: 28.7 us per loop

    $ bin/isympy -q -C
    IPython console for SymPy 0.7.1 (Python 2.7.1-64-bit) (ground types: gmpy, cache: off)

    In [1]: %timeit sin(2*pi);
    100 loops, best of 3: 2.75 ms per loop

The difference between using and not using cache is two orders of magnitude.

Naming convention of trigonometric inverses
-------------------------------------------

SymPy uses different names than most computer algebra systems for some
of the commonly used elementary functions. In particular, the inverse
trigonometric and hyperbolic functions use Python's naming convention,
so we have :func:`asin`, :func:`asinh`, :func:`acos` and so on, instead
of :func:`arcsin`, :func:`arcsinh`, :func:`arccos`, etc.

Container indices start at zero
-------------------------------

It should be obvious for people using Python, even for beginners, that when
indexing containers like ``list`` or ``tuple``, indexes start at zero, not
one::

    >>> L = symbols('x:5')
    >>> L
    (x₀, x₁, x₂, x₃, x₄)

    >>> L[0]
    x₀
    >>> L[1]
    x₁

This is a common thing in general purpose programming languages. However,
most symbolic mathematics systems, especially those which invent their own
mathematical programming language, use `1`--based indexing, sometimes reserving
the `0`--th index for special purpose (e.g. head of expressions in Mathematica).
