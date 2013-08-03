Basic usage
===========================

In interactive code examples that follow, it will be assumed that
all items in the ``mpmath`` namespace have been imported::

    >>> from mpmath import *

Importing everything can be convenient, especially when using mpmath interactively, but be
careful when mixing mpmath with other libraries! To avoid inadvertently overriding
other functions or objects, explicitly import only the needed objects, or use
the ``mpmath.`` or ``mp.`` namespaces::

    from mpmath import sin, cos
    sin(1), cos(1)

    import mpmath
    mpmath.sin(1), mpmath.cos(1)

    from mpmath import mp    # mp context object -- to be explained
    mp.sin(1), mp.cos(1)


Number types
------------

Mpmath provides the following numerical types:

    +------------+----------------+
    | Class      | Description    |
    +============+================+
    | ``mpf``    | Real float     |
    +------------+----------------+
    | ``mpc``    | Complex float  |
    +------------+----------------+
    | ``matrix`` | Matrix         |
    +------------+----------------+

The following section will provide a very short introduction to the types ``mpf`` and ``mpc``. Intervals and matrices are described further in the documentation chapters on interval arithmetic and matrices / linear algebra.

The ``mpf`` type is analogous to Python's built-in ``float``. It holds a real number or one of the special values ``inf`` (positive infinity), ``-inf`` (negative infinity) and ``nan`` (not-a-number, indicating an indeterminate result). You can create ``mpf`` instances from strings, integers, floats, and other ``mpf`` instances:

    >>> mpf(4)
    mpf('4.0')
    >>> mpf(2.5)
    mpf('2.5')
    >>> mpf("1.25e6")
    mpf('1250000.0')
    >>> mpf(mpf(2))
    mpf('2.0')
    >>> mpf("inf")
    mpf('+inf')

The ``mpc`` type represents a complex number in rectangular form as a pair of ``mpf`` instances. It can be constructed from a Python ``complex``, a real number, or a pair of real numbers:

    >>> mpc(2,3)
    mpc(real='2.0', imag='3.0')
    >>> mpc(complex(2,3)).imag
    mpf('3.0')

You can mix ``mpf`` and ``mpc`` instances with each other and with Python numbers:

    >>> mpf(3) + 2*mpf('2.5') + 1.0
    mpf('9.0')
    >>> mp.dps = 15      # Set precision (see below)
    >>> mpc(1j)**0.5
    mpc(real='0.70710678118654757', imag='0.70710678118654757')


Setting the precision
---------------------

Mpmath uses a global working precision; it does not keep track of the precision or accuracy of individual numbers. Performing an arithmetic operation or calling ``mpf()`` rounds the result to the current working precision. The working precision is controlled by a context object called ``mp``, which has the following default state:

    >>> print mp
    Mpmath settings:
      mp.prec = 53                [default: 53]
      mp.dps = 15                 [default: 15]
      mp.trap_complex = False     [default: False]

The term **prec** denotes the binary precision (measured in bits) while **dps** (short for *decimal places*) is the decimal precision. Binary and decimal precision are related roughly according to the formula ``prec = 3.33*dps``. For example, it takes a precision of roughly 333 bits to hold an approximation of pi that is accurate to 100 decimal places (actually slightly more than 333 bits is used).

Changing either precision property of the ``mp`` object automatically updates the other; usually you just want to change the ``dps`` value:

    >>> mp.dps = 100
    >>> mp.dps
    100
    >>> mp.prec
    336

When the precision has been set, all ``mpf`` operations are carried out at that precision::

    >>> mp.dps = 50
    >>> mpf(1) / 6
    mpf('0.16666666666666666666666666666666666666666666666666656')
    >>> mp.dps = 25
    >>> mpf(2) ** mpf('0.5')
    mpf('1.414213562373095048801688713')

The precision of complex arithmetic is also controlled by the ``mp`` object:

    >>> mp.dps = 10
    >>> mpc(1,2) / 3
    mpc(real='0.3333333333321', imag='0.6666666666642')

There is no restriction on the magnitude of numbers. An ``mpf`` can for example hold an approximation of a large Mersenne prime:

    >>> mp.dps = 15
    >>> print mpf(2)**32582657 - 1
    1.24575026015369e+9808357

Or why not 1 googolplex:

    >>> print mpf(10) ** (10**100)  # doctest:+ELLIPSIS
    1.0e+100000000000000000000000000000000000000000000000000...

The (binary) exponent is stored exactly and is independent of the precision.

Temporarily changing the precision
..................................

It is often useful to change the precision during only part of a calculation. A way to temporarily increase the precision and then restore it is as follows:

    >>> mp.prec += 2
    >>> # do_something()
    >>> mp.prec -= 2

In Python 2.5, the ``with`` statement along with the mpmath functions ``workprec``, ``workdps``, ``extraprec`` and ``extradps`` can be used to temporarily change precision in a more safe manner:

    >>> from __future__ import with_statement
    >>> with workdps(20):  # doctest: +SKIP
    ...     print mpf(1)/7
    ...     with extradps(10):
    ...         print mpf(1)/7
    ...
    0.14285714285714285714
    0.142857142857142857142857142857
    >>> mp.dps
    15

The ``with`` statement ensures that the precision gets reset when exiting the block, even in the case that an exception is raised. (The effect of the ``with`` statement can be emulated in Python 2.4 by using a ``try/finally`` block.)

The ``workprec`` family of functions can also be used as function decorators:

    >>> @workdps(6)
    ... def f():
    ...     return mpf(1)/3
    ...
    >>> f()
    mpf('0.33333331346511841')


Some functions accept the ``prec`` and ``dps`` keyword arguments and this will override the global working precision. Note that this will not affect the precision at which the result is printed, so to get all digits, you must either use increase precision afterward when printing or use ``nstr``/``nprint``:

    >>> mp.dps = 15
    >>> print exp(1)
    2.71828182845905
    >>> print exp(1, dps=50)      # Extra digits won't be printed
    2.71828182845905
    >>> nprint(exp(1, dps=50), 50)
    2.7182818284590452353602874713526624977572470937

Finally, instead of using the global context object ``mp``, you can create custom contexts and work with methods of those instances instead of global functions. The working precision will be local to each context object:

    >>> mp2 = mp.clone()
    >>> mp.dps = 10
    >>> mp2.dps = 20
    >>> print mp.mpf(1) / 3
    0.3333333333
    >>> print mp2.mpf(1) / 3
    0.33333333333333333333

**Note**: the ability to create multiple contexts is a new feature that is only partially implemented. Not all mpmath functions are yet available as context-local methods. In the present version, you are likely to encounter bugs if you try mixing different contexts.

Providing correct input
-----------------------

Note that when creating a new ``mpf``, the value will at most be as accurate as the input. *Be careful when mixing mpmath numbers with Python floats*. When working at high precision, fractional ``mpf`` values should be created from strings or integers:

    >>> mp.dps = 30
    >>> mpf(10.9)   # bad
    mpf('10.9000000000000003552713678800501')
    >>> mpf('10.9')  # good
    mpf('10.8999999999999999999999999999997')
    >>> mpf(109) / mpf(10)   # also good
    mpf('10.8999999999999999999999999999997')
    >>> mp.dps = 15

(Binary fractions such as 0.5, 1.5, 0.75, 0.125, etc, are generally safe as input, however, since those can be represented exactly by Python floats.)

Printing
--------

By default, the ``repr()`` of a number includes its type signature. This way ``eval`` can be used to recreate a number from its string representation:

    >>> eval(repr(mpf(2.5)))
    mpf('2.5')

Prettier output can be obtained by using ``str()`` or ``print``, which hide the ``mpf`` and ``mpc`` signatures and also suppress rounding artifacts in the last few digits:

    >>> mpf("3.14159")
    mpf('3.1415899999999999')
    >>> print mpf("3.14159")
    3.14159
    >>> print mpc(1j)**0.5
    (0.707106781186548 + 0.707106781186548j)

Setting the ``mp.pretty`` option will use the ``str()``-style output for ``repr()`` as well:

    >>> mp.pretty = True
    >>> mpf(0.6)
    0.6
    >>> mp.pretty = False
    >>> mpf(0.6)
    mpf('0.59999999999999998')

The number of digits with which numbers are printed by default is determined by the working precision. To specify the number of digits to show without changing the working precision, use :func:`mpmath.nstr` and :func:`mpmath.nprint`:

    >>> a = mpf(1) / 6
    >>> a
    mpf('0.16666666666666666')
    >>> nstr(a, 8)
    '0.16666667'
    >>> nprint(a, 8)
    0.16666667
    >>> nstr(a, 50)
    '0.16666666666666665741480812812369549646973609924316'
