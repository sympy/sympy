Term rewriting
==============

Term rewriting is a very general class of functionalities which are used to
convert expressions of one type in terms of expressions of different kind. For
example expanding, combining and converting expressions apply to term
rewriting, and also simplification routines can be included here. Currently
SymPy has several functions and basic built-in methods for performing various
types of rewriting.

Expanding
---------

The simplest rewrite rule is expanding expressions into a _sparse_ form.
Expanding has several flavors and include expanding complex valued expressions,
arithmetic expand of products and powers but also expanding functions in terms
of more general functions is possible. Below are listed all currently available
expand rules.

Expanding of arithmetic expressions involving products and powers:
    >>> from sympy import *
    >>> x, y, z = symbols('x,y,z')
    >>> ((x + y)*(x - y)).expand(basic=True)
    x**2 - y**2
    >>> ((x + y + z)**2).expand(basic=True)
    x**2 + 2*x*y + 2*x*z + y**2 + 2*y*z + z**2

Arithmetic expand is done by default in ``expand()`` so the keyword ``basic`` can
be omitted. However you can set ``basic=False`` to avoid this type of expand if
you use rules described below. This give complete control on what is done with
the expression.

Another type of expand rule is expanding complex valued expressions and putting
them into a normal form. For this ``complex`` keyword is used. Note that it will
always perform arithmetic expand to obtain the desired normal form:

    >>> (x + I*y).expand(complex=True)
    re(x) + I*re(y) + I*im(x) - im(y)

    >>> sin(x + I*y).expand(complex=True)
    sin(re(x) - im(y))*cosh(re(y) + im(x)) + I*cos(re(x) - im(y))*sinh(re(y) + im(x))

Note also that the same behavior can be obtained by using ``as_real_imag()``
method. However it will return a tuple containing the real part in the first
place and the imaginary part in the other. This can be also done in a two step
process by using ``collect`` function:

    >>> (x + I*y).as_real_imag()
    (re(x) - im(y), re(y) + im(x))

    >>> collect((x + I*y).expand(complex=True), I, evaluate=False)
    {1: re(x) - im(y), I: re(y) + im(x)}

There is also possibility for expanding expressions in terms of expressions of
different kind. This is very general type of expanding and usually you would
use ``rewrite()`` to do specific type of rewrite::

    >>> GoldenRatio.expand(func=True)
    1/2 + sqrt(5)/2

Common Subexpression Detection and Collection
---------------------------------------------

.. module:: sympy.simplify.cse_main

Before evaluating a large expression, it is often useful to identify common
subexpressions, collect them and evaluate them at once. This is implemented
in the ``cse`` function. Examples::

    >>> from sympy import cse, sqrt, sin, pprint
    >>> from sympy.abc import x

    >>> pprint(cse(sqrt(sin(x))), use_unicode=True)
    ⎛    ⎡  ________⎤⎞
    ⎝[], ⎣╲╱ sin(x) ⎦⎠

    >>> pprint(cse(sqrt(sin(x)+5)*sqrt(sin(x)+4)), use_unicode=True)
    ⎛                ⎡  ________   ________⎤⎞
    ⎝[(x₀, sin(x))], ⎣╲╱ x₀ + 4 ⋅╲╱ x₀ + 5 ⎦⎠

    >>> pprint(cse(sqrt(sin(x+1) + 5 + cos(y))*sqrt(sin(x+1) + 4 + cos(y))),
    ...     use_unicode=True)
    ⎛                             ⎡  ________   ________⎤⎞
    ⎝[(x₀, sin(x + 1) + cos(y))], ⎣╲╱ x₀ + 4 ⋅╲╱ x₀ + 5 ⎦⎠

    >>> pprint(cse((x-y)*(z-y) + sqrt((x-y)*(z-y))), use_unicode=True)
    ⎛                                     ⎡  ____     ⎤⎞
    ⎝[(x₀, -y), (x₁, (x + x₀)⋅(x₀ + z))], ⎣╲╱ x₁  + x₁⎦⎠

Optimizations to be performed before and after common subexpressions
elimination can be passed in the``optimizations`` optional argument. A set of
predefined basic optimizations can be applied by passing
``optimizations='basic'``::

    >>> pprint(cse((x-y)*(z-y) + sqrt((x-y)*(z-y)), optimizations='basic'),
    ...     use_unicode=True)
    ⎛                          ⎡  ____     ⎤⎞
    ⎝[(x₀, -(x - y)⋅(y - z))], ⎣╲╱ x₀  + x₀⎦⎠

However, these optimizations can be very slow for large expressions. Moreover,
if speed is a concern, one can pass the option ``order='none'``. Order of
terms will then be dependent on hashing algorithm implementation, but speed
will be greatly improved.

More information:

.. autofunction:noindex: cse
