.. _polys-wester:

==============================
Examples from Wester's Article
==============================

Introduction
============

In this tutorial we present examples from Wester's article concerning
comparison and critique of mathematical abilities of several computer
algebra systems (see [Wester1999]_). All the examples are related to
polynomial and algebraic computations and SymPy specific remarks were
added to all of them.

Examples
========

All examples in this tutorial are computable, so one can just copy and
paste them into a Python shell and do something useful with them. All
computations were done using the following setup::

    >>> from sympy import *

    >>> init_printing(use_unicode=True, wrap_line=False, no_global=True)

    >>> var('x,y,z,s,c,n')
    (x, y, z, s, c, n)

Simple univariate polynomial factorization
------------------------------------------

To obtain a factorization of a polynomial use :func:`factor` function.
By default :func:`factor` returns the result in unevaluated form, so the
content of the input polynomial is left unexpanded, as in the following
example::

    >>> factor(6*x - 10)
    2⋅(3⋅x - 5)

To achieve the same effect in a more systematic way use :func:`primitive`
function, which returns the content and the primitive part of the input
polynomial::

    >>> primitive(6*x - 10)
    (2, 3⋅x - 5)

.. note::

    The content and the primitive part can be computed only over a ring. To
    simplify coefficients of a polynomial over a field use :func:`monic`.

Univariate GCD, resultant and factorization
-------------------------------------------

Consider univariate polynomials ``f``, ``g`` and ``h`` over integers::

    >>> f = 64*x**34 - 21*x**47 - 126*x**8 - 46*x**5 - 16*x**60 - 81
    >>> g = 72*x**60 - 25*x**25 - 19*x**23 - 22*x**39 - 83*x**52 + 54*x**10 + 81
    >>> h = 34*x**19 - 25*x**16 + 70*x**7 + 20*x**3 - 91*x - 86

We can compute the greatest common divisor (GCD) of two polynomials using
:func:`gcd` function::

    >>> gcd(f, g)
    1

We see that ``f`` and ``g`` have no common factors. However, ``f*h`` and ``g*h``
have an obvious factor ``h``::

    >>> gcd(expand(f*h), expand(g*h)) - h
    0

The same can be verified using the resultant of univariate polynomials::

    >>> resultant(expand(f*h), expand(g*h))
    0

Factorization of large univariate polynomials (of degree 120 in this case) over
integers is also possible::

    >>> factor(expand(f*g))
     ⎛    60       47       34        8       5     ⎞ ⎛    60       52     39       25       23       10     ⎞
    -⎝16⋅x   + 21⋅x   - 64⋅x   + 126⋅x  + 46⋅x  + 81⎠⋅⎝72⋅x   - 83⋅x - 22⋅x   - 25⋅x   - 19⋅x   + 54⋅x   + 81⎠

Multivariate GCD and factorization
----------------------------------

What can be done in univariate case, can be also done for multivariate
polynomials. Consider the following polynomials ``f``, ``g`` and ``h``
in `\mathbb{Z}[x,y,z]`::

    >>> f = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    >>> g = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    >>> h = 11*x**12*y**7*z**13 - 23*x**2*y**8*z**10 + 47*x**17*y**5*z**8

As previously, we can verify that ``f`` and ``g`` have no common factors::

    >>> gcd(f, g)
    1

However, ``f*h`` and ``g*h`` have an obvious factor ``h``::

    >>> gcd(expand(f*h), expand(g*h)) - h
    0

Multivariate factorization of large polynomials is also possible::

    >>> factor(expand(f*g))
        7   ⎛   9  9  3       7  6       5    12       7⎞ ⎛   22       17  5  8      15  9  2         19  8    ⎞
    -2⋅y ⋅z⋅⎝6⋅x ⋅y ⋅z  + 10⋅x ⋅z  + 17⋅x ⋅y⋅z   + 40⋅y ⎠⋅⎝3⋅x   + 47⋅x  ⋅y ⋅z  - 6⋅x  ⋅y ⋅z  - 24⋅x⋅y  ⋅z  - 5⎠

Support for symbols in exponents
--------------------------------

Polynomial manipulation functions provided by :mod:`sympy.polys` are mostly
used with integer exponents. However, it's perfectly valid to compute with
symbolic exponents, e.g.::

    >>> gcd(2*x**(n + 4) - x**(n + 2), 4*x**(n + 1) + 3*x**n)
     n
    x

Testing if polynomials have common zeros
----------------------------------------

To test if two polynomials have a root in common we can use :func:`resultant`
function. The theory says that the resultant of two polynomials vanishes if
there is a common zero of those polynomials. For example::

    >>> resultant(3*x**4 + 3*x**3 + x**2 - x - 2, x**3 - 3*x**2 + x + 5)
    0

We can visualize this fact by factoring the polynomials::

    >>> factor(3*x**4 + 3*x**3 + x**2 - x - 2)
            ⎛   3        ⎞
    (x + 1)⋅⎝3⋅x  + x - 2⎠

    >>> factor(x**3 - 3*x**2 + x + 5)
            ⎛ 2          ⎞
    (x + 1)⋅⎝x  - 4⋅x + 5⎠

In both cases we obtained the factor `x + 1` which tells us that the common
root is `x = -1`.

Normalizing simple rational functions
-------------------------------------

To remove common factors from the numerator and the denominator of a rational
function the elegant way, use :func:`cancel` function. For example::

    >>> cancel((x**2 - 4)/(x**2 + 4*x + 4))
    x - 2
    ─────
    x + 2

Expanding expressions and factoring back
----------------------------------------

One can work easily we expressions in both expanded and factored forms.
Consider a polynomial ``f`` in expanded form. We differentiate it and
factor the result back::

    >>> f = expand((x + 1)**20)

    >>> g = diff(f, x)

    >>> factor(g)
              19
    20⋅(x + 1)

The same can be achieved in factored form::

    >>> diff((x + 1)**20, x)
              19
    20⋅(x + 1)

Factoring in terms of cyclotomic polynomials
--------------------------------------------

SymPy can very efficiently decompose polynomials of the form `x^n \pm 1` in
terms of cyclotomic polynomials::

    >>> factor(x**15 - 1)
            ⎛ 2        ⎞ ⎛ 4    3    2        ⎞ ⎛ 8    7    5    4    3       ⎞
    (x - 1)⋅⎝x  + x + 1⎠⋅⎝x  + x  + x  + x + 1⎠⋅⎝x  - x  + x  - x  + x - x + 1⎠

The original Wester`s example was `x^{100} - 1`, but was truncated for
readability purpose. Note that this is not a big struggle for :func:`factor`
to decompose polynomials of degree 1000 or greater.

Univariate factoring over Gaussian numbers
------------------------------------------

Consider a univariate polynomial ``f`` with integer coefficients::

    >>> f = 4*x**4 + 8*x**3 + 77*x**2 + 18*x + 153

We want to obtain a factorization of ``f`` over Gaussian numbers. To do this
we use :func:`factor` as previously, but this time we set ``gaussian`` keyword
to ``True``::

    >>> factor(f, gaussian=True)
      ⎛    3⋅ⅈ⎞ ⎛    3⋅ⅈ⎞
    4⋅⎜x - ───⎟⋅⎜x + ───⎟⋅(x + 1 - 4⋅ⅈ)⋅(x + 1 + 4⋅ⅈ)
      ⎝     2 ⎠ ⎝     2 ⎠

As the result we got a splitting factorization of ``f`` with monic factors
(this is a general rule when computing in a field with SymPy). The ``gaussian``
keyword is useful for improving code readability, however the same result can
be computed using more general syntax::

    >>> factor(f, extension=I)
      ⎛    3⋅ⅈ⎞ ⎛    3⋅ⅈ⎞
    4⋅⎜x - ───⎟⋅⎜x + ───⎟⋅(x + 1 - 4⋅ⅈ)⋅(x + 1 + 4⋅ⅈ)
      ⎝     2 ⎠ ⎝     2 ⎠

Computing with automatic field extensions
-----------------------------------------

Consider two univariate polynomials ``f`` and ``g``::

    >>> f = x**3 + (sqrt(2) - 2)*x**2 - (2*sqrt(2) + 3)*x - 3*sqrt(2)
    >>> g = x**2 - 2

We would like to reduce degrees of the numerator and the denominator of a
rational function ``f/g``. Do do this we employ :func:`cancel` function::

    >>> cancel(f/g)
     3      2       2
    x  - 2⋅x  + √2⋅x  - 3⋅x - 2⋅√2⋅x - 3⋅√2
    ───────────────────────────────────────
                      2
                     x  - 2

Unfortunately nothing interesting happened. This is because by default SymPy
treats `\sqrt{2}` as a generator, obtaining a bivariate polynomial for the
numerator. To make :func:`cancel` recognize algebraic properties of `\sqrt{2}`,
one needs to use ``extension`` keyword::

    >>> cancel(f/g, extension=True)
     2
    x  - 2⋅x - 3
    ────────────
       x - √2

Setting ``extension=True`` tells :func:`cancel` to find minimal algebraic
number domain for the coefficients of ``f/g``. The automatically inferred
domain is `\mathbb{Q}(\sqrt{2})`. If one doesn't want to rely on automatic
inference, the same result can be obtained by setting the ``extension``
keyword with an explicit algebraic number::

    >>> cancel(f/g, extension=sqrt(2))
     2
    x  - 2⋅x - 3
    ────────────
       x - √2

Univariate factoring over various domains
-----------------------------------------

Consider a univariate polynomial ``f`` with integer coefficients::

    >>> f = x**4 - 3*x**2 + 1

With :mod:`sympy.polys` we can obtain factorizations of ``f`` over different
domains, which includes:

* rationals::

    >>> factor(f)
    ⎛ 2        ⎞ ⎛ 2        ⎞
    ⎝x  - x - 1⎠⋅⎝x  + x - 1⎠

* finite fields::

    >>> factor(f, modulus=5)
           2        2
    (x - 2) ⋅(x + 2)

* algebraic numbers::

    >>> alg = AlgebraicNumber((sqrt(5) - 1)/2, alias='alpha')

    >>> factor(f, extension=alg)
    (x - α)⋅(x + α)⋅(x - α - 1)⋅(x + α + 1)

Factoring polynomials into linear factors
-----------------------------------------

Currently SymPy can factor polynomials into irreducibles over various domains,
which can result in a splitting factorization (into linear factors). However,
there is currently no systematic way to infer a splitting field (algebraic
number field) automatically. In future the following syntax will be
implemented::

    >>> factor(x**3 + x**2 - 7, split=True)
    Traceback (most recent call last):
    ...
    NotImplementedError: 'split' option is not implemented yet

Note this is different from ``extension=True``, because the later only tells how
expression parsing should be done, not what should be the domain of computation.
One can simulate the ``split`` keyword for several classes of polynomials using
:func:`solve` function.

Advanced factoring over finite fields
-------------------------------------

Consider a univariate polynomial ``f`` with integer coefficients::

    >>> f = x**11 + x + 1

We can factor ``f`` over a large finite field `F_{65537}`::

    >>> factor(f, modulus=65537)
    ⎛ 2        ⎞ ⎛ 9    8    6    5    3    2    ⎞
    ⎝x  + x + 1⎠⋅⎝x  - x  + x  - x  + x  - x  + 1⎠

and expand the resulting factorization back::

    >>> expand(_)
     11
    x   + x + 1

obtaining polynomial ``f``. This was done using symmetric polynomial
representation over finite fields The same thing can be done using
non-symmetric representation::

    >>> factor(f, modulus=65537, symmetric=False)
    ⎛ 2        ⎞ ⎛ 9          8    6          5    3          2    ⎞
    ⎝x  + x + 1⎠⋅⎝x  + 65536⋅x  + x  + 65536⋅x  + x  + 65536⋅x  + 1⎠

As with symmetric representation we can expand the factorization
to get the input polynomial back. This time, however, we need to
truncate coefficients of the expanded polynomial modulo 65537::

    >>> trunc(expand(_), 65537)
     11
    x   + x + 1

Working with expressions as polynomials
---------------------------------------

Consider a multivariate polynomial ``f`` in `\mathbb{Z}[x,y,z]`::

    >>> f = expand((x - 2*y**2 + 3*z**3)**20)

We want to compute factorization of ``f``. To do this we use ``factor`` as
usually, however we note that the polynomial in consideration is already
in expanded form, so we can tell the factorization routine to skip
expanding ``f``::

    >>> factor(f, expand=False)
                     20
    ⎛       2      3⎞
    ⎝x - 2⋅y  + 3⋅z ⎠

The default in :mod:`sympy.polys` is to expand all expressions given as
arguments to polynomial manipulation functions and :class:`Poly` class.
If we know that expanding is unnecessary, then by setting ``expand=False``
we can save quite a lot of time for complicated inputs. This can be really
important when computing with expressions like::

    >>> g = expand((sin(x) - 2*cos(y)**2 + 3*tan(z)**3)**20)

    >>> factor(g, expand=False)
                                     20
    ⎛               2           3   ⎞
    ⎝-sin(x) + 2⋅cos (y) - 3⋅tan (z)⎠

Computing reduced Gröbner bases
-------------------------------

To compute a reduced Gröbner basis for a set of polynomials use
:func:`groebner` function. The function accepts various monomial
orderings, e.g.: ``lex``, ``grlex`` and ``grevlex``, or a user
defined one, via ``order`` keyword. The ``lex`` ordering is the
most interesting because it has elimination property, which means
that if the system of polynomial equations to :func:`groebner` is
zero-dimensional (has finite number of solutions) the last element
of the basis is a univariate polynomial. Consider the following example::

    >>> f = expand((1 - c**2)**5 * (1 - s**2)**5 * (c**2 + s**2)**10)

    >>> groebner([f, c**2 + s**2 - 1])
                 ⎛⎡ 2    2       20      18       16       14      12    10⎤                           ⎞
    GroebnerBasis⎝⎣c  + s  - 1, c   - 5⋅c   + 10⋅c   - 10⋅c   + 5⋅c   - c  ⎦, s, c, domain=ℤ, order=lex⎠

The result is an ordinary Python list, so we can easily apply a function to
all its elements, for example we can factor those elements::

    >>> list(map(factor, _))
    ⎡ 2    2       10        5        5⎤
    ⎣c  + s  - 1, c  ⋅(c - 1) ⋅(c + 1) ⎦

From the above we can easily find all solutions of the system of polynomial
equations. Or we can use :func:`solve` to achieve this in a more systematic
way::

    >>> solve([f, s**2 + c**2 - 1], c, s)
    [(-1, 0), (0, -1), (0, 1), (1, 0)]

Multivariate factoring over algebraic numbers
---------------------------------------------

Computing with multivariate polynomials over various domains is as simple as
in univariate case. For example consider the following factorization over
`\mathbb{Q}(\sqrt{-3})`::

    >>> factor(x**3 + y**3, extension=sqrt(-3))
            ⎛      ⎛  1   √3⋅ⅈ⎞⎞ ⎛      ⎛  1   √3⋅ⅈ⎞⎞
    (x + y)⋅⎜x + y⋅⎜- ─ - ────⎟⎟⋅⎜x + y⋅⎜- ─ + ────⎟⎟
            ⎝      ⎝  2    2  ⎠⎠ ⎝      ⎝  2    2  ⎠⎠

.. note:: Currently multivariate polynomials over finite fields aren't supported.

Partial fraction decomposition
------------------------------

Consider a univariate rational function ``f`` with integer coefficients::

    >>> f = (x**2 + 2*x + 3)/(x**3 + 4*x**2 + 5*x + 2)

To decompose ``f`` into partial fractions use :func:`apart` function::

    >>> apart(f)
      3       2        2
    ───── - ───── + ────────
    x + 2   x + 1          2
                    (x + 1)

To return from partial fractions to the rational function use
a composition of :func:`together` and :func:`cancel`::

    >>> cancel(together(_))
         2
        x  + 2⋅x + 3
    ───────────────────
     3      2
    x  + 4⋅x  + 5⋅x + 2

Literature
==========

.. [Wester1999] Michael J. Wester, A Critique of the Mathematical Abilities of
    CA Systems, 1999, `<http://www.math.unm.edu/~wester/cas/book/Wester.pdf>`_
