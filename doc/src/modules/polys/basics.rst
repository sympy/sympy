.. _polys-basics:

=================================
Basic functionality of the module
=================================

Introduction
============

This tutorial tries to give an overview of the functionality concerning
polynomials within SymPy. All code examples assume::

    >>> from sympy import *
    >>> x, y, z = symbols('x,y,z')
    >>> init_printing(use_unicode=False, wrap_line=False, no_global=True)

Basic concepts
==============

Polynomials
-----------

Given a family `(x_i)` of symbols, or other suitable objects,
expressions derived from them by repeated addition, subtraction
and multiplication are called *polynomial expressions in the
generators* `x_i`.

By the distributive law it is possible to do perform
multiplications before additions and subtractions.
The products of generators thus obtained are called
*monomials*. They are usually written in the form
`x_1^{\nu_1}x_2^{\nu_2}\cdots x_n^{\nu_n}` where the exponents `\nu_i`
are nonnegative integers. It is often convenient to write this briefly
as `x^\nu` where `x = (x_1, x_2, \ldots, x_n)` denotes the family of
generators and `\nu = (\nu_1, \nu_2, \ldots, \nu_n)` is the
family of exponents.

When all monomials having the same exponents are combined, the polynomial
becomes a sum of products `c_\nu x^\nu`, called the *terms* of the polynomial,
where the *coefficients* `c_\nu` are integers.
If some of the `x_i` are integers, they are incorporated
in the coefficients and not regarded as generators.

The set of all polynomials
in the generators `x_i` is a *ring*, i.e., the sums, differences and
products of its elements are again polynomials in the same generators.
This ring is denoted `\mathbb{Z}[x_1, x_2, \ldots, x_n]`, or
`\mathbb{Z}[(x_i)]`, and called
the *ring of polynomials in the* `x_i` *with integer coefficients*.

Polynomial expressions can be transformed into polynomials by the
method :func:`as_poly`::

    >>> e = (x + y)*(y - 2*z)
    >>> e.as_poly()
    Poly(x*y - 2*x*z + y**2 - 2*y*z, x, y, z, domain='ZZ')

More generally, it is possible to construct polynomials with
other types of coefficients.
If a polynomial expression contains numbers that are not integers,
they are regarded as coefficients and the coefficient ring is
extended accordingly. In particular, division by integers
leads to rational coefficients::

    >>> e = (3*x/2 + y)*(z - 1)
    >>> e.as_poly()
    Poly(3/2*x*z - 3/2*x + y*z - y, x, y, z, domain='QQ')

Symbolic numbers are considered generators unless they are explicitly
excluded, in which case they are adjoined to the coefficient ring::

    >>> e = (x + 2*pi)*y
    >>> e.as_poly()
    Poly(x*y + 2*y*pi, x, y, pi, domain='ZZ')
    >>> e.as_poly(x, y)
    Poly(x*y + 2*pi*y, x, y, domain='EX')

Here ``EX`` denotes the domain of all expressions. Since operations
in ``EX`` are slow, it is often advisable to define the coefficient
ring explicitly::

    >>> e = (x + 2*pi)*y
    >>> e.as_poly(domain=ZZ[pi])
    Poly(x*y + 2*pi*y, x, y, domain='ZZ[pi]')

Note that the ring `\mathbb{Z}[\pi][x, y]` of polynomials in `x` and `y`
with coefficients in `\mathbb{Z}[\pi]` is mathematically equivalent to
`\mathbb{Z}[\pi, x, y]`, only their implementations differ.

If an expression contains functions of the generators, other
than their positive integer powers, these are interpreted as new
generators::

    >>> e = x*sin(y) - y
    >>> e.as_poly()
    Poly(x*(sin(y)) - y, x, y, sin(y), domain='ZZ')

Since `y` and `\sin(y)` are algebraically independent they can both
appear as generators in a polynomial. However, *polynomial expressions
must not contain negative powers of generators*::

    >>> e = x - 1/x
    >>> e.as_poly()
    Poly(x - (1/x), x, 1/x, domain='ZZ')

It is important to realize that the generators `x` and `1/x = x^{-1}` are
treated as algebraically independent variables. In particular, their product
is not equal to 1. Hence *generators in denominators should be avoided even
if they raise no error*. Similar problems emerge with
rational powers of generators. So, for example, ``x`` and
``sqrt(x) = x**(1/2)`` are not recognized as algebraically dependent.

If there are algebraic numbers in an expression, it is possible to
adjoin them to the coefficient ring by setting the keyword ``extension``::

    >>> e = x + sqrt(2)
    >>> e.as_poly()
    Poly(x + (sqrt(2)), x, sqrt(2), domain='ZZ')
    >>> e.as_poly(extension=True)
    Poly(x + sqrt(2), x, domain='QQ<sqrt(2)>')

With the default setting ``extension=False`` `x` and `\sqrt 2` are
wrongly considered algebraically independent. With coefficients in
the extension field `\mathbb{Q}(\sqrt 2)` the square root is treated
properly as an algebraic number.


Basic functionality
===================

These functions provide different algorithms dealing with polynomials in the
form of SymPy expression, like symbols, sums etc.

Division
--------

The function :func:`div` provides division of polynomials with remainder.
That is, for polynomials ``f`` and ``g``, it computes ``q`` and ``r``, such
that `f = g \cdot q + r` and `\deg(r) < q`. For polynomials in one variables
with coefficients in a field, say, the rational numbers, ``q`` and ``r`` are
uniquely defined this way::

    >>> f = 5*x**2 + 10*x + 3
    >>> g = 2*x + 2

    >>> q, r = div(f, g, domain='QQ')
    >>> q
    5*x   5
    --- + -
     2    2
    >>> r
    -2
    >>> (q*g + r).expand()
       2
    5*x  + 10*x + 3

As you can see, ``q`` has a non-integer coefficient. If you want to do division
only in the ring of polynomials with integer coefficients, you can specify an
additional parameter::

    >>> q, r = div(f, g, domain='ZZ')
    >>> q
    0
    >>> r
       2
    5*x  + 10*x + 3

But be warned, that this ring is no longer Euclidean and that the degree of the
remainder doesn't need to be smaller than that of ``f``. Since 2 doesn't divide 5,
`2 x` doesn't divide `5 x^2`, even if the degree is smaller. But::

    >>> g = 5*x + 1

    >>> q, r = div(f, g, domain='ZZ')
    >>> q
    x
    >>> r
    9*x + 3
    >>> (q*g + r).expand()
       2
    5*x  + 10*x + 3

This also works for polynomials with multiple variables::

    >>> f = x*y + y*z
    >>> g = 3*x + 3*z

    >>> q, r = div(f, g, domain='QQ')
    >>> q
    y
    -
    3
    >>> r
    0

In the last examples, all of the three variables ``x``, ``y`` and ``z`` are
assumed to be variables of the polynomials. But if you have some unrelated
constant as coefficient, you can specify the variables explicitly::

    >>> a, b, c = symbols('a,b,c')
    >>> f = a*x**2 + b*x + c
    >>> g = 3*x + 2
    >>> q, r = div(f, g, domain='QQ')
    >>> q
    a*x   2*a   b
    --- - --- + -
     3     9    3

    >>> r
    4*a   2*b
    --- - --- + c
     9     3

GCD and LCM
-----------

With division, there is also the computation of the greatest common divisor and
the least common multiple.

When the polynomials have integer coefficients, the contents' gcd is also
considered::

    >>> f = (12*x + 12)*x
    >>> g = 16*x**2
    >>> gcd(f, g)
    4*x

But if the polynomials have rational coefficients, then the returned polynomial is
monic::

    >>> f = 3*x**2/2
    >>> g = 9*x/4
    >>> gcd(f, g)
    x

It also works with multiple variables. In this case, the variables are ordered
alphabetically, be default, which has influence on the leading coefficient::

    >>> f = x*y/2 + y**2
    >>> g = 3*x + 6*y

    >>> gcd(f, g)
    x + 2*y

The lcm is connected with the gcd and one can be computed using the other::

    >>> f = x*y**2 + x**2*y
    >>> g = x**2*y**2
    >>> gcd(f, g)
    x*y
    >>> lcm(f, g)
     3  2    2  3
    x *y  + x *y
    >>> (f*g).expand()
     4  3    3  4
    x *y  + x *y
    >>> (gcd(f, g, x, y)*lcm(f, g, x, y)).expand()
     4  3    3  4
    x *y  + x *y

Square-free factorization
-------------------------

The square-free factorization of a univariate polynomial is the product of all
factors (not necessarily irreducible) of degree 1, 2 etc.::

    >>> f = 2*x**2 + 5*x**3 + 4*x**4 + x**5

    >>> sqf_list(f)
    (1, [(x + 2, 1), (x, 2), (x + 1, 2)])

    >>> sqf(f)
     2        2
    x *(x + 1) *(x + 2)

Factorization
-------------

This function provides factorization of univariate and multivariate polynomials
with rational coefficients::

    >>> factor(x**4/2 + 5*x**3/12 - x**2/3)
     2
    x *(2*x - 1)*(3*x + 4)
    ----------------------
              12

    >>> factor(x**2 + 4*x*y + 4*y**2)
             2
    (x + 2*y)

Groebner bases
--------------

Buchberger's algorithm is implemented, supporting various monomial orders::

    >>> groebner([x**2 + 1, y**4*x + x**3], x, y, order='lex')
                 /[ 2       4    ]                            \
    GroebnerBasis\[x  + 1, y  - 1], x, y, domain=ZZ, order=lex/


    >>> groebner([x**2 + 1, y**4*x + x**3, x*y*z**3], x, y, z, order='grevlex')
                 /[ 4       3   2    ]                                   \
    GroebnerBasis\[y  - 1, z , x  + 1], x, y, z, domain=ZZ, order=grevlex/

Solving Equations
-----------------

We have (incomplete) methods to find the complex or even symbolic roots of
polynomials and to solve some systems of polynomial equations::

    >>> from sympy import roots, solve_poly_system

    >>> solve(x**3 + 2*x + 3, x)
               ____          ____
         1   \/ 11 *I  1   \/ 11 *I
    [-1, - - --------, - + --------]
         2      2      2      2

    >>> p = Symbol('p')
    >>> q = Symbol('q')

    >>> solve(x**2 + p*x + q, x)
              __________           __________
             /  2                 /  2
       p   \/  p  - 4*q     p   \/  p  - 4*q
    [- - - -------------, - - + -------------]
       2         2          2         2

    >>> solve_poly_system([y - x, x - 5], x, y)
    [(5, 5)]

    >>> solve_poly_system([y**2 - x**3 + 1, y*x], x, y)
                                       ___                 ___
                                 1   \/ 3 *I         1   \/ 3 *I
    [(0, -I), (0, I), (1, 0), (- - - -------, 0), (- - + -------, 0)]
                                 2      2            2      2
