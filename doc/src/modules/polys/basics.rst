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
    >>> init_printing(use_unicode=False)

Basic concepts
==============

Polynomials
-----------

Given a family `(x_i)` of symbols, or other suitable objects, including numbers,
expressions derived from them by repeated addition, subtraction and
multiplication are called *polynomial expressions in the generators* `x_i`.

By the distributive law it is possible to perform multiplications before
additions and subtractions. The products of generators thus obtained are called
*monomials*. They are usually written in the form `x_1^{\nu_1}x_2^{\nu_2}\cdots
x_n^{\nu_n}` where the exponents `\nu_i` are nonnegative integers. It is often
convenient to write this briefly as `x^\nu` where `x = (x_1, x_2, \ldots, x_n)`
denotes the family of generators and `\nu = (\nu_1, \nu_2, \ldots, \nu_n)` is
the family of exponents.

When all monomials having the same exponents are combined, the polynomial
expression becomes a sum of products `c_\nu x^\nu`, called the *terms* of the
polynomial, where the *coefficients* `c_\nu` are integers. If some of the `x_i`
are manifest numbers, they are incorporated in the coefficients and not regarded
as generators. Such coefficients are typically rational, real or complex
numbers. Some symbolic numbers, e.g., ``pi``, can be either coefficients or
generators.

A polynomial expression that is a sum of terms with different monomials is
uniquely determined by its family of coefficients `(c_\nu)`. Such an expression
is customarily called a *polynomial*, though, more properly, that name does
stand for the coefficient family once the generators are given. SymPy implements
polynomials by default as dictionaries with monomials as keys and coefficients
as values. Another implementation consists of nested lists of coefficients.

.. _polys-ring:

The set of all polynomials with integer coefficients in the generators `x_i` is
a *ring*, i.e., the sums, differences and products of its elements are again
polynomials in the same generators. This ring is denoted `\mathbb{Z}[x_1, x_2,
\ldots, x_n]`, or `\mathbb{Z}[(x_i)]`, and called the *ring of polynomials in
the* `x_i` *with integer coefficients*.

More generally, the coefficients of a polynomial can be elements of any
commutative ring `A`, and the corresponding polynomial ring is then denoted
`A[x_1, x_2, \dots, x_n]`. The ring `A` can also be a polynomial ring. In SymPy,
the coefficient ring is called the ``domain`` of the polynomial ring, and it can
be given as a keyword parameter. By default, it is determined by the
coefficients of the polynomial arguments.

Polynomial expressions can be transformed into polynomials by the method
:obj:`sympy.core.expr.Expr.as_poly`::

    >>> e = (x + y)*(y - 2*z)
    >>> e.as_poly()
    Poly(x*y - 2*x*z + y**2 - 2*y*z, x, y, z, domain='ZZ')

If a polynomial expression contains numbers that are not integers, they are
regarded as coefficients and the coefficient ring is extended accordingly. In
particular, division by integers leads to rational coefficients::

    >>> e = (3*x/2 + y)*(z - 1)
    >>> e.as_poly()
    Poly(3/2*x*z - 3/2*x + y*z - y, x, y, z, domain='QQ')

Symbolic numbers are considered generators unless they are explicitly excluded,
in which case they are adjoined to the coefficient ring::

    >>> e = (x + 2*pi)*y
    >>> e.as_poly()
    Poly(x*y + 2*y*pi, x, y, pi, domain='ZZ')
    >>> e.as_poly(x, y)
    Poly(x*y + 2*pi*y, x, y, domain='ZZ[pi]')

Alternatively, the coefficient domain can be specified by means of a keyword
argument::

    >>> e = (x + 2*pi)*y
    >>> e.as_poly(domain=ZZ[pi])
    Poly(x*y + 2*pi*y, x, y, domain='ZZ[pi]')

Note that the ring `\mathbb{Z}[\pi][x, y]` of polynomials in `x` and `y` with
coefficients in `\mathbb{Z}[\pi]` is mathematically equivalent to
`\mathbb{Z}[\pi, x, y]`, only their implementations differ.

If an expression contains functions of the generators, other than their positive
integer powers, these are interpreted as new generators::

    >>> e = x*sin(y) - y
    >>> e.as_poly()
    Poly(x*(sin(y)) - y, x, y, sin(y), domain='ZZ')

Since `y` and `\sin(y)` are algebraically independent they can both appear as
generators in a polynomial. However, *polynomial expressions must not contain
negative powers of generators*::

    >>> e = x - 1/x
    >>> e.as_poly()
    Poly(x - (1/x), x, 1/x, domain='ZZ')

It is important to realize that the generators `x` and `1/x = x^{-1}` are
treated as algebraically independent variables. In particular, their product is
not equal to 1. Hence *generators in denominators should be avoided even if they
raise no error in the current implementation*. This behavior is undesirable and
may change in the future. Similar problems emerge with rational powers of
generators. So, for example, `x` and `\sqrt x = x^{1/2}` are not recognized as
algebraically dependent.

If there are algebraic numbers in an expression, it is possible to adjoin them
to the coefficient ring by setting the keyword ``extension``::

    >>> e = x + sqrt(2)
    >>> e.as_poly()
    Poly(x + (sqrt(2)), x, sqrt(2), domain='ZZ')
    >>> e.as_poly(extension=True)
    Poly(x + sqrt(2), x, domain='QQ<sqrt(2)>')

With the default setting ``extension=False``, both `x` and `\sqrt 2` are
incorrectly considered algebraically independent variables. With coefficients in
the extension field `\mathbb{Q}(\sqrt 2)` the square root is treated properly as
an algebraic number. Setting ``extension=True`` whenever algebraic numbers are
involved is definitely recommended even though it is not forced in the current
implementation.

Divisibility
------------

The fourth rational operation, division, or inverted multiplication, is not
generally possible in rings. If `a` and `b` are two elements of a ring `A`, then
there may exist a third element `q` in `A` such that `a = bq`. In fact, there
may exist several such elements.

If also `a = bq'` for some `q'` in `A`, then `b(q - q') = 0`. Hence either `b`
or `q - q'` is zero, or they are both *zero divisors*, nonzero elements whose
product is zero.

Integral domains
````````````````
Commutative rings with no zero divisors are called *integral domains*. Most of
the commonly encountered rings, the ring of integers, fields, and polynomial
rings over integral domains are integral domains.

Assume now that `A` is an integral domain, and consider the set `P` of its
nonzero elements, which is closed under multiplication. If `a` and `b` are in
`P`, and there exists an element `q` in `P` such that `a = bq`, then `q` is
unique and called the *quotient*, `a/b`, of `a`  by `b`. Moreover, it is said
that

- `a` is *divisible* by `b`,

- `b` is a *divisor* of `a`,

- `a` is a *multiple* of `b`,

- `b` is a *factor* of `a`.

An element `a` of `P` is a divisor of `1` if and only if it is *invertible* in
`A`, with the inverse `a^{-1} = 1/a`. Such elements are called *units*. The
units of the ring of integers are `1` and `-1`. The invertible elements in a
polynomial ring over a field are the nonzero constant polynomials.

If two elements of `P`, `a` and `b`, are divisible by each other, then the
quotient `a/b` is invertible with inverse `b/a`, or equivalently, `b = ua` where
`u` is a unit. Such elements are said to be *associated* with, or *associates*
of, each other. The associates of an integer `n` are `n` and `-n`. In a
polynomial ring over a field the associates of a polynomial are its constant
multiples.

Each element of `P` is divisible by its associates and the units. An element is
*irreducible* if it has no other divisors and is not a unit. The irreducible
elements in the ring of integers are the prime numbers `p` and their opposites
`-p`. In a field, every nonzero element is invertible and there are no
irreducible elements.

Factorial domains
`````````````````
In the ring of integers, each nonzero element can be represented as a product of
irreducible elements and optionally a unit `\pm 1`. Moreover, any two such
products have the same number of irreducible factors which are associated with
each other in a suitable order. Integral domains having this property are called
*factorial*, or *unique factorization domains*. In addition to the ring of
integers, all polynomial rings over a field are factorial, and so are more
generally polynomial rings over any factorial domain. Fields are trivially
factorial since there are only units. The irreducible elements of a factorial
domain are usually called *primes*.

A family of integers has only a finite number of common divisors and the
greatest of them is divisible by all of them. More generally, given a family of
nonzero elements `(a_i)` in an integral domain, a common divisor `d` of the
elements is called a *greatest common divisor*, abbreviated *gcd*, of the family
if it is a multiple of all common divisors. A greatest common divisor, if it
exists, is not unique in general; all of its associates have the same property.
It is denoted by `d = \gcd(a_1,\ldots,a_n)` if there is no danger of confusion.
A *least common multiple*, or *lcm*, of a family `(a_i)` is defined analogously
as a common multiple `m` that divides all common multiples. It is denoted by `m
= \operatorname{lcm}(a_1,\dots,a_n)`.

In a factorial domain, greatest common divisors always exists. They can be
found, at least in principle, by factoring each element of a family into a
product of prime powers and an optional unit, and, for each prime, taking the
least power that appears in the factorizations. The product of these prime
powers is then a greatest common divisor. A least common multiple can be
obtained from the same factorizations as the product of the greatest powers for
each prime.

Euclidean domains
`````````````````
A practical algorithm for computing a greatest common divisor can be implemented
in *Euclidean domains*. They are integral domains that can be endowed with a
function `w` assigning a nonnegative integer to each nonzero element of the
domain and having the following property:

    if `a` and `b` are nonzero, there are `q` and `r` that satisfy the *division
    identity*

        `a = qb + r`

    such that either `r = 0` or `w(r) < w(b)`.


The ring of integers and all univariate polynomial rings over fields are
Euclidean domains with `w(a) = |a|` resp. `w(a) = \deg(a)`.

The division identity for integers is implemented in Python as the built-in
function ``divmod`` that can also be applied to SymPy Integers::

    >>> divmod(Integer(53), Integer(7))
    (7, 4)

For polynomials the division identity is given in SymPy by the function
:func:`~.div`::

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

The division identity can be used to determine the divisibility of elements in a
Euclidean domain. If `r = 0` in the division identity, then `a` is divisible by
`b`. Conversely, if `a = cb` for some element `c`, then `(c - q)b = r`. It
follows that `c = q` and `r = 0` if `w` has the additional property:

    if `a` and `b` are nonzero, then `w(ab) \ge w(b)`.

This is satisfied by the functions given above. (And it is always possible to
redefine `w(a)` by taking the minimum of the values `w(xa)` for `x \ne 0`.)

The principal application of the division identity is the efficient computation
of a greatest common divisor by means of the `Euclidean algorithm
<https://en.wikipedia.org/wiki/Euclidean_algorithm>`_. It applies to two
elements of a Euclidean domain. A gcd of several elements can be obtained by
iteration.

The function for computing the greatest common divisor of integers in SymPy is
currently :func:`~.igcd`::

    >>> igcd(2, 4)
    2
    >>> igcd(5, 10, 15)
    5

For univariate polynomials over a field the function has its common name
:func:`~.gcd`, and the returned polynomial is monic::

    >>> f = 4*x**2 - 1
    >>> g = 8*x**3 + 1
    >>> gcd(f, g, domain=QQ)
    x + 1/2


Divisibility of polynomials
```````````````````````````
The ring `A = \mathbb{Z}[x]` of univariate polynomials over the ring of integers
is not Euclidean but it is still factorial. To see this, consider the
divisibility in `A`.

Let `f` and `g` be two nonzero polynomials in `A`. If `f` is divisible by `g` in
`A`, then it is also divisible in the ring `B = \mathbb{Q}[x]` of polynomials
with rational coefficients. Since `B` is Euclidean, this can be determined by
means of the division identity.

Assume, conversely, that `f = gh` for some polynomial `h` in `B`. Then `f` is
divisible by `g` in `A` if and only if the coefficients of `h` are integers. To
find out when this is true it is necessary to consider the divisibility of the
coefficients.

For a polynomial `f` in `A`, let `c` be the greatest common divisor of its
coefficients. Then `f` is divisible by the constant polynomial `c` in `A`, and
the quotient `f/c= p` is a polynomial whose coefficients are integers that have
no common divisor apart from the units. Such polynomials are called *primitive*.
A polynomial with rational coefficients can also be written as `f = cp`, where
`c` is a rational number and `p` is a primitive polynomial. The constant `c` is
called the *content* of `f`, and `p` is its *primitive part*. These components
can be found by the method :obj:`sympy.core.expr.Expr.as_content_primitive`::

    >>> f = 6*x**2 - 3*x + 9
    >>> c, p = f.as_content_primitive()
    >>> c, p
           2
    (3, 2*x  - x + 3)
    >>> f = x**2/3 - x/2 + 1
    >>> c, p = f.as_content_primitive()
    >>> c, p
             2
    (1/6, 2*x  - 3*x + 6)

Let `f`, `f'` be polynomials with contents `c`, `c'` and primitive parts `p`,
`p'`. Then `ff' = (cc')(pp')` where the product `pp'` is primitive by `Gauss's
lemma <https://en.wikipedia.org/wiki/Gauss%27s_lemma_(polynomial)>`_. It follows
that

    the content of a product of polynomials is the product of their contents and
    the primitive part of the product is the product of the primitive parts.

Returning to the divisibility in the ring `\mathbb{Z}[x]`, assume that `f` and
`g` are two polynomials with integer coefficients such that the division
identity in `\mathbb{Q}[x]` yields the equality `f = gh` for some polynomial `h`
with rational coefficients. Then the content of `f` is equal to the content of
`g` multiplied by the content of `h`. As `h` has integer coefficients if and
only if its content is an integer, we get the following criterion:

    `f` is divisible by `g` in the ring `\mathbb{Z}[x]` if and only if

    i. `f` is divisible by `g` in `\mathbb{Q}[x]`, and
    ii. the content of `f` is divisible by the content of `g` in `\mathbb{Z}`.

If `f = cp` is irreducible in `\mathbb{Z}[x]`, then either `c` or `p` must be a
unit. If `p` is not a unit, it must be irreducible also in `\mathbb{Q}[x]`. For
if it is a product of two polynomials, it is also the product of their primitive
parts, and one of them must be a unit. Hence there are two kinds of irreducible
elements in `\mathbb{Z}[x]`:

i. prime numbers of `\mathbb{Z}`, and
ii. primitive polynomials that are irreducible in `\mathbb{Q}[x]`.

It follows that each polynomial in `\mathbb{Z}[x]` is a product of irreducible
elements. It suffices to factor its content and primitive part separately. These
products are essentially unique; hence `\mathbb{Z}[x]` is also factorial.

Another important consequence is that a greatest common divisor of two
polynomials in `\mathbb{Z}[x]` can be found efficiently by applying the
Euclidean algorithm separately to their contents and primitive parts in the
Euclidean domains `\mathbb{Z}` and `\mathbb{Q}[x]`. This is also implemented in
SymPy::

    >>> f = 4*x**2 - 1
    >>> g = 8*x**3 + 1
    >>> gcd(f, g)
    2*x + 1
    >>> gcd(6*f, 3*g)
    6*x + 3

Basic functionality
===================

These functions provide different algorithms dealing with polynomials in the
form of SymPy expression, like symbols, sums etc.

Division
--------

The function :func:`~.div` provides division of polynomials with remainder. That
is, for polynomials ``f`` and ``g``, it computes ``q`` and ``r``, such that `f =
g \cdot q + r` and `\deg(r) < \deg(q)`. For polynomials in one variables with
coefficients in a field, say, the rational numbers, ``q`` and ``r`` are uniquely
defined this way::

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
remainder doesn't need to be smaller than that of ``f``. Since 2 doesn't divide
5, `2 x` doesn't divide `5 x^2`, even if the degree is smaller. But::

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

But if the polynomials have rational coefficients, then the returned polynomial
is monic::

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
                       2
    (1, [(x + 2, 1), (x  + x, 2)])

    >>> sqf(f)
                    2
            / 2    \
    (x + 2)*\x  + x/

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
