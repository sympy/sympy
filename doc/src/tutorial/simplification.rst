================
 Simplification
================

To make this document easier to read, we are going to enable pretty printing.

    >>> from sympy import *
    >>> x, y, z = symbols('x y z')
    >>> init_printing(use_unicode=True)

``simplify``
============

Now that we know how to create Symbols, let's jump in and do some interesting
mathematics.  One of the most useful features of a symbolic manipulation
system is the ability to simplify mathematical expressions.  SymPy has dozens
of functions to perform various kinds of simplification.  There is also one
general function called ``simplify`` that attempts to apply all of these
functions in an intelligent way to arrive at the simplest form of an
expression.  Here are some examples

    >>> simplify(sin(x)**2 + cos(x)**2)
    1
    >>> simplify((x**3 + x**2 - x - 1)/(x**2 + 2*x + 1))
    x - 1
    >>> simplify(gamma(x)/gamma(x - 2))
    (x - 2)⋅(x - 1)

Here, ``gamma(x)`` is `\Gamma(x)`, the gamma function.  We see that
``simplify`` is capable of handling a large class of expressions.

But ``simplify`` has a pitfall.  It just applies all the major simplification
operations in SymPy, and uses heuristics to determine the simplest result. But
"simplest" is not a well-defined term.  For example, say we wanted to
"simplify" `x^2 + 2x + 1` into `(x + 1)^2`:

    >>> simplify(x**2 + 2*x + 1)
     2
    x  + 2⋅x + 1

We did not get what we want.  There is a function to perform this
simplification, called ``factor``, which will be discussed below.

Another pitfall to ``simplify``, which is that it can be unnecessarily slow,
since it tries many kinds of simplifications before picking the best one.  If
you already know exactly what kind of simplification you are after, it is
better to apply the specific simplification function(s) that apply those
simplifications.

Applying specific simplification functions instead of ``simplify`` also has
the advantage that specific functions have certain guarantees about the form
of their output.  These will be discussed with each function below.  For
example, ``factor``, when called on a polynomial with rational coefficients,
is guaranteed to factor the polynomial into irreducible factors.  ``simplify``
has no guarantees.  It is entirely heuristical, and, as we saw above, it may
even miss a possible type of simplification that SymPy is capable of doing.

``simplify`` is best when used interactively, when you just want to whittle
down an expression to a simpler form.  You may then choose to apply specific
functions once you see what ``simplify`` returns, to get a more precise
result.  It is also useful when you have no idea what form an expression will
take, and you need a catchall function to simplify it.

Polynomial/Rational Function Simplification
===========================================

``expand``
----------

``expand`` is one of the most common simplification functions in SymPy.
Although it has a lot of scopes, for now, we will consider its function in
expanding polynomial expressions. For example:

    >>> expand((x + 1)**2)
     2
    x  + 2⋅x + 1
    >>> expand((x + 2)*(x - 3))
     2
    x  - x - 6

Given a polynomial, ``expand`` will put it into a canonical form of a sum of
monomials.

``expand`` may not sound like a simplification function.  After all, by its
very name, it makes expressions bigger, not smaller.  Usually, this is the
case, but often, an expression will become smaller upon calling ``expand`` on
it due to cancellation

    >>> expand((x + 1)*(x - 2) - (x - 1)*x)
    -2

``factor``
----------

``factor`` takes a polynomial and factors it into irreducible factors over the
rational numbers.  For example:

    >>> factor(x**3 - x**2 + x - 1)
            ⎛ 2    ⎞
    (x - 1)⋅⎝x  + 1⎠
    >>> factor(x**2*z + 4*x*y*z + 4*y**2*z)
               2
    z⋅(x + 2⋅y)

For polynomials, ``factor`` is the opposite of ``expand``.  ``factor`` uses a
complete multivariate factorization algorithm over the rational numbers, which
means that each of the factors returned by ``factor`` is guaranteed to be
irreducible.

If you are interested in the factors themselves, ``factor_list`` returns a
more structured output.

    >>> factor_list(x**2*z + 4*x*y*z + 4*y**2*z)
    (1, [(z, 1), (x + 2⋅y, 2)])

Note that the input to ``factor`` and ``expand`` need not be polynomials in
the strict sense.  They will intelligently factor or expand any kind of
expression (though note that the factors may not be irreducible if the input
is no longer a polynomial over the rationals).

    >>> expand((cos(x) + sin(x))**2)
       2                           2
    sin (x) + 2⋅sin(x)⋅cos(x) + cos (x)
    >>> factor(cos(x)**2 + 2*cos(x)*sin(x) + sin(x)**2)
                     2
    (sin(x) + cos(x))

``collect``
-----------

``collect`` collects common powers of a term in an expression.  For example

    >>> a = x*y + x - 3 + 2*x**2 - z*x**2 + x**3
    >>> a
     3    2        2
    x  - x ⋅z + 2⋅x  + x⋅y + x - 3
    >>> collected_a = collect(a, x)
    >>> collected_a
     3    2
    x  + x ⋅(-z + 2) + x⋅(y + 1) - 3

``collect`` is particularly useful in conjunction with the ``.coeff`` method,
which will be discussed in more detail later.

    >>> collected_a.coeff(x, 2)
    -z + 2

``cancel``
----------

``cancel`` will take any rational function and put it into the standard
canonical form, `\frac{p}{q}`, where `p` and `q` are expanded polynomials with
no common factors, and the leading coefficients of `p` and `q` do not have
denominators (i.e., are integers).

    >>> cancel((x**2 + 2*x + 1)/(x**2 + x))
    x + 1
    ─────
      x

    >>> a = 1/x + (3*x/2 - 2)/(x - 4)
    >>> a
    3⋅x
    ─── - 2
     2        1
    ─────── + ─
     x - 4    x
    >>> cancel(a)
       2
    3⋅x  - 2⋅x - 8
    ──────────────
         2
      2⋅x  - 8⋅x

    >>> a = (x*y**2 - 2*x*y*z + x*z**2 + y**2 - 2*y*z + z**2)/(x**2 - 1)
    >>> a
       2                2    2            2
    x⋅y  - 2⋅x⋅y⋅z + x⋅z  + y  - 2⋅y⋅z + z
    ───────────────────────────────────────
                      2
                     x  - 1
    >>> cancel(a)
     2            2
    y  - 2⋅y⋅z + z
    ───────────────
         x - 1

Note that since ``factor`` will completely factorize both the numerator and
the denominator of an expression, it can also be used to do the same thing:

    >>> factor(a)
           2
    (y - z)
    ────────
     x - 1

However, if you are only interested in making sure that the expression is in
canceled form, ``cancel`` is more efficient than ``factor``.

``apart``
---------

``apart`` performs a partial fraction decomposition on a rational function.

    >>> a = (4*x**3 + 21*x**2 + 10*x + 12)/(x**4 + 5*x**3 + 5*x**2 + 4*x)
    >>> a
       3       2
    4⋅x  + 21⋅x  + 10⋅x + 12
    ────────────────────────
      4      3      2
     x  + 5⋅x  + 5⋅x  + 4⋅x
    >>> apart(a)
     2⋅x - 1       1     3
    ────────── - ───── + ─
     2           x + 4   x
    x  + x + 1
