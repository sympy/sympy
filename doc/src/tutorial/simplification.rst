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

    >>> expr = x*y + x - 3 + 2*x**2 - z*x**2 + x**3
    >>> expr
     3    2        2
    x  - x ⋅z + 2⋅x  + x⋅y + x - 3
    >>> collected_expr = collect(expr, x)
    >>> collected_expr
     3    2
    x  + x ⋅(-z + 2) + x⋅(y + 1) - 3

``collect`` is particularly useful in conjunction with the ``.coeff`` method,
which will be discussed in more detail later.

    >>> collected_expr.coeff(x, 2)
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

    >>> expr = 1/x + (3*x/2 - 2)/(x - 4)
    >>> expr
    3⋅x
    ─── - 2
     2        1
    ─────── + ─
     x - 4    x
    >>> cancel(expr)
       2
    3⋅x  - 2⋅x - 8
    ──────────────
         2
      2⋅x  - 8⋅x

    >>> expr = (x*y**2 - 2*x*y*z + x*z**2 + y**2 - 2*y*z + z**2)/(x**2 - 1)
    >>> expr
       2                2    2            2
    x⋅y  - 2⋅x⋅y⋅z + x⋅z  + y  - 2⋅y⋅z + z
    ───────────────────────────────────────
                      2
                     x  - 1
    >>> cancel(expr)
     2            2
    y  - 2⋅y⋅z + z
    ───────────────
         x - 1

Note that since ``factor`` will completely factorize both the numerator and
the denominator of an expression, it can also be used to do the same thing:

    >>> factor(expr)
           2
    (y - z)
    ────────
     x - 1

However, if you are only interested in making sure that the expression is in
canceled form, ``cancel`` is more efficient than ``factor``.

``apart``
---------

``apart`` performs a partial fraction decomposition on a rational function.

    >>> expr = (4*x**3 + 21*x**2 + 10*x + 12)/(x**4 + 5*x**3 + 5*x**2 + 4*x)
    >>> expr
       3       2
    4⋅x  + 21⋅x  + 10⋅x + 12
    ────────────────────────
      4      3      2
     x  + 5⋅x  + 5⋅x  + 4⋅x
    >>> apart(expr)
     2⋅x - 1       1     3
    ────────── - ───── + ─
     2           x + 4   x
    x  + x + 1

Trigonometric Simplification
============================

.. note::

   SymPy follows Python's naming conventions for inverse trigonometric
   functions, which is to append an ``a`` to the front of the function's
   name.  For example, the inverse cosine, or arc cosine, is called ``acos``.

   >>> acos(x)
   acos(x)
   >>> cos(acos(x))
   x
   >>> asin(1)
   π
   ─
   2

.. TODO: Can we actually do anything with inverse trig functions,
   simplification wise?

``trigsimp``
------------

To simplify expressions using trigonometric identities, use ``trigsimp``.

    >>> trigsimp(sin(x)**2 + cos(x)**2)
    1
    >>> trigsimp(sin(x)**4 - 2*cos(x)**2*sin(x)**2 + cos(x)**4)
    cos(4⋅x)   1
    ──────── + ─
       2       2
    >>> trigsimp(sin(x)*tan(x)/sec(x))
       2
    sin (x)

``trigsimp`` also works with hyperbolic trig functions.

    >>> trigsimp(cosh(x)**2 + sinh(x)**2)
    cosh(2⋅x)
    >>> trigsimp(sinh(x)/tanh(x))
    cosh(x)

Much like ``simplify``, ``trigsimp`` applies various trigonometric identities to
the input expression, and then uses a heuristic to return the "best" one.

``expand_trig``
---------------

To expand trigonometric functions, that is, apply the sum or double angle
identities, use ``expand_trig``.

    >>> expand_trig(sin(x + y))
    sin(x)⋅cos(y) + sin(y)⋅cos(x)
    >>> expand_trig(tan(2*x))
       2⋅tan(x)
    ─────────────
         2
    - tan (x) + 1

Because ``expand_trig`` tends to make trigonometric expressions larger, and
``trigsimp`` tends to make them smaller, these identities can be applied in
reverse using ``trigsimp``

.. TODO: It would be much better to teach individual trig rewriting functions
   here, but they don't exist yet.  See
   https://code.google.com/p/sympy/issues/detail?id=357.

Powers
======

There are three kinds of identities satisfied by exponents

1. `x^ax^b = x^{a + b}`
2. `x^ay^a = (xy)^a`
3. `(x^a)^b = x^{ab}`

Identity 1 is always true.

Identity 2 is not always true.  For example, if `x = y = -1` and `a = \frac{1}{2}`,
then `x^ay^a = \sqrt{-1}\sqrt{-1} = i\cdot i = -1`, whereas `(xy)^a =
\sqrt{-1\cdot-1} = \sqrt{1} = 1`.  However, identity 2 is true at least if `x`
and `y` are nonnegative and `a` is real (it may also be true under other
conditions as well).  A common consequence of the failure of identity 2 is
that `\sqrt{x}\sqrt{y} \neq \sqrt{xy}`.

Identity 3 is not always true.  For example, if `x = -1`, `a = 2`, and `b =
\frac{1}{2}`, then `(x^a)^b =  {\left ((-1)^2\right )}^{1/2} = \sqrt{1} = 1` and `x^{ab} =
(-1)^{2\cdot1/2} = (-1)^1 = -1`.  However, identity 3 is true at least if `x`
is positive or `b` is an integer (again, it may also hold in other cases as
well).  Two common consequences of the failure of identity 3 are that
`\sqrt{x^2}\neq x` and that `\sqrt{\frac{1}{x}} \neq \frac{1}{\sqrt{x}}`.

This is important to remember, because by default, SymPy will not perform
simplifications if they are not true in general.

In order to make SymPy perform simplifications involving identities that are
only true under certain assumptions, we need to put assumptions on our
Symbols.  We will undertake a full discussion of the assumptions system later,
but for now, all we need to know are the following

- By default, SymPy Symbols are assumed to be complex.  That is,
  simplifications will not be applied to an expression with a given Symbol
  unless it holds for all complex numbers.

- Symbols can be given different assumptions by passing the assumption to
  ``symbols``.  For the rest of this section, we will be assuming that ``x``,
  ``y`` are positive, and that ``a`` and ``b`` are real.  We will leave ``z``,
  ``t``, and ``c`` as arbitrary complex Symbols to demonstrate what happens in
  that case.

    >>> x, y = symbols('x y', positive=True)
    >>> a, b = symbols('a b', real=True)
    >>> z, t, c = symbols('z t c')

  .. TODO: Rewrite this using the new assumptions

- In SymPy, ``sqrt(x)`` is just a shortcut to ``x**Rational(1, 2)``.  They are
  exactly the same object.

    >>> sqrt(x) == x**Rational(1, 2)
    True

``powsimp``
-----------

``powsimp`` applies identities 1 and 2 from above, from left to right.


   >>> powsimp(x**a*x**b)
     a + b
    x
   >>> powsimp(x**a*y**a)
        a
   (x⋅y)

Notice that ``powsimp`` refuses to do the simplification if it is not valid.

    >>> powsimp(t**c*z**c)
     c  c
    t ⋅z

If you know that you want to apply this simplification, but you don't want to
mess with assumptions, you can pass the ``force=True`` flag.  This will force
the simplification to take place, regardless of assumptions.

    >>> powsimp(t**c*z**c, force=True)
         c
    (t⋅z)

Note that in some instances, in particular, when the exponents are integers or
rational numbers, and identity 2 holds, it will be applied automatically

   >>> (z*t)**2
     2  2
    t ⋅z
   >>> sqrt(x*y)
      ___   ___
    ╲╱ x ⋅╲╱ y

This means that it will be impossible to undo this identity with ``powsimp``,
because even if ``powsimp`` were to put the bases together, they would be
automatically split apart again.

   >>> powsimp(z**2*t**2)
     2  2
    t ⋅z
   >>> powsimp(sqrt(x)*sqrt(y))
      ___   ___
    ╲╱ x ⋅╲╱ y

``expand_power_exp``/``expand_power_base``
------------------------------------------

``expand_power_exp`` and ``expand_power_base`` apply identities 1 and 2 from
right to left, respectively.

    >>> expand_power_exp(x**(a + b))
     a  b
    x ⋅x

    >>> expand_power_base((x*y)**a)
     a  a
    x ⋅y

As with ``powsimp``, identity 2 is not applied if it is not valid.

    >>> expand_power_base((z*t)**c)
         c
    (t⋅z)

And as with ``powsimp``, you can force the expansion to happen without
fiddling with assumptions by using ``force=True``.

   >>> expand_power_base((z*t)**c, force=True)
     c  c
    t ⋅z

As with identity 2, identity 1 is applied automatically if the power is a
number, and hence cannot be undone with ``expand_power_exp``.

   >>> x**2*x**3
     5
    x
   >>> expand_power_exp(x**5)
     5
    x

``powdenest``
-------------

``powdenest`` applies identity 3, for left to right.

    >>> powdenest((x**a)**b)
     a⋅b
    x

As before, the identity is not applied if it is not true under the given
assumptions.

    >>> powdenest((z**a)**b)
        b
    ⎛ a⎞
    ⎝z ⎠

And as before, this can be manually overridden with ``force=True``.

    >>> powdenest((z**a)**b, force=True)
     a⋅b
    z

Exponentials and logarithms
===========================

.. note::

   In SymPy, as in Python and most programming languages, ``log`` is the
   natural logarithm, also known as ``ln``.  SymPy automatically provides an
   alias ``ln = log`` in case you forget this.

    >>> ln(x)
    log(x)

Logarithms have similar issues as powers.  There are two main identities

1. `\log{(xy)} = \log{(x)} + \log{(y)}`
2. `\log{(x^n)} = n\log{(x)}`

Neither identity is true for arbitrary complex `x` and `y`, due to the branch
cut in the complex plane for the complex logarithm.  However, sufficient
conditions for the identities to hold are if `x` and `y` are positive, and if
`n` is real.

    >>> x, y = symbols('x y', positive=True)
    >>> n = symbols('n', real=True)

As before, ``z`` and ``t`` will be Symbols with no additional assumptions.

Note that the identity `\log{\left (\frac{x}{y}\right )} = \log(x) - \log(y)`
is a special case of identities 1 and 2 by `\log{\left (\frac{x}{y}\right )}
=` `\log{\left (x\cdot\frac{1}{y}\right )} =` `\log(x) + \log{\left(
y^{-1}\right )} =` `\log(x) - \log(y)`, and thus it also holds if `x` and `y`
are positive.

We also see that `\log{\left( e^x \right)} = x` comes from `\log{\left ( e^x
\right)} = x\log(e) = x`, and thus holds when `x` is real (however, it can be
verified that it does not hold in general for complex `x`, for example,
`\log{\left (e^{x + 2\pi i}\right)} = \log{\left (e^x\right )} = x \neq x +
2\pi i`).

``expand_log``
--------------

To apply identities 1 and 2 from left to right, use ``expand_log``.  As
always, the identities will not be applied unless they are valid.

    >>> expand_log(log(x*y))
    log(x) + log(y)
    >>> expand_log(log(x/y))
    log(x) - log(y)
    >>> expand_log(log(x**2))
    2⋅log(x)
    >>> expand_log(log(x**n))
    n⋅log(x)
    >>> expand_log(log(z*t))
    log(t⋅z)

As with ``powsimp`` and ``powdenest``, ``expand_log`` has a ``force`` option
that can be used to ignore assumptions.

    >>> expand_log(log(z**2))
       ⎛ 2⎞
    log⎝z ⎠
    >>> expand_log(log(z**2), force=True)
    2⋅log(z)

``logcombine``
--------------

To apply identities 1 and 2 from right to left, use ``logcombine``.

    >>> logcombine(log(x) + log(y))
    log(x⋅y)
    >>> logcombine(n*log(x))
       ⎛ n⎞
    log⎝x ⎠
    >>> logcombine(n*log(z))
    n⋅log(z)

``logcombine`` also has a ``force`` option that can be used to ignore
assumptions.

    >>> logcombine(n*log(z), force=True)
       ⎛ n⎞
    log⎝z ⎠

Special Functions
=================

SymPy implements dozens of special functions, whose use ranges from
combinatorics to mathematical physics.

An extensive list of the special functions included with SymPy and their
documentation is at the :ref:`functions` page.
