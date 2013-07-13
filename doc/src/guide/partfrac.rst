
.. include:: ../definitions.def

Partial fraction decomposition
==============================

The partial fraction decomposition of a univariate rational function:

.. math::

    f(x) = \frac{p(x)}{q(x)}

where `p` and `q` are co-prime and `\deg(p) < \deg(q)`, is an expression
of the form:

.. math::

    \sum_{i=1}^k \sum_{j=1}^{n_i} \frac{a_{ij}(x)}{q_i^j(x)}

where `q_i` for `i=1 \ldots k` are factors (e.g. over rationals or Gaussian
rationals) of `q`:

.. math::

    q(x) = \prod_{i=1}^k q_i^{n_i}

If `p` and `q` aren't co-prime, we can use :func:`cancel` to remove common
factors and if `\deg(p) >= \deg(q)`, then :func:`div` can be used to extract
the polynomial part of `f` and reduce the degree of `p`.

Suppose we would like to compute partial fraction decomposition of::

    >>> from sympy import *
    >>> init_printing(use_unicode=True, no_global=True)

    >>> var('x')
    x

    >>> f = 1/(x**2*(x**2 + 1))
    >>> f
         1
    ───────────
     2 ⎛ 2    ⎞
    x ⋅⎝x  + 1⎠

This can be achieved with SymPy's built-in function :func:`apart`::

    >>> apart(f)
        1      1
    - ────── + ──
       2        2
      x  + 1   x

We can use :func:`together` to verify this result::

    >>> together(_)
         1
    ───────────
     2 ⎛ 2    ⎞
    x ⋅⎝x  + 1⎠

Now we would like to compute this decomposition step-by-step. The rational
function `f` is already in factored form and has two factors `x^2` and
`x^2 + 1`. If `f` was in expanded from, we could use :func:`factor` to
obtain the desired factorization::

    >>> numer(f)/expand(denom(f))
       1
    ───────
     4    2
    x  + x

    >>> factor(_)
         1
    ───────────
     2 ⎛ 2    ⎞
    x ⋅⎝x  + 1⎠

Based on the definition, the partial fraction expansion of `f` will be of the
following form:

.. math::

    \frac{A}{x} + \frac{B}{x^2} + \frac{C x + D}{x^2 + 1}

Let's do this with SymPy. We will use undetermined coefficients method to
solve this problem. Let's start by defining some symbols::

    >>> var('A:D')
    (A, B, C, D)

We use here the lexicographic syntax of :func:`var`. Next we can define three
rational functions::

    >>> p1 = A/x
    >>> p2 = B/x**2
    >>> p3 = (C*x + D)/(x**2 + 1)

    >>> p1, p2, p3
    ⎛A  B   C⋅x + D⎞
    ⎜─, ──, ───────⎟
    ⎜x   2    2    ⎟
    ⎝   x    x  + 1⎠

Let's add them together to get the desired form::

    >>> h = sum(_)
    >>> h
    A   B    C⋅x + D
    ─ + ── + ───────
    x    2     2
        x     x  + 1

The next step is to rewrite this expression as rational function in `x`::

    >>> together(h)
        ⎛ 2    ⎞     ⎛ 2    ⎞    2
    A⋅x⋅⎝x  + 1⎠ + B⋅⎝x  + 1⎠ + x ⋅(C⋅x + D)
    ────────────────────────────────────────
                   2 ⎛ 2    ⎞
                  x ⋅⎝x  + 1⎠

    >>> factor(_, x)
               3            2
    A⋅x + B + x ⋅(A + C) + x ⋅(B + D)
    ─────────────────────────────────
                2 ⎛ 2    ⎞
               x ⋅⎝x  + 1⎠

Let's now visually compare the last expression with `f`::

    >>> Eq(_, f)
               3            2
    A⋅x + B + x ⋅(A + C) + x ⋅(B + D)        1
    ───────────────────────────────── = ───────────
                2 ⎛ 2    ⎞               2 ⎛ 2    ⎞
               x ⋅⎝x  + 1⎠              x ⋅⎝x  + 1⎠

Our task boils down to finding `A`, `B`, `C` and `D`. We notice that
denominators are equal so we will proceed only with numerators::

    >>> eq = Eq(numer(_.lhs), numer(_.rhs))
    >>> eq
               3            2
    A⋅x + B + x ⋅(A + C) + x ⋅(B + D) = 1

To solve this equation, we use :func:`solve_undetermined_coeffs`::

    >>> solve_undetermined_coeffs(eq, [A, B, C, D], x)
    {A: 0, B: 1, C: 0, D: -1}

This gave us values for our parameters, which now can be put into the initial
expression::

    >>> h.subs(_)
        1      1
    - ────── + ──
       2        2
      x  + 1   x

This result is identical to the result we got from ``apart(f)``. Suppose
however, we would like to see how undetermined coefficients method works.
First we have to extract coefficients of `x` of both sides of the equation::

    >>> lhs, rhs = Poly(eq.lhs, x), Poly(eq.rhs, x)

    >>> lhs
    Poly((A + C)*x**3 + (B + D)*x**2 + A*x + B, x, domain='ZZ[A,B,C,D]')
    >>> rhs
    Poly(1, x, domain='ZZ')

Now we can use :func:`Poly.nth` to obtain coefficients of `x`::

    >>> [ Eq(lhs.nth(i), rhs.nth(i)) for i in range(4) ]
    [B = 1, A = 0, B + D = 0, A + C = 0]

Solving this system of linear equations gives the same solution set as
previously::

    >>> solve(_)
    {A: 0, B: 1, C: 0, D: -1}

    >>> f.subs(_)
         1
    ───────────
     2 ⎛ 2    ⎞
    x ⋅⎝x  + 1⎠

There are several other ways we can approach undetermined coefficients
method. For example we could use :func:`collect` for this::

    >>> collect(eq.lhs - eq.rhs, x, evaluate=False)
    ⎧                 2          3       ⎫
    ⎨1: B - 1, x: A, x : B + D, x : A + C⎬
    ⎩                                    ⎭

    >>> solve(_.values())
    {A: 0, B: 1, C: 0, D: -1}

Notice that even though the expressions were not :func:`Eq`'s, this still
worked. This is because SymPy assumes by default that expressions are
identically equal to 0, so ``solve(Eq(expr, 0))`` is the same as
``solve(expr)``.

This approach is even simpler than using :func:`Poly.nth`. Finally we use a
little trick with :class:`Symbol` and visually present solution to partial
fraction decomposition of `f`::

    >>> Eq(Symbol('apart')(f), f.subs(_))
         ⎛     1     ⎞        1
    apart⎜───────────⎟ = ───────────
         ⎜ 2 ⎛ 2    ⎞⎟    2 ⎛ 2    ⎞
         ⎝x ⋅⎝x  + 1⎠⎠   x ⋅⎝x  + 1⎠

Tasks
-----

1. Compute partial fraction decomposition of:

  * `\frac{3 x + 5}{(2 x + 1)^2}`
  * `\frac{3 x + 5}{(u x + v)^2}`
  * `\frac{(3 x + 5)^2}{(2 x + 1)^2}`

  (:ref:`solution <solution_partfrac_1>`)

2. Can you use :func:`Expr.coeff` in place of :func:`Poly.nth`?

  (:ref:`solution <solution_partfrac_2>`)
