
.. include:: ../definitions.def

Summing roots of polynomials
============================

Let's suppose we are given a univariate polynomial `f(z)` and a univariate
rational function `g(z)`, and we wish to compute:

.. math::

    g(r_1) + g(r_2) + \ldots + g(r_n)

where `r_i` for `i = 1 \ldots n` are the roots of `f` (i.e. `f(r_i) = 0`).

In theory this is a very simple task. We just have to compute roots of `f`,
using the :func:`roots` function, substitute those roots for `z` in `g` and add
resulting values together.

Let's consider the following polynomial and rational function::

    >>> f = z**5 + z + 3
    >>> f
     5
    z  + z + 3

    >>> g = 1/z
    >>> g
    1
    ─
    z

Following the trivial approach, let's compute the roots of `f`::

    >>> roots(f)
    {}

We got a very unfortunate result: no roots! By the fundamental theorem
of algebra we should get five, possibly complex, roots, including
multiplicities. Unfortunately, there is no way to express roots in terms
of radicals of some polynomials of degree five and higher. For certain
instances of polynomials of this kind it may be possible to compute
their roots (e.g. :func:`roots` recognizes cyclotomic polynomials of
high degree), but in general we will most likely be unlucky.

Instead, we could switch to numerical root finding algorithms and compute
approximations of roots of `f` and proceed with summation of roots. This
can be done by using :func:`nroots`::

    >>> R = nroots(f)

    >>> for ri, r in zip(numbered_symbols('r'), R):
    ...     pprint(Eq(ri, r))
    ...
    r₀ = -1.13299756588507
    r₁ = -0.47538075666955 - 1.12970172509541⋅ⅈ
    r₂ = -0.47538075666955 + 1.12970172509541⋅ⅈ
    r₃ = 1.04187953961208 - 0.822870338109958⋅ⅈ
    r₄ = 1.04187953961208 + 0.822870338109958⋅ⅈ

We can substitute those roots for `z` in `g` and add together::

    >>> sum([ g.subs(z, r) for r in R ]).evalf(chop=True)
    -0.333333333333332

It was necessary to evaluate this sum with :func:`evalf`, because otherwise
we would get an unsimplified result. The additional parameter ``chop=True`` was
necessary to remove a tiny and insignificant imaginary part. Next we can use
:func:`nsimplify` to get an exact result from numerical approximation::

    >>> nsimplify(_)
    -1/3

Is this result correct? The best way is to figure out a purely symbolic
method that doesn't require computing roots of `f`. In SymPy it possible
to represent a root of a univariate polynomial with rational coefficients
using :class:`RootOf`::

    >>> RootOf(f, 0)
          ⎛ 5           ⎞
    RootOf⎝z  + z + 3, 0⎠

    >>> _.evalf()
    -1.13299756588507

We can obtain all roots using list comprehensions::

    >>> R = [ RootOf(f, i) for i in xrange(degree(f)) ]

    >>> for r in R:
    ...     pprint(r)
    ...
          ⎛ 5           ⎞
    RootOf⎝z  + z + 3, 0⎠
          ⎛ 5           ⎞
    RootOf⎝z  + z + 3, 1⎠
          ⎛ 5           ⎞
    RootOf⎝z  + z + 3, 2⎠
          ⎛ 5           ⎞
    RootOf⎝z  + z + 3, 3⎠
          ⎛ 5           ⎞
    RootOf⎝z  + z + 3, 4⎠

Alternatively we can use ``Poly(f).all_roots()`` which gives the same
result, but is much faster when `f` is a composite polynomial, because
the preprocessing step in :class:`RootOf` is executed only once.

Unfortunately we can't get anywhere from here, because SymPy is not yet
capable of simplifying expressions with :class:`RootOf`::

    >>> G = sum([ g.subs(z, r) for r in R ])
    >>> isinstance(G, Add)
    True

    >>> _ = simplify(G)
    >>> isinstance(_, Add)
    True

We can, however, evaluate sums of :class:`RootOf`'s using :func:`evalf`::

    >>> G.evalf()
    -0.333333333333333

    >>> nsimplify(_)
    -1/3

which gave us the same result as before. The difference is that now numerical
approximations of roots of `f` were computed using a hybrid symbolic--numeric
method, where first disjoint isolating intervals (rectangles) where computed
for all roots of `f` and then a numerical root finding algorithm was used in
each interval.

Let's approach this problem differently, using a purely symbolic
approach. We know that a polynomial of degree `n` has exactly `n`
complex roots, counting multiplicities. In our case `f` has five roots::

    >>> R = var('r:5')
    >>> R
    (r₀, r₁, r₂, r₃, r₄)

Let's now substitute those "roots" for `z` in `g`::

    >>> [ g.subs(z, r) for r in R ]
    ⎡1   1   1   1   1 ⎤
    ⎢──, ──, ──, ──, ──⎥
    ⎣r₀  r₁  r₂  r₃  r₄⎦

and add those expressions together::

    >>> sum(_)
    1    1    1    1    1
    ── + ── + ── + ── + ──
    r₄   r₃   r₂   r₁   r₀

We got a sum of simple rational functions. The next step is to put those
rational functions over a common denominator::

    >>> G = together(_)
    >>> G
    r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄
    ───────────────────────────────────────────────────────────────────
                               r₀⋅r₁⋅r₂⋅r₃⋅r₄

We got very peculiar numerator and denominator, which are very specific
functions of roots of `f` (symmetric polynomials). Polynomials of this
kind can be generated using :func:`viete`::

    >>> V = viete(f, R, z)

    >>> for lhs, rhs in V:
    ....     pprint(Eq(lhs, rhs))
    ....
    r₀ + r₁ + r₂ + r₃ + r₄ = 0
    r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r₃⋅r₄ = 0
    r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r₃ + r₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄ = 0
    r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄ = 1
    r₀⋅r₁⋅r₂⋅r₃⋅r₄ = -3

Viete formulas show the relationship between roots of a polynomial and
its coefficients:

.. math::

    V_{i-1} = (-1)^i \frac{a_{n-i}}{a_n}

where `f(z)=a_nz^n + a_{n-1}z^{n-1} + \ldots + a_1z + a_0` and `i = 1 \ldots n`. To obtain the final
result it sufficient to take `V_3` and `V_4` and substitute in `G`::

    >>> numer(G).subs(*V[3])/denom(G).subs(*V[4])
    -1/3

Or we could simply use ``G.subs(V)``, but due to a bug in SymPy (`#2552 <http://code.google.com/p/sympy/issues/detail?id=2552>`_) this
doesn't work as expected, leaving the denominator unchanged.

We obtained the same result as before, just this time using purely symbolic
techniques. This simple procedure can be extended to form an algorithm for
solving the root summation problem in the general setup. SymPy implements this
algorithm in :class:`RootSum`::

    >>> RootSum(f, Lambda(z, g))
    -1/3

The choice of `g` allowed us to recognize Viete formulas very easily in
`G`, but is this the case also for more complicated rational functions?
Let's modify `g` a little::

    >>> g = 1/(z + 2)
      1
    ─────
    z + 2

Now let's repeat the procedure for the new `g`::

    >>> G = together(sum([ g.subs(z, r) for r in R ]))

    >>> p = expand(numer(G))
    >>> q = expand(denom(G))

    >>> p
    r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + 4⋅r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃⋅r₄ + 4⋅r₀⋅r₁⋅r₃ + 4⋅r₀⋅r₁⋅r₄ + 12⋅r₀⋅r₁ + r₀⋅r₂⋅r₃⋅r₄ + \
    4⋅r₀⋅r₂⋅r₃ + 4⋅r₀⋅r₂⋅r₄ + 12⋅r₀⋅r₂ + 4⋅r₀⋅r₃⋅r₄ + 12⋅r₀⋅r₃ + 12⋅r₀⋅r₄ + 32⋅r₀ + r₁⋅r₂⋅r₃⋅r₄ + 4⋅r₁⋅r₂⋅r₃ + \
    4⋅r₁⋅r₂⋅r₄ + 12⋅r₁⋅r₂ + 4⋅r₁⋅r₃⋅r₄ + 12⋅r₁⋅r₃ + 12⋅r₁⋅r₄ + 32⋅r₁ + 4⋅r₂⋅r₃⋅r₄ + 12⋅r₂⋅r₃ + 12⋅r₂⋅r₄ + 32⋅r₂ + \
    12⋅r₃⋅r₄ + 32⋅r₃ + 32⋅r₄ + 80

    >>> q
    r₀⋅r₁⋅r₂⋅r₃⋅r₄ + 2⋅r₀⋅r₁⋅r₂⋅r₃ + 2⋅r₀⋅r₁⋅r₂⋅r₄ + 4⋅r₀⋅r₁⋅r₂ + 2⋅r₀⋅r₁⋅r₃⋅r₄ + 4⋅r₀⋅r₁⋅r₃ + 4⋅r₀⋅r₁⋅r₄ + \
    8⋅r₀⋅r₁ + 2⋅r₀⋅r₂⋅r₃⋅r₄ + 4⋅r₀⋅r₂⋅r₃ + 4⋅r₀⋅r₂⋅r₄ + 8⋅r₀⋅r₂ + 4⋅r₀⋅r₃⋅r₄ + 8⋅r₀⋅r₃ + 8⋅r₀⋅r₄ + 16⋅r₀ + \
    2⋅r₁⋅r₂⋅r₃⋅r₄ + 4⋅r₁⋅r₂⋅r₃ + 4⋅r₁⋅r₂⋅r₄ + 8⋅r₁⋅r₂ + 4⋅r₁⋅r₃⋅r₄ + 8⋅r₁⋅r₃ + 8⋅r₁⋅r₄ + 16⋅r₁ + 4⋅r₂⋅r₃⋅r₄ + \
    8⋅r₂⋅r₃ + 8⋅r₂⋅r₄ + 16⋅r₂ + 8⋅r₃⋅r₄ + 16⋅r₃ + 16⋅r₄ + 32

This doesn't look that familiar anymore. Let's try to apply Viete formulas
to the numerator and denominator::

    >>> p.subs(V).has(*R)
    True
    >>> q.subs(V).has(*R)
    True

We weren't able to get rid of the symbolic roots of `f`. We can, however, try
to rewrite `p` and `q` as polynomials in elementary symmetric polynomials.
This procedure is called symmetric reduction, and an algorithm for this is
implemented in :func:`symmetrize`::

    >>> (P, Q), mapping = symmetrize((p, q), R, formal=True)

    >>> P
    (32⋅s₁ + 12⋅s₂ + 4⋅s₃ + s₄ + 80, 0)
    >>> Q
    (16⋅s₁ + 8⋅s₂ + 4⋅s₃ + 2⋅s₄ + s₅ + 32, 0)

    >>> for s, poly in mapping:
    ...     pprint(Eq(s, poly))
    ...
    s₁ = r₀ + r₁ + r₂ + r₃ + r₄
    s₂ = r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r₃⋅r₄
    s₃ = r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r₃ + r₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄
    s₄ = r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄
    s₅ = r₀⋅r₁⋅r₂⋅r₃⋅r₄

Here we performed the formal simultaneous symmetric reduction of the polynomials `p`
and `q`, obtaining their representation in terms of elementary symmetric
polynomials, non-symmetric remainders, and elementary symmetric polynomials.
Remainders are always zero for symmetric inputs.

We can zip this mapping and Viete formulas together, obtaining::

    >>> [ (s, c) for (s, _), (_, c) in zip(mapping, V) ]
    [(s₁, 0), (s₂, 0), (s₃, 0), (s₄, 1), (s₅, -3)]

Now we can take head of ``P`` and ``Q`` and perform substitution::

    >>> P[0].subs(_)/Q[0].subs(_)
    81
    ──
    31

Let's verify this result using :class:`RootSum`::

    >>> RootSum(f, Lambda(z, g))
    81
    ──
    31

The numerical approach also works in this case::

    >>> sum([ g.subs(z, r) for r in Poly(f).all_roots() ]).evalf()
    2.61290322580645

    >>> nsimplify(_)
    81
    ──
    31

Tasks
-----

1. Repeat this procedure for:

 * `f = z^5 + z + a` and `g = \frac{1}{z + 1}`
 * `f = z^5 + z + a` and `g = \frac{1}{z + b}`

 (:ref:`solution <solution_rootsum_1>`)

2. Can this or a similar procedure be used with other classes of expressions
   than rational functions? If so, what kind of expressions can be allowed?

 (:ref:`solution <solution_rootsum_2>`)
