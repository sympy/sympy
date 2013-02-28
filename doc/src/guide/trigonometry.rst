
.. include:: ../definitions.def

Deriving trigonometric identities
=================================

Let's assume that we need a formula for `\sin(a + b)` in terms of `\sin(a)`,
`\sin(b)`, `\cos(a)` and `\cos(b)`, but we don't remember it, nor do we
know how to get it easily with SymPy. We will derive this formula from
scratch using Taylor series expansions and a little symbolic manipulation.

Let's start with definition of symbols and the expression in consideration::

    >>> var('a,b')
    (a, b)

    >>> f = sin(a + b)
    >>> f
    sin(a + b)

Now let's expand `f` as a power series with respect to `b` around 0::

    >>> f.series(b, 0, 10)
                         2           3           4           5           6           7           8           9
                        b ⋅sin(a)   b ⋅cos(a)   b ⋅sin(a)   b ⋅cos(a)   b ⋅sin(a)   b ⋅cos(a)   b ⋅sin(a)   b ⋅cos(a)
    sin(a) + b⋅cos(a) - ───────── - ───────── + ───────── + ───────── - ───────── - ───────── + ───────── + ───────── + O(b**10)
                            2           6           24         120         720         5040       40320       362880

This isn't very readable but we can clearly see a pattern around `\sin(a)`
and `\cos(a)`. Let's collect terms with respect to those two expressions::

    >>> collect(_, [sin(a), cos(a)])
    ⎛   9       7      5    3    ⎞          ⎛   8      6    4    2    ⎞
    ⎜  b       b      b    b     ⎟          ⎜  b      b    b    b     ⎟
    ⎜────── - ──── + ─── - ── + b⎟⋅cos(a) + ⎜───── - ─── + ── - ── + 1⎟⋅sin(a) + O(b**10)
    ⎝362880   5040   120   6     ⎠          ⎝40320   720   24   2     ⎠

    >>> _.removeO()
    ⎛   8      6    4    2    ⎞          ⎛   9       7      5    3    ⎞
    ⎜  b      b    b    b     ⎟          ⎜  b       b      b    b     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅sin(a) + ⎜────── - ──── + ─── - ── + b⎟⋅cos(a)
    ⎝40320   720   24   2     ⎠          ⎝362880   5040   120   6     ⎠

    >>> g = _

We got two subexpression that look very familiar. Let's expand `\sin(b)`
in `b` around 0 and remove the order term::

    >>> sin(b).series(b, 0, 10)
         3     5     7       9
        b     b     b       b
    b - ── + ─── - ──── + ────── + O(b**10)
        6    120   5040   362880

    >>> _.removeO()
       9       7      5    3
      b       b      b    b
    ────── - ──── + ─── - ── + b
    362880   5040   120   6

This is clearly the second subexpression, so let's substitute it for
`\sin(b)`::

    >>> g.subs(_, sin(b))
    ⎛   8      6    4    2    ⎞
    ⎜  b      b    b    b     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅sin(a) + sin(b)⋅cos(a)
    ⎝40320   720   24   2     ⎠

    >>> h = _

Now let's repeat this procedure for `\cos(b)`::

    >>> cos(b).series(b, 0, 10)
         2    4     6      8
        b    b     b      b
    1 - ── + ── - ─── + ───── + O(b**10)
        2    24   720   40320

    >>> _.removeO()
       8      6    4    2
      b      b    b    b
    ───── - ─── + ── - ── + 1
    40320   720   24   2

    >>> h.subs(_, cos(b))
    sin(a)⋅cos(b) + sin(b)⋅cos(a)

This gave us a formula for `\sin(a + b)`::

    >>> Eq(f, _)
    sin(a + b) = sin(a)⋅cos(b) + sin(b)⋅cos(a)

There is, however, a much simpler way to get the same result::

    >>> Eq(f, sin(a + b).expand(trig=True))
    sin(a + b) = sin(a)⋅cos(b) + sin(b)⋅cos(a)

Tasks
-----

1. Repeat this procedure but expand with respect to `a` in the first step.

   (:ref:`solution <solution_trig_1>`)

2. Use this procedure to derive a formula for `\cos(a + b)`.

   (:ref:`solution <solution_trig_2>`)
