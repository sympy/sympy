===========================================
Hongguang Fu's Trigonometric Simplification
===========================================

.. automodule:: sympy.simplify.fu
.. currentmodule:: sympy.simplify.fu

Implementation of the trigsimp algorithm by Fu et al.

The idea behind the Fu algorithm is to use a sequence of rules
that students learn during their pre-calculus courses.
The rules are applied heuristically and it uses a greedy algorithm to
apply multiple rules simultaneously and choose the result with the least
leaf counts.

There are transform rules in which a single rule is applied to the
expression tree. The following are just mnemonic in nature; see the
docstrings for examples.

- :func:`TR0` - simplify expression
- :func:`TR1` - sec-csc to cos-sin
- :func:`TR2` - tan-cot to sin-cos ratio
- :func:`TR2i` - sin-cos ratio to tan
- :func:`TR3` - angle canonicalization
- :func:`TR4` - functions at special angles
- :func:`TR5` - powers of sin to powers of cos
- :func:`TR6` - powers of cos to powers of sin
- :func:`TR7` - reduce cos power (increase angle)
- :func:`TR8` - expand products of sin-cos to sums
- :func:`TR9` - contract sums of sin-cos to products
- :func:`TR10` - separate sin-cos arguments
- :func:`TR10i` - collect sin-cos arguments
- :func:`TR11` - reduce double angles
- :func:`TR12` - separate tan arguments
- :func:`TR12i` - collect tan arguments
- :func:`TR13` - expand product of tan-cot
- :func:`TRmorrie` - prod(cos(x*2**i), (i, 0, k - 1)) -> sin(2**k*x)/(2**k*sin(x))
- :func:`TR14` - factored powers of sin or cos to cos or sin power
- :func:`TR15` - negative powers of sin to cot power
- :func:`TR16` - negative powers of cos to tan power
- :func:`TR22` - tan-cot powers to negative powers of sec-csc functions
- :func:`TR111` - negative sin-cos-tan powers to csc-sec-cot

There are 4 combination transforms (CTR1 - CTR4) in which a sequence of
transformations are applied and the simplest expression is selected from
a few options.

Finally, there are the 2 rule lists (RL1 and RL2), which apply a
sequence of transformations and combined transformations, and the ``fu``
algorithm itself, which applies rules and rule lists and selects the
best expressions. There is also a function ``L`` which counts the number
of trigonometric functions that appear in the expression.

Other than TR0, re-writing of expressions is not done by the transformations.
e.g. TR10i finds pairs of terms in a sum that are in the form like
``cos(x)*cos(y) + sin(x)*sin(y)``. Such expression are targeted in a bottom-up
traversal of the expression, but no manipulation to make them appear is
attempted. For example,

Set-up for examples below:

    >>> from sympy.simplify.fu import fu, L, TR9, TR10i, TR11
    >>> from sympy import factor, sin, cos, powsimp
    >>> from sympy.abc import x, y, z, a
    >>> from time import time

    >>> eq = cos(x + y)/cos(x)
    >>> TR10i(eq.expand(trig=True))
    -sin(x)*sin(y)/cos(x) + cos(y)

If the expression is put in "normal" form (with a common denominator) then
the transformation is successful:

    >>> TR10i(_.normal())
    cos(x + y)/cos(x)

TR11's behavior is similar. It rewrites double angles as smaller angles but
doesn't do any simplification of the result.

    >>> TR11(sin(2)**a*cos(1)**(-a), 1)
    (2*sin(1)*cos(1))**a/cos(1)**a
    >>> powsimp(_)
    (2*sin(1))**a

The temptation is to try make these TR rules "smarter" but that should really
be done at a higher level; the TR rules should try maintain the "do one thing
well" principle.  There is one exception, however. In TR10i and TR9 terms are
recognized even when they are each multiplied by a common factor:

    >>> fu(a*cos(x)*cos(y) + a*sin(x)*sin(y))
    a*cos(x - y)

Factoring with ``factor_terms`` is used but it is "JIT"-like, being delayed
until it is deemed necessary. Furthermore, if the factoring does not
help with the simplification, it is not retained, so
``a*cos(x)*cos(y) + a*sin(x)*sin(z)`` does not become a factored
(but unsimplified in the trigonometric sense) expression:

    >>> fu(a*cos(x)*cos(y) + a*sin(x)*sin(z))
    a*sin(x)*sin(z) + a*cos(x)*cos(y)

In some cases factoring might be a good idea, but the user is left
to make that decision. For example:

    >>> expr=((15*sin(2*x) + 19*sin(x + y) + 17*sin(x + z) + 19*cos(x - z) +
    ... 25)*(20*sin(2*x) + 15*sin(x + y) + sin(y + z) + 14*cos(x - z) +
    ... 14*cos(y - z))*(9*sin(2*y) + 12*sin(y + z) + 10*cos(x - y) + 2*cos(y -
    ... z) + 18)).expand(trig=True).expand()

In the expanded state, there are nearly 1000 trig functions:

    >>> L(expr)
    932

If the expression were factored first, this would take time but the
resulting expression would be transformed very quickly:

    >>> def clock(f, n=2):
    ...    t=time(); f(); return round(time()-t, n)
    ...
    >>> clock(lambda: factor(expr))  # doctest: +SKIP
    0.86
    >>> clock(lambda: TR10i(expr), 3)  # doctest: +SKIP
    0.016

If the unexpanded expression is used, the transformation takes longer but
not as long as it took to factor it and then transform it:

    >>> clock(lambda: TR10i(expr), 2)  # doctest: +SKIP
    0.28

So neither expansion nor factoring is used in ``TR10i``: if the
expression is already factored (or partially factored) then expansion
with ``trig=True`` would destroy what is already known and take
longer; if the expression is expanded, factoring may take longer than
simply applying the transformation itself.

Although the algorithms should be canonical, always giving the same
result, they may not yield the best result. This, in general, is
the nature of simplification where searching all possible transformation
paths is very expensive. Here is a simple example. There are 6 terms
in the following sum:

    >>> expr = (sin(x)**2*cos(y)*cos(z) + sin(x)*sin(y)*cos(x)*cos(z) +
    ... sin(x)*sin(z)*cos(x)*cos(y) + sin(y)*sin(z)*cos(x)**2 + sin(y)*sin(z) +
    ... cos(y)*cos(z))
    >>> args = expr.args

Serendipitously, fu gives the best result:

    >>> fu(expr)
    3*cos(y - z)/2 - cos(2*x + y + z)/2

But if different terms were combined, a less-optimal result might be
obtained, requiring some additional work to get better simplification,
but still less than optimal. The following shows an alternative form
of ``expr`` that resists optimal simplification once a given step
is taken since it leads to a dead end:

    >>> TR9(-cos(x)**2*cos(y + z) + 3*cos(y - z)/2 +
    ...     cos(y + z)/2 + cos(-2*x + y + z)/4 - cos(2*x + y + z)/4)
    sin(2*x)*sin(y + z)/2 - cos(x)**2*cos(y + z) + 3*cos(y - z)/2 + cos(y + z)/2

Here is a smaller expression that exhibits the same behavior:

    >>> a = sin(x)*sin(z)*cos(x)*cos(y) + sin(x)*sin(y)*cos(x)*cos(z)
    >>> TR10i(a)
    sin(x)*sin(y + z)*cos(x)
    >>> newa = _
    >>> TR10i(expr - a)  # this combines two more of the remaining terms
    sin(x)**2*cos(y)*cos(z) + sin(y)*sin(z)*cos(x)**2 + cos(y - z)
    >>> TR10i(_ + newa) == _ + newa  # but now there is no more simplification
    True

Without getting lucky or trying all possible pairings of arguments, the
final result may be less than optimal and impossible to find without
better heuristics or brute force trial of all possibilities.

Rules
=====

.. autofunction:: TR0

.. autofunction:: TR1

.. autofunction:: TR2

.. autofunction:: TR2i

.. autofunction:: TR3

.. autofunction:: TR4

.. autofunction:: TR5

.. autofunction:: TR6

.. autofunction:: TR7

.. autofunction:: TR8

.. autofunction:: TR9

.. autofunction:: TR10

.. autofunction:: TR10i

.. autofunction:: TR11

.. autofunction:: TR12

.. autofunction:: TR12i

.. autofunction:: TR13

.. autofunction:: TRmorrie

.. autofunction:: TR14

.. autofunction:: TR15

.. autofunction:: TR16

.. autofunction:: TR111

.. autofunction:: TR22

.. autofunction:: TRpower

.. autofunction:: fu

Notes
=====

This work was started by Dimitar Vlahovski at the Technological School
"Electronic systems" (30.11.2011).

Beyond TR13, other rules are not from the original paper, but extended
in SymPy.

References
==========

.. [1] Fu, Hongguang, Xiuqin Zhong, and Zhenbing Zeng.
    "Automated and readable simplification of trigonometric expressions."
    Mathematical and computer modelling 44.11 (2006): 1169-1177.
    https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.657.2478&rep=rep1&type=pdf

.. [2] A formula sheet for trigonometric functions.
    http://www.sosmath.com/trig/Trig5/trig5/pdf/pdf.html
