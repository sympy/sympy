"""
Implementation of the trigsimp algorithm by Fu et al.

The idea behind the ``fu`` algorithm is to use a sequence of rules, applied
in what is heuristically known to be a smart order, to select a simpler
expression that is equivalent to the input.

There are 15 transform rules (TR0 - TR10, TR10i, TR11 - TR13) in which a
single rule is applied to the expression tree. The following are just
mnemonic in nature; see the docstrings for examples.

    TR0 - non-trig simplification
    TR1 - sec, csc -> 1/cos, 1/sin
    TR2 - tan, cot -> sin/cos and cos/sin
    TR3 - angle canonicalization
    TR4 - replace functions at special angles with values
    TR5 - sin**2(a) -> (1 - cos(a)**2) and sin(a)**4 to the square of the same
    TR6 - cos**2(a) -> (1 - sin(a)**2) and cos(a)**4 to the square of the same
    TR7 - cos(a)**2 -> cos(2*a)
    TR8 - f(a)*g(b) -> f(a+b) + g(a-b);             f and g are sin or cos
    TR9 - f(a)*g(b) <- f(a+b) + g(a-b);             f and g are sin or cos
    TR10 - f(a+b) -> f(a)*g(b) + g(a)*f(b);         f and g are sin or cos
    TR10i - f(a+b) <- f(a)*g(b) + g(a)*f(b);        f and g are sin or cos
    TR11 - f(2*m*a) -> 2*f(m*a)*g(m*a);             f and g are sin or cos
    TR12 - tan(a + b) -> (tan sum)/(tan product)
    TR13 - f(a)*f(b) -> 1 +/- (f(a)+f(b))*g(a+b);   f = tan or cot

There are 4 combination transforms (CTR1 - CTR4) in which a seqence of
transformations are applied and the simplest expression is selected from
a few options.

Finally, there are the 2 rule lists (RL1 and RL2), which apply a
sequence of transformations and combined transformations, and the ``fu``
algorithm itself, which applies rules and rule lists and selects the
best expressions. There is also a function ``L`` which counts the number
of trigonometric funcions that appear in the expression.

Other than TR0, re-writing of expressions is not done by the transformations.
e.g. TR10i finds pairs of terms in a sum that are in the form like
``cos(x)*cos(y) + sin(x)*sin(y)``. Such expression are targeted in a bottom-up
traversal of the expression, but no manipulation to make
them appear is attempted. For example,

    Set-up for examples below:

    >>> from sympy.simplify.fu import fu, L, TR9, TR10i
    >>> from sympy import factor, sin, cos
    >>> from sympy.abc import x, y, z, a
    >>> from time import time

>>> eq = cos(x + y)/cos(x)
>>> TR10i(eq.expand(trig=True)) == eq
True
>>> TR10i(eq.expand(trig=True).expand())
-sin(x)*sin(y)/cos(x) + cos(y)

If the expression is put in "normal" form (with a common denominator) then
the transformation is successful:

>>> TR10i(_.normal())
cos(x + y)/cos(x)

There is one exception, however. In TR10i and TR9 terms are recognized even
when they are each multiplied by a common factor:

>>> fu(a*cos(x)*cos(y) + a*sin(x)*sin(y))
a*cos(x - y)

Factoring with ``factor_terms`` is used but it it "JIT"-like, being delayed
until it is deemed necessary. Furthermore, if the factoring does not
help with the simplification, it is not retained, so
``a*cos(x)*cos(y) + a*sin(x)*sin(z)`` does not become the factored
(but unsimplified in the trigonometric sense) expression:

>>> fu(a*cos(x)*cos(y) + a*sin(x)*sin(z))
a*sin(x)*sin(z) + a*cos(x)*cos(y)

In some cases factoring might be a good idea, but the user is left
to make that decision. For example:

>>> expr=((15*sin(2*x) + 19*sin(x + y) + 17*sin(x + z) + 19*cos(x - z) +
... 25)*(20*sin(2*x) + 15*sin(x + y) + sin(y + z) + 14*cos(x - z) +
... 14*cos(y - z))*(9*sin(2*y) + 12*sin(y + z) + 10*cos(x - y) + 2*cos(y -
... z) + 18)).expand(trig=True).expand()

In the expanded state, there are nearly 1000 trig functions

>>> L(expr)
932

If the expression where factored first, this would take time but the
resulting expression would be ransformed very quickly:

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

So neither expansion nor fatoring is used in ``TR10i``: if the
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

Notes
=====

This work was started by Dimitar Vlahovski at the Technological School
"Electronic systems" (30.11.2011).

References
==========
http://rfdz.ph-noe.ac.at/fileadmin/Mathematik_Uploads/ACDCA/
DESTIME2006/DES_contribs/Fu/simplification.pdf
"""

from collections import defaultdict

from sympy.simplify.simplify import simplify, powsimp, ratsimp, combsimp
from sympy.core.sympify import sympify
from sympy.functions.elementary.trigonometric import cos, sin, tan, cot, sec, csc, sqrt
from sympy.functions.elementary.hyperbolic import cosh, sinh
from sympy.core.compatibility import ordered, combinations
from sympy.core.core import C
from sympy.core.mul import Mul
from sympy.core.function import expand_mul, count_ops
from sympy.core.add import Add
from sympy.core.symbol import Dummy
from sympy.core.exprtools import Factors
from sympy.core.rules import Transform
from sympy.core.basic import S
from sympy.core.numbers import Integer, pi
from sympy.rules import minimize, chain, debug
from sympy.rules.strat_pure import identity
from sympy import SYMPY_DEBUG
from sympy.ntheory.factor_ import perfect_power


def TR0(rv):
    """Simplification of rational polynomials, trying to simplify
    the expression, e.g. combine things like 3*x + 2*x, etc....
    """
    rv = simplify(rv, fu=True)  # disable calls to trigsimp
    rv = powsimp(rv)
    rv = ratsimp(rv)
    rv = combsimp(rv)
    rv = expand_mul(rv)
    return rv


def TR1(rv):
    """Replace sec, csc with 1/cos, 1/sin

    Examples
    ========

    >>> from sympy.simplify.fu import TR1, sec, csc
    >>> from sympy.abc import x
    >>> TR1(2*csc(x) + sec(x))
    1/cos(x) + 2/sin(x)
    """
    rv = bottom_up(rv, TR1)
    if rv.func is sec:
        a = rv.args[0]
        return S.One/cos(a)
    elif rv.func is csc:
        a = rv.args[0]
        return S.One/sin(a)
    return rv


def TR2(rv):
    """Replace tan and cot with sin/cos and cos/sin

    Examples
    ========

    >>> from sympy.simplify.fu import TR2
    >>> from sympy.abc import x
    >>> from sympy import tan, cot, sin, cos
    >>> TR2(tan(x))
    sin(x)/cos(x)
    >>> TR2(cot(x))
    cos(x)/sin(x)
    >>> TR2(tan(tan(x) - sin(x)/cos(x)))
    0

    """
    rv = bottom_up(rv, TR2)
    if rv.func is tan:
        a = rv.args[0]
        return sin(a)/cos(a)
    elif rv.func is cot:
        a = rv.args[0]
        return cos(a)/sin(a)
    return rv


def TR3(rv):
    """Induced formula: example sin(-a) = -sin(a)

    Examples
    ========

    >>> from sympy.simplify.fu import TR3
    >>> from sympy.abc import x, y
    >>> from sympy import pi
    >>> from sympy import cos
    >>> TR3(cos(y - x*(y - x)))
    cos(x*(x - y) + y)
    >>> cos(pi/2 + x)
    -sin(x)
    >>> cos(30*pi/2 + x)
    -cos(x)

    >>> from sympy.utilities.randtest import test_numerically
    >>> from sympy import cos, sin, tan, cot, csc, sec
    >>> for f in (cos, sin, tan, cot, csc, sec):
    ...   i = f(3*pi/7)
    ...   j = TR3(i)
    ...   assert test_numerically(i, j) and i.func != j.func
    ...
    """
    from sympy.simplify.simplify import signsimp

    # Negative argument (already automatic for funcs like sin(-x) -> -sin(x)
    # but more complicated expressions can use it, too). Also, trig angles
    # between pi/4 and pi/2 are not reduced to an angle between 0 and pi/4.
    rv = bottom_up(rv, TR3)
    if not isinstance(rv, C.TrigonometricFunction):
        return rv
    rv = rv.func(signsimp(rv.args[0]))
    if S.Pi/4 < rv.args[0] < S.Pi/2:
        fmap = {cos: sin, sin: cos, tan: cot, cot: tan, sec: csc, csc: sec}
        rv = fmap[rv.func](S.Pi/2 - rv.args[0])
    return rv
    #The following are automatically handled
    #Argument of type: pi/2 +/- angle
    #Argument of type: pi +/- angle
    #Argument of type : 2k*pi +/- angle


def TR4(rv):
    """Identify values of special angles.

        a=  0   pi/6        pi/4        pi/3        pi/2
    ----------------------------------------------------
    cos(a)  0   1/2         sqrt(2)/2   sqrt(3)/2   1
    sin(a)  1   sqrt(3)/2   sqrt(2)/2   1/2         0
    tan(a)  0   sqt(3)/3    1           sqrt(3)     --

    Examples
    ========

    >>> from sympy.simplify.fu import TR4
    >>> from sympy import pi
    >>> from sympy import cos, sin, tan, cot
    >>> for s in (0, pi/6, pi/4, pi/3, pi/2):
    ...    for f in (cos, sin, tan, cot):
    ...      print f(s),
    ...    print
    ...
    1 0 0 zoo
    sqrt(3)/2 1/2 sqrt(3)/3 sqrt(3)
    sqrt(2)/2 sqrt(2)/2 1 1
    1/2 sqrt(3)/2 sqrt(3) sqrt(3)/3
    0 1 zoo 0
    """
    # special values at 0, pi/6, pi/4, pi/3, pi/2 already handled
    return rv


def _TR56(rv, f, g, max, pow):
    """Helper for TR5 and TR6 to replace f**2 with 1 - g**2

    Options
    =======

    max :   controls size of exponent that can appear on f
            e.g. if max=4 then f**4 will be changed to (1 - g**2)**2
            and (if pow is False) then f**6 will be changed to (1 - g**2)**3
    pow :   controls whether the exponent must be a perfect power of 2
            e.g. if pow=True (and max >= 6) then f**6 will not be changed
            but f**8 will be changed to (1 - g**2)**4

    >>> from sympy.simplify.fu import _TR56 as T
    >>> from sympy.abc import x
    >>> from sympy import sin, cos
    >>> T(sin(x)**3, sin, cos, 4, False)
    sin(x)**3
    >>> T(sin(x)**6, sin, cos, 6, False)
    (-cos(x)**2 + 1)**3
    >>> T(sin(x)**6, sin, cos, 6, True)
    sin(x)**6
    >>> T(sin(x)**8, sin, cos, 10, True)
    (-cos(x)**2 + 1)**4
    """

    rv = bottom_up(rv, lambda x: _TR56(x, f, g, max, pow))

    # I'm not sure if this transformation should target all even powers
    # or only those expressible as powers of 2. Also, should it only
    # make the changes in powers that appear in sums -- making an isolated
    # change is not going to allow a simplification as far as I can tell.
    if not (rv.is_Pow and rv.base.func == f):
        return rv

    if rv.exp < 0:
        return rv  # or return 1/_TR56(1/rv, f, g, max)
    if rv.exp > max:
        return rv
    if rv.exp == 2:
        return 1 - g(rv.base.args[0])**2
    else:
        if rv.exp == 4:
            e = 2
        elif not pow:
            if rv.exp % 2:
                return rv
            e = rv.exp//2
        else:
            p = perfect_power(rv.exp)
            if not p:
                return rv
            e = rv.exp//2
        return (1 - g(rv.base.args[0])**2)**e

def TR5(rv):
    """Replacement of sin**2 with 1 - cos(x)**2.

    Examples
    ========

    >>> from sympy.simplify.fu import TR5
    >>> from sympy.abc import x
    >>> from sympy import sin
    >>> TR5(sin(x)**2)
    -cos(x)**2 + 1
    >>> TR5(sin(x)**-2)  # unchanged
    sin(x)**(-2)
    >>> TR5(sin(x)**4)
    (-cos(x)**2 + 1)**2
    """
    return _TR56(rv, sin, cos, max=4, pow=False)


def TR6(rv):
    """Replacement of cos**2 with 1 - sin(x)**2.

    Examples
    ========

    >>> from sympy.simplify.fu import TR6
    >>> from sympy.abc import x
    >>> from sympy import cos
    >>> TR6(cos(x)**2)
    -sin(x)**2 + 1
    >>> TR6(cos(x)**-2)  #unchanged
    cos(x)**(-2)
    >>> TR6(cos(x)**4)
    (-sin(x)**2 + 1)**2
    """
    return _TR56(rv, cos, sin, max=4, pow=False)


def TR7(rv):
    """Lowering the degree of cos(x)**2

    Examples
    ========

    >>> from sympy.simplify.fu import TR7
    >>> from sympy.abc import x
    >>> from sympy import cos
    >>> TR7(cos(x)**2)
    cos(2*x)/2 + 1/2
    >>> TR7(cos(x)**2 + 1)
    cos(2*x)/2 + 3/2

    """
    rv = bottom_up(rv, TR7)
    if not (rv.is_Pow and rv.base.func == cos and rv.exp == 2):
        return rv
    return (1 + cos(2*rv.base.args[0]))/2


def TR8(rv):
    """Converting products of ``cos`` and/or ``sin`` to a sum or
    difference of ``cos`` and or ``sin`` terms.

    Examples
    ========

    >>> from sympy.simplify.fu import TR8
    >>> from sympy import cos, sin
    >>> TR8(cos(2)*cos(3))
    cos(5)/2 + cos(1)/2
    >>> TR8(cos(2)*sin(3))
    sin(5)/2 + sin(1)/2
    >>> TR8(sin(2)*sin(3))
    -cos(5)/2 + cos(1)/2
    >>> TR8(cos(2)*cos(3)*cos(4)*cos(5))  # (cos(5)/2 + cos(1)/2)*(cos(9)/2 + cos(1)/2)
    cos(4)/4 + cos(10)/8 + cos(8)/8 + cos(14)/8 + cos(1)**2/4 + cos(6)/8
    >>> TR8(cos(2)*cos(3)*cos(4)*cos(5)*cos(6))  # (cos(1)/2 + cos(7)/2)*(cos(11)/2 + cos(1)/2)*cos(2)
    cos(10)/8 + cos(16)/16 + cos(4)/16 + cos(1)**2*cos(2)/4 + cos(2)/16 + cos(8)/8 + cos(14)/16 + cos(20)/16 + cos(12)/16 + cos(6)/8
    """
    rv = bottom_up(rv, TR8)
    if not rv.is_Mul:
        return rv

    args = {cos: [], sin: [], None: []}
    for a in ordered(rv.args):
        if a.func in (cos, sin):
            args[a.func].append(a.args[0])
        else:
            args[None].append(a)
    c = args[cos]
    s = args[sin]
    if not (c and s or len(c) > 1 or len(s) > 1):
        return rv
    args = args[None]
    len0 = len(args)
    n = min(len(c), len(s))
    for i in range(n):
        a1 = s.pop()
        a2 = c.pop()
        args.append((sin(a1 + a2) + sin(a1 - a2))/2)
    while len(c) > 1:
        a1 = c.pop()
        a2 = c.pop()
        args.append((cos(a1 + a2) + cos(a1 - a2))/2)
    if c:
        args.append(cos(c.pop()))
    while len(s) > 1:
        a1 = s.pop()
        a2 = s.pop()
        args.append((-cos(a1 + a2) + cos(a1 - a2))/2)
    if s:
        args.append(sin(s.pop()))
    rv = Mul(*args)
    if len(args) - len0 > 1:
        rv = TR8(expand_mul(rv))
    return rv


def TR9(rv):
    """Sum of ``cos`` or ``sin`` terms as a product of ``cos`` or ``sin``.

    Examples
    ========

    >>> from sympy.simplify.fu import TR9
    >>> from sympy import cos, sin
    >>> TR9(cos(1) + cos(2))
    2*cos(1/2)*cos(3/2)
    >>> TR9(cos(1) - cos(2))
    2*sin(1/2)*sin(3/2)
    >>> TR9(sin(1) - sin(2))
    -2*sin(1/2)*cos(3/2)
    >>> TR9(sin(1) + sin(2))
    2*sin(3/2)*cos(1/2)
    >>> TR9(cos(1) + 2*sin(1) + 2*sin(2))
    cos(1) + 4*sin(3/2)*cos(1/2)
    >>> TR9(cos(4) + cos(2) + 2*cos(1)*cos(3))
    4*cos(1)*cos(3)
    >>> TR9((cos(4) + cos(2))/cos(3)/2 + cos(3))
    2*cos(1)*cos(2)

    If no change is made by TR9, no re-arrangement of the
    expression will be made. For example, though factoring
    of common term is attempted, if the factored expression
    wasn't changed, the original expression will be returned:

    >>> TR9(cos(3) + cos(3)*cos(2))
    cos(3) + cos(2)*cos(3)

    >>> from sympy.abc import x
    >>> from sympy import Add, Mul
    >>> c = cos(x); s = sin(x)
    >>> for si in ((1,1),(1,-1),(-1,1),(-1,-1)):
    ...   for a in ((c, s), (s, c)):
    ...    args = zip(si, a)
    ...    ex = Add(*[Mul(*ai) for ai in args])
    ...    t = TR9(ex)
    ...    if (a[0].func == a[1].func and (ex - t.expand(trig=True) or t.is_Add)
    ...         or a[1].func != a[0].func and ex != t):
    ...        print 'fail', ex
    ...
    """
    rv = bottom_up(rv, TR9)
    if not rv.is_Add:
        return rv

    def do(rv, first=True):
        # cos(a)+/-cos(b) can be combined into a product of cosines and
        # sin(a)+/-sin(b) can be combined into a product of cosine and sine.
        #
        # If there are more than two args, the pairs which "work" will have
        # a gcd extractable and the remaining two terms will have the above
        # structure -- all pairs must be checked to find the ones that work.
        # args that don't have a common set of symbols are skipped since this
        # doesn't lead to a simpler formula and also has the arbitrariness of
        # combining, for example, the x and y term instead of the y and z term
        # in something like cos(x) + cos(y) + cos(z).


        if not rv.is_Add:
            return rv

        args = list(ordered(rv.args))
        if len(args) != 2:
            hit = False
            for i in range(len(args)):
                ai = args[i]
                if ai is None:
                    continue
                for j in range(i + 1, len(args)):
                    aj = args[j]
                    if aj is None:
                        continue
                    was = ai + aj
                    new = do(was)
                    if new != was:
                        args[i] = new  # update in place
                        args[j] = None
                        hit = True
                        break  # go to next i
            if hit:
                rv = Add(*filter(None, args))
                if rv.is_Add:
                    rv = do(rv)

            return rv

        # two-arg Add
        split = trig_split(*args)
        if not split:
            return rv
        gcd, n1, n2, a, b, iscos = split

        # application of rule if possible
        if iscos:
            if n1 == n2:
                return gcd*n1*2*cos((a + b)/2)*cos((a - b)/2)
            if n1 < 0:
                a, b = b, a
            return -2*gcd*sin((a + b)/2)*sin((a - b)/2)
        else:
            if n1 == n2:
                return gcd*n1*2*sin((a + b)/2)*cos((a - b)/2)
            if n1 < 0:
                a, b = b, a
            return 2*gcd*cos((a + b)/2)*sin((a - b)/2)
        return rv

    return process_common_addends(rv, do)  # DON'T sift by free symbols


def TR10(rv, first=True):
    """Separate sums in ``cos`` and ``sin``.

    Examples
    ========

    >>> from sympy.simplify.fu import TR10
    >>> from sympy.abc import a, b, c
    >>> from sympy import cos, sin
    >>> TR10(cos(a + b))
    -sin(a)*sin(b) + cos(a)*cos(b)
    >>> TR10(sin(a + b))
    sin(a)*cos(b) + sin(b)*cos(a)
    >>> TR10(sin(a + b + c))
    (-sin(a)*sin(b) + cos(a)*cos(b))*sin(c) + (sin(a)*cos(b) + sin(b)*cos(a))*cos(c)
    """
    rv = bottom_up(rv, TR10)
    if not rv.func in (cos, sin):
        return rv

    f = rv.func
    arg = rv.args[0]
    if arg.is_Add:
        if first:
            args = list(ordered(arg.args))
        else:
            args = list(arg.args)
        a = args.pop()
        b = Add._from_args(args)
        if b.is_Add:
            if f == sin:
                return sin(a)*TR10(cos(b), first=False) + cos(a)*TR10(sin(b), first=False)
            else:
                return cos(a)*TR10(cos(b), first=False) - sin(a)*TR10(sin(b), first=False)
        else:
            if f == sin:
                return sin(a)*cos(b) + cos(a)*sin(b)
            else:
                return cos(a)*cos(b) - sin(a)*sin(b)
    return rv


def TR10i(rv):
    """Sum of products to function of sum.

    Examples
    ========

    >>> from sympy.simplify.fu import TR10i
    >>> from sympy import cos, sin, pi, Add, Mul, sqrt, Symbol
    >>> from sympy.abc import x, y

    >>> TR10i(cos(1)*cos(3) + sin(1)*sin(3))
    cos(2)
    >>> TR10i(cos(1)*cos(3) - sin(1)*sin(3))
    cos(4)
    >>> TR10i(cos(1)*sin(3) - sin(1)*cos(3))
    sin(2)
    >>> TR10i(cos(1)*sin(3) + sin(1)*cos(3))
    sin(4)
    >>> TR10i(cos(1)*sin(3) + sin(1)*cos(3) + 7)
    sin(4) + 7
    >>> TR10i(cos(1)*sin(3) + sin(1)*cos(3) + cos(3))
    cos(3) + sin(4)
    >>> TR10i(2*cos(1)*sin(3) + 2*sin(1)*cos(3)+cos(3))
    2*sin(4) + cos(3)
    >>> TR10i(cos(2)*cos(3)+sin(2)*(cos(1)*sin(2)+cos(2)*sin(1)))
    cos(1)
    >>> eq = (cos(2)*cos(3)+sin(2)*(cos(1)*sin(2)+cos(2)*sin(1)))*cos(5) + sin(1)*sin(5)
    >>> TR10i(eq) == TR10i(eq.expand()) == cos(4)
    True
    >>> TR10i(sqrt(2)*cos(x)*x + sqrt(6)*sin(x)*x)
    2*sqrt(2)*x*sin(x + pi/6)
    >>> TR10i(cos(x)/sqrt(6) + sin(x)/sqrt(2) + cos(x)/sqrt(6)/3 + sin(x)/sqrt(2)/3)
    4*sqrt(6)*sin(x + pi/6)/9
    >>> TR10i(cos(x)/sqrt(6) + sin(x)/sqrt(2) + cos(y)/sqrt(6)/3 + sin(y)/sqrt(2)/3)
    sqrt(6)*sin(x + pi/6)/3 + sqrt(6)*sin(y + pi/6)/9
    >>> TR10i(cos(x) + sqrt(3)*sin(x) + 2*sqrt(3)*cos(x + pi/6))
    4*cos(x)
    >>> TR10i(cos(x) + sqrt(3)*sin(x) + 2*sqrt(3)*cos(x + pi/6) + 4*sin(x))
    4*sqrt(2)*sin(x + pi/4)

    >>> A = Symbol('A', commutative=False)
    >>> TR10i(sqrt(2)*cos(x)*A + sqrt(6)*sin(x)*A)
    2*sqrt(2)*sin(x + pi/6)*A


    >>> c = cos(x); s = sin(x); h = sin(y); r = cos(y)
    >>> for si in ((1,1),(1,-1),(-1,1),(-1,-1)):
    ...   for a in ((c*r, s*h), (c*h,s*r)): # explicit 2-args
    ...    args = zip(si, a)
    ...    ex = Add(*[Mul(*ai) for ai in args])
    ...    t = TR10i(ex)
    ...    if ex - t.expand(trig=True) or t.is_Add:
    ...        print 'fail', ex
    ...

    >>> c = cos(x); s = sin(x); h = sin(pi/6); r = cos(pi/6)
    >>> for si in ((1,1),(1,-1),(-1,1),(-1,-1)):
    ...   for a in ((c*r, s*h), (c*h,s*r)): # induced
    ...    args = zip(si, a)
    ...    ex = Add(*[Mul(*ai) for ai in args])
    ...    t = TR10i(ex)
    ...    if ex - t.expand(trig=True) or t.is_Add:
    ...     print 'fail', ex
    ...
    """
    global _ROOT2, _ROOT3, _invROOT3
    if _ROOT2 is None:
        _roots()

    rv = bottom_up(rv, TR10i)
    if not rv.is_Add:
        return rv

    def do(rv, first=True):
        # args which can be expressed as A*(cos(a)*cos(b) +/- sin(a)*sin(b)) or
        # B*(cos(a)*sin(b) +/- cos(b)*sin(a)) can be combined into A*f(a+/-b)
        # where f is either sin or cos.
        #
        # If there are more than two args, the pairs which "work" will have
        # a gcd extractable and the remaining two terms will have the above
        # structure -- all pairs must be checked to find the ones that work.

        if not rv.is_Add:
            return rv

        args = list(ordered(rv.args))
        if len(args) != 2:
            hit = False
            for i in range(len(args)):
                ai = args[i]
                if ai is None:
                    continue
                for j in range(i + 1, len(args)):
                    aj = args[j]
                    if aj is None:
                        continue
                    was = ai + aj
                    new = do(was)
                    if new != was:
                        args[i] = new  # update in place
                        args[j] = None
                        hit = True
                        break  # go to next i
            if hit:
                rv = Add(*filter(None, args))
                if rv.is_Add:
                    rv = do(rv)

            return rv

        # two-arg Add
        split = trig_split(*args, **dict(two=True))
        if not split:
            return rv
        gcd, n1, n2, a, b, same = split

        # identify and get c1 to be cos then apply rule if possible
        if same:  # coscos, sinsin
            gcd = n1*gcd
            if n1 == n2:
                return gcd*cos(a - b)
            return gcd*cos(a + b)
        else:  #cossin, cossin
            gcd = n1*gcd
            if n1 == n2:
                return gcd*sin(a + b)
            return gcd*sin(b - a)
        return rv

    rv = process_common_addends(rv, do, lambda x: tuple(ordered(x.free_symbols)))

    # need to check for induceable pairs in ratio of sqrt(3):1 that appeared
    # in different lists when sorting by coefficient
    while rv.is_Add:
        byrad = defaultdict(list)
        for a in rv.args:
            hit = 0
            if a.is_Mul:
                for ai in a.args:
                    if ai.is_Pow and ai.exp is S.Half and ai.base.is_Integer:
                        byrad[ai].append(a)
                        hit = 1
                        break
            if not hit:
                byrad[S.One].append(a)

        # no need to check all pairs -- just check for the onees
        # that have the right ratio
        args = []
        for a in byrad:
            for b in [_ROOT3*a, _invROOT3]:
                if b in byrad:
                    for i in range(len(byrad[a])):
                        if byrad[a][i] is None:
                            continue
                        for j in range(len(byrad[b])):
                            if byrad[b][j] is None:
                                continue
                            was = Add(byrad[a][i] + byrad[b][j])
                            new = do(was)
                            if new != was:
                                args.append(new)
                                byrad[a][i] = None
                                byrad[b][j] = None
                                break
        if args:
            rv = Add(*(args + [Add(*filter(None, v)) for v in byrad.values()]))
        else:
            rv = do(rv)  # final pass to resolve any new unduceable pairs
            break

    return rv


def TR11(rv):
    """Function of double angle to product.

    Examples
    ========

    >>> from sympy.simplify.fu import TR11
    >>> from sympy import cos, sin
    >>> from sympy.abc import x
    >>> TR11(sin(2*x))
    2*sin(x)*cos(x)
    >>> TR11(sin(4*x))
    4*(-sin(x)**2 + cos(x)**2)*sin(x)*cos(x)
    >>> TR11(sin(4*x/3))
    4*(-sin(x/3)**2 + cos(x/3)**2)*sin(x/3)*cos(x/3)

    >>> TR11(cos(2*x))
    -sin(x)**2 + cos(x)**2
    >>> TR11(cos(4*x))
    (-sin(x)**2 + cos(x)**2)**2 - 4*sin(x)**2*cos(x)**2

    If the arguments are simply integers, no change is made:

    >>> TR11(cos(2))
    cos(2)

    """
    rv = bottom_up(rv, TR11)
    if not (rv.func in (sin, cos) and not rv.args[0].is_Number):
        return rv

    c, m = rv.args[0].as_coeff_Mul()
    if c.p % 2 == 0:
        arg = c.p//2*m/c.q
        c = TR11(cos(arg))
        s = TR11(sin(arg))
        if rv.func == sin:
            rv = 2*s*c
        else:
            rv = c**2 - s**2
    return rv


def TR12(rv):
    """Separate sums in ``tan``.

    Examples
    ========

    >>> from sympy.simplify.fu import TR12
    >>> from sympy.abc import x, y
    >>> from sympy import tan
    >>> from sympy.simplify.fu import TR12
    >>> TR12(tan(x + y))
    (tan(x) + tan(y))/(-tan(x)*tan(y) + 1)
    """
    rv = bottom_up(rv, TR12)
    if not rv.func == tan:
        return rv

    arg = rv.args[0]  # should expand_mul be used?
    if arg.is_Add:
        a, b = arg.as_two_terms()
        if b.is_Add:
            tb = TR12(tan(b))
        else:
            tb = tan(b)
        return (tan(a) + tb)/(1 - tan(a)*tb)
    return rv


def TR13(rv):
    """Change products of ``tan`` or ``cot``.

    Examples
    ========

    >>> from sympy.simplify.fu import TR13
    >>> from sympy import tan, cot, cos
    >>> TR13(tan(3)*tan(2))
    -(tan(2) + tan(3))*cot(5) + 1
    >>> TR13(cot(3)*cot(2))
    1 + (cot(3) + cot(2))*cot(5)
    >>> TR13((1 + tan(3)*tan(2)/(1 - (tan(2) + tan(3))*cot(5)))*cos(3))
    2*cos(3)
    """
    rv = bottom_up(rv, TR13)
    if not rv.is_Mul:
        return rv

    args = {tan: [], cot: [], None: []}
    for a in ordered(Mul.make_args(rv)):
        if a.func in (tan, cot):
            args[a.func].append(a.args[0])
        else:
            args[None].append(a)
    t = args[tan]
    c = args[cot]
    if len(t) < 2 and len(c) < 2:
        return rv
    args = args[None]
    while len(t) > 1:
        t1 = t.pop()
        t2 = t.pop()
        args.append(1 - (tan(t1) + tan(t2))*cot(t1 + t2))
    if t:
        args.append(tan(t.pop()))
    while len(c) > 1:
        t1 = c.pop()
        t2 = c.pop()
        args.append(1 + (cot(t1) + cot(t2))*cot(t1 + t2))
    if c:
        args.append(cot(t.pop()))
    return Mul(*args)


def L(rv):
    """Return count of trigonometric functions in expression.

    Examples
    ========

    >>> from sympy.simplify.fu import L
    >>> from sympy.abc import x
    >>> from sympy import cos, sin
    >>> L(cos(x)+sin(x))
    2
    """
    return S(rv.count(C.TrigonometricFunction))

if SYMPY_DEBUG:
    TR0,TR1,TR2,TR3,TR4,TR5,TR6,TR7,TR8,TR9,TR10,TR11,TR12,TR13 = map(debug,
            (TR0,TR1,TR2,TR3,TR4,TR5,TR6,TR7,TR8,TR9,TR10,TR11,TR12,TR13))

_CTR1 = [TR5, TR0], [TR6, TR0], [identity]

_CTR2 = [TR11, TR5, TR0], [TR11, TR6, TR0], [TR11, TR0]

_CTR3 = [TR8, TR0], [TR8, TR10i, TR0], [identity]

_CTR4 = [TR4, TR10i], [identity]


def CTRstrat(lists):
    return minimize(*[chain(*list) for list in lists],
        **dict(objective=lambda x: (L(x), x.count_ops())))

CTR1, CTR2, CTR3, CTR4 = map(CTRstrat, (_CTR1, _CTR2, _CTR3, _CTR4))

_RL1 = [TR4, TR3, TR4, TR12, TR4, TR13, TR4, TR0]


# XXX it's a little unclear how this one is to be implemented
# see Fu paper of reference, page 7. What is the Union symbol refering to?
# The diagram shows all these as one chain of transformations, but the
# text refers to them being applied independently. Also, a break
# if L starts to increase has not been implemented.
_RL2 = [
    [TR4, TR3, TR10, TR4, TR3, TR11],
    [TR5, TR7, TR11, TR4],
    [CTR3, TR0, CTR1, TR9, CTR2, TR4, TR9, TR0, TR9, CTR4],
    [identity],
    ]


def RLstrat(rls):
    return chain(*rls)

RL1 = RLstrat(_RL1)
RL2 = CTRstrat(_RL2)


def fu(rv):
    """Attempt to simplify expression by using transformation rules given
    in the algorithm by Fu et al.

    Examples
    ========

    >>> from sympy.simplify.fu import fu
    >>> from sympy import cos, sin, tan, pi, S, sqrt
    >>> from sympy.abc import x, y, a, b

    >>> fu(sin(50)**2 + cos(50)**2 + sin(pi/6))
    3/2
    >>> fu(sqrt(6)*cos(x) + sqrt(2)*sin(x))
    2*sqrt(2)*sin(x + pi/3)

    CTR1 example

    >>> eq = sin(x)**4 - cos(y)**2 + sin(y)**2 + 2*cos(x)**2
    >>> fu(eq)
    cos(x)**4 - 2*cos(y)**2 + 2

    CTR2 example

    >>> fu(S.Half - cos(2*x)/2)
    sin(x)**2

    CTR3 example

    >>> fu(sin(a)*(cos(b) - sin(b)) + cos(a)*(sin(b) + cos(b)))
    sqrt(2)*sin(a + b + pi/4)

    CTR4 example

    >>> fu(sqrt(3)*cos(x)/2 + sin(x)/2)
    sin(x + pi/3)

    Example 1

    >>> fu(1-sin(2*x)**2/4-sin(y)**2-cos(x)**4)
    -cos(x)**2 + cos(y)**2

    Example 2

    >>> fu(cos(4*pi/9))
    sin(pi/18)
    >>> fu(cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9))
    1/16

    Example 3

    >>> fu(tan(7*pi/18)+tan(5*pi/18)-sqrt(3)*tan(5*pi/18)*tan(7*pi/18))
    -sqrt(3)

    References
    ==========
    http://rfdz.ph-noe.ac.at/fileadmin/Mathematik_Uploads/ACDCA/
    DESTIME2006/DES_contribs/Fu/simplification.pdf
    """
    was = rv
    rv = sympify(rv)
    rv = TR0(rv)
    rv = TR1(rv)
    if rv.has(tan, cot):
        rv1 = RL1(rv)
        if (_L(rv1) < _L(rv)):
            rv = rv1
        if rv.has(tan, cot):
            rv = TR2(rv)
    rv = TR0(rv)
    if rv.has(sin, cos):
        rv1 = RL2(rv)
        rv2 = TR8(rv1)
        rv = ordered([was, rv, rv1, rv2], keys=(L, count_ops), default=False).next()
    return rv

_L = lambda x: (L(x), x.count_ops())

def bottom_up(rv, F):
    """Apply ``F`` to all expressions in an expression tree from the
    bottom up.
    """
    if rv.args:
        args = tuple([F(a) for a in rv.args])
        if args != rv.args:
            rv = rv.func(*args)
    return rv


def process_common_addends(rv, do, key2=None, key1=True):
    """Apply ``do`` to addends of ``rv`` that (if key1=True) share at least
    a common absolute value of their coefficient and the value of key2 when
    applied to the argument. If key1 is False key2 must be supplied and will
    be the only key applied.
    """

    # collect by absolute value of coefficient and key2
    absc = defaultdict(list)
    if key1:
        for a in rv.args:
            c, a = a.as_coeff_Mul()
            if c < 0:
                c = -c
                a = -a  # put the sign on `a`
            absc[(c, key2(a) if key2 else 1)].append(a)
    elif key2:
        for a in rv.args:
            absc[(S.One, key2(a))].append(a)
    else:
        raise ValueError('must have at least one key')

    args = []
    hit = False
    for k in absc:
        v = absc[k]
        c, _ = k
        if len(v) > 1:
            e = Add(*v)
            new = do(e)
            if new != e:
                e = new
                hit = True
            args.append(c*e)
        else:
            args.append(c*v[0])
    if hit:
        rv = Add(*args)

    return rv


FU = dict(zip('TR0, TR1, TR2, TR3, TR4, TR5, TR6, TR7, TR8, TR9, TR10, TR10i, TR11, TR12, TR13, CTR1, CTR2, CTR3, CTR4, RL1, RL2, L'.split(', '),
              (TR0, TR1, TR2, TR3, TR4, TR5, TR6, TR7, TR8, TR9, TR10, TR10i, TR11, TR12, TR13, CTR1, CTR2, CTR3, CTR4, RL1, RL2, L)))


def _roots():
    global _ROOT2, _ROOT3, _invROOT3
    _ROOT2, _ROOT3 = sqrt(2), sqrt(3)
    _invROOT3 = 1/_ROOT3
_ROOT2 = None


def trig_split(a, b, two=False):
    """Return the gcd, s1, s2, a1, a2, bool where

    If two is False (default) then::
        a + b = gcd*(s1*f(a1) + s2*f(a2)) where f = cos if bool else sin
    else:
        if bool, a + b was +/- cos(a1)*cos(a2) +/- sin(a1)*sin(a2) and equals
            n1*gcd*cos(a - b) if n1 == n2 else
            n1*gcd*cos(a + b)
        else a + b was +/- cos(a1)*sin(a2) +/- sin(a1)*cos(a2) and equals
            n1*gcd*sin(a + b) if n1 = n2 else
            n1*gcd*sin(b - a)

    Examples
    ========
    >>> from sympy.simplify.fu import trig_split
    >>> from sympy.abc import x, y, z
    >>> from sympy import cos, sin, sqrt

    >>> trig_split(cos(x), cos(y))
    (1, 1, 1, x, y, True)
    >>> trig_split(2*cos(x), -2*cos(y))
    (2, 1, -1, x, y, True)
    >>> trig_split(cos(x)*sin(y), cos(y)*sin(y))
    (sin(y), 1, 1, x, y, True)

    >>> trig_split(cos(x), -sqrt(3)*sin(x), two=True)
    (2, 1, -1, x, pi/6, False)
    >>> trig_split(cos(x), sin(x), two=True)
    (sqrt(2), 1, 1, x, pi/4, False)
    >>> trig_split(cos(x), -sin(x), two=True)
    (sqrt(2), 1, -1, x, pi/4, False)
    >>> trig_split(sqrt(2)*cos(x), -sqrt(6)*sin(x), two=True)
    (2*sqrt(2), 1, -1, x, pi/6, False)
    >>> trig_split(-sqrt(6)*cos(x), -sqrt(2)*sin(x), two=True)
    (-2*sqrt(2), 1, 1, x, pi/3, False)
    >>> trig_split(cos(x)/sqrt(6), sin(x)/sqrt(2), two=True)
    (sqrt(6)/3, 1, 1, x, pi/6, False)
    >>> trig_split(-sqrt(6)*cos(x)*sin(y), -sqrt(2)*sin(x)*sin(y), two=True)
    (-2*sqrt(2)*sin(y), 1, 1, x, pi/3, False)

    >>> trig_split(cos(x), sin(x))
    >>> trig_split(cos(x), sin(z))
    >>> trig_split(2*cos(x), -sin(x))
    >>> trig_split(cos(x), -sqrt(3)*sin(x))
    >>> trig_split(cos(x)*cos(y), sin(x)*sin(z))
    >>> trig_split(cos(x)*cos(y), sin(x)*sin(y))
    >>> trig_split(-sqrt(6)*cos(x), sqrt(2)*sin(x)*sin(y), two=True)
    """
    global _ROOT2, _ROOT3, _invROOT3
    if _ROOT2 is None:
        _roots()

    a, b = [Factors(i) for i in (a, b)]
    ua, ub = a.normal(b)
    gcd = a.gcd(b).as_expr()
    if S.NegativeOne in ua.factors:
        ua = ua.quo(S.NegativeOne)
        n1 = -1
        n2 = 1
    elif S.NegativeOne in ub.factors:
        ub = ub.quo(S.NegativeOne)
        n1 = 1
        n2 = -1
    else:
        n1 = n2 = 1
    a, b = [i.as_expr() for i in (ua, ub)]

    def pow_cos_sin(a, two):
        """Return ``a`` as a tuple (r, c, s) such that ``a = (r or 1)*(c or 1)*(s or 1)``.

        Three arguments are returned (radical, c-factor, s-factor) as
        long as the conditions set by ``two`` are met; otherwise None is
        returned. If ``two`` is True there will be one or two non-None
        values in the tuple: c and s or c and r or s and r or s or c with c
        being a cosine function (if possible) else a sine, and s being a sine
        function (if possible) else oosine. If ``two`` is False then there
        will only be a c or s term in the tuple.

        ``two`` also require that either two cos and/or sin be present (with
        the condition that if the functions are the same the arguments are
        different or vice versa) or that a single cosine or a single sine
        be present with an optional radical.

        If the above conditions dictated by ``two`` are not met then None
        is returned.
        """
        c = s = None
        co = S.One
        if a.is_Mul:
            co, a = a.as_coeff_Mul()
            if len(a.args) > 2 or not two:
                return None
            if a.is_Mul:
                args = list(a.args)
            else:
                args = [a]
            a = args.pop(0)
            if a.func is cos:
                c = a
            elif a.func is sin:
                s = a
            elif a.is_Pow and a.exp is S.Half:  # autoeval doesn't allow -1/2
                co *= a
            else:
                return None
            if args:
                b = args[0]
                if b.func is cos:
                    if c:
                        s = b
                    else:
                        c = b
                elif b.func is sin:
                    if s:
                        c = b
                    else:
                        s = b
                elif b.is_Pow and b.exp is S.Half:
                    co *= b
                else:
                    return None
            return co if co is not S.One else None, c, s
        elif a.func is cos:
            c = a
        elif a.func is sin:
            s = a
        if c is None and s is None:
            return
        co = co if co is not S.One else None
        if not two and co:
            return None
        return co, c, s

    # get the parts
    m = pow_cos_sin(a, two)
    if m is None:
        return
    coa, ca, sa = m
    m = pow_cos_sin(b, two)
    if m is None:
        return
    cob, cb, sb = m

    # check them
    if (not ca) and cb or ca and ca.func is sin:
        coa, ca, sa, cob, cb, sb = cob, cb, sb, coa, ca, sa
        n1, n2 = n2, n1
    if two is False:  # need cos(x) and cos(y) or sin(x) and sin(y)
        if ca and sa or cb and sb:
            return
        c = ca or sa
        s = cb or sb
        if c.func is not s.func:
            return None
        return gcd, n1, n2, c.args[0], s.args[0], c.func is cos
    else:
        if not coa and not cob:
            if (ca and cb and sa and sb):
                if not ((ca.func is sa.func) is (cb.func is sb.func)):
                    return
                args = set([j.args for j in (ca, sa)])
                if not all(i.args in args for i in (cb, sb)):
                    return
                return gcd, n1, n2, ca.args[0], sa.args[0], ca.func is sa.func
        if ca and sa or cb and sb:
            return
        c = ca or sa
        s = cb or sb
        if c.args != s.args:
            return
        if c.func is s.func:
            return
        if not coa:
            coa = S.One
        if not cob:
            cob = S.One
        if coa is cob:
            gcd *= _ROOT2
            return gcd, n1, n2, c.args[0], pi/4, False
        elif coa/cob == _ROOT3:
            gcd *= 2*cob
            return gcd, n1, n2, c.args[0], pi/3, False
        elif coa/cob == _invROOT3:
            gcd *= 2*coa
            return gcd, n1, n2, c.args[0], pi/6, False
