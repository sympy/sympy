"""
Implementation of the trigsimp algorithm by Fu et al

The idea behind the ``fu`` algorithm is to use a sequence of rules, applied
in what is heuristically known to be a smart order, to select a simpler
expression that is equivalent to the input.

There are 15 transform rules (TR0 - TR10, TR10i, TR11 - TR13) in which a
single rule is applied to the expression tree. The following are not
mnemonic in nature; see the docstrings for examples.

    TR0 - non-trig simplification
    TR1 - sec, csc -> 1/cos, 1/sin
    TR2 - tan, cot -> sin/cos and cos/sin
    TR3 - sign simplification of arguments
    TR4 - special angle values
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
traversal of the expression, but if they don't appear no manipulation to make
them appear is attempted, being left instead to the user. For example,

    Set-up for examples below:

    >>> from sympy.simplify.fu import fu, L, TR9, TR10i
    >>> from sympy import factor, sin, cos
    >>> from sympy.abc import x, y, z, A
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

>>> fu(A*cos(x)*cos(y) + A*sin(x)*sin(y))
A*cos(x - y)

Factoring with ``factor_terms`` is used but it it "JIT"-like, being delayed
until it is deemed necessary. Furthermore, if the factoring does not
help with the simplification, it is not retained, so
``A*cos(x)*cos(y) + A*sin(x)*sin(z)`` does not become the factored
(but unsimplified in the trigonometric sense) expression:

>>> fu(A*cos(x)*cos(y) + A*sin(x)*sin(z))
A*sin(x)*sin(z) + A*cos(x)*cos(y)

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
resulting expression would be factored very quickly:

>>> def clock(f, n=2):
...    t=time(); f(); return round(time()-t, n)
...
>>> clock(lambda: factor(expr))  # doctest: +SKIP
0.86
>>> clock(lambda: TR10i(expr), 3)  # doctest: +SKIP
0.016

If the unexpanded expression is used, the time takes longer but
not as long as the combined factoring and simplification times:

>>> clock(lambda: TR10i(expr), 2)  # doctest: +SKIP
0.28

So neither expansion nor fatoring is used in ``TR10i``: if the
expression is already factored (or partially factored) then expansion
with ``trig=True`` would destroy what is already known and take
longer; if the expression is expanded, factoring may be more time
consuming than the algorithm itself.

Although the algorithms should be canonical, always giving the same
result, they may not yield the best result. This, in general, is
the nature of simplification where searching all possible transformation
paths is very expensive. Here is a simple example. There are 6 terms
in the following sum:

>>> expr = (sin(x)**2*cos(y)*cos(z) + sin(x)*sin(y)*cos(x)*cos(z) +
... sin(x)*sin(z)*cos(x)*cos(y) + sin(y)*sin(z)*cos(x)**2 + sin(y)*sin(z) +
... cos(y)*cos(z))
>>> args = expr.args

Although ``TR10i`` currently returns the best result, though the full
Fu algoirthm gives a different result (which can be simplified a little
more with single application of TR9):

>>> TR10i(expr)
sin(x + y)*sin(x + z) + cos(y - z)
>>> fu(expr)
-cos(x)**2*cos(y + z) + 3*cos(y - z)/2 + cos(y + z)/2 + cos(-2*x + y + z)/4 - cos(2*x + y + z)/4
>>> TR9(_)
sin(2*x)*sin(y + z)/2 - cos(x)**2*cos(y + z) + 3*cos(y - z)/2 + cos(y + z)/2

The order in which terms are combined at each step can lead to a dead-end
short of the optimal expression:

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
from sympy.functions.elementary.trigonometric import cos, sin, tan, cot
from sympy.functions.elementary.hyperbolic import cosh, sinh
from sympy.core.compatibility import ordered
from sympy.core.core import C
from sympy.core.mul import Mul
from sympy.core.function import expand_mul, count_ops
from sympy.core.add import Add
from sympy.core.symbol import Wild
from sympy.core.exprtools import Factors
from sympy.core.rules import Transform
from sympy.core.basic import S
from sympy.core.numbers import Integer
from sympy.rules import minimize, chain, debug
from sympy.rules.strat_pure import identity
from sympy import SYMPY_DEBUG


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
    """Replace sec, csc with 1/cos, 1/sin"""
    try:
        rv = bottom_up(rv, TR1)
        if rv.func == sec:
            a = rv.args[0]
            return S.One/cos(a)
        elif rv.func == csc:
            a = rv.args[0]
            return S.One/sin(a)
    except:
        pass
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
    if rv.func == tan:
        a = rv.args[0]
        return sin(a)/cos(a)
    elif rv.func == cot:
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
        fmap = {cos: sin, sin: cos, tan: cot, cot: tan}
        try:
            fmap[sec] = csc
            fmap[csc] = sec
        except:
            pass
        rv = fmap[rv.func](S.Pi/2 - rv.args[0])
    return rv
    #The following are automatically handled
    #Argument of type: pi/2 +/- angle
    #Argument of type: pi +/- angle
    #Argument of type : 2k*pi +/- angle


def TR4(rv):
    """Identify values of special angles.

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


def TR5(rv):
    """Replacement of sin**2 or sin**4 factors.

    Examples
    ========

    >>> from sympy.simplify.fu import TR5
    >>> from sympy.abc import x
    >>> from sympy import sin
    >>> TR5(sin(x)**2)
    -cos(x)**2 + 1
    >>> TR5(sin(x)**4)
    (-cos(x)**2 + 1)**2
    >>> TR5(sin(x)**6)
    sin(x)**6
    """
    # XXX should this do all even powers? all those that are powers of 2?
    rv = bottom_up(rv, TR5)
    if not (rv.is_Pow and rv.base.func == sin and rv.exp in (2, 4)):
        return rv
    return (1 - cos(rv.base.args[0])**2)**(rv.exp//2)


def TR6(rv):
    """Replacement of cos**2 or cos**4 factors.

    Examples
    ========

    >>> from sympy.simplify.fu import TR6
    >>> from sympy.abc import x
    >>> from sympy import cos
    >>> TR6(cos(x)**4)
    (-sin(x)**2 + 1)**2
    >>> TR6(cos(x)**2)
    -sin(x)**2 + 1
    >>> TR6(cos(x)**6)
    cos(x)**6
    """
    rv = bottom_up(rv, TR6)
    if not (rv.is_Pow and rv.base.func == cos and rv.exp in (2, 4)):
        return rv
    return (1 - sin(rv.base.args[0])**2)**(rv.exp//2)


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
    for a in ordered(Mul.make_args(rv)):
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

        def ok(a):
            # return ``s, f`` if ``a`` is ``s*f(a)`` where ``f`` is ``cos`` or ``sin``,
            # else return None.
            if not (a.func in (cos, sin) or (a.is_Mul and (-a).func in (cos, sin))):
                return
            n = 1
            if a.is_Mul:
                n = -1
                a = -a
            return n, a

        def ok2(rv):
            # see if a gcd can be pulled out of the two args and if so, return both ``ok`` values
            # or None if either ``ok`` fails
            a, b = [Factors(i) for i in rv.args]
            A = a
            a, b = a.normal(b)
            gcd = A.div(a)[0]
            if gcd.is_one:
                return  # there was no gcd
            o1 = ok(a.as_expr())
            if not o1:
                return
            o2 = ok(b.as_expr())
            if not o2:
                return
            return gcd.as_expr(), o1, o2

        if len(rv.args) != 2:
            args = list(ordered(rv.args))
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
        # arg 1
        gcd = S.One
        o1 = ok(rv.args[0])
        if o1:
            n1, f = o1
        # arg 2
        o2 = ok(rv.args[1])
        if o2:
            if not o1:
                return rv  # both must succeed
            n2, g = o2
        elif o1:
            return rv  # both must succeed
        else:  # both failed
            if first:
                o = ok2(rv)
                if not o:
                    return rv
            else:
                return rv
            gcd, o1, o2 = o
            n1, f = o1
            n2, g = o2

        # application of rule if possible
        if f.func == g.func:
            a, b = f.args[0], g.args[0]
            if f.func == cos:
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


def TR10(rv):
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
        args = list(ordered(arg.args))
        a = args.pop()
        b = Add(*args)
        if b.is_Add:
            if f == sin:
                return sin(a)*TR10(cos(b)) + cos(a)*TR10(sin(b))
            else:
                return cos(a)*TR10(cos(b)) - sin(a)*TR10(sin(b))
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
    >>> from sympy import cos, sin
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
    """
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

        def ok(a):
            # return ``s, f, g`` if ``a`` is ``s*f(a)*g(b)`` where ``f`` and
            # ``g`` are ``cos`` and/or ``sin`` else return None; ``c`` will be
            # cos unless both ``f`` and ``g`` are ``sin``.
            if not a.is_Mul:
                return
            if a.args[0].is_Integer:
                if a.args[0] is not S.NegativeOne:  # gcd already removed
                    return
                a = -a
                n = -1
            else:
                n = 1
            if not a.is_Mul or len(a.args) != 2:
                return
            c, s = a.args
            if c.func == sin:
                c, s = c, s
            if all(f.func in (cos, sin) for f in (c, s)):
                return n, c, s

        def ok2(rv):
            # see if a gcd can be pulled out of the two args and if so, return both ``ok`` values
            # or None if either ``ok`` fails
            a, b = [Factors(i) for i in rv.args]
            A = a
            a, b = a.normal(b)
            gcd = A.div(a)[0]
            if gcd.is_one:
                return  # there was no gcd
            o1 = ok(a.as_expr())
            if not o1:
                return
            o2 = ok(b.as_expr())
            if not o2:
                return
            return gcd.as_expr(), o1, o2

        if len(rv.args) != 2:
            args = list(ordered(rv.args))
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
        # arg 1
        gcd = S.One
        o1 = ok(rv.args[0])
        if o1:
            n1, c1, s1 = o1
        # arg 2
        o2 = ok(rv.args[1])
        if o2:
            if o1 is None:  # both have to succeed
                return rv
            n2, c2, s2 = o2
        elif o1:
            return rv  # both have to succeed
        else:  # both failed
            if first:
                o = ok2(rv)
                if not o:
                    return rv
            else:
                return rv
            gcd, o1, o2 = o
            n1, c1, s1 = o1
            n2, c2, s2 = o2

        # identify and get c1 to be cos then apply rule if possible
        if c1.func == s1.func and c2.func == s2.func and c1.func != c2.func:
            if c1.func == sin:
                c1, c2 = c2, c1
                s1, s2 = s2, s1
            a, b = c1.args[0], s1.args[0]
            aa, bb = c2.args[0], s2.args[0]
            if not (a == aa and b == bb):
                return rv
            if n1 == n2:
                return gcd*cos(a - b)
            return gcd*cos(a + b)
        elif c1.func != s1.func and c2.func != s2.func:
            if c1.func == sin:
                c1, s1 = s1, c1
            if c2.func == sin:
                c2, s2 = s2, c2
            a, b = c1.args[0], s1.args[0]
            aa, bb = c2.args[0], s2.args[0]
            if not (a == bb and b == aa):
                return rv
            if n1 == n2:
                return gcd*sin(a + b)
            return gcd*n1*sin(b - a)
        return rv

    return process_common_addends(rv, do, lambda x: tuple(ordered(x.free_symbols)))


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

_CTR4 = [TR4, TR0], [identity]


def CTRstrat(lists):
    return minimize(*[chain(*list) for list in lists], objective=L)

CTR1, CTR2, CTR3, CTR4 = map(CTRstrat, (_CTR1, _CTR2, _CTR3, _CTR4))

_RL1 = [TR4, TR3, TR4, TR12, TR4, TR13, TR4, TR0]

_RL2 = [TR4, TR3, TR10, TR4, TR3, TR11, TR5, TR7, TR11, TR4, CTR3, TR0, CTR1,
        TR9, CTR2, TR4, TR9, TR0, TR9, CTR4]


def RLstrat(rls):
    return chain(*rls)

RL1, RL2 = map(RLstrat, (_RL1, _RL2))


def fu(rv):
    """Attempt to simplify expression by using transformation rules given
    in the Fu et al algorithm.

    Examples
    ========

    >>> from sympy.simplify.fu import fu
    >>> from sympy import cos, sin, pi
    >>> a = sin(50)**2 + cos(50)**2 + sin(pi/6)
    >>> fu(a)
    3/2

    >>> a = sin(3)**4 - cos(2)**2 + sin(2)**2 + 2*cos(3)**2
    >>> fu(a)  # -8*cos(1)**4 + cos(3)**4 + 8*cos(1)**2
    cos(6)**2/4 + cos(6)/2 - cos(4) + 5/4

    References
    ==========
    http://rfdz.ph-noe.ac.at/fileadmin/Mathematik_Uploads/ACDCA/
    DESTIME2006/DES_contribs/Fu/simplification.pdf
    """
    rv = sympify(rv)
    rv = TR0(rv)
    rv = TR1(rv)
    if rv.has(tan, cot):
        rv1 = RL1(rv)
        if (L(rv1) < L(rv)):
            rv = rv1
    if rv.has(tan, cot):
        rv = TR2(rv)
    rv = TR0(rv)
    if rv.has(sin, cos):
        rv1 = RL2(rv)
        rv2 = TR8(rv1)
        rv = ordered([rv, rv1, rv2], keys=(L, count_ops), default=False).next()
    return rv


def bottom_up(rv, F):
    """Apply ``F`` to all expressions in an expression tree from the
    bottom up.
    """
    if rv.args:
        args = tuple([F(a) for a in rv.args])
        if args != rv.args:
            rv = rv.func(*args)
    return rv


def process_common_addends(rv, do, key2=None):
    """Apply ``do`` to addends of ``rv`` that share at least a common
    absolute value of their coefficient and the value of key2 when applied
    to the argument.
    """
    # collect by absolute value of coefficient and key2
    absc = defaultdict(list)
    for a in rv.args:
        c, a = a.as_coeff_Mul()
        if c < 0:
            c = -c
            a = -a  # put the sign on `a`
        absc[(c, key2(a) if key2 else 1)].append(a)

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
