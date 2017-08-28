"""
Limits
======

Implemented according to the PhD thesis
http://www.cybertester.com/data/gruntz.pdf, which contains very thorough
descriptions of the algorithm including many examples.  We summarize here
the gist of it.

All functions are sorted according to how rapidly varying they are at
infinity using the following rules. Any two functions f and g can be
compared using the properties of L:

L=lim  log|f(x)| / log|g(x)|           (for x -> oo)

We define >, < ~ according to::

    1. f > g .... L=+-oo

        we say that:
        - f is greater than any power of g
        - f is more rapidly varying than g
        - f goes to infinity/zero faster than g

    2. f < g .... L=0

        we say that:
        - f is lower than any power of g

    3. f ~ g .... L!=0, +-oo

        we say that:
        - both f and g are bounded from above and below by suitable integral
          powers of the other

Examples
========
::
    2 < x < exp(x) < exp(x**2) < exp(exp(x))
    2 ~ 3 ~ -5
    x ~ x**2 ~ x**3 ~ 1/x ~ x**m ~ -x
    exp(x) ~ exp(-x) ~ exp(2x) ~ exp(x)**2 ~ exp(x+exp(-x))
    f ~ 1/f

So we can divide all the functions into comparability classes (x and x^2
belong to one class, exp(x) and exp(-x) belong to some other class). In
principle, we could compare any two functions, but in our algorithm, we
don't compare anything below the class 2~3~-5 (for example log(x) is
below this), so we set 2~3~-5 as the lowest comparability class.

Given the function f, we find the list of most rapidly varying (mrv set)
subexpressions of it. This list belongs to the same comparability class.
Let's say it is {exp(x), exp(2x)}. Using the rule f ~ 1/f we find an
element "w" (either from the list or a new one) from the same
comparability class which goes to zero at infinity. In our example we
set w=exp(-x) (but we could also set w=exp(-2x) or w=exp(-3x) ...). We
rewrite the mrv set using w, in our case {1/w, 1/w^2}, and substitute it
into f. Then we expand f into a series in w::

    f = c0*w^e0 + c1*w^e1 + ... + O(w^en),       where e0<e1<...<en, c0!=0

but for x->oo, lim f = lim c0*w^e0, because all the other terms go to zero,
because w goes to zero faster than the ci and ei. So::

    for e0>0, lim f = 0
    for e0<0, lim f = +-oo   (the sign depends on the sign of c0)
    for e0=0, lim f = lim c0

We need to recursively compute limits at several places of the algorithm, but
as is shown in the PhD thesis, it always finishes.

Important functions from the implementation:

compare(a, b, x) compares "a" and "b" by computing the limit L.
mrv(e, x) returns list of most rapidly varying (mrv) subexpressions of "e"
rewrite(e, Omega, x, wsym) rewrites "e" in terms of w
leadterm(f, x) returns the lowest power term in the series of f
mrv_leadterm(e, x) returns the lead term (c0, e0) for e
limitinf(e, x) computes lim e  (for x->oo)
limit(e, z, z0) computes any limit by converting it to the case x->oo

All the functions are really simple and straightforward except
rewrite(), which is the most difficult/complex part of the algorithm.
When the algorithm fails, the bugs are usually in the series expansion
(i.e. in SymPy) or in rewrite.

This code is almost exact rewrite of the Maple code inside the Gruntz
thesis.

Debugging
---------

Because the gruntz algorithm is highly recursive, it's difficult to
figure out what went wrong inside a debugger. Instead, turn on nice
debug prints by defining the environment variable SYMPY_DEBUG. For
example:

[user@localhost]: SYMPY_DEBUG=True ./bin/isympy

In [1]: limit(sin(x)/x, x, 0)
limitinf(_x*sin(1/_x), _x) = 1
+-mrv_leadterm(_x*sin(1/_x), _x) = (1, 0)
| +-mrv(_x*sin(1/_x), _x) = set([_x])
| | +-mrv(_x, _x) = set([_x])
| | +-mrv(sin(1/_x), _x) = set([_x])
| |   +-mrv(1/_x, _x) = set([_x])
| |     +-mrv(_x, _x) = set([_x])
| +-mrv_leadterm(exp(_x)*sin(exp(-_x)), _x, set([exp(_x)])) = (1, 0)
|   +-rewrite(exp(_x)*sin(exp(-_x)), set([exp(_x)]), _x, _w) = (1/_w*sin(_w), -_x)
|     +-sign(_x, _x) = 1
|     +-mrv_leadterm(1, _x) = (1, 0)
+-sign(0, _x) = 0
+-limitinf(1, _x) = 1

And check manually which line is wrong. Then go to the source code and
debug this function to figure out the exact problem.

"""
from __future__ import print_function, division

from sympy.core import S, oo, Dummy, Mul, Add, evaluate
from sympy.core.compatibility import default_sort_key
from sympy.functions import log, exp, sign as sgn
from sympy.series.order import Order
from sympy.simplify.powsimp import powsimp, powdenest
from sympy.core.cache import cacheit

from sympy.core.compatibility import reduce

from sympy.utilities.timeutils import timethis
timeit = timethis('gruntz')

from sympy.utilities.misc import debug_decorator as debug


def compare(a, b, x):
    r"""Determine order relation between two functons.

    Returns
    =======

    {1, 0, -1}
        Respectively, if `a(x) \succ b(x)`, `a(x) \asymp b(x)`
        or `b(x) \succ a(x)`.
    """
    # The log(exp(...)) must always be simplified here for termination.
    la = a.exp if a.is_Pow and a.base is S.Exp1 else log(a)
    lb = b.exp if b.is_Pow and b.base is S.Exp1 else log(b)

    c = limitinf(la/lb, x)
    if c.is_zero:
        return -1
    elif c.is_infinite:
        return 1
    else:
        return 0


def mrv(e, x):
    """Calculate the MRV set of expression."""
    if not e.has(x):
        return set()
    elif e == x:
        return {x}
    elif e.is_Mul or e.is_Add:
        a, b = e.as_two_terms()
        return mrv_max(mrv(a, x), mrv(b, x), x)
    elif e.is_Pow:
        if e.base is S.Exp1:
            if e.exp == x:
                return {e}
            elif any(a.is_infinite for a in Mul.make_args(limitinf(e.exp, x))):
                return mrv_max({e}, mrv(e.exp, x), x)
            else:
                return mrv(e.exp, x)
        else:
            assert not e.exp.has(x)
            return mrv(e.base, x)
    elif e.func is log:
        return mrv(e.args[0], x)
    elif e.is_Function:
        return reduce(lambda a, b: mrv_max(a, b, x), [mrv(a, x) for a in e.args])
    raise NotImplementedError("Don't know how to calculate the mrv of '%s'" % e)


def mrv_max(f, g, x):
    """Computes the maximum of two MRV sets."""
    if not f:
        return g
    elif not g:
        return f
    elif f & g:
        return f | g

    c = compare(list(f)[0], list(g)[0], x)
    if c > 0:
        return f
    elif c < 0:
        return g
    else:
        return f | g


@cacheit
@timeit
def sign(e, x):
    r"""
    Determine a sign of an expression at infinity.

    Returns
    =======

    {1, 0, -1}
        One or minus one, if `e > 0` or `e < 0` for `x` sufficiently
        large and zero if `e` is *constantly* zero for `x\to\infty`.

    The result of this function is currently undefined if `e` changes sign
    arbitarily often at infinity (e.g. `sin(x)`).

    Note that this returns zero only if e is *constantly* zero
    for x sufficiently large. [If e is constant, of course, this is just
    the same thing as the sign of e.]
    """
    if not e.has(x):
        return sgn(e).simplify()
    elif e == x:
        return 1
    elif e.is_Mul:
        a, b = e.as_two_terms()
        return sign(a, x)*sign(b, x)
    elif e.is_Pow:
        s = sign(e.base, x)
        if s == 1:
            return 1

    # if all else fails, do it the hard way
    c0, e0 = mrv_leadterm(e, x)
    return sign(c0, x)


@debug
@timeit
@cacheit
def limitinf(e, x):
    """Compute limit of the expression at the infinity."""
    # TODO restore: assert x.is_real and x.is_positive

    # Rewrite e in terms of tractable functions only:
    e = e.rewrite('tractable', deep=True)

    if not e.has(x):
        # This is a bit of a heuristic for nice results.  We always rewrite
        # tractable functions in terms of familiar intractable ones.
        # TODO: It might be nicer to rewrite the exactly to what they were
        # initially, but that would take some work to implement.
        return e.rewrite('intractable', deep=True)

    c0, e0 = mrv_leadterm(e, x)
    sig = sign(e0, x)
    if sig == 1:
        return S.Zero
    elif sig == -1:
        s = sign(c0, x)
        assert s != S.Zero
        return s*oo
    elif sig == 0:
        return limitinf(c0, x)  # e0=0: lim f = lim c0


@cacheit
def mrv_leadterm(e, x):
    """Compute the leading term of the series.

    Returns
    =======

    tuple
        The leading term `c_0 w^{e_0}` of the series of `e` in terms
        of the most rapidly varying subexpression `w` in form of
        the pair ``(c0, e0)`` of Expr.
    """
    if not e.has(x):
        return (e, S.Zero)

    e = e.replace(lambda f: f.is_Pow and f.base != S.Exp1 and f.exp.has(x),
                  lambda f: exp(log(f.base)*f.exp))
    e = e.replace(lambda f: f.is_Mul and sum(a.is_Pow for a in f.args) > 1,
                  lambda f: Mul(exp(Add(*[a.exp for a in f.args if a.is_Pow and a.base is S.Exp1])),
                                *[a for a in f.args if not a.is_Pow or a.base is not S.Exp1]))

    # The positive dummy, w, is used here so log(w*2) etc. will expand.
    # TODO: For limits of complex functions, the algorithm would have to
    # be improved, or just find limits of Re and Im components separately.
    w = Dummy("w", real=True, positive=True)
    e, logw = rewrite(e, x, w)

    lt = e.compute_leading_term(w, logx=logw)
    return lt.as_coeff_exponent(w)


def rewrite(e, x, w):
    """Rewrites expression in terms of the most rapidly varying subexpression.

    Parameters
    ==========

    e : Expr
        an expression
    x : Symbol
        variable of the `e`
    w : Symbol
        The symbol which is going to be used for substitution in place
        of the most rapidly varying in `x` subexpression.

    Returns
    =======

    tuple
        A pair: rewritten (in `w`) expression and `\log(w)`.
    """
    Omega = mrv(e, x)
    if not Omega:
        return e, None  # e really does not depend on x

    assert all(e.has(t) for t in Omega)

    if x in Omega:
        # Moving up in the asymptotical scale (exponentiate e and Omega):
        with evaluate(False):
            e = e.xreplace({x: exp(x)})
            Omega = {s.xreplace({x: exp(x)}) for s in Omega}

    # Use default_sort_key as a last resort to get deterministic output.
    Omega = sorted(Omega, key=lambda a: (-len(mrv(a, x)), default_sort_key(a)))

    for g in Omega:
        sig = sign(g.exp, x)
        if sig not in (1, -1):
            raise NotImplementedError('Result depends on the sign of %s' % sig)

    if sig == 1:
        w = 1/w  # if g goes to oo, substitute 1/w

    # Rewrite and substitute subexpressions in the Omega.
    for a in Omega:
        c = limitinf(a.exp/g.exp, x)
        b = exp(a.exp - c*g.exp)*w**c  # exponential must never be expanded here
        with evaluate(False):
            e = e.xreplace({a: b})

    return e, -sig*g.exp


def gruntz(e, z, z0, dir="+"):
    """
    Compute the limit of e(z) at the point z0 using the Gruntz algorithm.

    z0 can be any expression, including oo and -oo.

    For dir="+" (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-). For infinite z0
    (oo or -oo), the dir argument doesn't matter.

    This algorithm is fully described in the module docstring in the gruntz.py
    file. It relies heavily on the series expansion. Most frequently, gruntz()
    is only used if the faster limit() function (which uses heuristics) fails.
    """
    if not z.is_Symbol:
        raise NotImplementedError("Second argument must be a Symbol")

    # convert all limits to the limit z->oo; sign of z is handled in limitinf
    if z0 == oo:
        r = limitinf(e, z)
    elif z0 == -oo:
        r = limitinf(e.subs(z, -z), z)
    else:
        if str(dir) == "-":
            e0 = e.subs(z, z0 - 1/z)
        elif str(dir) == "+":
            e0 = e.subs(z, z0 + 1/z)
        else:
            raise NotImplementedError("dir must be '+' or '-'")
        r = limitinf(e0, z)

    # This is a bit of a heuristic for nice results... we always rewrite
    # tractable functions in terms of familiar intractable ones.
    # It might be nicer to rewrite the exactly to what they were initially,
    # but that would take some work to implement.
    return r.rewrite('intractable', deep=True)
