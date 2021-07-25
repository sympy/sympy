from sympy.core import S, Symbol, Add, sympify, Expr, PoleError, Mul
from sympy.core.exprtools import factor_terms
from sympy.functions.combinatorial.factorials import factorial
from sympy.functions.special.gamma_functions import gamma
from sympy.polys import PolynomialError, factor
from sympy.series.order import Order
from sympy.simplify.ratsimp import ratsimp
from sympy.simplify.simplify import together
from .gruntz import gruntz

def limit(e, z, z0, dir="+"):
    """Computes the limit of ``e(z)`` at the point ``z0``.

    Parameters
    ==========

    e : expression, the limit of which is to be taken

    z : symbol representing the variable in the limit.
        Other symbols are treated as constants. Multivariate limits
        are not supported.

    z0 : the value toward which ``z`` tends. Can be any expression,
        including ``oo`` and ``-oo``.

    dir : string, optional (default: "+")
        The limit is bi-directional if ``dir="+-"``, from the right
        (z->z0+) if ``dir="+"``, and from the left (z->z0-) if
        ``dir="-"``. For infinite ``z0`` (``oo`` or ``-oo``), the ``dir``
        argument is determined from the direction of the infinity
        (i.e., ``dir="-"`` for ``oo``).

    Examples
    ========

    >>> from sympy import limit, sin, oo
    >>> from sympy.abc import x
    >>> limit(sin(x)/x, x, 0)
    1
    >>> limit(1/x, x, 0) # default dir='+'
    oo
    >>> limit(1/x, x, 0, dir="-")
    -oo
    >>> limit(1/x, x, 0, dir='+-')
    zoo
    >>> limit(1/x, x, oo)
    0

    Notes
    =====

    First we try some heuristics for easy and frequent cases like "x", "1/x",
    "x**2" and similar, so that it's fast. For all other cases, we use the
    Gruntz algorithm (see the gruntz() function).

    See Also
    ========

     limit_seq : returns the limit of a sequence.
    """

    return Limit(e, z, z0, dir).doit(deep=False)


def heuristics(e, z, z0, dir):
    """Computes the limit of an expression term-wise.
    Parameters are the same as for the ``limit`` function.
    Works with the arguments of expression ``e`` one by one, computing
    the limit of each and then combining the results. This approach
    works only for simple limits, but it is fast.
    """

    from sympy.calculus.util import AccumBounds
    rv = None
    if abs(z0) is S.Infinity:
        rv = limit(e.subs(z, 1/z), z, S.Zero, "+" if z0 is S.Infinity else "-")
        if isinstance(rv, Limit):
            return
    elif e.is_Mul or e.is_Add or e.is_Pow or e.is_Function:
        r = []
        for a in e.args:
            l = limit(a, z, z0, dir)
            if l.has(S.Infinity) and l.is_finite is None:
                if isinstance(e, Add):
                    m = factor_terms(e)
                    if not isinstance(m, Mul): # try together
                        m = together(m)
                    if not isinstance(m, Mul): # try factor if the previous methods failed
                        m = factor(e)
                    if isinstance(m, Mul):
                        return heuristics(m, z, z0, dir)
                    return
                return
            elif isinstance(l, Limit):
                return
            elif l is S.NaN:
                return
            else:
                r.append(l)
        if r:
            rv = e.func(*r)
            if rv is S.NaN and e.is_Mul and any(isinstance(rr, AccumBounds) for rr in r):
                r2 = []
                e2 = []
                for ii in range(len(r)):
                    if isinstance(r[ii], AccumBounds):
                        r2.append(r[ii])
                    else:
                        e2.append(e.args[ii])

                if len(e2) > 0:
                    e3 = Mul(*e2).simplify()
                    l = limit(e3, z, z0, dir)
                    rv = l * Mul(*r2)

            if rv is S.NaN:
                try:
                    rat_e = ratsimp(e)
                except PolynomialError:
                    return
                if rat_e is S.NaN or rat_e == e:
                    return
                return limit(rat_e, z, z0, dir)
    return rv


class Limit(Expr):
    """Represents an unevaluated limit.

    Examples
    ========

    >>> from sympy import Limit, sin
    >>> from sympy.abc import x
    >>> Limit(sin(x)/x, x, 0)
    Limit(sin(x)/x, x, 0)
    >>> Limit(1/x, x, 0, dir="-")
    Limit(1/x, x, 0, dir='-')

    """

    def __new__(cls, e, z, z0, dir="+"):
        e = sympify(e)
        z = sympify(z)
        z0 = sympify(z0)

        if z0 is S.Infinity:
            dir = "-"
        elif z0 is S.NegativeInfinity:
            dir = "+"

        if(z0.has(z)):
            raise NotImplementedError("Limits approaching a variable point are"
                    " not supported (%s -> %s)" % (z, z0))
        if isinstance(dir, str):
            dir = Symbol(dir)
        elif not isinstance(dir, Symbol):
            raise TypeError("direction must be of type basestring or "
                    "Symbol, not %s" % type(dir))
        if str(dir) not in ('+', '-', '+-'):
            raise ValueError("direction must be one of '+', '-' "
                    "or '+-', not %s" % dir)

        obj = Expr.__new__(cls)
        obj._args = (e, z, z0, dir)
        return obj


    @property
    def free_symbols(self):
        e = self.args[0]
        isyms = e.free_symbols
        isyms.difference_update(self.args[1].free_symbols)
        isyms.update(self.args[2].free_symbols)
        return isyms


    def pow_heuristics(self):
        from sympy import exp, log
        expr, z, z0, _ = self.args
        b, e = expr.base, expr.exp
        if not b.has(z):
            res = limit(e*log(b), z, z0)
            return exp(res)

        ex_lim = limit(e, z, z0)
        base_lim = limit(b, z, z0)

        if base_lim is S.One:
            if ex_lim in (S.Infinity, S.NegativeInfinity):
                res = limit(e*(b - 1), z, z0)
                return exp(res)
        if base_lim is S.NegativeInfinity and ex_lim is S.Infinity:
            return S.ComplexInfinity


    def doit(self, **hints):
        """Evaluates the limit.

        Parameters
        ==========

        deep : bool, optional (default: True)
            Invoke the ``doit`` method of the expressions involved before
            taking the limit.

        hints : optional keyword arguments
            To be passed to ``doit`` methods; only used if deep is True.
        """
        from sympy import Abs, sign

        e, z, z0, dir = self.args

        if z0 is S.ComplexInfinity:
            raise NotImplementedError("Limits at complex "
                                    "infinity are not implemented")

        if hints.get('deep', True):
            e = e.doit(**hints)
            z = z.doit(**hints)
            z0 = z0.doit(**hints)

        if e == z:
            return z0

        if not e.has(z):
            return e

        if z0 is S.NaN:
            return S.NaN

        if e.has(S.Infinity, S.NegativeInfinity, S.ComplexInfinity, S.NaN):
            return self

        if e.is_Order:
            return Order(limit(e.expr, z, z0), *e.args[1:])

        cdir = 0
        if str(dir) == "+":
            cdir = 1
        elif str(dir) == "-":
            cdir = -1

        def set_signs(expr):
            if not expr.args:
                return expr
            newargs = tuple(set_signs(arg) for arg in expr.args)
            if newargs != expr.args:
                expr = expr.func(*newargs)
            abs_flag = isinstance(expr, Abs)
            sign_flag = isinstance(expr, sign)
            if abs_flag or sign_flag:
                sig = limit(expr.args[0], z, z0, dir)
                if sig.is_zero:
                    sig = limit(1/expr.args[0], z, z0, dir)
                if sig.is_extended_real:
                    if (sig < 0) == True:
                        return -expr.args[0] if abs_flag else S.NegativeOne
                    elif (sig > 0) == True:
                        return expr.args[0] if abs_flag else S.One
            return expr

        e = set_signs(e)

        if e.is_meromorphic(z, z0):
            if abs(z0) is S.Infinity:
                newe = e.subs(z, -1/z)
            else:
                newe = e.subs(z, z + z0)
            try:
                coeff, ex = newe.leadterm(z, cdir=cdir)
            except ValueError:
                pass
            else:
                if ex > 0:
                    return S.Zero
                elif ex == 0:
                    return coeff
                if str(dir) == "+" or not(int(ex) & 1):
                    return S.Infinity*sign(coeff)
                elif str(dir) == "-":
                    return S.NegativeInfinity*sign(coeff)
                else:
                    return S.ComplexInfinity

        if abs(z0) is S.Infinity:
            if e.is_Mul:
                e = factor_terms(e)
            newe = e.subs(z, 1/z)
            # cdir changes sign as oo- should become 0+
            cdir = -cdir
        else:
            newe = e.subs(z, z + z0)
        try:
            coeff, ex = newe.leadterm(z, cdir=cdir)
        except (ValueError, NotImplementedError, PoleError):
            # The NotImplementedError catching is for custom functions
            if e.is_Pow:
                r = self.pow_heuristics()
                if r is not None:
                    return r
        else:
            if coeff.has(S.Infinity, S.NegativeInfinity, S.ComplexInfinity):
                return self
            if not coeff.has(z):
                if ex.is_positive:
                    return S.Zero
                elif ex == 0:
                    return coeff
                elif ex.is_negative:
                    if ex.is_integer:
                        if str(dir) == "-" or ex.is_even:
                            return S.Infinity*sign(coeff)
                        elif str(dir) == "+":
                            return S.NegativeInfinity*sign(coeff)
                        else:
                            return S.ComplexInfinity
                    else:
                        if str(dir) == "+":
                            return S.Infinity*sign(coeff)
                        elif str(dir) == "-":
                            return S.NegativeInfinity*sign(coeff)*S.NegativeOne**(S.One + ex)
                        else:
                            return S.ComplexInfinity

        # gruntz fails on factorials but works with the gamma function
        # If no factorial term is present, e should remain unchanged.
        # factorial is defined to be zero for negative inputs (which
        # differs from gamma) so only rewrite for positive z0.
        if z0.is_extended_positive:
            e = e.rewrite(factorial, gamma)

        l = None

        try:
            if str(dir) == '+-':
                r = gruntz(e, z, z0, '+')
                l = gruntz(e, z, z0, '-')
                if l != r:
                    raise ValueError("The limit does not exist since "
                            "left hand limit = %s and right hand limit = %s"
                            % (l, r))
            else:
                r = gruntz(e, z, z0, dir)
            if r is S.NaN or l is S.NaN:
                raise PoleError()
        except (PoleError, ValueError):
            if l is not None:
                raise
            r = heuristics(e, z, z0, dir)
            if r is None:
                return self

        return r
