from sympy.core import (Basic, Expr, S, C, Symbol, Wild, Add, sympify, diff,
                        oo, Tuple, Interval)

from sympy.core.symbol import Dummy
from sympy.core.compatibility import is_sequence
from sympy.integrals.trigonometry import trigintegrate
from sympy.integrals.deltafunctions import deltaintegrate
from sympy.integrals.rationaltools import ratint
from sympy.integrals.risch import heurisch
from sympy.integrals.meijerint import meijerint_definite, meijerint_indefinite
from sympy.utilities import xthreaded, flatten
from sympy.utilities.misc import filldedent
from sympy.polys import Poly, PolynomialError
from sympy.solvers.solvers import solve, posify
from sympy.functions import Piecewise, sqrt, sign
from sympy.geometry import Curve
from sympy.functions.elementary.piecewise import piecewise_fold
from sympy.series import limit

def _process_limits(*symbols):
    """Convert the symbols-related limits into proper limits,
    storing them as Tuple(symbol, lower, upper). The sign of
    the function is also returned when the upper limit is missing
    so (x, 1, None) becomes (x, None, 1) and the sign is changed.
    """
    limits = []
    sign = 1
    for V in symbols:
        if isinstance(V, Symbol):
            limits.append(Tuple(V))
            continue
        elif is_sequence(V, Tuple):
            V = sympify(flatten(V))
            if V[0].is_Symbol:
                newsymbol = V[0]
                if len(V) == 2 and isinstance(V[1], Interval):
                    V[1:] = [V[1].start, V[1].end]

                if len(V) == 3:
                    if V[1] is None and V[2] is not None:
                        nlim = [V[2]]
                    elif V[1] is not None and V[2] is None:
                        sign *= -1
                        nlim = [V[1]]
                    elif V[1] is None and V[2] is None:
                        nlim = []
                    else:
                        nlim = V[1:]
                    limits.append(Tuple(newsymbol, *nlim ))
                    continue
                elif len(V) == 1 or (len(V) == 2 and V[1] is None):
                    limits.append(Tuple(newsymbol))
                    continue
                elif len(V) == 2:
                    limits.append(Tuple(newsymbol, V[1]))
                    continue

        raise ValueError('Invalid limits given: %s' % str(symbols))

    return limits, sign

class Integral(Expr):
    """Represents unevaluated integral."""

    __slots__ = ['is_commutative']

    def __new__(cls, function, *symbols, **assumptions):
        """Create an unevaluated integral.

        Arguments are an integrand followed by one or more limits.

        If no limits are given and there is only one free symbol in the
        expression, that symbol will be used, otherwise an error will be
        raised.

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> Integral(x)
        Integral(x, x)
        >>> Integral(y)
        Integral(y, y)

        When limits are provided, they are interpreted as follows (using
        ``x`` as though it were the variable of integration):

            (x,) or x - indefinite integral
            (x, a) - "evaluate at" integral
            (x, a, b) - definite integral

        Although the same integral will be obtained from an indefinite
        integral and an "evaluate at" integral when ``a == x``, they
        respond differently to substitution:

        >>> i = Integral(x, x)
        >>> at = Integral(x, (x, x))
        >>> i.doit() == at.doit()
        True
        >>> i.subs(x, 1)
        Integral(1, x)
        >>> at.subs(x, 1)
        Integral(x, (x, 1))

        The ``as_dummy`` method can be used to see which symbols cannot be
        targeted by subs: those with a preppended underscore cannot be
        changed with ``subs``. (Also, the integration variables themselves --
        the first element of a limit -- can never be changed by subs.)

        >>> i.as_dummy()
        Integral(x, x)
        >>> at.as_dummy()
        Integral(_x, (_x, x))

        """

        # Any embedded piecewise functions need to be brought out to the
        # top level so that integration can go into piecewise mode at the
        # earliest possible moment.
        function = piecewise_fold(sympify(function))

        if function is S.NaN:
            return S.NaN

        if symbols:
            limits, sign = _process_limits(*symbols)
        else:
            # no symbols provided -- let's compute full anti-derivative
            free = function.free_symbols
            if len(free) != 1:
                raise ValueError("specify variables of integration for %s" % function)
            limits, sign = [Tuple(s) for s in free], 1

        while isinstance(function, Integral):
            # denest the integrand
            limits = list(function.limits) + limits
            function = function.function

        obj = Expr.__new__(cls, **assumptions)
        arglist = [sign*function]
        arglist.extend(limits)
        obj._args = tuple(arglist)
        obj.is_commutative = function.is_commutative # limits already checked

        return obj

    def __getnewargs__(self):
        return (self.function,) + tuple([tuple(xab) for xab in self.limits])

    @property
    def function(self):
        """Return the function to be integrated.

        Examples
        ========

        >>> from sympy import Integral
        >>> from sympy.abc import x
        >>> Integral(x**2, (x,)).function
        x**2

        See Also
        ========

        limits, variables, free_symbols
        """
        return self._args[0]

    @property
    def limits(self):
        """Return the limits of integration.

        Examples
        ========

        >>> from sympy import Integral
        >>> from sympy.abc import x, i
        >>> Integral(x**i, (i, 1, 3)).limits
        ((i, 1, 3),)

        See Also
        ========

        function, variables, free_symbols
        """
        return self._args[1:]

    @property
    def variables(self):
        """Return a list of the integration variables.

        Examples
        ========

        >>> from sympy import Integral
        >>> from sympy.abc import x, i
        >>> Integral(x**i, (i, 1, 3)).variables
        [i]

        See Also
        ========

        function, limits, free_symbols
        as_dummy : Replace integration variables with dummy ones
        transform : Perform mapping on the integration variable
        """
        return [l[0] for l in self.limits]

    @property
    def free_symbols(self):
        """
        This method returns the symbols that will exist when the
        integral is evaluated. This is useful if one is trying to
        determine whether an integral depends on a certain
        symbol or not.

        Examples
        ========

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> Integral(x, (x, y, 1)).free_symbols
        set([y])

        See Also
        ========

        function, limits, variables
        """
        function, limits = self.function, self.limits
        if function.is_zero:
            return set()
        isyms = function.free_symbols
        for xab in limits:
            if len(xab) == 1:
                isyms.add(xab[0])
                continue
            # take out the target symbol
            if xab[0] in isyms:
                isyms.remove(xab[0])
            if len(xab) == 3 and xab[1] == xab[2]:
                # if two limits are the same the integral is 0
                # and there are no symbols
                return set()
            # add in the new symbols
            for i in xab[1:]:
                isyms.update(i.free_symbols)
        return isyms

    @property
    def is_zero(self):
        """Since Integral doesn't autosimplify it it useful to see if
        it would simplify to zero or not in a trivial manner, i.e. when
        the function is 0 or two limits of a definite integral are the same.

        This is a very naive and quick test, not intended to check for special
        patterns like Integral(sin(m*x)*cos(n*x), (x, 0, 2*pi)) == 0.

        Examples
        ========

        >>> from sympy import Integral
        >>> from sympy.abc import x, y, z
        >>> Integral(1, (x, 1, 1)).is_zero
        True
        >>> Integral(0, (x, y, z)).is_zero
        True
        >>> Integral(1, (x, 1, 2)).is_zero
        False

        See Also
        ========

        is_number
        """
        if (self.function.is_zero or
            any(len(xab) == 3 and xab[1] == xab[2] for xab in self.limits)):
            return True
        if not self.free_symbols and self.function.is_number:
            # the integrand is a number and the limits are numerical
            return False

    @property
    def is_number(self):
        """
        Return True if the Integral will result in a number, else False.

        sympy considers anything that will result in a number to have
        is_number == True.

        >>> from sympy import log
        >>> log(2).is_number
        True

        Integrals are a special case since they contain symbols that can
        be replaced with numbers. Whether the integral can be done or not is
        another issue. But answering whether the final result is a number is
        not difficult.

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> Integral(x).is_number
        False
        >>> Integral(x, y).is_number
        False
        >>> Integral(x, (y, 1, x)).is_number
        False
        >>> Integral(x, (y, 1, 2)).is_number
        False
        >>> Integral(x, (y, 1, 1)).is_number
        True
        >>> Integral(x, (x, 1, 2)).is_number
        True
        >>> Integral(x*y, (x, 1, 2), (y, 1, 3)).is_number
        True
        >>> Integral(1, x, (x, 1, 2)).is_number
        True

        See Also
        ========

        is_zero
        """

        integrand, limits = self.function, self.limits
        isyms = integrand.atoms(Symbol)
        for xab in limits:
            if len(xab) == 1:
                isyms.add(xab[0])
                continue # it may be removed later
            elif len(xab) == 3 and xab[1] == xab[2]: # XXX naive equality test
                return True # integral collapsed
            if xab[0] in isyms:
                # take it out of the symbols since it will be replace
                # with whatever the limits of the integral are
                isyms.remove(xab[0])
            # add in the new symbols
            for i in xab[1:]:
                isyms.update(i.free_symbols)
        # if there are no surviving symbols then the result is a number
        return len(isyms) == 0

    def as_dummy(self):
        """
        Replace instances of the integration variables with their dummy
        counterparts to make clear what are dummy variables and what
        are real-world symbols in an Integral.

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> Integral(x, (x, x, y), (y, x, y)).as_dummy()
        Integral(_x, (_x, x, _y), (_y, x, y))

        The "integral at" limit that has a length of 1 is not treated as
        though the integration symbol is a dummy, but the explicit form
        of length 2 does treat the integration variable as a dummy.

        >>> Integral(x, x).as_dummy()
        Integral(x, x)
        >>> Integral(x, (x, x)).as_dummy()
        Integral(_x, (_x, x))

        If there were no dummies in the original expression, then the
        output of this function will show which symbols cannot be
        changed by subs(), those with an underscore prefix.

        See Also
        ========

        variables : Lists the integration variables
        transform : Perform mapping on the integration variable
        """
        reps = {}
        f = self.function
        limits = list(self.limits)
        for i in xrange(-1, -len(limits) - 1, -1):
            xab = list(limits[i])
            if len(xab) == 1:
                continue
            x = xab[0]
            xab[0] = x.as_dummy()
            for j in range(1, len(xab)):
                xab[j] = xab[j].subs(reps)
            reps[x] = xab[0]
            limits[i] = xab
        f = f.subs(reps)
        return Integral(f, *limits)

    def transform(self, x, u, inverse=False):
        r"""
        Performs a change of variables from `x` to `u` using the relationship
        given by `x` and `u` which will define the transformations `f` and `F`
        (which are inverses of each other) as follows:

        1) If `x` is a Symbol (which is a variable of integration) then `u`
           will be interpreted as some function, f(u), with inverse F(u).
           This, in effect, just makes the substitution of x with f(x).

        2) If `u` is a Symbol then `x` will be interpreted as some function,
           F(x), with inverse f(u). This is commonly referred to as
           u-substitution.

        The `inverse` option will reverse `x` and `u`. It is a deprecated option
        since `x` and `u` can just be passed in reverse order.

        Once f and F have been identified, the transformation is made as
        follows:

        .. math:: \int_a^b x \mathrm{d}x \rightarrow \int_{F(a)}^{F(b)} f(x)
                  \frac{\mathrm{d}}{\mathrm{d}x}

        where `F(x)` is the inverse of `f(x)` and the limits and integrand have
        been corrected so as to retain the same value after integration.

        Notes
        =====

        The mappings, F(x) or f(u), must lead to a unique integral. Linear
        or rational linear expression, `2*x`, `1/x` and `sqrt(x)`, will
        always work; quadratic expressions like `x**2 - 1` are acceptable
        as long as the resulting integrand does not depend on the sign of
        the solutions (see examples).

        The integral will be returned unchanged if `x` is not a variable of
        integration.

        `x` must be (or contain) only one of of the integration variables. If
        `u` has more than one free symbol then it should be sent as a tuple
        (`u`, `uvar`) where `uvar` identifies which variable is replacing
        the integration variable.
        XXX can it contain another integration variable?

        Examples
        ========

        >>> from sympy.abc import a, b, c, d, x, u, y
        >>> from sympy import Integral, S, cos, sqrt

        >>> i = Integral(x*cos(x**2 - 1), (x, 0, 1))

        transform can change the variable of integration

        >>> i.transform(x, u)
        Integral(u*cos(u**2 - 1), (u, 0, 1))

        transform can perform u-substitution as long as a unique
        integrand is obtained:

        >>> i.transform(x**2 - 1, u)
        Integral(cos(u)/2, (u, -1, 0))

        This attempt fails because x = +/-sqrt(u + 1) and the
        sign does not cancel out of the integrand:

        >>> Integral(cos(x**2 - 1), (x, 0, 1)).transform(x**2 - 1, u)
        Traceback (most recent call last):
        ...
        ValueError:
        The mapping between F(x) and f(u) did not give a unique integrand.

        transform can do a substitution. Here, the previous
        result is transformed back into the original expression
        using "u-substitution":

        >>> ui = _
        >>> _.transform(sqrt(u + 1), x) == i
        True

        We can accomplish the same with a regular substitution:

        >>> ui.transform(u, x**2 - 1) == i
        True

        If the `x` does not contain a symbol of integration then
        the integral will be returned unchanged. Integral `i` does
        not have an integration variable `a` so no change is made:

        >>> i.transform(a, x) == i
        True

        When `u` has more than one free symbol the symbol that is
        replacing `x` must be identified by passing `u` as a tuple:

        >>> Integral(x, (x, 0, 1)).transform(x, (u + a, u))
        Integral(a + u, (u, -a, -a + 1))
        >>> Integral(x, (x, 0, 1)).transform(x, (u + a, a))
        Integral(a + u, (a, -u, -u + 1))

        See Also
        ========

        variables : Lists the integration variables
        as_dummy : Replace integration variables with dummy ones
        """

        if inverse:
            # when this is removed, update the docstring
            from sympy.utilities.exceptions import SymPyDeprecationWarning
            SymPyDeprecationWarning(
            feature="transform(x, f(x), inverse=True)",
            useinstead="transform(f(x), x)",
            issue=3380, deprecated_since_version="0.7.2",
            ).warn()
            # in the old style x and u contained the same variable so
            # don't worry about using the old-style feature with the
            # new style input...but it will still work:
            # i.transform(x, u).transform(x, u, inverse=True) -> i
            x, u = u, x

        d = Dummy('d')

        xfree = x.free_symbols.intersection(self.variables)
        if len(xfree) > 1:
            raise ValueError('F(x) can only contain one of: %s' % self.variables)
        xvar = xfree.pop() if xfree else d

        if xvar not in self.variables:
            return self

        u = sympify(u)
        if isinstance(u, Expr):
            ufree = u.free_symbols
            if len(ufree) != 1:
                raise ValueError(filldedent('''
                When f(u) has more than one free symbol, the one replacing x
                must be identified: pass f(u) as (f(u), u)'''))
            uvar = ufree.pop()
        else:
            u, uvar = u
            if uvar not in u.free_symbols:
                raise ValueError(filldedent('''
                Expecting a tuple (expr, symbol) where symbol identified
                a free symbol in expr, but symbol is not in expr's free
                symbols.'''))
            if not isinstance(uvar, Symbol):
                raise ValueError(filldedent('''
                Expecting a tuple (expr, symbol) but didn't get
                a symbol; got %s''' % uvar))

        if x.is_Symbol and u.is_Symbol:
            return self.xreplace({x: u})

        if not x.is_Symbol and not u.is_Symbol:
            raise ValueError('either x or u must be a symbol')

        if uvar == xvar:
            return self.transform(x, u.subs(uvar, d)).xreplace({d: uvar})

        if uvar in self.limits:
            raise ValueError(filldedent('''
            u must contain the same variable as in x
            or a variable that is not already an integration variable'''))

        if not x.is_Symbol:
            F = [x.subs(xvar, d)]
            soln = solve(u - x, xvar, check=False)
            if not soln:
                raise ValueError('no solution for solve(F(x) - f(u), x)')
            f = [fi.subs(uvar, d) for fi in soln]
        else:
            f = [u.subs(uvar, d)]
            pdiff, reps = posify(u - x)
            puvar = uvar.subs([(v, k) for k, v in reps.iteritems()])
            soln = [s.subs(reps) for s in solve(pdiff, puvar)]
            if not soln:
                raise ValueError('no solution for solve(F(x) - f(u), u)')
            F = [fi.subs(xvar, d) for fi in soln]

        newfuncs = set([(self.function.subs(xvar, fi)*fi.diff(d)
                        ).subs(d, uvar) for fi in f])
        if len(newfuncs) > 1:
            raise ValueError(filldedent('''
            The mapping between F(x) and f(u) did not give
            a unique integrand.'''))
        newfunc = newfuncs.pop()

        def _calc_limit_1(F, a, b):
            """
            replace d with a, using subs if possible, otherwise limit
            where sign of b is considered
            """
            wok = F.subs(d, a)
            if wok is S.NaN or wok.is_bounded is False and a.is_bounded:
                return limit(sign(b)*F, d, a)
            return wok

        def _calc_limit(a, b):
            """
            replace d with a, using subs if possible, otherwise limit
            where sign of b is considered
            """
            avals = list(set([_calc_limit_1(Fi, a, b) for Fi in F]))
            if len(avals) > 1:
                raise ValueError(filldedent('''
                The mapping between F(x) and f(u) did not
                give a unique limit.'''))
            return avals[0]

        newlimits = []
        for xab in self.limits:
            sym = xab[0]
            if sym == xvar:
                if len(xab) == 3:
                    a, b = xab[1:]
                    a, b = _calc_limit(a, b), _calc_limit(b, a)
                    if a > b:
                        a, b = b, a
                        newfunc = -newfunc
                    newlimits.append((uvar, a, b))
                elif len(xab) == 2:
                    a = _calc_limit(xab[1], 1)
                    newlimits.append((uvar, a))
                else:
                    newlimits.append(uvar)
            else:
                newlimits.append(xab)

        return Integral(newfunc, *newlimits)

    def doit(self, **hints):
        """
        Perform the integration using any hints given.

        Examples
        ========

        >>> from sympy import Integral
        >>> from sympy.abc import x, i
        >>> Integral(x**i, (i, 1, 3)).doit()
        x**3/log(x) - x/log(x)

        See Also
        ========

        sympy.integrals.trigonometry.trigintegrate
        sympy.integrals.risch.heurisch
        sympy.integrals.rationaltools.ratint
        as_sum : Approximate the integral using a sum
        """
        if not hints.get('integrals', True):
            return self

        deep = hints.get('deep', True)
        meijerg = hints.get('meijerg', None)
        conds = hints.get('conds', 'piecewise')
        if conds not in ['separate', 'piecewise', 'none']:
            raise ValueError('conds must be one of "separate", "piecewise", ' \
                             '"none", got: %s' % conds)

        # check for the trivial case of equal upper and lower limits
        if self.is_zero:
            return S.Zero

        # now compute and check the function
        function = self.function
        if deep:
            function = function.doit(**hints)

        if function.is_zero:
            return S.Zero

        # There is no trivial answer, so continue

        undone_limits = []
        ulj = set() # free symbols of any undone limits' upper and lower limits
        for xab in self.limits:
            # compute uli, the free symbols in the
            # Upper and Lower limits of limit I
            if len(xab) == 1:
                uli = set(xab[:1])
            elif len(xab) == 2:
                uli = xab[1].free_symbols
            elif len(xab) == 3:
                uli = xab[1].free_symbols.union(xab[2].free_symbols)
            # this integral can be done as long as there is no blocking
            # limit that has been undone. An undone limit is blocking if
            # it contains an integration variable that is in this limit's
            # upper or lower free symbols or vice versa
            if xab[0] in ulj or any(v[0] in uli for v in undone_limits):
                undone_limits.append(xab)
                ulj.update(uli)
                continue

            # There are a number of tradeoffs in using the meijer g method.
            # It can sometimes be a lot faster than other methods, and
            # sometimes slower. And there are certain types of integrals for
            # which it is more likely to work than others.
            # These heuristics are incorporated in deciding what integration
            # methods to try, in what order.
            # See the integrate() docstring for details.
            def try_meijerg(function, xab):
                ret = None
                if len(xab) == 3 and meijerg is not False:
                    x, a, b = xab
                    try:
                        res = meijerint_definite(function, x, a, b)
                    except NotImplementedError:
                        from sympy.integrals.meijerint import _debug
                        _debug('NotImplementedError from meijerint_definite')
                        res = None
                    if res is not None:
                        f, cond = res
                        if conds == 'piecewise':
                            ret = Piecewise((f, cond),
                                          (Integral(function, (x, a, b)), True))
                        elif conds == 'separate':
                            if len(self.limits) != 1:
                                raise ValueError('conds=separate not supported in ' \
                                                 'multiple integrals')
                            ret = f, cond
                        else:
                            ret = f
                return ret

            meijerg1 = meijerg
            if len(xab) == 3 and xab[1].is_real and xab[2].is_real \
               and not function.is_Poly and \
               (xab[1].has(oo, -oo) or xab[2].has(oo, -oo)):
                ret = try_meijerg(function, xab)
                if ret is not None:
                    function = ret
                    continue
                else:
                    meijerg1 = False

            # If the special meijerg code did not succeed finding a definite
            # integral, then the code using meijerint_indefinite will not either
            # (it might find an antiderivative, but the answer is likely to be
            #  nonsensical).
            # Thus if we are requested to only use meijer g-function methods,
            # we give up at this stage. Otherwise we just disable g-function
            # methods.
            if meijerg1 is False and meijerg is True:
                antideriv = None
            else:
                antideriv = self._eval_integral(function, xab[0], meijerg1)
                if antideriv is None and meijerg1 is True:
                    ret = try_meijerg(function, xab)
                    if ret is not None:
                        function = ret
                        continue

            if antideriv is None:
                undone_limits.append(xab)
            else:
                if len(xab) == 1:
                    function = antideriv
                else:
                    if len(xab) == 3:
                        x, a, b = xab
                    if len(xab) == 2:
                        x, b = xab
                        a = None

                    if deep:
                        if isinstance(a, Basic):
                            a = a.doit(**hints)
                        if isinstance(b, Basic):
                            b = b.doit(**hints)

                    if antideriv.is_Poly:
                        gens = list(antideriv.gens)
                        gens.remove(x)

                        antideriv = antideriv.as_expr()

                        function = antideriv._eval_interval(x, a, b)
                        function = Poly(function, *gens)
                    else:
                        try:
                            function = antideriv._eval_interval(x, a, b)
                        except NotImplementedError:
                            # This can happen if _eval_interval depends in a
                            # complicated way on limits that cannot be computed
                            undone_limits.append(xab)

        if undone_limits:
            return self.func(*([function] + undone_limits))
        return function

    def _eval_derivative(self, sym):
        """Evaluate the derivative of the current Integral object by
        differentiating under the integral sign [1], using the Fundamental
        Theorem of Calculus [2] when possible.

        Whenever an Integral is encountered that is equivalent to zero or
        has an integrand that is independent of the variable of integration
        those integrals are performed. All others are returned as Integral
        instances which can be resolved with doit() (provided they are integrable).

        References:
           [1] http://en.wikipedia.org/wiki/Differentiation_under_the_integral_sign
           [2] http://en.wikipedia.org/wiki/Fundamental_theorem_of_calculus

        Examples
        ========

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> i = Integral(x + y, y, (y, 1, x))
        >>> i.diff(x)
        Integral(x + y, (y, x)) + Integral(1, y, (y, 1, x))
        >>> i.doit().diff(x) == i.diff(x).doit()
        True
        >>> i.diff(y)
        0

        The previous must be true since there is no y in the evaluated integral:
        >>> i.free_symbols
        set([x])
        >>> i.doit()
        2*x**3/3 - x/2 - 1/6

        """

        # differentiate under the integral sign; we do not
        # check for regularity conditions (TODO), see issue 1116

        # get limits and the function
        f, limits = self.function, list(self.limits)

        # the order matters if variables of integration appear in the limits
        # so work our way in from the outside to the inside.
        limit = limits.pop(-1)
        if len(limit) == 3:
            x, a, b = limit
        elif len(limit) == 2:
            x, b = limit
            a = None
        else:
            a = b = None
            x = limit[0]

        if limits: # f is the argument to an integral
            f = Integral(f, *tuple(limits))

        # assemble the pieces
        def _do(f, ab):
            dab_dsym = diff(ab, sym)
            if not dab_dsym:
                return S.Zero
            if isinstance(f, Integral):
                limits = [(x, x) if (len(l) == 1 and l[0] == x) else l
                          for l in f.limits]
                f = Integral(f.function, *limits)
            return f.subs(x, ab)*dab_dsym
        rv = 0
        if b is not None:
            rv += _do(f, b)
        if a is not None:
            rv -= _do(f, a)
        if len(limit) == 1 and sym == x:
            # the dummy variable *is* also the real-world variable
            arg = f
            rv += arg
        else:
            # the dummy variable might match sym but it's
            # only a dummy and the actual variable is determined
            # by the limits, so mask off the variable of integration
            # while differentiating
            u = Dummy('u')
            arg = f.subs(x, u).diff(sym).subs(u, x)
            rv += Integral(arg, Tuple(x, a, b))
        return rv

    def _eval_integral(self, f, x, meijerg=None):
        """Calculate the anti-derivative to the function f(x).

        This is a powerful function that should in theory be able to integrate
        everything that can be integrated. If you find something, that it
        doesn't, it is easy to implement it.

        (1) Simple heuristics (based on pattern matching and integral table):

         - most frequently used functions (e.g. polynomials)
         - functions non-integrable by any of the following algorithms (e.g.
           exp(-x**2))

        (2) Integration of rational functions:

         (a) using apart() - apart() is full partial fraction decomposition
             procedure based on Bronstein-Salvy algorithm. It gives formal
             decomposition with no polynomial factorization at all (so it's
             fast and gives the most general results). However it needs an
             implementation of the RootsOf class.
         (b) using Trager's algorithm - possibly faster than (a) but needs
             implementation :)

        (3) Whichever implementation of pmInt (Mateusz, Kirill's or a
            combination of both).

          - this way we can handle efficiently huge class of elementary and
            special functions

        (4) Recursive Risch algorithm as described in Bronstein's integration
            tutorial.

          - this way we can handle those integrable functions for which (3)
            fails

        (5) Powerful heuristics based mostly on user defined rules.

         - handle complicated, rarely used cases
        """

        # if it is a poly(x) then let the polynomial integrate itself (fast)
        #
        # It is important to make this check first, otherwise the other code
        # will return a sympy expression instead of a Polynomial.
        #
        # see Polynomial for details.
        if isinstance(f, Poly) and not meijerg:
            return f.integrate(x)

        # Piecewise antiderivatives need to call special integrate.
        if f.func is Piecewise:
            return f._eval_integral(x)

        # let's cut it short if `f` does not depend on `x`
        if not f.has(x):
            return f*x

        # try to convert to poly(x) and then integrate if successful (fast)
        poly = f.as_poly(x)

        if poly is not None and not meijerg:
            return poly.integrate().as_expr()

        # since Integral(f=g1+g2+...) == Integral(g1) + Integral(g2) + ...
        # we are going to handle Add terms separately,
        # if `f` is not Add -- we only have one term
        parts = []
        args = Add.make_args(f)
        for g in args:
            coeff, g = g.as_independent(x)

            # g(x) = const
            if g is S.One and not meijerg:
                parts.append(coeff*x)
                continue

            # g(x) = expr + O(x**n)
            order_term = g.getO()

            if order_term is not None:
                h = self._eval_integral(g.removeO(), x)

                if h is not None:
                    h_order_expr = self._eval_integral(order_term.expr, x)

                    if h_order_expr is not None:
                        h_order_term = order_term.func(h_order_expr, *order_term.variables)
                        parts.append(coeff*(h + h_order_term))
                        continue

                # NOTE: if there is O(x**n) and we fail to integrate then there is
                # no point in trying other methods because they will fail anyway.
                return None

            #               c
            # g(x) = (a*x+b)
            if g.is_Pow and not g.exp.has(x) and not meijerg:
                a = Wild('a', exclude=[x])
                b = Wild('b', exclude=[x])

                M = g.base.match(a*x + b)

                if M is not None:
                    if g.exp == -1:
                        h = C.log(g.base)
                    else:
                        h = g.base**(g.exp + 1) / (g.exp + 1)

                    parts.append(coeff * h / M[a])
                    continue

            #        poly(x)
            # g(x) = -------
            #        poly(x)
            if g.is_rational_function(x) and not meijerg:
                parts.append(coeff * ratint(g, x))
                continue

            if not meijerg:
                # g(x) = Mul(trig)
                h = trigintegrate(g, x)
                if h is not None:
                    parts.append(coeff * h)
                    continue

                # g(x) has at least a DiracDelta term
                h = deltaintegrate(g, x)
                if h is not None:
                    parts.append(coeff * h)
                    continue

            if not meijerg:
                # fall back to the more general algorithm
                try:
                    h = heurisch(g, x, hints=[])
                except PolynomialError:
                    # XXX: this exception means there is a bug in the
                    # implementation of heuristic Risch integration
                    # algorithm.
                    h = None
            else:
                h = None

            if meijerg is not False and h is None:
                # rewrite using G functions
                try:
                    h = meijerint_indefinite(g, x)
                except NotImplementedError:
                    from sympy.integrals.meijerint import _debug
                    _debug('NotImplementedError from meijerint_definite')
                    res = None
                if h is not None:
                    parts.append(coeff * h)
                    continue

            # if we failed maybe it was because we had
            # a product that could have been expanded,
            # so let's try an expansion of the whole
            # thing before giving up; we don't try this
            # out the outset because there are things
            # that cannot be solved unless they are
            # NOT expanded e.g., x**x*(1+log(x)). There
            # should probably be a checker somewhere in this
            # routine to look for such cases and try to do
            # collection on the expressions if they are already
            # in an expanded form
            if not h and len(args) == 1:
                f = f.expand(mul=True, deep=False)
                if f.is_Add:
                    return self._eval_integral(f, x, meijerg)


            if h is not None:
                parts.append(coeff * h)
            else:
                return None

        return Add(*parts)

    def _eval_lseries(self, x):
        for term in self.function.lseries(x):
            yield integrate(term, *self.limits)

    def _eval_nseries(self, x, n, logx):
        terms, order = self.function.nseries(x, n=n, logx=logx).as_coeff_add(C.Order)
        return integrate(terms, *self.limits) + Add(*order)*x

    def _eval_subs(self, old, new):
        """
        Substitute old with new in the integrand and the limits, but don't
        change anything that is (or corresponds to) a dummy variable of
        integration.

        The normal substitution semantics -- traversing all arguments looking
        for matching patterns -- should not be applied to the Integrals since
        changing the integration variables should also entail a change in the
        integration limits (which should be done with the transform method). So
        this method just makes changes in the integrand and the limits.

        Not all instances of a given variable are conceptually the same: the
        first argument of the limit tuple with length greater than 1 and any
        corresponding variable in the integrand are dummy variables while
        every other symbol is a symbol that will be unchanged when the integral
        is evaluated. For example, the dummy variables for ``i`` can be seen
        as symbols with a preppended underscore:

        >>> from sympy import Integral
        >>> from sympy.abc import a, b, c, x, y
        >>> i = Integral(a + x, (a, a, b))
        >>> i.as_dummy()
        Integral(_a + x, (_a, a, b))

        If you want to change the lower limit to 1 there is no reason to
        prohibit this since it is not conceptually related to the integration
        variable, _a. Nor is there reason to disallow changing the b to 1.

        If a second limit were added, however, as in:

        >>> i = Integral(x + a, (a, a, b), (b, 1, 2))

        the dummy variables become:

        >>> i.as_dummy()
        Integral(_a + x, (_a, a, _b), (_b, 1, 2))

        Note that the ``b`` of the first limit is now a dummy variable since
        ``b`` is a dummy variable in the second limit.

        The "evaluate at" form of an integral allows some flexibility in how
        the integral will be treated by subs: if there is no second argument,
        none of the symbols matching the integration symbol are considered to
        be dummy variables, but if an explicit expression is given for a limit
        then the usual interpretation of the integration symbol as a dummy
        symbol applies:

        >>> Integral(x).as_dummy() # implicit integration wrt x
        Integral(x, x)
        >>> Integral(x, x).as_dummy()
        Integral(x, x)
        >>> _.subs(x, 1)
        Integral(1, x)
        >>> i = Integral(x, (x, x))
        >>> i.as_dummy()
        Integral(_x, (_x, x))
        >>> i.subs(x, 1)
        Integral(x, (x, 1))

        Summary: no variable of the integrand or limit can be the target of
        substitution if it appears as a variable of integration in a limit
        positioned to the right of it. The only exception is for a variable
        that defines an indefinite integral limit (a single symbol): that
        symbol *can* be replaced in the integrand.

        >>> i = Integral(a + x, (a, a, 3), (b, x, c))
        >>> i.free_symbols # only these can be changed
        set([a, c, x])
        >>> i.subs(a, c) # note that the variable of integration is unchanged
        Integral(a + x, (a, c, 3), (b, x, c))
        >>> i.subs(a + x, b) == i # there is no x + a, only x + <a>
        True
        >>> i.subs(x, y - c)
        Integral(a - c + y, (a, a, 3), (b, -c + y, c))
        """
        integrand, limits = self.function, self.limits
        old_atoms = old.free_symbols
        limits = list(limits)

        dummies = set()
        for i in xrange(-1, -len(limits) - 1, -1):
            xab = limits[i]
            if len(xab) == 1:
                continue
            if not dummies.intersection(old_atoms):
                limits[i] = Tuple(xab[0],
                                  *[l._subs(old, new) for l in xab[1:]])
            dummies.add(xab[0])
        if not dummies.intersection(old_atoms):
            integrand = integrand.subs(old, new)
        return Integral(integrand, *limits)

    def as_sum(self, n, method="midpoint"):
        """
        Approximates the integral by a sum.

        method ... one of: left, right, midpoint

        This is basically just the rectangle method [1], the only difference is
        where the function value is taken in each interval.

        [1] http://en.wikipedia.org/wiki/Rectangle_method

        **method = midpoint**:

        Uses the n-order midpoint rule to evaluate the integral.

        Midpoint rule uses rectangles approximation for the given area (e.g.
        definite integral) of the function with heights equal to the point on
        the curve exactly in the middle of each interval (thus midpoint
        method). See [1] for more information.

        Examples
        ========

        >>> from sympy import sqrt
        >>> from sympy.abc import x
        >>> from sympy.integrals import Integral
        >>> e = Integral(sqrt(x**3+1), (x, 2, 10))
        >>> e
        Integral(sqrt(x**3 + 1), (x, 2, 10))
        >>> e.as_sum(4, method="midpoint")
        4*sqrt(7) + 6*sqrt(14) + 4*sqrt(86) + 2*sqrt(730)
        >>> e.as_sum(4, method="midpoint").n()
        124.164447891310
        >>> e.n()
        124.616199194723

        **method=left**:

        Uses the n-order rectangle rule to evaluate the integral, at each
        interval the function value is taken at the left hand side of the
        interval.

        Examples
        ========

        >>> from sympy import sqrt
        >>> from sympy.abc import x
        >>> e = Integral(sqrt(x**3+1), (x, 2, 10))
        >>> e
        Integral(sqrt(x**3 + 1), (x, 2, 10))
        >>> e.as_sum(4, method="left")
        6 + 2*sqrt(65) + 2*sqrt(217) + 6*sqrt(57)
        >>> e.as_sum(4, method="left").n()
        96.8853618335341
        >>> e.n()
        124.616199194723

        See Also
        ========

        Integral.doit : Perform the integration using any hints
        """

        limits = self.limits
        if len(limits) > 1:
            raise NotImplementedError("Multidimensional midpoint rule not implemented yet")
        else:
            limit = limits[0]
        if n <= 0:
            raise ValueError("n must be > 0")
        if n == oo:
            raise NotImplementedError("Infinite summation not yet implemented")
        sym, lower_limit, upper_limit = limit
        dx = (upper_limit - lower_limit)/n
        result = 0.
        for i in range(n):
            if method == "midpoint":
                xi = lower_limit + i*dx + dx/2
            elif method == "left":
                xi = lower_limit + i*dx
            elif method == "right":
                xi = lower_limit + i*dx + dx
            else:
                raise NotImplementedError("Unknown method %s" % method)
            result += self.function.subs(sym, xi)
        return result*dx


@xthreaded
def integrate(*args, **kwargs):
    """integrate(f, var, ...)

    Compute definite or indefinite integral of one or more variables
    using Risch-Norman algorithm and table lookup. This procedure is
    able to handle elementary algebraic and transcendental functions
    and also a huge class of special functions, including Airy,
    Bessel, Whittaker and Lambert.

    var can be:

    - a symbol                   -- indefinite integration
    - a tuple (symbol, a)        -- indefinite integration with result
                                    given with `a` replacing `symbol`
    - a tuple (symbol, a, b)     -- definite integration

    Several variables can be specified, in which case the result is
    multiple integration. (If var is omitted and the integrand is
    univariate, the indefinite integral in that variable will be performed.)

    Indefinite integrals are returned without terms that are independent
    of the integration variables. (see examples)

    Definite improper integrals often entail delicate convergence
    conditions. Pass conds='piecewise', 'separate' or 'none' to have
    these returned, respectively, as a Piecewise function, as a separate
    result (i.e. result will be a tuple), or not at all (default is
    'piecewise').

    **Strategy**

    SymPy uses various approaches to integration. One method is to find
    an antiderivative for the integrand, and then use the fundamental
    theorem of calculus. Various functions are implemented to integrate
    polynomial, rational and trigonometric functions, and integrands
    containing DiracDelta terms. There is also a (very successful,
    albeit somewhat slow) general implementation of the heuristic risch
    algorithm. See the docstring of Integral._eval_integral() for more
    details on computing the antiderivative using algebraic methods.

    Another family of strategies comes from re-writing the integrand in
    terms of so-called Meijer G-functions. Indefinite integrals of a
    single G-function can always be computed, and the definite integral
    of a product of two G-functions can be computed from zero to
    infinity. Various strategies are implemented to rewrite integrands
    as G-functions, and use this information to compute integrals (see
    the ``meijerint`` module).

    In general, the algebraic methods work best for computing
    antiderivatives of (possibly complicated) combinations of elementary
    functions. The G-function methods work best for computing definite
    integrals from zero to infinity of moderately complicated
    combinations of special functions, or indefinite integrals of very
    simple combinations of special functions.

    The strategy employed by the integration code is as follows:

    - If computing a definite integral, and both limits are real,
      and at least one limit is +- oo, try the G-function method of
      definite integration first.

    - Try to find an antiderivative, using all available methods, ordered
      by performance (that is try fastest method first, slowest last; in
      particular polynomial integration is tried first, meijer
      g-functions second to last, and heuristic risch last).

    - If still not successful, try G-functions irrespective of the
      limits.

    The option meijerg=True, False, None can be used to, respectively:
    always use G-function methods and no others, never use G-function
    methods, or use all available methods (in order as described above).
    It defaults to None.

    Examples
    ========

    >>> from sympy import integrate, log, exp, oo
    >>> from sympy.abc import a, x, y

    >>> integrate(x*y, x)
    x**2*y/2

    >>> integrate(log(x), x)
    x*log(x) - x

    >>> integrate(log(x), (x, 1, a))
    a*log(a) - a + 1

    >>> integrate(x)
    x**2/2

    Terms that are independent of x are dropped by indefinite integration:

    >>> from sympy import sqrt
    >>> integrate(sqrt(1 + x), (x, 0, x))
    2*(x + 1)**(3/2)/3 - 2/3
    >>> integrate(sqrt(1 + x), x)
    2*(x + 1)**(3/2)/3

    >>> integrate(x*y)
    Traceback (most recent call last):
    ...
    ValueError: specify integration variables to integrate x*y

    Note that ``integrate(x)`` syntax is meant only for convenience
    in interactive sessions and should be avoided in library code.

    >>> integrate(x**a*exp(-x), (x, 0, oo)) # same as conds='piecewise'
    Piecewise((gamma(a + 1), -re(a) < 1), (Integral(x**a*exp(-x), (x, 0, oo)), True))

    >>> integrate(x**a*exp(-x), (x, 0, oo), conds='none')
    gamma(a + 1)

    >>> integrate(x**a*exp(-x), (x, 0, oo), conds='separate')
    (gamma(a + 1), -re(a) < 1)

    See Also
    ========

    Integral, Integral.doit
    """
    meijerg = kwargs.pop('meijerg', None)
    conds = kwargs.pop('conds', 'piecewise')
    integral = Integral(*args, **kwargs)

    if isinstance(integral, Integral):
        return integral.doit(deep = False, meijerg = meijerg, conds = conds)
    else:
        return integral

@xthreaded
def line_integrate(field, curve, vars):
    """line_integrate(field, Curve, variables)

    Compute the line integral.

    Examples
    ========

    >>> from sympy import Curve, line_integrate, E, ln
    >>> from sympy.abc import x, y, t
    >>> C = Curve([E**t + 1, E**t - 1], (t, 0, ln(2)))
    >>> line_integrate(x + y, C, [x, y])
    3*sqrt(2)

    See Also
    ========

    integrate, Integral
    """
    F = sympify(field)
    if not F:
        raise ValueError("Expecting function specifying field as first argument.")
    if not isinstance(curve, Curve):
        raise ValueError("Expecting Curve entity as second argument.")
    if not is_sequence(vars):
        raise ValueError("Expecting ordered iterable for variables.")
    if len(curve.functions) != len(vars):
        raise ValueError("Field variable size does not match curve dimension.")

    if curve.parameter in vars:
        raise ValueError("Curve parameter clashes with field parameters.")

    # Calculate derivatives for line parameter functions
    # F(r) -> F(r(t)) and finally F(r(t)*r'(t))
    Ft = F
    dldt = 0
    for i, var in enumerate(vars):
        _f = curve.functions[i]
        _dn = diff(_f, curve.parameter)
        # ...arc length
        dldt = dldt + (_dn * _dn)
        Ft = Ft.subs(var, _f)
    Ft = Ft * sqrt(dldt)

    integral = Integral(Ft, curve.limits).doit(deep = False)
    return integral
