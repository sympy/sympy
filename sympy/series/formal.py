"""Formal Power Series"""

from __future__ import print_function, division

from sympy import oo, zoo, nan
from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.function import Derivative, Function
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.core.symbol import Wild, Dummy, symbols, Symbol
from sympy.core.relational import Eq
from sympy.sets.sets import Interval
from sympy.functions.combinatorial.factorials import binomial, factorial, rf
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.series.sequences import sequence
from sympy.series.series_class import SeriesBase
from sympy.series.order import Order


def rational_algorithm(f, x, k, order=4, full=False):
    """Rational algorithm for computing
    formula of coefficients of Formal Power Series
    of a function.

    Applicable when f(x) or some derivative of f(x)
    is a rational function in x.

    ``rational_algorithm`` uses :func:`apart` function for partial fraction
    decomposition. ``apart`` by default uses 'undetermined coefficients method'.
    By setting full=True, 'Bronstein's algorithm' can be used instead.

    Looks for derivative of a function upto 4'th order(by default).
    This can be overrided using order option.

    returns a Tuple of (formula, series independent terms, order) if successful
    otherwise None.

    Examples
    ========

    >>> from sympy import log, atan, I
    >>> from sympy.series.formal import rational_algorithm as ra
    >>> from sympy.abc import x, k

    >>> ra(1 / (1 - x), x, k)
    (1, 0, 0)
    >>> ra(log(1 + x), x, k)
    (-(-1)**(-k)/k, 0, 1)

    >>> ra(atan(x), x, k, full=True)
    ((-I*(-I)**(-k)/2 + I*I**(-k)/2)/k, 0, 1)

    Notes
    =====

    By setting full=True, range of admissible functions increases.
    This option should be used carefully as it can signifcantly
    slow down the computation as ``doit`` is performed on the :class:`RootSum`
    object returned by the ``apart`` function. Use full=False whenever
    possible.

    See Also
    ========

    sympy.polys.partfrac.apart

    References
    ==========

    .. [1] Formal Power Series - Dominik Gruntz, Wolfram Koepf
    .. [2] Power Series in Computer Algebra - Wolfram Koepf
    """

    from sympy.polys import RootSum, apart

    diff = f
    ds = []  # list of diff

    for i in range(order + 1):
        if i:
            diff = diff.diff(x)

        if diff.is_rational_function(x):
            coeff, sep = S.Zero, S.Zero

            terms = apart(diff, x, full=full)
            if terms.has(RootSum):
                terms = terms.doit()

            for t in Add.make_args(terms):
                num, den = t.as_numer_denom()
                if not den.has(x):
                    sep += t
                else:
                    if isinstance(den, Mul):
                        # m*(n*x - a)**j -> (n*x - a)**j
                        ind = den.as_independent(x)
                        den = ind[1]
                        num /= ind[0]

                    # (n*x - a)**j -> (x - b)
                    den, j = den.as_base_exp()
                    a, xterm = den.as_coeff_add(x)

                    # term -> m/x**n
                    if not a:
                        sep += t
                        continue

                    xc = xterm[0].coeff(x)
                    a /= -xc
                    num /= xc**j

                    ak = ((-1)**j * num * binomial(j + k - 1, k)
                          .rewrite(factorial) / a**(j + k))
                    coeff += ak

            # Hacky, better way?
            if coeff is S.Zero:
                return None
            if (coeff.has(x) or coeff.has(zoo) or coeff.has(oo) or
                    coeff.has(nan)):
                return None

            from sympy.integrals import integrate
            for j in range(i):
                coeff = (coeff / (k + j + 1))
                sep = integrate(sep, x)
                c = ds.pop().limit(x, 0) - sep.limit(x, 0)
                if c.is_finite is False:
                    return None
                sep += c
            return (coeff.subs(k, k - i), sep, i)

        else:
            ds.append(diff)

    return None


def rational_independent(terms, x):
    """returns a list of all the rationally independent terms

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.series.formal import rational_independent
    >>> from sympy.abc import x

    >>> rational_independent([cos(x), sin(x)], x)
    [cos(x), sin(x)]
    >>> rational_independent([x**2, sin(x), x*sin(x), x**3], x)
    [x**3 + x**2, x*sin(x) + sin(x)]
    """
    if not terms:
        return []

    ind = terms[0:1]

    for t in terms[1:]:
        dep = False
        n = t.as_independent(x)[1]
        for i in range(len(ind)):
            d = ind[i].as_independent(x)[1]
            q = (n / d).cancel()
            if q.is_rational_function(x):
                ind[i] += t
                dep = True
                break
        if not dep:
            ind.append(t)
    return ind


def simpleDE(f, x, g, order=4):
    """
    Computes a simple DE of the form

    .. math::
        f^k(x) + \sum\limits\_{j=0}^{k-1} A_j f^j(x) = 0

    where :math:`A_j` should be rational function in x.

    By default DE is tried uptill order 4. If it is not
    found ``None`` is returned. By increasing order, higher
    order DE's can be found.

    Examples
    ========

    >>> from sympy import Function, exp, ln
    >>> from sympy.series.formal import simpleDE
    >>> from sympy.abc import x
    >>> f = Function('f')

    >>> simpleDE(exp(x), x, f)
    -f(x) + Derivative(f(x), x)
    >>> simpleDE(ln(1 + x), x, f)
    (x + 1)*Derivative(f(x), x, x) + Derivative(f(x), x)
    """
    from sympy.solvers import solve
    a = symbols('a:%d' % (order))

    def _makeDE(k):
        eq = f.diff(x, k) + Add(*[a[i]*f.diff(x, i) for i in range(0, k)])
        DE = g(x).diff(x, k) + Add(*[a[i]*g(x).diff(x, i) for i in range(0, k)])
        return eq, DE

    eq, DE = _makeDE(order)

    for k in range(1, order + 1):
        eq, DE = _makeDE(k)
        eq = eq.expand()
        terms = eq.as_ordered_terms()
        ind = rational_independent(terms, x)
        if len(ind) == k:
            sol = solve(ind, a, dict=True)
            if sol:
                for k, v in sol[0].items():
                    DE = DE.subs(k, v)
                DE = DE.as_numer_denom()[0]
                DE = DE.factor().as_coeff_mul(Derivative)[1][0]
                return DE.collect(Derivative(g(x)))


def exp_re(DE, r, k):
    """
    Converts a DE with constant coefficients (explike)
    into a RE

    substitutes :math:`f^j(x) \to r(k + j)`

    Normalises so that lowest term is always r(k).

    Examples
    ========

    >>> from sympy import Function, Derivative
    >>> from sympy.series.formal import exp_re
    >>> from sympy.abc import x, k
    >>> f, r = Function('f'), Function('r')

    >>> exp_re(-f(x) + Derivative(f(x)), r, k)
    -r(k) + r(k + 1)
    >>> exp_re(Derivative(f(x), x) + Derivative(f(x), x, x), r, k)
    r(k) + r(k + 1)

    See Also
    ========

    sympy.series.formal.hyper_re
    """
    RE = S.Zero

    g = DE.atoms(Function).pop()

    mini = None
    for t in Add.make_args(DE):
        coeff, d = t.as_independent(g)
        if isinstance(d, Derivative):
            j = len(d.args[1:])
        else:
            j = 0
        if mini is None or j < mini:
            mini = j
        RE += coeff * r(k + j)
    if mini:
        RE = RE.subs(k, k - mini)
    return RE


def hyper_re(DE, r, k):
    """
    Converts a DE into a RE

    substitutes :math:`x^l f^j(x) \to (k + 1 - l)_j . a_{k + j - l}`

    Normalises so that lowest term is always r(k).

    Examples
    ========

    >>> from sympy import Function, Derivative
    >>> from sympy.series.formal import hyper_re
    >>> from sympy.abc import x, k
    >>> f, r = Function('f'), Function('r')

    >>> hyper_re(-f(x) + Derivative(f(x)), r, k)
    (k + 1)*r(k + 1) - r(k)
    >>> hyper_re(-x*f(x) + Derivative(f(x), x, x), r, k)
    (k + 2)*(k + 3)*r(k + 3) - r(k)

    See Also
    ========

    sympy.series.formal.exp_re
    """
    RE = S.Zero

    g = DE.atoms(Function).pop()
    x = g.atoms(Symbol).pop()

    mini = None
    for t in Add.make_args(DE.expand()):
        coeff, d = t.as_independent(g)
        c, v = coeff.as_independent(x)
        l = v.as_coeff_exponent(x)[1]
        if isinstance(d, Derivative):
            j = len(d.args[1:])
        else:
            j = 0
        RE += c * rf(k + 1 - l, j) * r(k + j - l)
        if mini is None or j - l < mini:
            mini = j - l

    RE = RE.subs(k, k - mini)

    m = Wild('m')
    return RE.collect(r(k + m))


def rsolve_hypergeometric(f, x, P, Q, k, m):
    """
    Attempts to solve RE of the form

    Q(k)*a(k + m) - P(k)*a(k)

    Transformations that preserve Hypergeometric type:
        a. x**n*f(x): b(k + m) = R(k - n)*b(k)
        b. f(A*x): b(k + m) = A**m*R(k)*b(k)
        c. f(x**n): b(k + n*m) = R(k/n)*b(k)
        d. f(x**(1/m)): b(k + 1) = R(k*m)*b(k)
        e. f'(x): b(k + m) = ((k + m + 1)/(k + 1))*R(k + 1)*b(k)

    Some of these transformations have been used to solve
    the RE

    returns a Tuple of (formula, series independent terms, order) if successful
    otherwise None.

    Examples
    ========

    >>> from sympy import exp, ln, S
    >>> from sympy.series.formal import rsolve_hypergeometric as rh
    >>> from sympy.abc import x, k

    >>> rh(exp(x), x, -S.One, (k + 1), k, 1)
    (Piecewise((1/(factorial(k)), Eq(Mod(k, 1), 0)), (0, True)), 1, 1)

    >>> rh(ln(1 + x), x, k**2, k*(k + 1), k, 1)
    (Piecewise(((-1)**(k - 1)*factorial(k - 1)/RisingFactorial(2, k - 1),
     Eq(Mod(k, 1), 0)), (0, True)), x, 2)

    References
    ==========

    .. [1] Formal Power Series - Dominik Gruntz, Wolfram Koepf
    .. [2] Power Series in Computer Algebra - Wolfram Koepf
    """
    from sympy.polys import lcm, roots

    # tranformation - c
    proots, qroots = roots(P, k), roots(Q, k)
    sol = dict(proots)
    sol.update(qroots)
    scale = lcm([r.as_numer_denom()[1] for r, t in sol.items()
                 if r.is_rational])

    f = f.subs(x, x**scale)
    m *= scale
    P = P.subs(k, k / scale)
    Q = Q.subs(k, k / scale)

    # transformation - a
    qroots = roots(Q, k)
    k_min = Min(*qroots.keys())
    shift = k_min + m
    f *= x**(-shift)
    P = P.subs(k, k + shift)
    Q = Q.subs(k, k + shift)

    if (x*f).limit(x, 0) is not S.Zero:
        return None

    qroots = roots(Q, k)
    k_max = Max(*qroots.keys())
    ind = S.Zero
    if k_max + m >= 0:
        for i in range(k_max + m + 1):
            r = f.diff(x, i).limit(x, 0) / factorial(i)
            ind += r*x**(i + shift)
        ind = ind.subs(x, x**(1/scale))
        s, e = k_max + 1, m + k_max + 1
    else:
        s, e = 0, m

    sol = []
    for i in range(s, e):
        res = S.One

        r = f.diff(x, i).limit(x, 0) / factorial(i)
        if r is S.Zero:
            continue
        elif r.is_finite is False:
            return None

        res *= r

        p = P.subs(k, m*k + i)
        q = Q.subs(k, m*k + i)
        c1 = p.subs(k, 1/k).leadterm(k)[0]
        c2 = q.subs(k, 1/k).leadterm(k)[0]
        res *= (-c1 / c2)**k

        for r, mul in roots(p, k).items():
            res *= rf(-r, k)**mul
        for r, mul in roots(q, k).items():
            res /= rf(-r, k)**mul

        t_p = (m*k + i + shift) / scale
        j, mk = t_p.as_coeff_Add()
        c = mk.coeff(k)

        res = res.subs(k, (k - j) / c)

        sol.append((res, Eq(k % c, j % m)))

    sol.append((S.Zero, True))

    if k_max + m >= 0:
        return (Piecewise(*sol), ind, k_max + m + 1 + shift)
    else:
        return (Piecewise(*sol), S.Zero, S.Zero)


def solve_re(f, x, RE, g, k):
    """
    Solves the RE

    If The RE is of the form Q(k)*a(k + m) - P(k)*a(k),
    uses :func:`rsolve_hypergeometric`,
    otherwise fallsback to :func:`rsolve`

    returns a Tuple of (formula, series independent terms, order) if successful
    otherwise None.

    Examples
    ========

    >>> from sympy import exp, ln
    >>> from sympy.series.formal import solve_re
    >>> from sympy.abc import x, k, f

    >>> solve_re(exp(x), x, (k+1)*f(k+1) - f(k), f, k)
    (Piecewise((1/(factorial(k)), Eq(Mod(k, 1), 0)), (0, True)), 1, 1)

    >>> solve_re(ln(1 + x), x, k*(k + 1)*f(k + 1) + k**2*f(k), f, k)
    (Piecewise(((-1)**(k - 1)*factorial(k - 1)/RisingFactorial(2, k - 1),
     Eq(Mod(k, 1), 0)), (0, True)), x, 2)

    See Also
    ========

    sympy.series.formal.rsolve_hypergeometric
    sympy.solvers.recurr.rsolve

    """
    m = Wild('m')
    RE = RE.collect(g(k + m))
    terms = Add.make_args(RE)

    if len(terms) == 2:
        gs = list(RE.atoms(Function))
        P, Q = map(RE.coeff, gs)
        m = gs[1].args[0] - gs[0].args[0]
        if m < 0:
            P, Q = Q, P
            m = abs(m)
        return rsolve_hypergeometric(f, x, P, Q, k, m)

    init = {}
    for i in range(len(terms)):
        if i:
            f = f.diff(x)
        init[g(k).subs(k, i)] = f.subs(x, 0) / factorial(i)

    from sympy.solvers import rsolve
    sol = rsolve(RE, g(k), init)

    if sol:
        return (sol, S.Zero, S.Zero)


def hyper_algorithm(f, x, k, order=4):
    """Hypergeometric algorithm

    Steps:
        * Compute a simpleDE
        * Convert the DE into RE
        * Solve The RE

    Examples
    ========

    >>> from sympy import exp, ln
    >>> from sympy.series.formal import hyper_algorithm

    >>> from sympy.abc import x, k

    >>> hyper_algorithm(exp(x), x, k)
    (Piecewise((1/(factorial(k)), Eq(Mod(k, 1), 0)), (0, True)), 1, 1)

    >>> hyper_algorithm(ln(1 + x), x, k)
    (Piecewise(((-1)**(k - 1)*factorial(k - 1)/RisingFactorial(2, k - 1),
     Eq(Mod(k, 1), 0)), (0, True)), x, 2)

    See Also
    ========

    sympy.series.formal.simpleDE
    sympy.series.formal.hyper_re
    sympy.series.formal.solve_re
    """
    g = Function('g')

    DE = simpleDE(f, x, g, order)

    if DE is None:
        return None

    RE = hyper_re(DE, g, k)
    sol = solve_re(f, x, RE, g, k)

    if sol is None:
        return None

    return sol


def compute_fps(f, x, x0=0, dir=1, hyper=True, order=4, rational=True,
                full=False):
    """Computes the formula for Formal Power Series of a function

    returns (sequence of coefficients, sequence of x, independent terms)

    Tries to compute the formula by applying the following techniques
    (in order)

    * rational_algorithm
    * Hypergeomitric algorithm (TODO)

    Arguments
    =========

    * x0 = series expansion around ``x = x0``(Default = 0).

    * dir = For ``dir=1`` (default) the series is calculated from the right and
    for ``dir=-1`` the series from the left. For smooth functions this
    flag will not alter the results.

    * hyper = set hyper=False, for not using hypergeometric algorithm.

    * order = Order of the derivative of f, till which algorithms
    are run.

    * rational = set rational=False to skip rational algorithm

    * full = full by default is False. Try setting full to True
    to increase the range of rational algorithm. See :func:`rational_algorithm`
    for details.

    See Also
    ========

    sympy.series.rational_algorithm
    """
    if x0 in [S.Infinity, -S.Infinity]:
        dir = {S.Infinity: S.One, S.NegativeInfinity: -S.One}[x0]
        temp = f.subs(x, 1/x)
        result = compute_fps(temp, x, 0, dir, hyper, order, rational, full)
        if result is None:
            return None
        return result[0], result[1].subs(x, 1/x), result[2].subs(x, 1/x)
    elif x0 or dir == -S.One:
        if dir == -S.One:
            rep = -x + x0
            rep2 = -x
            rep2b = x0
        else:
            rep = x + x0
            rep2 = x
            rep2b = -x0
        temp = f.subs(x, rep)
        result = compute_fps(temp, x, 0, S.One, hyper, order, rational, full)
        if result is None:
            return None
        return (result[0], result[1].subs(x, rep2 + rep2b),
                result[2].subs(x, rep2 + rep2b))

    #  Break instances of Add
    #  this allows application of different
    #  algorithms on different terms increasing the
    #  range of admissible functions.
    if isinstance(f, Add):
        result = False
        ak = sequence(S.Zero, (0, oo))
        ind, xk = S.Zero, None
        for t in Add.make_args(f):
            res = compute_fps(t, x, 0, S.One, hyper, order, rational, full)
            if res:
                if not result:
                    result = True
                    xk = res[1]
                if res[0].start > ak.start:
                    seq = ak
                    s, f = ak.start, res[0].start
                else:
                    seq = res[0]
                    s, f = res[0].start, ak.start
                save = Add(*[z[0]*z[1] for z in zip(seq[0:(f - s)], xk[s:f])])
                ak += res[0]
                ind += res[2] + save
            else:
                ind += t
        if result:
            return ak, xk, ind
        return None

    if f.is_polynomial(x):
        return None

    result = None

    # from here on it's x0=0 and dir=1 handling
    k = Dummy('k')
    if rational:
        result = rational_algorithm(f, x, k, order, full)

    if result is None and hyper:
        result = hyper_algorithm(f, x, k, order)

    if result is None:
        return None

    ak = sequence(result[0], (k, result[2], oo))
    xk = sequence(x**k, (k, 0, oo))
    ind = result[1]

    return ak, xk, ind


class FormalPowerSeries(SeriesBase):
    """Represents Formal Power Series

    No computation is performed.
    This class should only to be used to represent
    a series. No checks are performed.

    For computing a series use :func:`fps`.

    See Also
    ========

    sympy.series.formal.fps
    sympy.series.formal.compute_fps
    """
    def __new__(cls, *args):
        args = map(sympify, args)
        return Expr.__new__(cls, *args)

    @property
    def function(self):
        return self.args[0]

    @property
    def x(self):
        return self.args[1]

    @property
    def x0(self):
        return self.args[2]

    @property
    def ak(self):
        return self.args[4][0]

    @property
    def xk(self):
        return self.args[4][1]

    @property
    def ind(self):
        return self.args[4][2]

    @property
    def interval(self):
        return Interval(0, oo)

    @property
    def start(self):
        return self.interval.inf

    @property
    def stop(self):
        return self.interval.sup

    @property
    def length(self):
        return oo

    @property
    def infinite(self):
        """returns an infinite representation of the series"""
        from sympy.concrete import Sum

        ak, xk = self.ak, self.xk
        k = ak.variables[0]

        return self.ind + Sum(ak.formula * xk.formula, (k, ak.start, ak.stop))

    def polynomial(self, n=6):
        """truncated series as polynomial.

        returns series sexpansion of f upto order ``O(x**n)``
        as a polynomial (without ``O`` term)
        """
        terms = []
        for i, t in enumerate(self):
            if i >= n:
                break
            if t is not S.Zero:
                terms.append(t)

        return Add(*terms)

    def truncate(self, n=6):
        """truncated series.

        returns truncated series expansion of f upto
        order ``O(x**n)``

        if n is ``None``, returns an infinite iterator
        """
        if n is None:
            return iter(self)

        x, x0 = self.x, self.x0
        pt_xk = self.xk.coeff(n)
        if x0 is S.NegativeInfinity:
            x0 = S.Infinity

        return self.polynomial(n) + Order(pt_xk, (x, x0))

    def _eval_term(self, pt):
        pt_xk = self.xk.coeff(pt)

        def _get_xterm(t):
            if self.x0 in [S.Infinity, -S.Infinity]:
                t = t.as_numer_denom()[1]
            else:
                t = t.as_numer_denom()[0]
            return t.as_independent(self.x)[1]

        ind = S.Zero
        if self.ind:
            for t in Add.make_args(self.ind):
                xterm = _get_xterm(t)
                if pt_xk == 1 and xterm == 1:
                    ind += t
                elif xterm != 1 and pt_xk != 1:
                    j = xterm.as_base_exp()[1]
                    if j == pt:
                        ind += t

        try:
            pt_ak = self.ak.coeff(pt).simplify()  # TODO: Thoughts?
        except IndexError:
            pt_ak = S.Zero

        term = (pt_ak * pt_xk) + ind

        return term.collect(self.x)

    def _eval_subs(self, old, new):
        x = self.x
        if old.has(x):
            return self

    def _eval_as_leading_term(self, x):
        for t in self:
            if t is not S.Zero:
                return t


def fps(f, x=None, x0=0, dir=1, hyper=True, order=4, rational=True, full=False):
    """Generates Formal Power Series of f

    returns the formal series expansion of ``f`` around ``x = x0``
    with respect to ``x`` in the form of a ``FormalPowerSeries`` object

    Formal Power Series is represented using an explicit formula
    computed using different algorithms.

    See :func:`compute_fps` for the more details regarding the computation
    of formula.

    Arguments
    =========

    * x = if ``x=None`` and ``f`` is univariate, the univariate
    symbols will be supplied, otherwise an error will be raised.

    * x0 = series expansion around ``x = x0``(Default = 0).

    * dir = For ``dir=1`` (default) the series is calculated from the right and
    for ``dir=-1`` the series from the left. For smooth functions this
    flag will not alter the results.

    * hyper = set hyper=False, for not using hypergeometric algorithm.

    * order = Order of the derivative of f, till which algorithms
    are run.

    * rational = set rational=False to skip rational algorithm

    * full = full by default is False. Try setting full to True
    to increase the range of rational algorithm. See :func:`rational_algorithm`
    for details.

    Examples
    ========

    >>> from sympy import fps, O, ln, atan
    >>> from sympy.abc import x

    Rational Functions

    >>> fps(ln(1 + x)).truncate()
    x - x**2/2 + x**3/3 - x**4/4 + x**5/5 + O(x**6)

    >>> fps(atan(x), full=True).truncate()
    x - x**3/3 + x**5/5 + O(x**6)

    See Also
    ========

    sympy.series.formal.FormalPowerSeries
    sympy.series.formal.compute_fps
    """
    f = sympify(f)

    if x is None:
        free = f.free_symbols
        if len(free) == 1:
            x = free.pop()
        elif not free:
            return f
        else:
            raise NotImplementedError("multivariate formal power series")
    else:
        x = sympify(x)

    if not f.has(x):
        return f

    x0 = sympify(x0)
    dir = sympify(dir)

    if dir not in [S.One, -S.One]:
        raise ValueError("Dir must be 1 or -1")

    result = compute_fps(f, x, x0, dir, hyper, order, rational, full)

    if result is None:
        return f

    ak, xk, ind = result

    return FormalPowerSeries(f, x, x0, dir, (ak, xk, ind))
