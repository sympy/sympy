"""Formal Power Series"""

from __future__ import print_function, division

from sympy import oo, zoo, nan
from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.core.symbol import Wild, Dummy, symbols
from sympy.sets.sets import Interval
from sympy.functions.combinatorial.factorials import binomial, factorial
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
    ds = [] # list of diff

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
            if coeff.has(x) or coeff.has(zoo) or coeff.has(oo) or coeff.has(nan):
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


def compute_fps(f, x, x0=0, dir=1, hyper=True, order=4, rational=True, full=False):
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

    result = None

    # from here on it's x0=0 and dir=1 handling
    if rational:
        k = Dummy('k')
        result = rational_algorithm(f, x, k, order, full)

    if result is None and hyper:
        return None # TODO

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
            pt_ak = self.ak.coeff(pt).simplify() # TODO: Thoughts?
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
    f= sympify(f)

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
