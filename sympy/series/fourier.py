"""Fourier Series"""

from __future__ import print_function, division

from sympy import pi, oo
from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.sympify import sympify
from sympy.core.singleton import S
from sympy.core.symbol import Dummy, Symbol
from sympy.core.compatibility import is_sequence
from sympy.core.containers import Tuple
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.sets.sets import Interval
from sympy.series.series_class import SeriesBase
from sympy.series.sequences import SeqFormula


def fourier_cos_seq(func, limits, n):
    """Returns the cos sequence in a fourier series"""
    from sympy.integrals import integrate
    x, L = limits[0], limits[2] - limits[1]
    cos_term = cos(2*n*pi*x / L)
    return SeqFormula(2 * cos_term * integrate(func * cos_term, limits)
                        / L, (n, 0, oo))


def fourier_sin_seq(func, limits, n):
    """Returns the sin sequence in a fourier series"""
    from sympy.integrals import integrate
    x, L = limits[0], limits[2] - limits[1]
    sin_term = sin(2*n*pi*x / L)
    return SeqFormula(2 * sin_term * integrate(func * sin_term, limits)
                        / L, (n, 1, oo))


def _process_limits(func, limits):
    """
    limits should be of the form (x, start, stop)
    x should be a symbol. Both start and stop should be bounded

    * If x is not given, x is determined from func
    * if limits is None. limit of the form (x, -pi, pi) is returned

    Examples
    ========

    >>> from sympy import pi
    >>> from sympy.series.fourier import _process_limits as pari
    >>> from sympy.abc import x
    >>> pari(x**2, (x, -2, 2))
    (x, -2, 2)
    >>> pari(x**2, (-2, 2))
    (x, -2, 2)
    >>> pari(x**2, None)
    (x, -pi, pi)
    """
    def _find_x(func):
        free = func.free_symbols
        if len(func.free_symbols) == 1:
            return free.pop()
        elif len(func.free_symbols) == 0:
            return Dummy('k')
        else:
            raise ValueError(
                " specify dummy variables for %s. If the function contains"
                " more than one free symbol, a dummy variable should be"
                " supplied explicitly e.g., FourierSeries(m*n**2, (n, -pi, pi))"
                % func)

    x, start, stop = None, None, None
    if limits is None:
        x, start, stop = _find_x(func), -pi, pi
    if is_sequence(limits, Tuple):
        if len(limits) == 3:
            x, start, stop = limits
        elif len(limits) == 2:
            x = _find_x(func)
            start, stop = limits

    if not isinstance(x, Symbol) or start is None or stop is None:
        raise ValueError('Invalid limits given: %s' % str(limits))

    unbounded = [S.NegativeInfinity, S.Infinity]
    if start in unbounded or stop in unbounded:
            raise ValueError("Both the start and end value"
                             " should be bounded")

    return sympify((x, start, stop))


class FourierSeries(SeriesBase):
    r"""Represents fourier sine/cosine series

    Examples
    ========

    >>> from sympy import FourierSeries, pi, cos
    >>> from sympy.abc import x

    >>> s = FourierSeries(x**2, (x, -pi, pi))
    >>> s.as_series(n=3)
    -4*cos(x) + cos(2*x) + pi**2/3

    Shifting

    >>> s.shift(1)
    FourierSeries(x**2 + 1, (x, -pi, pi))
    >>> s.shiftx(1)
    FourierSeries((x + 1)**2, (x, -pi, pi))

    Scaling

    >>> s.scale(2)
    FourierSeries(2*x**2, (x, -pi, pi))
    >>> s.scalex(2)
    FourierSeries(4*x**2, (x, -pi, pi))

    Differentiating

    >>> s.diff(x)
    FourierSeries(2*x, (x, -pi, pi))

    Notes
    =====

    Computing a fourier series can be slow
    due to the integration required in computing
    an, bn.

    It is faster to compute fourier series of a function
    by using shifting and scaling on an already
    computed fourier series rather than computing
    again.

    eg. If Fourier series of x**2 is known
    fourier series of x**2 - 1 can be found by shifting by 1

    References
    =========

    .. [1] mathworld.wolfram.com/FourierSeries.html
    """
    def __new__(cls, func, limits=None, **kwargs):
        func = sympify(func)

        limits = _process_limits(func, limits)
        x, lower, upper = limits

        if x not in func.free_symbols:
            return func

        kan = kwargs.pop('an', None)
        kbn = kwargs.pop('bn', None)

        if kan is None and kbn is None:
            n = Dummy('n')
            neg_func = func.subs(x, -x)
            if func == neg_func:
                an = fourier_cos_seq(func, limits, n)
                bn = SeqFormula(0, (1, oo))
            elif func == -neg_func:
                an = SeqFormula(0, (1, oo))
                bn = fourier_sin_seq(func, limits, n)
            else:
                an = fourier_cos_seq(func, limits, n)
                bn = fourier_sin_seq(func, limits, n)
        else:
            an, bn = kan, kbn

        a0 = kwargs.pop('a0', None)
        if a0 is None:
            a0 = an[0] / 2


        obj = Expr.__new__(cls, func, limits)
        obj.a0, obj.an, obj.bn = a0, an, bn
        return obj

    @property
    def function(self):
        return self.args[0]

    @property
    def x(self):
        return self.args[1][0]

    @property
    def period(self):
        return (self.args[1][1], self.args[1][2])

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

    def as_series(self, n=3):
        """"returns the first n terms(non-zero) of the series
        if n is none returns an iterator"""
        if n is None:
            return iter(self)

        i, j = 0, 0
        terms = []
        while(i < n):
                t = self.term(j)
                if t != S.Zero:
                    terms.append(t)
                    i += 1
                j += 1

        return Add(*terms)

    def shift(self, s):
        """
        Shift the function by a
        term independent of x

        f(x) -> f(x) + s

        This is fast, if fourier series of f(x) is already
        computed.

        Examples
        ========

        >>> from sympy import FourierSeries, pi
        >>> from sympy.abc import x
        >>> s = FourierSeries(x**2, (x, -pi, pi))
        >>> s.shift(1)
        FourierSeries(x**2 + 1, (x, -pi, pi))
        """
        s, x = sympify(s), self.x

        if x in s.free_symbols:
            raise ValueError("'%s' should be independent of %s" %(s, x))

        an = self.an
        bn = self.bn
        a0 = self.a0 + s
        sfunc = self.args[0] + s

        return FourierSeries(sfunc, self.args[1], a0=a0, an=an, bn=bn)

    def shiftx(self, s):
        """
        Shift x by a
        term independent of x

        f(x) -> f(x + s)

        This is fast, if fourier series of f(x) is already
        computed.

        Examples
        ========

        >>> from sympy import FourierSeries, pi
        >>> from sympy.abc import x
        >>> s = FourierSeries(x**2, (x, -pi, pi))
        >>> s.shiftx(1)
        FourierSeries((x + 1)**2, (x, -pi, pi))
        """
        s, x = sympify(s), self.x

        if x in s.free_symbols:
            raise ValueError("'%s' should be independent of %s" %(s, x))

        an = self.an.subs(x, x + s)
        bn = self.bn.subs(x, x + s)
        a0 = self.a0
        sfunc = self.args[0].subs(x, x + s)

        return FourierSeries(sfunc, self.args[1], a0=a0, an=an, bn=bn)

    def scale(self, s):
        """
        Scale function by a
        term independent of x

        f(x) -> s * f(x)

        This is fast, if fourier series of f(x) is already
        computed.

        Examples
        ========

        >>> from sympy import FourierSeries, pi
        >>> from sympy.abc import x
        >>> s = FourierSeries(x**2, (x, -pi, pi))
        >>> s.scale(2)
        FourierSeries(2*x**2, (x, -pi, pi))
        """
        s, x = sympify(s), self.x

        if x in s.free_symbols:
            raise ValueError("'%s' should be independent of %s" %(s, x))

        an = self.an.coeff_mul(s)
        bn = self.bn.coeff_mul(s)
        a0 = self.a0 * s
        sfunc = self.args[0] * s

        return FourierSeries(sfunc, self.args[1], a0=a0, an=an, bn=bn)

    def scalex(self, s):
        """
        Scale x by a
        term independent of x

        f(x) -> f(s*x)

        This is fast, if fourier series of f(x) is already
        computed.

        Examples
        ========

        >>> from sympy import FourierSeries, pi
        >>> from sympy.abc import x
        >>> s = FourierSeries(x**2, (x, -pi, pi))
        >>> s.scalex(2)
        FourierSeries(4*x**2, (x, -pi, pi))
        """
        s, x = sympify(s), self.x

        if x in s.free_symbols:
            raise ValueError("'%s' should be independent of %s" %(s, x))

        an = self.an.subs(x, x * s)
        bn = self.bn.subs(x, x * s)
        a0 = self.a0
        sfunc = self.args[0].subs(x, x * s)

        return FourierSeries(sfunc, self.args[1], a0=a0, an=an, bn=bn)

    def _eval_derivative(self, x):
        an, bn = self.an, self.bn
        van, vbn = an.variables[0], bn.variables[0]

        fan = an.formula.diff(x)
        fbn = bn.formula.diff(x)

        an = SeqFormula(fan, an.args[1])
        bn = SeqFormula(fbn, bn.args[1])

        dfunc = self.args[0].diff(x)

        return FourierSeries(dfunc, self.args[1], an=an, bn=bn)

    def _eval_as_leading_term(self, x):
        for t in self:
            if t != S.Zero:
                return t

    def _eval_term(self, pt):
        if pt == 0:
            return self.a0
        return self.an.coeff(pt) + self.bn.coeff(pt)

    def __neg__(self):
        return self.scale(-1)

    def __add__(self, other):
        if isinstance(other, FourierSeries):
            if self.period != other.period:
                raise ValueError("Both the series should have same periods")

            x, y = self.x, other.x
            function = self.function + other.function.subs(y, x)
            an = self.an + other.an
            bn = self.bn + other.bn
            a0 = self.a0 + other.a0
            return FourierSeries(function, self.args[1], a0=a0, an=an, bn=bn)

        return Add(self, other)

    def __sub__(self, other):
        return self.__add__(-other)
