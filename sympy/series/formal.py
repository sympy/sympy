from __future__ import print_function, division

from sympy import C, S, collect, Function, Add, Mul, simplify, cancel, sympify, Basic
from sympy.core.containers import Stream, Tuple
from sympy.core.symbol import Symbol, symbols, Dummy, Wild
from sympy.core.sympify import sympify
from sympy.core.relational import Eq
from sympy.functions import sign
from sympy.functions.elementary.piecewise import Piecewise
from sympy.concrete import Sum
from sympy.core import Add, Pow
from sympy.core.expr import Expr
from sympy.polys import roots, lcm, degree
from sympy.polys.partfrac import apart
from sympy.solvers import solve, rsolve


def FormalSeries(x=None, x0=0, dir='+', *args, **kwargs):
    """
    Finds generator for given function and returns a ``Lazyseries`` object.

    Kwargs
    ======

    function
    sequence
    generator

    Examples
    ========

    >>> from sympy.series.formal import FormalSeries
    >>> from sympy.abc import x
    >>> FormalSeries(x, function=1/(1-x))
    Lazyseries(...)
    >>> FormalSeries(x, sequence=(1, 2, 3))
    Lazyseries(...)

    """
    k = Dummy('k', integer=True)
    generator = kwargs.pop("generator", None)
    if generator:
        c, e, k = generator

    sequence = kwargs.pop('sequence', None)
    if sequence:
        c, e = findgen_seq(sequence, k)

    function = kwargs.pop('function', None)
    if function is not None:
        function = sympify(function)

        if len(dir) != 1 or dir not in '+-':
            raise ValueError("Dir must be '+' or '-'")

        if x0 in [S.Infinity, S.NegativeInfinity]:
            dir = {S.Infinity: '+', S.NegativeInfinity: '-'}[x0]
            return FormalSeries(x, dir=dir, function=function.subs(x, 1/x)).subs(x, 1/x)

        if x0 or dir == '-':
            if dir == '-':
                rep = -x + x0
                rep2 = -x + x0
            else:
                rep = x + x0
                rep2 = x - x0
            return FormalSeries(x, function=function.subs(x, rep)).subs(x, rep2)

        if x.is_positive is x.is_negative is None:
            xpos = Dummy('x', positive=True, bounded=True)
            return FormalSeries(xpos, function=function.subs(x, xpos)).subs(xpos, x)

        generator = findgen(x, function, k)
        return sum(FormalSeries(x, generator=(c, e, k)) for c, e in generator)

    if not e.has(k):
        def gen():
            yield c*x**e
            raise StopIteration
    else:
        def gen():
            n = 0
            while True:
                nc, ne = c.subs(k, n), e.subs(k, n)
                yield nc*x**ne
                n = n + 1
    return Lazyseries(x, Stream(gen()))


def findgen_seq(sequence, k):
    """
    Returns generator for a sequence in terms of k

    Examples
    ========

    >>> from sympy.series.formal import findgen_seq
    >>> from sympy.abc import x, k
    >>> findgen_seq((1, 2, 3), k)
    (Piecewise((1, Mod(k, 3) == 0), (2, Mod(k, 3) == 1), (3, Mod(k, 3) == 2)), k)
    >>> findgen_seq([1], k)
    (Piecewise((1, Mod(k, 1) == 0)), k)

    """
    l = len(sequence)
    cond = [(sequence[i], Eq(k%l, i)) for i in range(l)]
    return (Piecewise(*cond), k)


def findgen_rational(x, function, k):
    """
    Returns the generator for a given rational function in x in terms of k

    Examples
    ========

    >>> from sympy.series.formal import findgen_rational
    >>> from sympy.abc import x, k
    >>> findgen_rational(x, 1/(1-x), k)
    [(1, k)]
    >>> findgen_rational(x, 1/(1+x), k)
    [(-(-1)**(-k - 1), k)]

    """
    gen = []
    m = Wild('m', exclude=(0,))
    terms = Add.make_args(apart(function))
    for t in terms:
        n, d = t.as_numer_denom()
        if not d.has(x):
            c, e = t.as_coeff_exponent(x)
            gen += [(c, e)]
        elif not d.match(x+m):
            raise ValueError('Expected denominator of type "x + a", got: %s' % d)
        else:
            d, j = d.as_base_exp()
            a = -d.as_coeff_add()[0]
            c, e = (-1)**j * n * C.binomial(j+k-1, k).rewrite(C.factorial) / a**(j+k), k
            gen += [(c, e)]
    return gen


def rational_independent(x, terms):
    """
    Returns list of rationally independent terms

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.series.formal import rational_independent
    >>> from sympy.abc import x
    >>> rational_independent(x, [sin(x), cos(x)])
    [sin(x), cos(x)]
    >>> rational_independent(x, [sin(x), cos(x), x*sin(x)])
    [x*sin(x) + sin(x), cos(x)]

    """
    ind = terms[0:1]
    for t in terms[1:]:
        dep = False
        for i in range(len(ind)):
            n = t.as_independent(x)[1]
            d = ind[i].as_independent(x)[1]
            q = simplify(n/d)
            if q.is_rational_function():
                ind[i] = ind[i] + t
                dep = True
                break
        if not dep:
            ind.append(t)
    return ind


def simpleDE(x, function, f):
    """
    Converts a function into a simple differential equation.

    A simple differential equation is a linear homogeneous differential equation
    with rational coefficients.

    Examples:
    ========

    >>> from sympy import sin, exp, Function
    >>> from sympy.series.formal import simpleDE
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> simpleDE(x, exp(x), f)
    -f(x) + Derivative(f(x), x)
    >>> simpleDE(x, sin(x), f)
    f(x) + Derivative(f(x), x, x)

    """
    a = symbols('a:4')

    makeDE = lambda k: (function.diff(x, k) + \
            Add(*[a[i]*function.diff(x, i) for i in range(0, k)]), f(x).diff(x, k) +
            Add(*[a[i]*f(x).diff(x, i) for i in range(0, k)]))

    # Solve for case k=1
    eq, DE = makeDE(1)
    sol = solve(eq, a[0])
    if sol and sol[0].is_rational_function():
        sol[0] = cancel(sol[0])
        DE = DE.subs(a[0], sol[0])
        DE = simplify(DE)
        return DE.as_numer_denom()[0]

    # Solve for cases from k=2
    # Search upto 4th order DE
    # Good enough for most functions
    for k in range(2, 5):
        eq, DE = makeDE(k)
        eq = eq.expand()
        terms = [t for t in eq.as_ordered_terms()]
        ind = rational_independent(x, terms)
        if len(ind) == k:
            sol = solve(ind, a, dict=True)
            if sol:
                for key, val in sol[0].iteritems():
                    val = cancel(val)
                    DE = DE.subs(key, val)
            DE = simplify(DE)
            return DE.as_numer_denom()[0]

    raise NotImplementedError('Cannot find simple differential equation for %s' % function)


def DEtoRE(DE, r, k):
    """
    Converts a differential equation into a recurrence equation

    Parameters:
        r: Function in which RE whill be expressed
        n: Argument of recurr function

    Examples
    ========

    >>> from sympy import Derivative, Function
    >>> from sympy.series.formal import DEtoRE
    >>> from sympy.abc import x, k
    >>> f, r = Function('f'), Function('r')
    >>> DE = -f(x) + Derivative(f(x), x)
    >>> DEtoRE(DE, r, k)
    (k + 1)*r(k + 1) - r(k)

    """
    RE = S.Zero
    DE = DE.expand()

    f = DE.atoms(Function).pop()
    x = f.atoms(Symbol).pop()

    # Minimum argument of ``r`` in the RE
    # This converts recurrences like
    #   r(k+1) - r(k-1)
    # to
    #   r(k+2) - r(k)
    m = None

    for t in Add.make_args(DE):
        coeff, d = t.as_independent(f)
        c, v = coeff.as_independent(x)
        l = v.as_coeff_exponent(x)[1]
        if d.is_Derivative:
            j = len(d.args) - 1
        else:
            j = 0
        RE += c * C.RisingFactorial(k+1-l, j) * r(k+j-l)
        if not m or k+j-l < m:
            m = k+j-l
    sol = solve(m, k)
    if sol and sol[0] > 0:
        RE = RE.subs(k, k+sol[0])
    return RE


def rsolve_hypergeometric(P, Q, m, k, function):
    """
    Solves recurrence relation of hypergeometric type

    Q(k)*r(k + m) = P(k)*r(k)

    """
    # Following transformation preserve hypergeometric type
    #   (a) x**n*f(x)   ->  r(k + m) = R(k - n) * b(k)
    #   (b) f(A*x)      ->  r(k + m) = A**m * R(k) * b(k)
    #   (c) f(x**n)     ->  r(k + n*m) = R(k/n) * b(k)

    # Transformation (b)
    x = function.atoms(Symbol).pop()

    sol = solve(P, k) + solve(Q, k)
    scale = lcm([s.as_numer_denom()[1] for s in sol])
    m = m*scale
    function = function.subs(x, x**scale)
    P = P.subs(k, k/scale)
    Q = Q.subs(k, k/scale)

    # Transformation (a)
    k_min = min(solve(P, k) + solve(Q, k))
    shift = k_min + m
    function = function * x**(-shift)
    P = P.subs(k, k + shift)
    Q = Q.subs(k, k + shift)

    gen = []
    for j in range(m):
        init = function.diff(x, j).subs(x, 0) / C.factorial(j)
        if init == 0:
            continue
        p = P.subs(k, m*k + j)
        q = Q.subs(k, m*k + j)
        c1 = p.subs(k, 1/k).leadterm(k)[0]
        c2 = q.subs(k, 1/k).leadterm(k)[0]
        c = -c1/c2

        res = S.One
        for sol, mul in roots(p, k).iteritems():
            res *= C.RisingFactorial(-sol, k)**mul
        for sol, mul in roots(q, k).iteritems():
            res /= C.RisingFactorial(-sol, k)**mul
        res *= c**k
        gen += [(res, (m*k + j + shift)/scale)]
    return gen


def solveRE(RE, r, k, function):
    """
    Returns generator for a given recurrence equation

    Examples
    ========

    >>> from sympy import exp, Function
    >>> from sympy.series.formal import solveRE
    >>> from sympy.abc import x, k
    >>> r = Function('r')
    >>> RE = (k+1)*r(k+1) - r(k)
    >>> solveRE(RE, r, k, exp(x))
    [(1/(factorial(k)), k)]

    """
    RE = RE.expand().collect(r(k).func(Wild('m')))
    terms = Add.make_args(RE)

    if len(terms) == 2:
        f = list(RE.atoms(Function))
        P, Q = RE.coeff(f[0]), RE.coeff(f[1])
        m = f[1].args[0] - f[0].args[0]
        if m < 0:
            P, Q = Q, P
            m = -m
        return rsolve_hypergeometric(P, Q, m, k, function)

    x = function.atoms(Symbol).pop()
    key = [r(i) for i in range(0, len(terms)-1)]
    val = [function.diff(x, i).subs(x, 0)/C.factorial(i) for i in range(0, len(terms)-1)]
    init = dict(zip(key, val))

    sol = rsolve(RE, r(k), init)
    if sol:
        return [(sol, k)]

    raise NotImplementedError('Cannot solve %s' % RE)


def findgen(x, function, k):
    """
    Returns the generator for a given function in x

    The algorithm has three parts:
        (1) Convert function to a simple differential equation.
        (2) Convert differential equation to recurrence equation.
        (3) Solve the recurrence equation to get the generator.

    Examples
    ========

    >>> from sympy import sin, exp
    >>> from sympy.series.formal import findgen
    >>> from sympy.abc import x, k
    >>> findgen(x, sin(x), k)
    [((-1/4)**k/(RisingFactorial(3/2, k)*factorial(k)), 2*k + 1)]
    >>> findgen(x, exp(x), k)
    [(1/(factorial(k)), k)]

    References
    ==========

    .. [1] Formal Power Series - Dominik Gruntz, Wolfram Koepf
    .. [2] Power Series in Computer Algebra - Wolfram Koepf

    See Also
    ========

    simplDE, DEtoRE, solveRE
    """
    # Search upto 4th order DE
    # Good enough for most functions
    for order in range(0, 5):
        diff = function.diff(x, order)
        if diff.is_rational_function():
            try:
                gen = findgen_rational(x, diff, k)  # Integrate order times
                integral = lambda e, k: e if k == 0 else integral(e.integrate(x, conds="none"), k - 1)
                return [integral(c*x**e, order).as_coeff_exponent(x) for c, e in gen]
            except ValueError:
                break

    f = Function('f')
    r = Function('r')

    DE = simpleDE(x, function, f)
    RE = DEtoRE(DE, r, k)
    gen = solveRE(RE, r, k, function)

    return gen


from heapq import merge
from itertools import imap, islice, takewhile

class Lazyseries(Expr):
    """
    Represents an infinite series
    Use ``FormalSeries`` to create a Lazyseries object

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.series.formal import FormalSeries
    >>> from sympy.abc import x
    >>> s = FormalSeries(x, function=sin(x))
    >>> s.as_series(n=6)
    x - x**3/6 + x**5/120 + O(x**6)
    >>> t = FormalSeries(x, function=cos(x))
    >>> t.as_series(n=6)
    1 - x**2/2 + x**4/24 + O(x**6)
    >>> u = s * t
    >>> u.as_series(n=6)
    x - 2*x**3/3 + 2*x**5/15 + O(x**6)

    """
    def __init__(self, x, gen):
        self.x = x
        self.sym = x.atoms(Symbol).pop()
        self.gen = gen

    def __iter__(self):
        return self.gen.__iter__()

    def __getitem__(self, index):
        return self.gen.__getitem__(index)

    def as_series(self, n=6):
        sym = self.sym
        it = takewhile(lambda t: degree(t, sym) < n, self.gen)
        s = [t for t in it]
        return Add(*s) + C.Order(sym**n, sym)

    @property
    def is_empty(self):
        return self.gen.is_empty

    def inverse(self):
        return self.__class__(self.sym, Stream(Lazyseries.series_inverse(self.sym, self.gen)))

    @staticmethod
    def series_inverse(sym, self):
        def gen():
            yield S.One
            raise StopIteration
        for t in Lazyseries.series_div(sym, Stream(gen()), self):
            yield t

    def compose(self, other):
        if self.sym != other.sym:
            raise ValueError('Cannot compose series of different variables')
        return self.__class__(self.sym, Stream(Lazyseries.series_compose(self.sym, self.gen, other.gen)))

    @staticmethod
    def series_compose(sym, self, other):
        if self.is_empty:
            raise StopIteration
        yield self[0]
        self = Stream(imap(lambda t: t/sym, self.shift(1)))
        it = Lazyseries.series_mul(sym, other, Stream(Lazyseries.series_compose(sym, self, other)))
        for t in it:
            yield t

    def __add__(self, other):
        if self.sym != other.sym:
            raise ValueError('Cannot add series of different variables')
        sym = self.sym
        return self.__class__(sym, Stream(Lazyseries.series_add(sym, self.gen, other.gen)))

    @staticmethod
    def series_add(sym, *iterables):
        """ Custom add for series to produce terms in sorted order """
        # Wrap the original term around Term to define __lt__ method for merge
        class Term():
            def __init__(self, term):
                self.term = term
            def __lt__(self, other):
                return degree(self.term, sym) < degree(other.term, sym)
        input_iter = [imap(lambda t: Term(t), it) for it in iterables]
        return imap(lambda t: t.term, merge(*input_iter))

    def __mul__(self, other):
        if not isinstance(other, Lazyseries):
            return self.__class__(self.x, Stream(imap(lambda t: t*other, self.gen)))
        if self.sym != other.sym:
            raise ValueError('Cannot multiply series of different variables')
        return self.__class__(self.sym, Stream(Lazyseries.series_mul(self.sym, self.gen, other.gen)))

    @staticmethod
    def series_mul(sym, self, other):
        try:
            a, b = self[0], other[0]
        except IndexError:
            raise StopIteration
        yield a * b
        self, other = self.shift(1), other.shift(1)
        it = Stream(Lazyseries.series_add(sym,
                Stream(imap(lambda t: t*b, self)),
                Stream(imap(lambda t: t*a, other)),
                Lazyseries.series_mul(sym, self, other)))
        for t in it:
            yield t

    def __truediv__(self, other):
        if not isinstance(other, Lazyseries):
            return self.__class__(self.x, Stream(imap(lambda t: t/other, self.gen)))
        if self.sym != other.sym:
            raise ValueError('Cannot multiply series of different variables')
        return self.__class__(self.sym, Stream(Lazyseries.series_div(self.sym, self.gen, other.gen)))
    __div__ = __truediv__

    @staticmethod
    def series_div(sym, self, other):
        if self.is_empty:
            raise StopIteration
        if other.is_empty:
            raise ZeroDivisionError
        a, b = self[0], other[0]
        self, other = self.shift(1), other.shift(1)
        def gen():
            yield a/b
            q = imap(lambda t: t*S.NegativeOne/b, Lazyseries.series_mul(sym, other, Stream(gen())))
            for t in Lazyseries.series_add(sym, self, Stream(q)):
                yield t
        for t in gen():
            yield t
