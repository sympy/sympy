from __future__ import print_function, division

from sympy import C, S, collect, Function, Add, Mul, simplify, cancel, sympify, Basic
from sympy.core.symbol import Symbol, symbols, Dummy, Wild
from sympy.core.sympify import sympify
from sympy.core.relational import Eq
from sympy.core.compatibility import integer_types
from sympy.functions import sign
from sympy.functions.elementary.piecewise import Piecewise
from sympy.concrete import Sum, summation
from sympy.core import Add, Pow
from sympy.core.expr import Expr
from sympy.polys import roots, lcm, degree
from sympy.polys.partfrac import apart
from sympy.solvers import solve, rsolve


def findgen_seq(sequence):
    """
    Returns generator for a sequence in terms of k

    Examples
    ========

    >>> from sympy.series.formal import findgen_seq
    >>> from sympy.abc import x, k
    >>> f = findgen_seq((1, 2, 3))
    >>> f(k)
    (Piecewise((1, Mod(k, 3) == 0), (2, Mod(k, 3) == 1), (3, Mod(k, 3) == 2)), k)
    >>> f = findgen_seq([1])
    >>> f(k)
    (Piecewise((1, Mod(k, 1) == 0)), k)

    """
    l = len(sequence)
    def gen(k):
        cond = [(sequence[i], Eq(k%l, i)) for i in range(l)]
        return (Piecewise(*cond), k)
    return gen


def findgen_rational(function, x=None):
    """
    Returns the generator for a given rational function in x in terms of k

    Examples
    ========

    >>> from sympy.series.formal import findgen_rational
    >>> from sympy.abc import x, k
    >>> f = findgen_rational(1/(1-x), x)
    >>> f(k)
    [(1, k)]
    >>> f = findgen_rational(1/(1+x), x)
    >>> f(k)
    [(-(-1)**(-k - 1), k)]

    """
    if x is None:
        syms = function.atoms(C.Symbol)
        if len(syms) > 1:
            raise ValueError('x must be given for multivariate functions')
        x = syms.pop()

    def gen(k):
        g = []
        m = Wild('m', exclude=(0,))
        terms = Add.make_args(apart(function))
        for t in terms:
            n, d = t.as_numer_denom()
            if not d.has(x):
                c, e = t.as_coeff_exponent(x)
                g += [(c, e)]
            elif not d.match(x+m):
                raise ValueError('Expected denominator of type "x + a", got: %s' % d)
            else:
                d, j = d.as_base_exp()
                a = -d.as_coeff_add()[0]
                c, e = (-1)**j * n * C.binomial(j+k-1, k).rewrite(C.factorial) / a**(j+k), k
                g += [(c, e)]
        return g
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


def simpleDE(x, function, f, order=4):
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
    a = symbols('a:%d' % order)

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
    for k in range(2, order+1):
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


def DEtoRE(DE, r):
    """
    Converts a differential equation into a recurrence equation

    Parameters:
        DE: Function in which RE whill be expressed
        r: Argument of recurr function

    Examples
    ========

    >>> from sympy import Derivative, Function
    >>> from sympy.series.formal import DEtoRE
    >>> from sympy.abc import x, k
    >>> f, r = Function('f'), Function('r')
    >>> DE = -f(x) + Derivative(f(x), x)
    >>> DEtoRE(DE, r(k))
    (k + 1)*r(k + 1) - r(k)

    """
    RE = S.Zero
    DE = DE.expand()

    f = DE.atoms(Function).pop()
    x = f.atoms(Symbol).pop()
    k = r.args[0]

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
        RE += c * C.RisingFactorial(k+1-l, j) * r.subs(k, k+j-l)
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


def solveRE(RE, r, function):
    """
    Returns generator for a given recurrence equation

    Examples
    ========

    >>> from sympy import exp, Function
    >>> from sympy.series.formal import solveRE
    >>> from sympy.abc import x, k
    >>> r = Function('r')
    >>> RE = (k+1)*r(k+1) - r(k)
    >>> solveRE(RE, r(k), exp(x))
    [(1/(factorial(k)), k)]

    """
    RE = RE.expand().collect(r.func(Wild('m')))
    k = r.args[0]
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
    key = [r.subs(k, i) for i in range(0, len(terms)-1)]
    val = [function.diff(x, i).subs(x, 0)/C.factorial(i) for i in range(0, len(terms)-1)]
    init = dict(zip(key, val))

    sol = rsolve(RE, r, init)
    if sol:
        return [(sol, k)]

    raise NotImplementedError('Cannot solve %s' % RE)


def findgen(function, x=None):
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
    >>> f = findgen(sin(x), x)
    >>> f(k)
    [((-1/4)**k/(RisingFactorial(3/2, k)*factorial(k)), 2*k + 1)]
    >>> f = findgen(exp(x), x)
    >>> f(k)
    [(1/(factorial(k)), k)]

    References
    ==========

    .. [1] Formal Power Series - Dominik Gruntz, Wolfram Koepf
    .. [2] Power Series in Computer Algebra - Wolfram Koepf

    See Also
    ========

    simplDE, DEtoRE, solveRE
    """
    if x is None:
        syms = function.atoms(C.Symbol)
        if len(syms) > 1:
            raise ValueError('x must be given for multivariate functions')
        x = syms.pop()

    # Search upto 4th order DE
    # Good enough for most functions
    for order in range(0, 5):
        diff = function.diff(x, order)
        if diff.is_rational_function():
            try:
                gen = findgen_rational(diff, x)
                integral = lambda e, k: e if k == 0 else integral(e.integrate(x, conds="none"), k - 1)
                return lambda k: [integral(c*x**e, order).as_coeff_exponent(x) for c, e in gen(k)]
            except ValueError:
                break

    f = Function('f')
    r = Function('r')

    def gen(k):
        DE = simpleDE(x, function, f)
        RE = DEtoRE(DE, r(k))
        return solveRE(RE, r(k), function)

    return gen


class FormalSeries(Expr):
    """
    Represents an infinite series

    Kwargs
    ======

    function
    sequence
    generator

    Examples
    ========

    >>> from sympy import sin
    >>> from sympy.series.formal import FormalSeries
    >>> from sympy.abc import x
    >>> s = FormalSeries(x, function=sin(x))
    >>> s.as_series()
    x - x**3/6 + x**5/120 + O(x**6)
    >>> s = FormalSeries(x, sequence=(1, 2, 3))
    >>> s.as_series()
    1 + 2*x + 3*x**2 + x**3 + 2*x**4 + 3*x**5 + O(x**6)

    """
    def __new__(cls, x, *args, **kwargs):
        k = Dummy('k', integer=True)

        generator = kwargs.pop("generator", None)
        if generator:
            gen, k = generator

        sequence = kwargs.pop('sequence', None)
        if sequence:
            c, e = findgen_seq(sequence)(k)
            gen = c*x**e

        function = kwargs.pop('function', None)
        if function is not None:
            function = sympify(function)
            generator = findgen(function, x)
            gen = sum(c*x**e for c, e in generator(k))

        obj = Expr.__new__(cls, x, gen, k)
        obj.x, obj.gen, obj.k = x, gen, k
        return obj

    def __getitem__(self, index):
        k = self.k
        if isinstance(index, integer_types):
            if index < 0:
                raise ValueError("Argument must be greater than 0")
            terms = Add.make_args(self.gen)
            val = S.Zero
            for t in terms:
                if not t.has(self.k):
                    if index == 0:
                        val += t.doit()
                else:
                    val += t.subs(k, index).doit()
            return val
        elif isinstance(index, slice):
            if index.step == 0:
                raise ValueError("Step must not be 0")
            return [self[i] for i in xrange(index.start, index.stop, index.step or 1)]

    @property
    def free_symbols(self):
        return self.x.free_symbols

    def _eval_derivative(self, x):
        gen = self.gen.diff(x)
        return FormalSeries(self.x, generator=(gen, self.k))

    def _eval_as_leading_term(self, x):
        if x == self.x:
            return self[0]
        else:
            return self

    def __add__(self, other):
        if not isinstance(other, FormalSeries):
            gen = self.gen + other
            return FormalSeries(self.x, generator=(gen, self.k))
        if self.x != other.x:
            raise ValueError('Cannot add series of different variables')
        gen = self.gen + other.gen.subs(other.k, self.k)
        return FormalSeries(self.x, generator=(gen, self.k))

    def __sub__(self, other):
        if not isinstance(other, FormalSeries):
            gen = self.gen - other
            return FormalSeries(self.x, generator=(gen, self.k))
        if self.x != other.x:
            raise ValueError('Cannot subtract series of different variables')
        gen = self.gen - other.gen.subs(other.k, self.k)
        return FormalSeries(self.x, generator=(gen, self.k))

    def __mul__(self, other):
        if not isinstance(other, FormalSeries):
            gen = self.gen * other
            return FormalSeries(self.x, generator=(gen, self.k))
        if self.x != other.x:
            raise ValueError('Cannot multiply series of different variables')
        j = Dummy('j', integer=True)
        gen = self.gen.subs(self.k, j) * other.gen.subs(other.k, self.k - j)
        return FormalSeries(self.x, generator=(Sum(gen, (j, 0, self.k)), self.k))

    def as_series(self, n=6):
        x = self.x
        s = []
        i = 0
        if not self.gen.has(self.k):
            return self[0]
        while True:
            if self[i].has(self.x) and (self[i] + C.Order(x**n)).is_Order:
                break
            s.append(self[i])
            i += 1
        return Add(*s) + C.Order(x**n)

