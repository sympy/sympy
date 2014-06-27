from __future__ import print_function, division

from sympy import C, S, collect, Function, Add, Mul, simplify, cancel, sympify
from sympy.core.symbol import Symbol, symbols, Dummy, Wild
from sympy.core.sympify import sympify
from sympy.core.relational import Eq
from sympy.functions import sign
from sympy.functions.elementary.piecewise import Piecewise
from sympy.concrete import Sum
from sympy.core import Add, Pow
from sympy.core.expr import Expr
from sympy.polys import roots, lcm
from sympy.polys.partfrac import apart
from sympy.solvers import solve, rsolve


def FormalSeries(x=None, n=0, *args, **kwargs):
    """
    Finds generator for given function and returns a ``Stream`` class.

    ``n`` represents the starting term for the infinite series

    Kwargs:
        function:
        sequence:
        generator:

    Examples
    ========

    >>> from sympy.series.formal import FormalSeries
    >>> from sympy.abc import x
    >>> FormalSeries(x, function=1/(1-x))
    Stream(...)
    >>> FormalSeries(x, sequence=(1, 2, 3))
    Stream(...)

    """
    k = Dummy('k')
    generator = kwargs.pop("generator", None)
    if generator:
        c, e, k = generator

    sequence = kwargs.pop('sequence', None)
    if sequence:
        c, e = findgen_seq(x, sequence, k)

    function = kwargs.pop('function', None)
    if function:
        generator = findgen(x, function, k)
        return sum([FormalSeries(x, generator=(c, e, k)) for c, e in generator])

    # If generator is constant like x**2 or 1
    # Then all except first terms are zero
    if not e.has(k):
        general_term = c*x**e
        if n > 0:
            nth_term = (S.Zero, S.Infinity)
        else:
            nth_term = (c, e)
    else:
        general_term = Sum(c*x**e, (k, S.Zero, S.Infinity))
        nth_term = (c.subs(k, n), e.subs(k, n))
    compute_tail = lambda: FormalSeries(x, n=n+1, generator=(c, e, k))

    return Stream(x, nth_term, general_term, compute_tail)


def findgen_seq(x, sequence, k):
    """
    Returns generator for a sequence in terms of k

    Examples
    ========

    >>> findgen_seq(x, (1, 2, 3), k)
    (Piecewise((1, Eq(k%3, 0)), (2, Eq(k%3, 1)), (3, Eq(k%3, 2))), k, k)
    >>> findgen_seq(x, [1], k)
    (1, k, k)

    """
    l = len(sequence)
    cond = [(sequence[i], Eq(k%l, i)) for i in range(l)]
    return (Piecewise(*cond), k)


def findgen_rational(x, function, k):
    """
    Returns the generator for a given rational function in x in terms of k

    Examples
    ========

    >>> findgen_rational(x, 1/(1-x), k)
    (1, k, k)
    >>> findgen_rational(x, 1/(1+x), k)
    (-(-1)**(-k - 1), k, k)

    """
    gen = []
    m = Wild('m')
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

    >>> rational_independent(x, [sin(x), cos(x)])
    [sin(x), cos(x)]
    >>> rational_independent(x, [sin(x), cos(x), x*sin(x)])
    [x*sin(x) + sin(x), cos(x)]

    """
    # XXX: Not sure if this is most efficient way
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


def simpleDE(x, function, func):
    """
    Converts a function into a simple differential equation.

    A simple differential equation is a linear homogeneous differential equation
    with rational coefficients.

    Examples:
    ========

    >>> f = Function('f')
    >>> simpleDE(x, exp(x), f(x))
    -f(x) + Derivative(f(x), x, x)
    >>> simpleDE(x, sin(x), f(x))
    f(x) + Derivative(f(x), x, x)

    """
    a = symbols('a:4')

    makeDE = lambda k: (function.diff(x, k) + \
            Add(*[a[i]*function.diff(x, i) for i in range(0, k)]), func(x).diff(x, k) +
            Add(*[a[i]*func(x).diff(x, i) for i in range(0, k)]))

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
        recurr: Function in which RE whill be expressed
        n: Argument of recurr function

    Examples
    ========

    >>> r = Function('r')
    >>> DE = -f(x) + Derivative(f(x), x)
    >>> DEtoRE(DE, r, k)
    (k + 1)*r(k + 1) - r(k)

    """
    RE = S.Zero
    DE = DE.expand()

    f = DE.atoms(Function).pop()
    x = DE.atoms(Symbol).pop()

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

    >>> RE = (k+1)*r(k+1) - r(k)
    >>> solveRE(x, RE, r, k, exp(x))
    (1/k!, k)

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
    >>> findgen(x, sin(x), k)


    References
    ==========

    .. [1] Formal Power Series - Dominik Gruntz, Wolfram Koepf
    .. [2] Power Series in Computer Algebra - Wolfram Koepf

    See Also
    ========

    simplDE, DEtoRE, solveRE
    """
    function = sympify(function)

    if x.is_positive is not True:
        xpos = Dummy('x', positive=True)
        function = function.subs(x, xpos)
        return findgen(xpos, function, k)

    # Search upto 4th order DE
    # Good enough for most functions
    for order in range(0, 5):
        diff = function.diff(x, order)
        if diff.is_rational_function():
            try:
                return findgen_rational(x, diff, k)  # Integrate order times
            except ValueError:
                break

    f = Function('f')
    r = Function('r')

    DE = simpleDE(x, function, f)
    RE = DEtoRE(DE, r, k)
    gen = solveRE(RE, r, k, function)

    return gen

class EmptyStream():
    """ Represents an empty stream """
    pass


class Stream(Expr):
    """
    Stream class to represent infinite series.

    Use ``FormalSeries`` to generate a stream object representing infinite series of a function.

    Attributes:
        head: first term of the series represented as a tuple (c, e)
        compute_tail: function to generate rest of the series, returns a Stream class
        tail: stream object representing rest of the series generated by compute_tail

    References
    ==========

    .. [1] Gruntz' Thesis pp. 100 Secion 7.1.2

    """
    def __init__(self, x, head, general_term, compute_tail=lambda: EmptyStream()):
        self.x = x
        self.head = head
        self.general_term = general_term
        self._compute_tail = compute_tail

    def __add__(self, other):
        if self is EmptyStream():
            return other
        if other is EmptyStream():
            return self
        if self.x != other.x:
            raise ValueError('Expected streams of same variable, got: %s, %s' % (self.x, other.x))

        c1, e1 = self.head
        c2, e2 = other.head

        gen1, gen2 = self.general_term, other.general_term
        if not gen1 or not gen2:
            gen = None
        else:
            gen = gen1 + gen2

        s = sign(e2 - e1)
        if s == 1:
            return Stream(self.x, self.head, gen, lambda: self.tail + other)
        if s == -1:
            return Stream(self.x, other.head, gen, lambda: self + other.tail)
        if s == 0:
            return Stream(self.x, (c1 + c2, e1), gen, lambda: self.tail + other.tail)
        raise NotImplementedError('Result depends on sign of %s' % s)

    def __sub__(self, other):
        if self is None:
            return other.scale(-1)
        if other is None:
            return self
        if self.x != other.x:
            raise ValueError('Expected streams of same variable, got: %s, %s' % (self.x, other.x))

        c1, e1 = self.head
        c2, e2 = other.head

        gen1, gen2 = self.general_term, other.general_term
        if not gen1 or not gen2:
            gen = None
        else:
            gen = gen1 - gen2

        s = sign(e2 - e1)
        if s == 1:
            return Stream(self.x, self.head, gen, lambda: self.tail - other)
        if s == -1:
            return Stream(self.x, other.head, gen, lambda: self - other.tail)
        if s == 0:
            return Stream(self.x, (c1 - c2, e1), gen, lambda: self.tail - other.tail)
        raise NotImplementedError('Result depends on sign of %s' % s)

    def __mul__(self, other):
        if self is EmptyStream() or other is EmptyStream():
            return EmptyStream()
        if self.x != other.x:
            raise ValueError('Expected streams of same variable, got: %s, %s' % (self.x, other.x))

        c1, e1 = self.head
        c2, e2 = other.head
        gen = None  # TODO: cauchy product formula

        compute_tail = lambda: other.tail.shift(e1).scale(c1) + self.tail.shift(e2).scale(c2) + self.tail*other.tail
        return Stream(self.x, (c1*c2, e1+e2), gen, compute_tail)

    def __truediv__(self, other):
        if self is EmptyStream():
            return EmptyStream()
        if other is EmptyStream():
            raise ZeroDivisionError('Division by zero')
        if self.x != other.x:
            raise ValueError('Expected streams of same variable, got: %s, %s' % (self.x, other.x))

        c1, e1 = self.head
        c2, e2 = other.head
        gen = None

        q = lambda: Stream(self.x, (c1/c2, e1-e2), gen, lambda: (self.tail + (other.tail * q()).scale(-1)).scale(S.One/c2).shift(-e2))
        return q()

    def compose(self, other):
        c1, e1 = self.head
        c2, e2 = other.head
        gen = None
        return Stream(self.x, self.head, gen, lambda: other * self.tail.shift(-1).compose(other))

    def invert(self):
        c1, e1 = self.head
        gen = None
        r = lambda: Stream(self.x, (S.One/c1, -e1), gen, lambda: self.tail.scale(-1) * r())
        return r()

    def scale(self, n):
        c1, e1 = self.head

        if not self.general_term:
            gen = None
        else:
            t, limits = self.general_term.args
            t = t*n
            gen = Sum(t, limits)

        return Stream(self.x, (c1*n, e1), gen, lambda: self.tail.scale(n))

    def shift(self, n):
        x = self.x
        c1, e1 = self.head

        if not self.general_term:
            gen = None
        else:
            t, limits = self.general_term.args
            t = t*x**n
            gen = Sum(t, limits)

        return Stream(self.x, (c1, e1+n), gen, lambda: self.tail.shift(n))

    @property
    def tail(self):
        if self._compute_tail is not None:
            self._tail = self._compute_tail()
            self._compute_tail = None
        return self._tail

    @property
    def term(self):
        coeff, expo = self.head
        x = self.x
        return coeff * x**expo

    def as_summation(self):
        if self.general_term:
            return self.general_term
        raise ValueError

    def as_series(self, n=6):
        x = self.x
        s = self.term
        it = self.tail
        # Go for 9 more terms until you hit the limit
        # Or you get an empty stream
        for i in range(n + 9):
            c, e = it.head
            if it == EmptyStream() or e > n:
                break
            s += it.term
            it = it.tail
        return s + C.Order(x**n, x)
