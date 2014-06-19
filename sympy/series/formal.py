from __future__ import print_function, division

from sympy import C, S, symbols, Symbol, collect, Function, Add, Mul, simplify
from sympy.core.sympify import sympify
from sympy.core.relational import Eq
from sympy.functions import sign
from sympy.functions.elementary.piecewise import Piecewise
from sympy.concrete import Sum
from sympy.core import Add, Pow
from sympy.core.expr import Expr
from sympy.polys import lcm
from sympy.polys.partfrac import apart
from sympy.solvers import solve, rsolve

def FormalSeries(x=None, n=0, *args, **kwargs):
    """
    Finds generator for given function and returns a ``Stream`` class.

    Generator is a tuple - (c, e, k)
        c: coefficient of the kth term
        e: exponent of the kth term
        k: symbol used for representing kth term

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
    k = C.Dummy('k', integer=True)
    generator = kwargs.pop("generator", None)

    sequence = kwargs.pop('sequence', None)
    if sequence:
        generator = findgen_seq(x, sequence, k)

    function = kwargs.pop('function', None)
    if function:
        generator = findgen(x, function, k)

    c, e, k = generator
    nth_term = (c.subs(k, n), e.subs(k, n))
    compute_tail = lambda: FormalSeries(x, n=n+1, generator=generator)

    return Stream(x, nth_term, compute_tail)


def findgen_seq(x, sequence, k):
    """ Returns generator for a sequence in terms of k """
    l = len(sequence)
    cond = [(sequence[i], Eq(k%l, i)) for i in range(l)]
    return (Piecewise(*cond), k, k)


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
    coeff = S.Zero
    terms = apart(function).as_ordered_terms()
    for t in terms:
        c, d = t.as_numer_denom()
        d, j = d.as_base_exp()
        a = -d.as_coeff_add()[0]
        coeff += (-1)**j * c * C.binomial(j+k-1, k).rewrite(C.factorial) / a**(j+k)
    return (coeff, k, k)


def rational_independent(x, terms):
    """ Returns list of rationally independent terms """
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
        return DE.subs(a[0], sol[0])

    for k in range(2, 5):
        eq, DE = makeDE(k)
        eq = eq.expand()
        terms = [t for t in eq.as_ordered_terms()]
        ind = rational_independent(x, terms)
        if len(ind) == k:
            sol = solve(ind, a, dict=True)
            if sol:
                for key, val in sol[0].iteritems():
                    DE = DE.subs(key, val)
            DE = DE.expand()
            DE = simplify(DE)
            return DE.as_numer_denom()[0]

    raise NotImplementedError('Cannot find simple differential equation for %s' % function)


def DEtoRE(x, DE, recurr, n):
    """
    Converts a differential equation into a recurrence equation

    Parameters:
        recurr: Function in which RE whill be expressed
        n: Argument of recurr function

    Examples
    ========

    >>> recurr = Function('r')
    >>> DE = -f(x) + Derivative(f(x), x)
    >>> DEtoRE(x, DE, f, r, n)
    (n + 1)*r(n + 1) - r(n)

    """
    RE = S.Zero
    DE = DE.expand()

    for t in DE.as_ordered_terms():
        a = t.as_ordered_factors()
        a = [Mul(*a[:-1]), a[-1]]
        if len(a) == 1:
            l = S.Zero
            # Length of a Derivate args determine its order
            # XXX: Is there a better method?
            if a[0].is_Derivative:
                j = len(a[0].args) - 1
            else:
                j = 0
            RE += C.RisingFactorial(n+1-l, j) * recurr(n+j-l)
        else:
            c, v = a[0].as_independent(x)
            if v.has(x):
                l = v.as_base_exp()[1]
            else:
                l = 0
            if a[1].is_Derivative:
                j = len(a[1].args) - 1
            else:
                j = 0
            RE += c * C.RisingFactorial(n+1-l, j) * recurr(n+j-l)
    return RE


def solveRE(x, RE, recurr, n, function):
    """
    Returns generator for a given recurrence equation
    """
    keys = [recurr(i) for i in range(0, 5)]
    values = [function.diff(x, i).subs(x, 0)/C.factorial(i) for i in range(0, 5)]
    init = dict(zip(keys, values))
    sol = rsolve(RE, recurr(n), init)
    return sol


def findgen(x, function, k):
    """
    Returns the generator for a given function in x
    """
    for order in range(0, 5):
        diff = function.diff(x, order)
        if diff.is_rational_function():
            return findgen_rational(x, function)  # Integrate k times

    func = Function('f')
    recurr = Function('r')
    n = C.Dummy('n', integer=True)

    DE = simpleDE(x, function, func)
    RE = DEtoRE(x, DE, recurr, n)
    return solveRE(x, RE, recurr, n, function)


class Stream(Expr):
    """
    Stream class to represent infinite series.

    Use ``FormalSeries`` to generate a stream object representing infinite series of a function.

    Attributes:
        head: first term of the series represented as a tuple (c, e)
        compute_tail: function to generate rest of the series, returns a Stream class
        tail: stream object representing rest of the series generated by compute_tail
    """
    def __init__(self, x, head, compute_tail=lambda: None):
        self.head = head
        self._compute_tail = compute_tail
        self.x = x

    def __add__(self, other):
        if self is None:
            return other
        if other is None:
            return self
        if self.x != other.x:
            raise ValueError
        c1, e1 = self.head
        c2, e2 = other.head
        s = sign(e2 - e1)
        if s == 1:
            return Stream(self.x, self.head, lambda: self.tail + other)
        if s == -1:
            return Stream(self.x, other.head, lambda: self + other.tail)
        if s == 0:
            return Stream(self.x, (c1 + c2, e1), lambda: self.tail + other.tail)
        raise NotImplementedError('Result depends on sign of %s' % s)

    def __sub__(self, other):
        if self is None:
            return other.scale(-1)
        if other is None:
            return self
        if self.x != other.x:
            raise ValueError
        c1, e1 = self.head
        c2, e2 = other.head
        s = sign(e2 - e1)
        if s == 1:
            return Stream(self.x, self.head, lambda: self.tail - other)
        if s == -1:
            return Stream(self.x, other.head, lambda: self - other.tail)
        if s == 0:
            return Stream(self.x, (c1 - c2, e1), lambda: self.tail - other.tail)
        raise NotImplementedError('Result depends on sign of %s' % s)

    def __mul__(self, other):
        if self is None or other is None:
            return None
        if self.x != other.x:
            raise ValueError
        c1, e1 = self.head
        c2, e2 = other.head
        compute_tail = lambda: other.tail.shift(e1).scale(c1) + self.tail.shift(e2).scale(c2) + self.tail*other.tail
        return Stream(self.x, (c1*c2, e1+e2), compute_tail)

    def __truediv__(self, other):
        if self is None:
            return None
        if other is None:
            raise ZeroDivisionError('Division by zero')
        if self.x != other.x:
            raise ValueError
        c1, e1 = self.head
        c2, e2 = other.head
        q = lambda: Stream(self.x, (c1/c2, e1-e2), lambda: (self.tail + (other.tail * q()).scale(-1)).scale(S.One/c2).shift(-e2))
        return q()

    def compose(self, other):
        c1, e1 = self.head
        c2, e2 = other.head
        return Stream(self.x, self.head, lambda: other.tail * self.tail.shift(-1).compose(other))

    def invert(self):
        c1, e1 = self.head
        r = lambda: Stream(self.x, (S.One/c1, -e1), lambda: self.tail.scale(-1) * r())
        return r()

    def scale(self, n):
        c1, e1 = self.head
        return Stream(self.x, (c1*n, e1), lambda: self.tail.scale(n))

    def shift(self, n):
        c1, e1 = self.head
        return Stream(self.x, (c1, e1+n), lambda: self.tail.shift(n))

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

    def as_series(self, n=6):
        x = self.x
        it = self.tail
        s = self.term
        while it != None and n > 0:
            s = s + it.term
            it = it.tail
            n = n - 1
        return s
