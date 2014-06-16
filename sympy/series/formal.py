from __future__ import print_function, division

from sympy import C, S, symbols, Symbol, collect, Function, Add
from sympy.core.sympify import sympify
from sympy.core.relational import Eq
from sympy.functions import sign
from sympy.functions.elementary.piecewise import Piecewise
from sympy.concrete import Sum
from sympy.core import Add, Pow
from sympy.core.expr import Expr
from sympy.polys.partfrac import apart
from sympy.solvers import solve, rsolve

import gc, traceback

def FormalSeries(x=None, n=0, *args, **kwargs):
    generator = kwargs.pop("generator", None)

    formula = kwargs.pop("formula", None)
    if formula:
        k = Symbol('k', integer=True)
        coeff = formula[0]
        generator = ((coeff, k), k)

    sequence = kwargs.pop('sequence', None)
    if sequence:
        k = Symbol('k', integer=True)
        l = len(sequence)
        cond = [(sequence[i], Eq(k%l, i)) for i in range(l)]
        generator = ((Piecewise(*cond), k), k)

    function = kwargs.pop('function', None)
    if function:
        generator = findGen(x, function)

    def compute_tail():
        return FormalSeries(x, n=n+1, generator=generator)

    coeff, expo = generator[0]
    k = generator[1]
    nth_term = (coeff.subs(k, n), expo.subs(k, n))

    return Stream(x, nth_term, compute_tail)


def findGenR(x, function):
    """
    Finds the generator for a given rational function in x
    """
    coeff = S.Zero
    k = Symbol('k', integer=True)
    terms = apart(function).as_ordered_terms()
    for t in terms:
        c, d = t.as_numer_denom()
        d, j = d.as_base_exp()
        a = -d.as_coeff_add()[0]
        coeff += (-1)**j * c * C.binomial(j+k-1, k).rewrite(C.factorial) / a**(j+k)
    return ((coeff, k), k)


def simpleDE(x, function):
    """
    Converts a function into a simple differential equation
    """
    def makeDE(x, function, func, k, a):
        eq = function.diff(x, k)
        DE = func.diff(x, k)
        for i in range(0, k):
            eq += a[i] * function.diff(x, i)
            DE += a[i] * func.diff(x, i)
        return eq, DE

    def independent(terms, k):
        terms = list(set(terms))  # Remove repeating terms
        if len(terms) == k:
            return True
        return False

    a = symbols('a:4')
    func = Function('f')(x)

    # Solve for case k=1
    eq, DE = makeDE(x, function, func, S.One, a)
    sol = solve(eq, a[0])
    if sol and sol[0].is_rational_function():
        return DE.subs(a[0], sol[0]), func

    for k in range(2, 5):
        eq, DE = makeDE(x, function, func, k, a)
        terms = [t.as_independent(x)[1] for t in eq.as_ordered_terms()]
        if independent(terms, k):
            eq = collect(eq, [function.diff(x, i) for i in range(0, k+1)])
            coeff = [t.as_independent(x)[0] for t in eq.as_ordered_terms()]
            sol = solve(coeff, a, dict=True)
            if sol:
                for key, val in sol[0].iteritems():
                    DE = DE.subs(key, val)
            return DE, func

    return DE, func


def DEtoRE(x, DE, func):
    def findOrder(x, derivative, func):
        order = S.Zero
        while derivative != func:
            order += 1
            derivative = derivative.integrate(x)
        return order

    k = Symbol('k', integer=True)
    r = Function('r')
    RE = S.Zero

    DE = DE.expand()
    for t in DE.as_ordered_terms():
        a = t.as_ordered_factors()
        if len(a) == 1:
            l = S.Zero
            j = findOrder(x, a[0], func)
            RE += C.RisingFactorial(k+1-l, j) * r(k+j-l)
        else:
            l = a[0].as_base_exp()[1]
            j = findOrder(a[1])
            RE += C.RisingFactorial(k+1-l, j) * r(k+j-l)
    return RE, r, k


def solveRE(x, RE, func, k):
    sol = rsolve(RE, func(k))
    print (RE, sol)


def findGen(x, function):
    """
    Finds the generator for a given function in x
    """
    for k in range(0, 5):
        diff = function.diff(x, k)
        if diff.is_rational_function():
            return findGenR(x, function)  # Integrate k times

    DE, func = simpleDE(x, function)
    RE, func, k = DEtoRE(x, DE, func)
    generator = solveRE(x, RE, func, k)

    k = Symbol('k', integer=True)
    return ((S.One/C.factorial(k), k), k)


class Stream(Expr):
    def __init__(self, x, head, compute_tail=lambda: None):
        self.head = head
        self._compute_tail = compute_tail
        self.x = x

    def __add__(self, other):
        s = sign(other.head[1] - self.head[1])
        if s == 1:
            return Stream(self.x, self.head, lambda: self.tail + other)
        if s == -1:
            return Stream(self.x, other.head, lambda: self + other.tail)
        h = (self.head[0] + other.head[0], self.head[1])
        return Stream(self.x, h, lambda: self.tail + other.tail)

    def __sub__(self, other):
        s = sign(other.head[1] - self.head[1])
        if s == 1:
            return Stream(self.x, self.head, lambda: self.tail - other)
        if s == -1:
            return Stream(self.x, other.head, lambda: self - other.tail)
        h = (self.head[0] - other.head[0], self.head[1])
        return Stream(self.x, h, lambda: self.tail - other.tail)

    def __mul__(self, other):
        c1, e1 = self.head
        c2, e2 = other.head
        coeff = c1 * c2
        expo = e1 + e2
        return Stream(self.x, (coeff, expo), lambda: other.tail.shift(e1).scale(c1) + self.tail.shift(e2).scale(c2) + self.tail * other.tail)

    def __div__(self, other):
        c1, e1 = self.head
        c2, e2 = other.head
        # q is the fixed-point operator
        def q():
            return Stream(self.x, (c1/c2, e1-e2), lambda: (self.tail + (other.tail * q()).scale(-1)).scale(S.One/c2).shift(-e2))
        return q()

    def scale(self, n):
        coeff = self.head[0] * n
        expo = self.head[1]
        return Stream(self.x, (coeff, expo), lambda: self.tail.scale(n))

    def shift(self, n):
        coeff = self.head[0]
        expo = self.head[1] + n
        return Stream(self.x, (coeff, expo), lambda: self.tail.shift(n))

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
            traceback.print_stack()
            s = s + it.term
            it = it.tail
            n = n - 1
        return s
