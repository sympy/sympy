from sympy.core.symbol import IntConst, Symbol, symbols
from sympy.core.relational import Eq
from sympy.core.function import Function
from sympy.core import I
from sympy.functions import exp, sin, cos, sinh, cosh
from sympy.core.multidimensional import vectorize
from sympy.utilities.pytest import XFAIL, raises
from sympy.core.compatibility import permutations


@vectorize(0)
def IntConst_to_Symbol(expr):
    """Copied from test_ode."""
    symbols_starting_with_C = [c for c in expr.free_symbols
                                  if c.name.startswith('C')
                                  and not isinstance(c, IntConst)]
    if symbols_starting_with_C:
        raise ValueError('integrations constant is not IntConst instance')
    old_constants = list(expr.atoms(IntConst))
    new_constants = [Symbol('C'+c.name) for c in old_constants]
    return expr.subs(zip(old_constants, new_constants))

import sympy
def dsolve(*args, **kwargs):
    """Copied (and simplified) from test_ode."""
    ret = sympy.dsolve(*args, **kwargs)
    return IntConst_to_Symbol(ret)

func = Function('f')
gunc = Function('g')
hunc = Function('h')
x = Symbol('x')
f = func(x)
f_ = f.diff(x)
f__ = f_.diff(x)
g = gunc(x)
g_ = g.diff(x)
g__ = g_.diff(x)
h = hunc(x)
h_ = h.diff(x)
h__ = h_.diff(x)

C1, C2, C3, C4, C5, C6 = symbols('C1:7')


def test_linear_with_init_cond():
    sys = [f_+f]
    sol = dsolve(sys, [f])
    assert sol[0] in [Eq(f, C1*exp(-x))]

    sys = [f_+f, func(0)-2]
    sol = dsolve(sys, [f])
    assert sol[0] in [Eq(f, 2*exp(-x))]

    sys = [f_+f, func(2)-3]
    sol = dsolve(sys, [f])
    assert sol[0] in [Eq(f, 3*exp(2)*exp(-x))]

    sys = [f__+f]
    sol = set(dsolve(sys, [f]))
    assert sol in [
        set([Eq(f, C1*sin(x) + C2*cos(x))]),
        set([Eq(f, C1*cos(x) + C2*sin(x))]),
        ]

    sys = [g-f_, f+g_]
    sol = set(dsolve(sys, [f, g]))
    assert sol in [
        set([Eq(g, -I*C1*exp(-I*x) + I*C2*exp(I*x)),
             Eq(f, C1*exp(-I*x) + C2*exp(I*x))]),
        set([Eq(g, C1*exp(-I*x) + C2*exp(I*x)),
             Eq(f, I*C1*exp(-I*x) - I*C2*exp(I*x))]),
        set([Eq(g, I*C1*exp(I*x) - I*C2*exp(-I*x)),
             Eq(f, C1*exp(I*x) + C2*exp(-I*x))]),
        set([Eq(g, C1*exp(I*x) + C2*exp(-I*x)),
             Eq(f, -I*C1*exp(I*x) + I*C2*exp(-I*x))]),
        ]

    sys = [g-f_, f+g_, func(0)-1, gunc(0)]
    sol = set(dsolve(sys, [f, g]))
    assert sol == set([Eq(g, I*exp(I*x)/2 - I*exp(-I*x)/2),
                       Eq(f, exp(I*x)/2 + exp(-I*x)/2)])


def test_initial_conditions_using_Subs():
    sys = [f__+f, func(0)-1, f_.subs(x,0)]
    sol = dsolve(sys, [f])
    assert sol[0] == Eq(f, cos(x))


def test_already_triangular():
    # Testing for [ 1 1 ]
    #             [ 0 1 ]
    # This test is separate from test_not_diagonal, because systems that are
    # triagonal are never sent for diagonalization.
    # See also: test_nonlinear_iteratively_solvable
    sys = [f+g-f_, g-g_]
    sol = set(dsolve(sys, [f,g]))
    possible_sols = [
            set([Eq(g(x), C1*exp(x)), Eq(f(x), (C2 + C1*x)*exp(x))]),
            set([Eq(g(x), C2*exp(x)), Eq(f(x), (C1 + C2*x)*exp(x))]),
            ]
    assert sol in possible_sols


def test_not_diagonal():
    # Testing for [ 3 1 ]
    #             [-1 1 ]
    sys = [3*f+g-f_, -f+g-g_]
    #TODO sol = set(dsolve(sys, [f,g]))
    #TODO sol_f = (C1(1+sqrt(3)) + exp(sqrt(3)*x)*((-1+sqrt(3))*C1 - C2) + C2)/2/sqrt(3)/exp((-1+sqrt(3))*x/2)
    #TODO sol_g = exp(x/2)*(sqrt(3)*C2*cosh(sqrt(3)*x/2) + (-2*C1+C2)*sinh(sqrt(3)*x/2))/sqrt(3)
    #TODO assert sol == set([sol_f, sol_g])
    raises(NotImplementedError, lambda : dsolve(sys, [f, g]))


@XFAIL
def test_not_homogeneous():
    # Testing for [ 3 1 ]
    #             [-1 1 ]
    sys = [3*f+g-f_+sin(x), -f+g-g_-cos(x)]
    raises(NotImplementedError, lambda : dsolve(sys, [f, g]))


def test_nonlinear_separable():
    sys = [f_+f**2, g_+h, h_-g]
    sol = set(dsolve(sys, [f, g, h]))
    sol_categories = [
        set([Eq(g, I*C1*exp(I*x) - I*C2*exp(-I*x)),
             Eq(f, 1/(C3 + x)),
             Eq(h, C1*exp(I*x) + C2*exp(-I*x))]),
        set([Eq(g, -I*C1*exp(-I*x) + I*C2*exp(I*x)),
             Eq(f, 1/(C3 + x)),
             Eq(h, C1*exp(-I*x) + C2*exp(I*x))]),
        set([Eq(g, C1*exp(I*x) + C2*exp(-I*x)),
             Eq(f, 1/(C3 + x)),
             Eq(h, -I*C1*exp(I*x) + I*C2*exp(-I*x))]),
        set([Eq(g, C1*exp(-I*x) + C2*exp(I*x)),
             Eq(f, 1/(C3 + x)),
             Eq(h, I*C1*exp(-I*x) - I*C2*exp(I*x))]),
        ]
    perms = list(permutations([C1, C2, C3]))
    all_sols =  [set([eq.subs(zip(perms[0], p), simultaneous=True)
                       for eq in s])
                       for p in perms for s in sol_categories]
    assert sol in all_sols


def test_nonlinear_iteratively_solvable():
    sys = [f_+f**2, g_*f-1]
    sol = set(dsolve(sys, [f,g]))
    possible_sols = [
            set([Eq(g(x), C1 + C2*x + x**2/2), Eq(f(x), 1/(C2 + x))]),
            set([Eq(g(x), C2 + C1*x + x**2/2), Eq(f(x), 1/(C1 + x))]),
            ]
    assert sol in possible_sols

    sys = [f_+f**2, g_*f-1, h_*f-g]
    sol = set(dsolve(sys, [f,g,h]))
    sol_categories = [
            set([Eq(g(x), C1 + C2*x + x**2/2),
                 Eq(f(x), 1/(C2 + x)),
                 Eq(h(x), C1*C2*x + C2*x**3/2 + C3 + x**4/8 + x**2*(C1/2 + C2**2/2))]),
            ]
    perms = list(permutations([C1, C2, C3]))
    all_sols =  [set([eq.subs(zip(perms[0], p), simultaneous=True)
                       for eq in s])
                       for p in perms for s in sol_categories]
    assert sol in all_sols


@XFAIL
def test_nonlinear_not_separable():
    sys = [f_*g_+f, f_+g]
    raises(NotImplementedError, lambda : dsolve(sys, [f, g]))


@XFAIL
def test_overdetermined_consistent():
    sys = [f_+f, f__+f]
    sol = dsolve(sys, [f])
    assert sol[0] == Eq(f, 0)

    sys = [f_+g, g_+f, g_-f]
    sol = set(dsolve(sys, [f, g]))
    assert sol == set([Eq(f, 0), Eq(g, 0)])

    sys = [f_+f+1, f__+f+1, f+1]
    sol = dsolve(sys, [f])
    assert sol[0] == Eq(f, -1)


def test_overdetermined_inconsistent():
    sys = [f_+h, f_+f+g, h_+g]
    raises(ValueError, lambda : dsolve(sys, [f]))

    sys = [f+1, f__+f]
    raises(ValueError, lambda : dsolve(sys, [f]))


@XFAIL
def test_implicit_solution():
    eq = f/x*cos(f/x) - (x/f*sin(f/x) + cos(f/x))*f_
    sol_f = Eq(f*sin(f/x), C1)
    sol = set(dsolve([eq, f+g_], [f,g]))
    assert sol == set([sol_f, f+g_])


@XFAIL
def test_multiple_solutions():
    eq = f*f_-1
    sol = dsolve([eq, f_-g], [f,g])
    sol_f = [Eq(f, -sqrt(2)*sqrt(C1 + x)),     Eq(f, sqrt(2)*sqrt(C1 + x))]
    sol_g = [Eq(g, -sqrt(2)/(2*sqrt(C1 + x))), Eq(g, sqrt(2)/(2*sqrt(C1 + x)))]
    assert sol == [set(sol_f[0], sol_g[0]), set(sol_f[1], sol_g[1])]


def test_inconsistent_init_conds():
    sys = [f__+f, func(0)-1, f.subs(x,0)]
    raises(ValueError, lambda : dsolve(sys, [f]))
