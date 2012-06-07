# TODO this file is staged but no real tests are added
from sympy import ode_system, Symbol, Function

func = Function('f')
gunc = Function('g')
x = Symbol('x')
f = func(x)
f_ = f.diff(x)
f__ = f_.diff(x)
g = gunc(x)
g_ = g.diff(x)
g__ = g_.diff(x)


def test_linear_with_init_cond():
    sys = [f_+f]
    sol = ode_system(sys, [f])

    sys = [f_+f, func(0)-2]
    sol = ode_system(sys, [f])

    sys = [f_+f, func(2)-3]
    sol = ode_system(sys, [f])

    sys = [f__+f]
    sol = ode_system(sys, [f])

    sys = [g-f_, f+g_] # The same as [f__+f] just substitute g_ = f__
    sol = ode_system(sys, [f, g])

    sys = [g-f_, f+g_, func(0)-1, gunc(0)]
    sol = ode_system(sys, [f, g])
