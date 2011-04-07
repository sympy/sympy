"""
Tests for solving systems of differential equations
"""
from sympy.core import S
from sympy.solvers.sysode import ode_1st_linear_system
from sympy.core.relational import Equality, Eq
from sympy.core.symbol import Symbol

from sympy.functions import exp, sin, cos, log
from sympy.abc import x,y,z

t = Symbol('t')
C0 = Symbol('C0')
C1 = Symbol('C1')
C2 = Symbol('C2')
C3 = Symbol('C3')
C4 = Symbol('C4')

def test_ode_1st_linear_system():
    sys_soln1 = ode_1st_linear_system([x,2*y,3*z], t, x, y, z)
    assert sys_soln1 == [Eq(x,C0*exp(t)), Eq(y,C1*exp(2*t)), Eq(z,C2*exp(3*t))]
    sys_soln2 = ode_1st_linear_system([x + sin(t),y + cos(t)], t, x, y)
    assert sys_soln2 == [Eq(x,(C0 - exp(-t)*sin(t)/2 - cos(t)*exp(-t)/2)*exp(t)),\
                         Eq(y,(C1 + exp(-t)*sin(t)/2 - cos(t)*exp(-t)/2)*exp(t))]
