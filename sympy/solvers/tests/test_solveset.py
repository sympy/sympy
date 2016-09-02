from sympy import (
    Abs, Dummy, Eq, Gt, Function, Mod,
    LambertW, Piecewise, Poly, Rational, S, Symbol, Matrix,
    asin, acos, acsc, asec, atan, atanh, cos, csc, erf, erfinv, erfc, erfcinv,
    exp, log, pi, sin, sinh, sec, sqrt, symbols,
    tan, tanh, atan2, arg,
    Lambda, imageset, cot, acot, I, EmptySet, Union, E, Interval, Intersection,
    oo)

from sympy.core.function import nfloat
from sympy.core.relational import Unequality as Ne
from sympy.functions.elementary.complexes import im, re
from sympy.functions.elementary.hyperbolic import HyperbolicFunction
from sympy.functions.elementary.trigonometric import TrigonometricFunction

from sympy.polys.rootoftools import CRootOf

from sympy.sets import (FiniteSet, ConditionSet, Complement, ImageSet)

from sympy.utilities.pytest import XFAIL, raises, skip, slow, SKIP
from sympy.utilities.randtest import verify_numerically as tn
from sympy.physics.units import cm
from sympy.core.containers import Dict

from sympy.solvers.solveset import (
    solveset_real, domain_check, solveset_complex, linear_eq_to_matrix,
    linsolve, _is_function_class_equation, invert_real, invert_complex,
    solveset, solve_decomposition, substitution, nonlinsolve, solvify)

a = Symbol('a', real=True)
b = Symbol('b', real=True)
c = Symbol('c', real=True)
x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
q = Symbol('q', real=True)
m = Symbol('m', real=True)
n = Symbol('n', real=True)





def test_issue_11534():
    # eq and eq2 should give the same solution as a Complement
    eq = -y + x/sqrt(-x**2 + 1)
    eq2 = -y**2 + x**2/(-x**2 + 1)
    soln = Complement(FiniteSet(-y/sqrt(y**2 + 1), y/sqrt(y**2 + 1)), FiniteSet(-1, 1))
    assert solveset(eq, x, S.Reals) == soln
    assert solveset(eq2, x, S.Reals) == soln
