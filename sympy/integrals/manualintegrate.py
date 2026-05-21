"""Integration method that emulates by-hand techniques.

This module also provides functionality to get the steps used to evaluate a
particular integral, in the ``integral_steps`` function. This will return
nested ``Rule`` s representing the integration rules used.

Each ``Rule`` class represents a (maybe parametrized) integration rule, e.g.
``SinRule`` for integrating ``sin(x)`` and ``ReciprocalSqrtQuadraticRule``
for integrating ``1/sqrt(a+b*x+c*x**2)``. The ``eval`` method returns the
integration result.

The ``manualintegrate`` function computes the integral by calling ``eval``
on the rule returned by ``integral_steps``.

The integrator can be extended with new heuristics and evaluation
techniques. To do so, extend the ``Rule`` class, implement ``eval`` method,
then write a function that accepts an ``IntegralInfo`` object and returns
either a ``Rule`` instance or ``None``. If the new technique requires a new
match, add the key and call to the antiderivative function to integral_steps.
To enable simple substitutions, add the match to find_substitutions.

"""

from __future__ import annotations
from typing import NamedTuple, Callable, Sequence, TYPE_CHECKING
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Mapping

from sympy.core.add import Add
from sympy.core.cache import cacheit
from sympy.core.containers import Dict
from sympy.core.function import Derivative
from sympy.core.logic import fuzzy_not
from sympy.core.mul import Mul
from sympy.core.numbers import Integer, Number, E
from sympy.core.power import Pow
from sympy.core.relational import Eq, Ne
from sympy.core.singleton import S
from sympy.core.sorting import ordered
from sympy.core.symbol import Dummy, Symbol, Wild
from sympy.core.exprtools import factor_terms
from sympy.core.function import WildFunction
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.hyperbolic import (HyperbolicFunction, csch,
    cosh, coth, sech, sinh, tanh, asinh)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import (TrigonometricFunction,
    cos, sin, tan, cot, csc, sec, acos, asin, atan, acot, acsc, asec)
from sympy.functions.special.delta_functions import Heaviside, DiracDelta
from sympy.functions.special.error_functions import (erf, erfc, erfi, fresnelc,
    fresnels, Ci, Chi, Si, Shi, Ei, li)
from sympy.functions.special.gamma_functions import uppergamma
from sympy.functions.special.elliptic_integrals import elliptic_e, elliptic_f
from sympy.functions.special.polynomials import (chebyshevt, chebyshevu,
    legendre, hermite, laguerre, assoc_laguerre, gegenbauer, jacobi,
    OrthogonalPolynomial)
from sympy.functions.special.zeta_functions import polylog
from .integrals import Integral
from sympy.logic.boolalg import And, Boolean
from sympy.ntheory.factor_ import primefactors
from sympy.polys.polytools import degree, factor_list, lcm_list, gcd_list, Poly
from sympy.simplify.radsimp import fraction
from sympy.simplify.simplify import simplify
from sympy.simplify.powsimp import powsimp
from sympy.solvers.solvers import solve
from sympy.strategies.core import switch, do_one, null_safe, condition
from sympy.utilities.iterables import iterable
from sympy.utilities.misc import debug

if TYPE_CHECKING:
    from sympy.core.expr import Expr
