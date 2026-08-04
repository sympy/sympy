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
either a ``Rule`` instance or ``None``. If the rule needs subproblems
solved, write it as a generator function that yields the ``IntegralInfo``
of each subproblem and receives the resulting ``Rule`` back at the yield
expression; the ``IntegrationSolver`` driving the rule performs the
recursive calls. If the new technique requires a new match, add the key
and call to the antiderivative function to the rule strategy in
``IntegrationSolver``. To enable simple substitutions, add the match to
find_substitutions.

"""

from __future__ import annotations
from typing import NamedTuple, Callable, Sequence, TYPE_CHECKING, cast
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Mapping
from functools import wraps
from inspect import signature
from types import GeneratorType

from sympy import SYMPY_DEBUG
from sympy.core.add import Add
from sympy.core.cache import cacheit
from sympy.core.basic import Basic
from sympy.core.containers import Dict
from sympy.core.function import Derivative, Function, Lambda, expand_trig
from sympy.core.logic import fuzzy_not
from sympy.core.mul import Mul
from sympy.core.numbers import Integer, Number, E, Rational, pi
from sympy.core.power import Pow
from sympy.core.relational import Eq, Ne
from sympy.core.singleton import S
from sympy.core.sorting import ordered
from sympy.core.traversal import preorder_traversal
from sympy.core.symbol import Dummy, Symbol, Wild
from sympy.core.exprtools import factor_terms
from sympy.core.function import WildFunction
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.hyperbolic import (HyperbolicFunction, csch,
    cosh, coth, sech, sinh, tanh, asinh, acosh, atanh, acoth, asech, acsch)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import (TrigonometricFunction,
    cos, sin, tan, cot, csc, sec, acos, asin, atan, acot, acsc, asec)
from sympy.functions.elementary.integers import ceiling, floor, frac
from sympy.core.mod import Mod
from sympy.concrete.summations import Sum
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
from .rationaltools import ratint
from sympy.logic.boolalg import And, Boolean
from sympy.ntheory.factor_ import primefactors
from sympy.simplify.fu import TR1, TR2
from sympy.polys.polytools import degree, factor_list, lcm_list, gcd_list, Poly
from sympy.polys.polyerrors import PolynomialError
from sympy.polys.rootoftools import RootSum
from sympy.simplify.radsimp import fraction
from sympy.simplify.simplify import simplify
from sympy.simplify.powsimp import powsimp
from sympy.solvers.solvers import solve
from sympy.strategies.core import switch, do_one, null_safe, condition
from sympy.utilities.iterables import iterable
from sympy.utilities.misc import debug

if TYPE_CHECKING:
    from sympy.core.expr import Expr

def _if_zero_implies_zero(P, Q):
    """
    Check if expression P = 0 implies Q = 0.

    Returns True if P is not zero or if substituting every irreducible
    factor of the numerator of P in the numerator of Q makes Q = 0.
    """
    if P.is_zero is False:
        return True
    if P.is_zero is True:
        return Q.is_zero
    num_p, _ = P.as_numer_denom()
    num_q, _ = Q.as_numer_denom()
    factors_P = {f for f, p in factor_list(num_p)[1]}
    # use factor() to help find substitutions (eg. (a**2 - 1) is zero if (a + 1) = 0)
    factored_num_q = num_q.factor()
    for factor in factors_P:
        if factored_num_q.subs(factor, 0) != 0:
            return False
    return True

class Rule(ABC):

    __slots__ = ('integrand', 'variable')

    integrand: Expr
    variable: Symbol

    def __init__(self, integrand: Expr, variable: Symbol):
        self.integrand = integrand
        self.variable = variable

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return all(
            getattr(self, attr) == getattr(other, attr) for attr in self._get_slots()
        )

    def __repr__(self) -> str:
        parts = [f"{self.__class__.__name__}("]
        for i, name in enumerate(self._get_slots()):
            if i:
                parts.append(", ")
            parts.append(f"{name}={getattr(self, name)!r}")
        parts.append(")")
        return "".join(parts)

    def _get_slots(self):
        seen = {()}
        slots = []
        for cls in self.__class__.__mro__[::-1]:
            key = tuple(getattr(cls, "__slots__", ()))
            if key not in seen:
                seen.add(key)
                slots.extend(key)
        return slots

    @abstractmethod
    def eval(self) -> Expr:
        pass

    @abstractmethod
    def contains_dont_know(self) -> bool:
        pass


class AtomicRule(Rule, ABC):
    """A simple rule that does not depend on other rules"""

    __slots__ = ()

    if SYMPY_DEBUG:
        def __init_subclass__(cls, **kwargs):
            super().__init_subclass__(**kwargs)
            if "eval" not in cls.__dict__:
                return
            original_eval = cls.eval

            @wraps(original_eval)
            def wrapped_eval(self, *args, **kwargs):
                sig = signature(type(self).__init__)
                params = ', '.join(
                    f"{name}={getattr(self, name)!r}"
                    for name in sig.parameters
                    if name != 'self' and hasattr(self, name)
                )
                debug(f"Rule calling {type(self).__name__}({params})")
                return original_eval(self, *args, **kwargs)
            cls.eval = wrapped_eval

    def contains_dont_know(self) -> bool:
        return False


class ConstantRule(AtomicRule):
    """integrate(a, x)  ->  a*x"""

    __slots__ = ()

    def eval(self) -> Expr:
        return self.integrand * self.variable


class ConstantTimesRule(Rule):
    """integrate(a*f(x), x)  ->  a*integrate(f(x), x)"""

    __slots__ = ('constant', 'other', 'substep')

    constant: Expr
    other: Expr
    substep: Rule

    def __init__(
        self,
        integrand: Expr,
        variable: Symbol,
        constant: Expr,
        other: Expr,
        substep: Rule,
    ) -> None:
        super().__init__(integrand, variable)
        self.constant = constant
        self.other = other
        self.substep = substep

    def eval(self) -> Expr:
        return self.constant * self.substep.eval()

    def contains_dont_know(self) -> bool:
        return self.substep.contains_dont_know()


class PowerRule(AtomicRule):
    """integrate(x**a, x)"""

    __slots__ = ("base", "exp")

    base: Expr
    exp: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, base: Expr, exp: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.base = base
        self.exp = exp

    def eval(self) -> Expr:
        return Piecewise(
            ((self.base**(self.exp + 1))/(self.exp + 1), Ne(self.exp, -1)),
            (log(self.base), True),
        )


class NestedPowRule(AtomicRule):
    """integrate((x**a)**b, x)"""

    __slots__ = ("base", "exp")

    base: Expr
    exp: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, base: Expr, exp: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.base = base
        self.exp = exp

    def eval(self) -> Expr:
        m = self.base * self.integrand
        return Piecewise((m / (self.exp + 1), Ne(self.exp, -1)),
                         (m * log(self.base), True))


class AddRule(Rule):
    """integrate(f(x) + g(x), x) -> integrate(f(x), x) + integrate(g(x), x)"""

    __slots__ = ("substeps",)

    substeps: list[Rule]

    def __init__(self, integrand: Expr, variable: Symbol, substeps: list[Rule]) -> None:
        super().__init__(integrand, variable)
        self.substeps = substeps

    def eval(self) -> Expr:
        return Add(*(substep.eval() for substep in self.substeps))

    def contains_dont_know(self) -> bool:
        return any(substep.contains_dont_know() for substep in self.substeps)


class URule(Rule):
    """integrate(f(g(x))*g'(x), x) -> integrate(f(u), u), u = g(x)"""

    __slots__ = ("u_var", "u_func", "substep")

    u_var: Symbol
    u_func: Expr
    substep: Rule

    def __init__(
        self,
        integrand: Expr,
        variable: Symbol,
        u_var: Symbol,
        u_func: Expr,
        substep: Rule,
    ) -> None:
        super().__init__(integrand, variable)
        self.u_var = u_var
        self.u_func = u_func
        self.substep = substep

    def eval(self) -> Expr:
        result = self.substep.eval()
        if self.u_func.is_Pow:
            base, exp_ = self.u_func.as_base_exp()
            if exp_ == -1:
                # avoid needless -log(1/x) from substitution
                result = result.subs(log(self.u_var), -log(base))
        return result.subs(self.u_var, self.u_func)

    def contains_dont_know(self) -> bool:
        return self.substep.contains_dont_know()


class PartsRule(Rule):
    """integrate(u(x)*v'(x), x) -> u(x)*v(x) - integrate(u'(x)*v(x), x)"""

    __slots__ = ("u", "dv", "v_step", "second_step")

    u: Symbol
    dv: Expr
    v_step: Rule
    second_step: Rule | None  # None when is a substep of CyclicPartsRule

    def __init__(
        self,
        integrand: Expr,
        variable: Symbol,
        u: Symbol,
        dv: Expr,
        v_step: Rule,
        second_step: Rule | None = None,
    ) -> None:
        super().__init__(integrand, variable)
        self.u = u
        self.dv = dv
        self.v_step = v_step
        self.second_step = second_step

    def eval(self) -> Expr:
        assert self.second_step is not None
        v = self.v_step.eval()
        return self.u * v - self.second_step.eval()

    def contains_dont_know(self) -> bool:
        return self.v_step.contains_dont_know() or (
            self.second_step is not None and self.second_step.contains_dont_know())


class CyclicPartsRule(Rule):
    """Apply PartsRule multiple times to integrate exp(x)*sin(x)"""

    __slots__ = ("parts_rules", "coefficient")

    parts_rules: list[PartsRule]
    coefficient: Expr

    def __init__(
        self,
        integrand: Expr,
        variable: Symbol,
        parts_rules: list[PartsRule],
        coefficient: Expr,
    ) -> None:
        super().__init__(integrand, variable)
        self.parts_rules = parts_rules
        self.coefficient = coefficient

    def eval(self) -> Expr:
        result = []
        sign = 1
        for rule in self.parts_rules:
            result.append(sign * rule.u * rule.v_step.eval())
            sign *= -1
        return Add(*result) / (1 - self.coefficient)

    def contains_dont_know(self) -> bool:
        return any(substep.contains_dont_know() for substep in self.parts_rules)


class TrigRule(AtomicRule, ABC):
    __slots__ = ()


class SinRule(TrigRule):
    """integrate(sin(x), x) -> -cos(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return -cos(self.variable)


class CosRule(TrigRule):
    """integrate(cos(x), x) -> sin(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return sin(self.variable)


class SecTanRule(TrigRule):
    """integrate(sec(x)*tan(x), x) -> sec(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return sec(self.variable)


class CscCotRule(TrigRule):
    """integrate(csc(x)*cot(x), x) -> -csc(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return -csc(self.variable)


class Sec2Rule(TrigRule):
    """integrate(sec(x)**2, x) -> tan(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return tan(self.variable)


class Csc2Rule(TrigRule):
    """integrate(csc(x)**2, x) -> -cot(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return -cot(self.variable)


class HyperbolicRule(AtomicRule, ABC):
    __slots__ = ()


class SinhRule(HyperbolicRule):
    """integrate(sinh(x), x) -> cosh(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return cosh(self.variable)


class CoshRule(HyperbolicRule):
    """integrate(cosh(x), x) -> sinh(x)"""

    __slots__ = ()

    def eval(self):
        return sinh(self.variable)


class ExpRule(AtomicRule):
    """integrate(a**x, x) -> a**x/ln(a)"""

    __slots__ = ("base", "exp")

    base: Expr
    exp: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, base: Expr, exp: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.base = base
        self.exp = exp

    def eval(self) -> Expr:
        return self.integrand / log(self.base)


class ReciprocalRule(AtomicRule):
    """integrate(1/x, x) -> ln(x)"""

    __slots__ = ("base",)

    base: Expr

    def __init__(self, integrand: Expr, variable: Symbol, base: Expr) -> None:
        super().__init__(integrand, variable)
        self.base = base

    def eval(self) -> Expr:
        return log(self.base)


class ArcsinRule(AtomicRule):
    """integrate(1/sqrt(1-x**2), x) -> asin(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return asin(self.variable)


class ArcsinhRule(AtomicRule):
    """integrate(1/sqrt(1+x**2), x) -> asin(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        return asinh(self.variable)


class ReciprocalSqrtQuadraticRule(AtomicRule):
    """integrate(1/sqrt(a+b*x+c*x**2), x) -> log(2*sqrt(c)*sqrt(a+b*x+c*x**2)+b+2*c*x)/sqrt(c)"""

    __slots__ = ("a", "b", "c")

    a: Expr
    b: Expr
    c: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, a: Expr, b: Expr, c: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c

    def eval(self) -> Expr:
        a, b, c, x = self.a, self.b, self.c, self.variable
        return log(2*sqrt(c)*sqrt(a+b*x+c*x**2)+b+2*c*x)/sqrt(c)


class SqrtQuadraticDenomRule(AtomicRule):
    """integrate(poly(x)/sqrt(a+b*x+c*x**2), x)"""

    __slots__ = ("a", "b", "c", "coeffs")

    a: Expr
    b: Expr
    c: Expr
    coeffs: list[Expr]

    def __init__(
        self,
        integrand: Expr,
        variable: Symbol,
        a: Expr,
        b: Expr,
        c: Expr,
        coeffs: list[Expr],
    ) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c
        self.coeffs = coeffs

    def eval(self) -> Expr:
        a, b, c, coeffs, x = self.a, self.b, self.c, self.coeffs.copy(), self.variable
        # Integrate poly/sqrt(a+b*x+c*x**2) using recursion.
        # coeffs are coefficients of the polynomial.
        # Let I_n = x**n/sqrt(a+b*x+c*x**2), then
        # I_n = A * x**(n-1)*sqrt(a+b*x+c*x**2) - B * I_{n-1} - C * I_{n-2}
        # where A = 1/(n*c), B = (2*n-1)*b/(2*n*c), C = (n-1)*a/(n*c)
        # See https://github.com/sympy/sympy/pull/23608 for proof.
        result_coeffs = []
        coeffs = coeffs.copy()
        for i in range(len(coeffs)-2):
            n = len(coeffs)-1-i
            coeff = coeffs[i]/(c*n)
            result_coeffs.append(coeff)
            coeffs[i+1] -= (2*n-1)*b/2*coeff
            coeffs[i+2] -= (n-1)*a*coeff
        d, e = coeffs[-1], coeffs[-2]
        s = sqrt(a+b*x+c*x**2)
        constant = d-b*e/(2*c)
        if constant == 0:
            I0 = 0
        else:
            gen = inverse_trig_rule(IntegralInfo(1/s, x), degenerate=False)
            step = IntegrationSolver().run_generator(gen)
            I0 = constant*step.eval()
        return Add(*(result_coeffs[i]*x**(len(coeffs)-2-i)
                     for i in range(len(result_coeffs))), e/c)*s + I0


class SqrtQuadraticRule(AtomicRule):
    """integrate(sqrt(a+b*x+c*x**2), x)"""

    __slots__ = ("a", "b", "c")

    a: Expr
    b: Expr
    c: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, a: Expr, b: Expr, c: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c

    def eval(self) -> Expr:
        gen = sqrt_quadratic_rule(IntegralInfo(self.integrand, self.variable), degenerate=False)
        step = IntegrationSolver().run_generator(gen)
        return step.eval()


class RatintRule(AtomicRule):
    """Integrate a rational function using ``ratint`` as a fallback."""

    __slots__ = ()

    def eval(self) -> Expr:
        return ratint(self.integrand, self.variable)


class AlternativeRule(Rule):
    """Multiple ways to do integration."""

    __slots__ = ("alternatives",)

    alternatives: list[Rule]

    def __init__(
        self, integrand: Expr, variable: Symbol, alternatives: list[Rule]
    ) -> None:
        super().__init__(integrand, variable)
        self.alternatives = alternatives

    def eval(self) -> Expr:
        return self.alternatives[0].eval()

    def contains_dont_know(self) -> bool:
        return any(substep.contains_dont_know() for substep in self.alternatives)


class DontKnowRule(Rule):
    """Leave the integral as is."""

    __slots__ = ()

    def eval(self) -> Expr:
        return Integral(self.integrand, self.variable)

    def contains_dont_know(self) -> bool:
        return True


class DerivativeRule(AtomicRule):
    """integrate(f'(x), x) -> f(x)"""

    __slots__ = ()

    def eval(self) -> Expr:
        assert isinstance(self.integrand, Derivative)
        variable_count = list(self.integrand.variable_count)
        for i, (var, count) in enumerate(variable_count):
            if var == self.variable:
                variable_count[i] = (var, count - 1)
                break
        return Derivative(self.integrand.expr, *variable_count)


class RewriteRule(Rule):
    """Rewrite integrand to another form that is easier to handle."""

    __slots__ = ("rewritten", "substep")

    rewritten: Expr
    substep: Rule

    def __init__(
        self, integrand: Expr, variable: Symbol, rewritten: Expr, substep: Rule
    ) -> None:
        super().__init__(integrand, variable)
        self.rewritten = rewritten
        self.substep = substep

    def eval(self) -> Expr:
        return self.substep.eval()

    def contains_dont_know(self) -> bool:
        return self.substep.contains_dont_know()


class CompleteSquareRule(RewriteRule):
    """Rewrite a+b*x+c*x**2 to a-b**2/(4*c) + c*(x+b/(2*c))**2"""
    __slots__ = ()


class PiecewiseRule(Rule):

    __slots__ = ("subfunctions",)

    subfunctions: Sequence[tuple[Rule, bool | Boolean]]

    def __init__(
        self,
        integrand: Expr,
        variable: Symbol,
        subfunctions: Sequence[tuple[Rule, bool | Boolean]],
    ) -> None:
        super().__init__(integrand, variable)
        self.subfunctions = subfunctions

    def eval(self) -> Expr:
        return Piecewise(*[(substep.eval(), cond)
                           for substep, cond in self.subfunctions])

    def contains_dont_know(self) -> bool:
        return any(substep.contains_dont_know() for substep, _ in self.subfunctions)


class HeavisideRule(Rule):

    __slots__ = ("harg", "ibnd", "substep")

    harg: Expr
    ibnd: Expr
    substep: Rule

    def __init__(
        self,
        integrand: Expr,
        variable: Symbol,
        harg: Expr,
        ibnd: Expr,
        substep: Rule,
    ) -> None:
        super().__init__(integrand, variable)
        self.harg = harg
        self.ibnd = ibnd
        self.substep = substep

    def eval(self) -> Expr:
        # If we are integrating over x and the integrand has the form
        #       Heaviside(m*x+b)*g(x) == Heaviside(harg)*g(symbol)
        # then there needs to be continuity at -b/m == ibnd,
        # so we subtract the appropriate term.
        result = self.substep.eval()
        return Heaviside(self.harg) * (result - result.subs(self.variable, self.ibnd))

    def contains_dont_know(self) -> bool:
        return self.substep.contains_dont_know()


class DiracDeltaRule(AtomicRule):

    __slots__ = ("n", "a", "b")

    n: Expr
    a: Expr
    b: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, n: Expr, a: Expr, b: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.n = n
        self.a = a
        self.b = b

    def eval(self) -> Expr:
        n, a, b, x = self.n, self.a, self.b, self.variable
        if n == 0:
            return Heaviside(a+b*x)/b
        return DiracDelta(a+b*x, n-1)/b


class TrigSubstitutionRule(Rule):

    __slots__ = ("theta", "func", "rewritten", "substep", "restriction")

    theta: Expr
    func: Expr
    rewritten: Expr
    substep: Rule
    restriction: bool | Boolean

    def __init__(
        self,
        integrand: Expr,
        variable: Symbol,
        theta: Expr,
        func: Expr,
        rewritten: Expr,
        substep: Rule,
        restriction: bool | Boolean,
    ) -> None:
        super().__init__(integrand, variable)
        self.theta = theta
        self.func = func
        self.rewritten = rewritten
        self.substep = substep
        self.restriction = restriction

    def eval(self) -> Expr:
        theta, func, x = self.theta, self.func, self.variable
        func = func.subs(sec(theta), 1/cos(theta))
        func = func.subs(csc(theta), 1/sin(theta))
        func = func.subs(cot(theta), 1/tan(theta))

        trig_function = list(func.find(TrigonometricFunction))
        assert len(trig_function) == 1
        trig_function = trig_function[0]
        relation = solve(x - func, trig_function)
        assert len(relation) == 1
        numer, denom = fraction(relation[0])

        if isinstance(trig_function, sin):
            opposite = numer
            hypotenuse = denom
            adjacent = sqrt(denom**2 - numer**2)
            inverse = asin(relation[0])
        elif isinstance(trig_function, cos):
            adjacent = numer
            hypotenuse = denom
            opposite = sqrt(denom**2 - numer**2)
            inverse = acos(relation[0])
        else:  # tan
            opposite = numer
            adjacent = denom
            hypotenuse = sqrt(denom**2 + numer**2)
            inverse = atan(relation[0])

        substitution = [
            (sin(theta), opposite/hypotenuse),
            (cos(theta), adjacent/hypotenuse),
            (tan(theta), opposite/adjacent),
            (theta, inverse)
        ]
        return Piecewise(
                (self.substep.eval().subs(substitution).trigsimp(), self.restriction) # type: ignore
        )

    def contains_dont_know(self) -> bool:
        return self.substep.contains_dont_know()


class ArctanRule(AtomicRule):
    """integrate(a/(b*x**2+c), x) -> a/b / sqrt(c/b) * atan(x/sqrt(c/b))"""

    __slots__ = ("a", "b", "c")

    a: Expr
    b: Expr
    c: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, a: Expr, b: Expr, c: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c

    def eval(self) -> Expr:
        a, b, c, x = self.a, self.b, self.c, self.variable
        return a/b / sqrt(c/b) * atan(x/sqrt(c/b))


class OrthogonalPolyRule(AtomicRule, ABC):

    __slots__ = ("n",)

    n: Expr

    def __init__(self, integrand: Expr, variable: Symbol, n: Expr) -> None:
        super().__init__(integrand, variable)
        self.n = n


class JacobiRule(OrthogonalPolyRule):

    __slots__ = ("a", "b")

    a: Expr
    b: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, n: Expr, a: Expr, b: Expr
    ) -> None:
        super().__init__(integrand, variable, n)
        self.a = a
        self.b = b

    def eval(self) -> Expr:
        n, a, b, x = self.n, self.a, self.b, self.variable
        return Piecewise(
            (2*jacobi(n + 1, a - 1, b - 1, x)/(n + a + b), Ne(n + a + b, 0)),
            (x, Eq(n, 0)),
            ((a + b + 2)*x**2/4 + (a - b)*x/2, Eq(n, 1)))


class GegenbauerRule(OrthogonalPolyRule):

    __slots__ = ("a",)

    a: Expr

    def __init__(self, integrand: Expr, variable: Symbol, n: Expr, a: Expr) -> None:
        super().__init__(integrand, variable, n)
        self.a = a

    def eval(self) -> Expr:
        n, a, x = self.n, self.a, self.variable
        return Piecewise(
            (gegenbauer(n + 1, a - 1, x)/(2*(a - 1)), Ne(a, 1)),
            (chebyshevt(n + 1, x)/(n + 1), Ne(n, -1)),
            (S.Zero, True))


class ChebyshevTRule(OrthogonalPolyRule):

    __slots__ = ()

    def eval(self) -> Expr:
        n, x = self.n, self.variable
        return Piecewise(
            ((chebyshevt(n + 1, x)/(n + 1) -
              chebyshevt(n - 1, x)/(n - 1))/2, Ne(Abs(n), 1)),
            (x**2/2, True))


class ChebyshevURule(OrthogonalPolyRule):

    __slots__ = ()

    def eval(self) -> Expr:
        n, x = self.n, self.variable
        return Piecewise(
            (chebyshevt(n + 1, x)/(n + 1), Ne(n, -1)),
            (S.Zero, True))


class LegendreRule(OrthogonalPolyRule):

    __slots__ = ()

    def eval(self) -> Expr:
        n, x = self.n, self.variable
        return(legendre(n + 1, x) - legendre(n - 1, x))/(2*n + 1)


class HermiteRule(OrthogonalPolyRule):

    __slots__ = ()

    def eval(self) -> Expr:
        n, x = self.n, self.variable
        return hermite(n + 1, x)/(2*(n + 1))


class LaguerreRule(OrthogonalPolyRule):

    __slots__ = ()

    def eval(self) -> Expr:
        n, x = self.n, self.variable
        return laguerre(n, x) - laguerre(n + 1, x)


class AssocLaguerreRule(OrthogonalPolyRule):

    __slots__ = ("a",)

    a: Expr

    def __init__(self, integrand: Expr, variable: Symbol, n: Expr, a: Expr) -> None:
        super().__init__(integrand, variable, n)
        self.a = a

    def eval(self) -> Expr:
        return -assoc_laguerre(self.n + 1, self.a - 1, self.variable)


class IRule(AtomicRule, ABC):

    __slots__ = ("a", "b")

    a: Expr
    b: Expr

    def __init__(self, integrand: Expr, variable: Symbol, a: Expr, b: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b


class CiRule(IRule):

    __slots__ = ()

    def eval(self) -> Expr:
        a, b, x = self.a, self.b, self.variable
        return cos(b)*Ci(a*x) - sin(b)*Si(a*x)


class ChiRule(IRule):

    __slots__ = ()

    def eval(self) -> Expr:
        a, b, x = self.a, self.b, self.variable
        return cosh(b)*Chi(a*x) + sinh(b)*Shi(a*x)


class EiRule(IRule):

    __slots__ = ()

    def eval(self) -> Expr:
        a, b, x = self.a, self.b, self.variable
        return exp(b)*Ei(a*x)


class SiRule(IRule):

    __slots__ = ()

    def eval(self) -> Expr:
        a, b, x = self.a, self.b, self.variable
        return sin(b)*Ci(a*x) + cos(b)*Si(a*x)


class ShiRule(IRule):

    __slots__ = ()

    def eval(self) -> Expr:
        a, b, x = self.a, self.b, self.variable
        return sinh(b)*Chi(a*x) + cosh(b)*Shi(a*x)


class LiRule(IRule):

    __slots__ = ()

    def eval(self) -> Expr:
        a, b, x = self.a, self.b, self.variable
        return li(a*x + b)/a


class ErfRule(AtomicRule):

    __slots__ = ("a", "b", "c")

    a: Expr
    b: Expr
    c: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, a: Expr, b: Expr, c: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c

    def eval(self) -> Expr:
        a, b, c, x = self.a, self.b, self.c, self.variable
        if a.is_extended_real:
            return Piecewise(
                (sqrt(S.Pi)/sqrt(-a)/2 * exp(c - b**2/(4*a)) *
                    erf((-2*a*x - b)/(2*sqrt(-a))), a < 0),
                (sqrt(S.Pi)/sqrt(a)/2 * exp(c - b**2/(4*a)) *
                    erfi((2*a*x + b)/(2*sqrt(a))), True))
        return sqrt(S.Pi)/sqrt(a)/2 * exp(c - b**2/(4*a)) * \
                erfi((2*a*x + b)/(2*sqrt(a)))


class FresnelCRule(AtomicRule):

    __slots__ = ("a", "b", "c")

    a: Expr
    b: Expr
    c: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, a: Expr, b: Expr, c: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c

    def eval(self) -> Expr:
        a, b, c, x = self.a, self.b, self.c, self.variable
        return sqrt(S.Pi)/sqrt(2*a) * (
            cos(b**2/(4*a) - c)*fresnelc((2*a*x + b)/sqrt(2*a*S.Pi)) +
            sin(b**2/(4*a) - c)*fresnels((2*a*x + b)/sqrt(2*a*S.Pi)))


class FresnelSRule(AtomicRule):

    __slots__ = ("a", "b", "c")

    a: Expr
    b: Expr
    c: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, a: Expr, b: Expr, c: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c

    def eval(self) -> Expr:
        a, b, c, x = self.a, self.b, self.c, self.variable
        return sqrt(S.Pi)/sqrt(2*a) * (
            cos(b**2/(4*a) - c)*fresnels((2*a*x + b)/sqrt(2*a*S.Pi)) -
            sin(b**2/(4*a) - c)*fresnelc((2*a*x + b)/sqrt(2*a*S.Pi)))


class PolylogRule(AtomicRule):

    __slots__ = ("a", "b")

    a: Expr
    b: Expr

    def __init__(self, integrand: Expr, variable: Symbol, a: Expr, b: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b

    def eval(self) -> Expr:
        return polylog(self.b + 1, self.a * self.variable)


class UpperGammaRule(AtomicRule):

    __slots__ = ("a", "e")

    a: Expr
    e: Expr

    def __init__(self, integrand: Expr, variable: Symbol, a: Expr, e: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.e = e

    def eval(self) -> Expr:
        a, e, x = self.a, self.e, self.variable
        return x**e * (-a*x)**(-e) * uppergamma(e + 1, -a*x)/a


class EllipticFRule(AtomicRule):

    __slots__ = ("a", "d")

    a: Expr
    d: Expr

    def __init__(self, integrand: Expr, variable: Symbol, a: Expr, d: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.d = d

    def eval(self) -> Expr:
        return elliptic_f(self.variable, self.d/self.a)/sqrt(self.a)


class EllipticERule(AtomicRule):

    __slots__ = ("a", "d")

    a: Expr
    d: Expr

    def __init__(self, integrand: Expr, variable: Symbol, a: Expr, d: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.d = d

    def eval(self) -> Expr:
        return elliptic_e(self.variable, self.d/self.a)*sqrt(self.a)


class IntegralInfo(NamedTuple):
    integrand: Expr
    symbol: Symbol


class PartsUCheck(NamedTuple):
    """Request, yielded by a generator rule, that the solver record one more
    use of ``u`` as the part to be differentiated in integration by parts
    and report whether picking it once more is still allowed.

    ``u_key`` is the chosen ``u`` with the integration variable replaced by
    a common dummy, so that uses are counted across subproblems.
    """
    u_key: Expr


class BranchQuery(NamedTuple):
    """Request, yielded by a generator rule, for the value of the solver's
    ``branch`` option: ``True`` if the rule should collect every applicable
    alternative into an :class:`AlternativeRule`, ``False`` if the first
    workable candidate is enough.
    """


def manual_diff(f, symbol):
    """Derivative of f in form expected by find_substitutions

    SymPy's derivatives for some trig functions (like cot) are not in a form
    that works well with finding substitutions; this replaces the
    derivatives for those particular forms with something that works better.

    """
    if f.args:
        arg = f.args[0]
        if isinstance(f, tan):
            return arg.diff(symbol) * sec(arg)**2
        elif isinstance(f, cot):
            return -arg.diff(symbol) * csc(arg)**2
        elif isinstance(f, sec):
            return arg.diff(symbol) * sec(arg) * tan(arg)
        elif isinstance(f, csc):
            return -arg.diff(symbol) * csc(arg) * cot(arg)
        elif isinstance(f, Add):
            return sum(manual_diff(arg, symbol) for arg in f.args)
        elif isinstance(f, Mul):
            if len(f.args) == 2 and isinstance(f.args[0], Number):
                return f.args[0] * manual_diff(f.args[1], symbol)
    return f.diff(symbol)

def manual_subs(expr, *args):
    """
    A wrapper for `expr.subs(*args)` with additional logic for substitution
    of invertible functions.
    """
    if len(args) == 1:
        sequence = args[0]
        if isinstance(sequence, (Dict, Mapping)):
            sequence = sequence.items()
        elif not iterable(sequence):
            raise ValueError("Expected an iterable of (old, new) pairs")
    elif len(args) == 2:
        sequence = [args]
    else:
        raise ValueError("subs accepts either 1 or 2 arguments")

    new_subs = []
    for old, new in sequence:
        if isinstance(old, log):
            # If log(x) = y, then exp(a*log(x)) = exp(a*y)
            # that is, x**a = exp(a*y). Replace nontrivial powers of x
            # before subs turns them into `exp(y)**a`, but
            # do not replace x itself yet, to avoid `log(exp(y))`.
            x0 = old.args[0]
            expr = expr.replace(lambda x: x.is_Pow and x.base == x0,
                lambda x: exp(x.exp*new))
            new_subs.append((x0, exp(new)))

    return expr.subs(list(sequence) + new_subs)

# Method based on that on SIN, described in "Symbolic Integration: The
# Stormy Decade"

inverse_trig_functions = (atan, asin, acos, acot, acsc, asec)


def find_substitutions(integrand, symbol, u_var):
    results = []

    def test_subterm(u, u_diff):
        if u_diff == 0:
            return False
        substituted = integrand / u_diff
        debug("substituted: {}, u: {}, u_var: {}".format(substituted, u, u_var))
        substituted = manual_subs(substituted, u, u_var).cancel()

        if substituted.has_free(symbol):
            return False
        # avoid increasing the degree of a rational function
        if integrand.is_rational_function(symbol) and substituted.is_rational_function(u_var):
            deg_before = max(degree(t, symbol) for t in integrand.as_numer_denom())
            deg_after = max(degree(t, u_var) for t in substituted.as_numer_denom())
            if deg_after > deg_before:
                return False
        return substituted.as_independent(u_var, as_Add=False)

    def exp_subterms(term: Expr):
        linear_coeffs = []
        terms = []
        n = Wild('n', properties=[lambda n: n.is_Integer])
        for exp_ in term.find(exp):
            arg = exp_.args[0]
            if symbol not in arg.free_symbols:
                continue
            match = arg.match(n*symbol)
            if match:
                linear_coeffs.append(match[n])
            else:
                terms.append(exp_)
        if linear_coeffs:
            terms.append(exp(gcd_list(linear_coeffs)*symbol))
        return terms

    def possible_subterms(term):
        if isinstance(term, (TrigonometricFunction, HyperbolicFunction,
                             *inverse_trig_functions,
                             exp, log, Heaviside)):
            return [term.args[0]]
        elif isinstance(term, (chebyshevt, chebyshevu,
                        legendre, hermite, laguerre)):
            return [term.args[1]]
        elif isinstance(term, (gegenbauer, assoc_laguerre)):
            return [term.args[2]]
        elif isinstance(term, jacobi):
            return [term.args[3]]
        elif isinstance(term, Mul):
            r = []
            for u in term.args:
                r.append(u)
                r.extend(possible_subterms(u))
            return r
        elif isinstance(term, Pow):
            r = [arg for arg in term.args if arg.has(symbol)]
            if term.exp.is_Integer:
                r.extend([term.base**d for d in primefactors(term.exp)
                    if 1 < d < abs(term.args[1])])
                if term.base.is_Add:
                    r.extend([t for t in possible_subterms(term.base)
                        if t.is_Pow])
            return r
        elif isinstance(term, Add):
            r = []
            for arg in term.args:
                r.append(arg)
                r.extend(possible_subterms(arg))
            return r
        return []

    for u in list(dict.fromkeys(possible_subterms(integrand) + exp_subterms(integrand))):
        if u == symbol:
            continue
        u_diff = manual_diff(u, symbol)
        new_integrand = test_subterm(u, u_diff)
        if new_integrand is not False:
            constant, new_integrand = new_integrand
            if new_integrand == integrand.subs(symbol, u_var):
                continue
            substitution = (u, constant, new_integrand)
            if substitution not in results:
                results.append(substitution)

    return results

def rewriter(condition, rewrite):
    """Strategy that rewrites an integrand."""
    def _rewriter(integral):
        integrand, symbol = integral
        debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewrite, symbol))
        if condition(*integral):
            rewritten = rewrite(*integral)
            if rewritten != integrand:
                substep = yield IntegralInfo(rewritten, symbol)
                if not isinstance(substep, DontKnowRule) and substep:
                    return RewriteRule(integrand, symbol, rewritten, substep)
    return _rewriter

def proxy_rewriter(condition, rewrite):
    """Strategy that rewrites an integrand based on some other criteria."""
    def _proxy_rewriter(criteria):
        criteria, integral = criteria
        integrand, symbol = integral
        debug("Integral: {} is rewritten with {} on symbol: {} and criteria: {}".format(integrand, rewrite, symbol, criteria))
        args = criteria + list(integral)
        if condition(*args):
            rewritten = rewrite(*args)
            if rewritten != integrand:
                substep = yield IntegralInfo(rewritten, symbol)
                return RewriteRule(integrand, symbol, rewritten, substep)
    return _proxy_rewriter

def multiplexer(conditions):
    """Apply the rule that matches the condition, else None"""
    def multiplexer_rl(expr):
        for key, rule in conditions.items():
            if key(expr):
                return rule(expr)
    return multiplexer_rl

def alternatives(*rules, branch=False):
    """Strategy that makes an AlternativeRule out of multiple possible results.

    With ``branch=False`` no ``AlternativeRule`` is built: the first rule
    whose result is free of ``DontKnowRule`` is returned directly, falling
    back to the first partial result if no rule solves the integral
    completely.
    """
    def _alternatives(integral):
        alts = []
        first_partial = None
        count = 0
        if branch:
            debug("List of Alternative Rules")
        for rule in rules:
            count = count + 1
            if branch:
                debug("Rule {}: {}".format(count, rule))

            result = rule(integral)
            if (not result or isinstance(result, DontKnowRule) or
                    result == integral):
                continue

            if not branch:
                if not result.contains_dont_know():
                    return result
                if first_partial is None:
                    first_partial = result
                continue

            if result not in alts:
                alts.append(result)

        if not branch:
            return first_partial

        if len(alts) == 1:
            return alts[0]
        elif alts:
            doable = [rule for rule in alts if not rule.contains_dont_know()]
            if doable:
                return AlternativeRule(*integral, doable)
            else:
                return AlternativeRule(*integral, alts)
    return _alternatives

def constant_rule(integral):
    return ConstantRule(*integral)

def power_rule(integral):
    integrand, symbol = integral
    base, expt = integrand.as_base_exp()

    if symbol not in expt.free_symbols and isinstance(base, Symbol):
        if simplify(expt + 1) == 0:
            return ReciprocalRule(integrand, symbol, base)
        return PowerRule(integrand, symbol, base, expt)
    elif symbol not in base.free_symbols and isinstance(expt, Symbol):
        rule = ExpRule(integrand, symbol, base, expt)

        if fuzzy_not(log(base).is_zero):
            return rule
        elif log(base).is_zero:
            return ConstantRule(1, symbol)

        return PiecewiseRule(integrand, symbol, [
            (rule, Ne(log(base), 0)),
            (ConstantRule(1, symbol), True)
        ])

def exp_rule(integral):
    integrand, symbol = integral
    if isinstance(integrand.args[0], Symbol):
        return ExpRule(integrand, symbol, E, integrand.args[0])


def combine_power_rule(integral):
    """
    Strategy that simplifies the exponent of a power.
    exp(a*x**2) * exp(b*x) -> exp((a*x**2 + b*x))
    For example, this is useful for the ErfRule.
    """
    integrand, symbol = integral
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    k = Wild('k', exclude=[symbol])
    rest = Wild('rest')

    match = integrand.match(rest * k**(a*symbol**2) * k**(b*symbol))

    if not match:
        return

    simplified = powsimp(integrand, combine='exp')

    if simplified != integrand:
        steps = yield IntegralInfo(simplified, symbol)
        return RewriteRule(integrand, symbol, simplified, steps)


def orthogonal_poly_rule(integral):
    orthogonal_poly_classes = {
        jacobi: JacobiRule,
        gegenbauer: GegenbauerRule,
        chebyshevt: ChebyshevTRule,
        chebyshevu: ChebyshevURule,
        legendre: LegendreRule,
        hermite: HermiteRule,
        laguerre: LaguerreRule,
        assoc_laguerre: AssocLaguerreRule
        }
    orthogonal_poly_var_index = {
        jacobi: 3,
        gegenbauer: 2,
        assoc_laguerre: 2
        }
    integrand, symbol = integral
    for klass in orthogonal_poly_classes:
        if isinstance(integrand, klass):
            var_index = orthogonal_poly_var_index.get(klass, 1)
            if (integrand.args[var_index] is symbol and not
                any(v.has(symbol) for v in integrand.args[:var_index])):
                    return orthogonal_poly_classes[klass](integrand, symbol, *integrand.args[:var_index])


_special_function_patterns: list[tuple[type, Expr, Callable | None, tuple]] = []
_wilds: list[Wild] = []
_symbol = Dummy('x')


def special_function_rule(integral):
    integrand, symbol = integral
    if not _special_function_patterns:
        a = Wild('a', exclude=[_symbol], properties=[lambda x: not x.is_zero])
        b = Wild('b', exclude=[_symbol])
        c = Wild('c', exclude=[_symbol])
        d = Wild('d', exclude=[_symbol], properties=[lambda x: not x.is_zero])
        e = Wild('e', exclude=[_symbol], properties=[
            lambda x: not (x.is_nonnegative and x.is_integer)])
        _wilds.extend((a, b, c, d, e))
        # patterns consist of a SymPy class, a wildcard expr, an optional
        # condition coded as a lambda (when Wild properties are not enough),
        # followed by an applicable rule
        linear_pattern = a*_symbol + b
        quadratic_pattern = a*_symbol**2 + b*_symbol + c
        _special_function_patterns.extend((
            (Mul, exp(linear_pattern, evaluate=False)/_symbol, None, EiRule),
            (Mul, cos(linear_pattern, evaluate=False)/_symbol, None, CiRule),
            (Mul, cosh(linear_pattern, evaluate=False)/_symbol, None, ChiRule),
            (Mul, sin(linear_pattern, evaluate=False)/_symbol, None, SiRule),
            (Mul, sinh(linear_pattern, evaluate=False)/_symbol, None, ShiRule),
            (Pow, 1/log(linear_pattern, evaluate=False), None, LiRule),
            (exp, exp(quadratic_pattern, evaluate=False), None, ErfRule),
            (sin, sin(quadratic_pattern, evaluate=False), None, FresnelSRule),
            (cos, cos(quadratic_pattern, evaluate=False), None, FresnelCRule),
            (Mul, _symbol**e*exp(a*_symbol, evaluate=False), None, UpperGammaRule),
            (Mul, polylog(b, a*_symbol, evaluate=False)/_symbol, None, PolylogRule),
            (Pow, 1/sqrt(a - d*sin(_symbol, evaluate=False)**2),
                lambda a, d: a != d, EllipticFRule),
            (Pow, sqrt(a - d*sin(_symbol, evaluate=False)**2),
                lambda a, d: a != d, EllipticERule),
        ))
    _integrand = integrand.subs(symbol, _symbol)
    for type_, pattern, constraint, rule in _special_function_patterns:
        if isinstance(_integrand, type_):
            match = _integrand.match(pattern)
            if match:
                wild_vals = tuple(match.get(w) for w in _wilds
                                  if match.get(w) is not None)
                if constraint is None or constraint(*wild_vals):
                    return rule(integrand, symbol, *wild_vals)


def _add_degenerate_step(generic_cond, generic_step: Rule, degenerate_step: Rule | None) -> Rule:
    if degenerate_step is None:
        return generic_step
    if isinstance(generic_step, PiecewiseRule):
        subfunctions = [(substep, (cond & generic_cond).simplify())
                        for substep, cond in generic_step.subfunctions]
    else:
        subfunctions = [(generic_step, generic_cond)]
    if isinstance(degenerate_step, PiecewiseRule):
        subfunctions += degenerate_step.subfunctions
    else:
        subfunctions.append((degenerate_step, S.true))
    return PiecewiseRule(generic_step.integrand, generic_step.variable, subfunctions)


def nested_pow_rule(integral: IntegralInfo):
    # nested (c*(a+b*x)**d)**e
    integrand, x = integral

    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    pattern = a_+b_*x
    generic_cond: Boolean = S.true

    class NoMatch(Exception):
        pass

    def _get_base_exp(expr: Expr) -> tuple[Expr, Expr]:
        if not expr.has_free(x):
            return S.One, S.Zero
        if expr.is_Mul:
            _, terms = expr.as_coeff_mul()
            if not terms:
                return S.One, S.Zero
            results = [_get_base_exp(term) for term in terms]
            bases = {b for b, _ in results}
            bases.discard(S.One)
            if len(bases) == 1:
                return bases.pop(), Add(*(e for _, e in results))
            raise NoMatch
        if expr.is_Pow:
            b, e = expr.base, expr.exp  # type: ignore
            if e.has_free(x):
                raise NoMatch
            base_, sub_exp = _get_base_exp(b)
            return base_, sub_exp * e
        match = expr.match(pattern)
        if match:
            a, b = match[a_], match[b_]
            base_ = x + a/b
            nonlocal generic_cond
            generic_cond = Ne(b, 0)
            return base_, S.One
        raise NoMatch

    try:
        base, exp_ = _get_base_exp(integrand)
    except NoMatch:
        return
    if generic_cond is S.true:
        degenerate_step = None
    else:
        # equivalent with subs(b, 0) but no need to find b
        degenerate_step = ConstantRule(integrand.subs(x, 0), x)
    generic_step = NestedPowRule(integrand, x, base, exp_)
    return _add_degenerate_step(generic_cond, generic_step, degenerate_step)


def inverse_trig_rule(integral: IntegralInfo, degenerate=True):
    """
    Set degenerate=False on recursive call where coefficient of quadratic term
    is assumed non-zero.
    """
    integrand, symbol = integral
    base, exp = integrand.as_base_exp()
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    c = Wild('c', exclude=[symbol, 0])
    match = base.match(a + b*symbol + c*symbol**2)

    if not match:
        return

    def make_inverse_trig(RuleClass, a, sign_a, c, sign_c, h) -> Rule:
        u_var = Dummy("u")
        rewritten = 1/sqrt(sign_a*a + sign_c*c*(symbol-h)**2)  # a>0, c>0
        quadratic_base = sqrt(c/a)*(symbol-h)
        constant = 1/sqrt(c)
        u_func = None
        if quadratic_base is not symbol:
            u_func = quadratic_base
            quadratic_base = u_var
        standard_form = 1/sqrt(sign_a + sign_c*quadratic_base**2)
        substep = RuleClass(standard_form, quadratic_base)
        if constant != 1:
            substep = ConstantTimesRule(constant*standard_form, symbol, constant, standard_form, substep)
        if u_func is not None:
            substep = URule(rewritten, symbol, u_var, u_func, substep)
        if h != 0:
            substep = CompleteSquareRule(integrand, symbol, rewritten, substep)
        return substep

    a, b, c = [match.get(i, S.Zero) for i in (a, b, c)]
    generic_cond = Ne(c, 0)
    if not degenerate or generic_cond is S.true:
        degenerate_step = None
    elif b.is_zero:
        degenerate_step = ConstantRule(a ** exp, symbol)
    else:
        degenerate_step = yield from sqrt_fractional_linear_rule(IntegralInfo((a + b * symbol) ** exp, symbol))

    if simplify(2*exp + 1) == 0:
        h, k = -b/(2*c), a - b**2/(4*c)  # rewrite base to k + c*(symbol-h)**2
        non_square_cond = Ne(k, 0)
        square_step = None
        if non_square_cond is not S.true:
            square_step = NestedPowRule(1/sqrt(c*(symbol-h)**2), symbol, symbol-h, S.NegativeOne)
        if non_square_cond is S.false:
            return square_step
        generic_step = ReciprocalSqrtQuadraticRule(integrand, symbol, a, b, c)
        step = _add_degenerate_step(non_square_cond, generic_step, square_step)
        if k.is_real and c.is_real:
            # list of ((rule, base_exp, a, sign_a, b, sign_b), condition)
            rules: list[tuple[Rule, Boolean]] = []
            for args, cond in (  # don't apply ArccoshRule to x**2-1
                ((ArcsinRule, k, 1, -c, -1, h), And(k > 0, c < 0)),  # 1-x**2
                ((ArcsinhRule, k, 1, c, 1, h), And(k > 0, c > 0)),  # 1+x**2
            ):
                if cond is S.true:
                    return make_inverse_trig(*args)
                if cond is not S.false:
                    rules.append((make_inverse_trig(*args), cond))
            if rules:
                if not k.is_positive:  # conditions are not thorough, need fall back rule
                    rules.append((generic_step, S.true))
                step = PiecewiseRule(integrand, symbol, rules)
            else:
                step = generic_step
        return _add_degenerate_step(generic_cond, step, degenerate_step)
    if exp == S.Half:
        step = SqrtQuadraticRule(integrand, symbol, a, b, c)
        return _add_degenerate_step(generic_cond, step, degenerate_step)


def add_rule(integral):
    integrand, symbol = integral
    results = []
    for g in integrand.as_ordered_terms():
        results.append((yield IntegralInfo(g, symbol)))
    return None if None in results else AddRule(integrand, symbol, results)


def mul_rule(integral: IntegralInfo):
    integrand, symbol = integral

    # Constant times function case
    coeff, f = integrand.as_independent(symbol)
    if coeff != 1:
        next_step = yield IntegralInfo(f, symbol)
        if next_step is not None:
            return ConstantTimesRule(integrand, symbol, coeff, f, next_step)


special_error_functions = (erf, erfc, erfi, fresnelc, fresnels, Ci, Chi, Si, Shi, Ei, li)


def _parts_rule_gen(integrand, symbol):
    # LIATE rule:
    # log, inverse trig, algebraic, trigonometric, exponential
    def pull_out_algebraic(integrand):
        integrand = integrand.cancel().together()
        # iterating over Piecewise args would not work here
        algebraic = ([] if isinstance(integrand, Piecewise) or not integrand.is_Mul
            else [arg for arg in integrand.args if arg.is_algebraic_expr(symbol)])
        if algebraic:
            u = Mul(*algebraic)
            dv = (integrand / u).cancel()
            return u, dv

    def pull_out_u(*functions) -> Callable[[Expr], tuple[Expr, Expr] | None]:
        def pull_out_u_rl(integrand: Expr) -> tuple[Expr, Expr] | None:
            if any(integrand.has(f) for f in functions):
                args = [arg for arg in integrand.args
                        if any(isinstance(arg, cls) for cls in functions)]
                if args:
                    u = Mul(*args) # type: ignore
                    dv = integrand / u
                    return u, dv
            return None

        return pull_out_u_rl

    def pull_out_dv(*functions) -> Callable[[Expr], tuple[Expr, Expr] | None]:
        # Prefer forms that are easier to integrate using special functions
        # x*exp(-x**2) instead of exp(-x**2) -> erf
        # x*sin(x**2) instead of sin(x**2) -> Fresnel
        def pull_out_dv_rl(integrand: Expr) -> tuple[Expr, Expr] | None:
            power = integrand.as_powers_dict().get(symbol)
            if (
                isinstance(power, Integer) and power >= 2 and
                integrand.has(*functions)
            ):
                for target in integrand.args:
                    if not any(isinstance(target, cls) for cls in functions):
                        continue
                    inner = target.args[0]
                    if (
                        inner.is_polynomial(symbol) and  # type: ignore
                        degree(inner, symbol) == 2
                    ):
                        dv = target * symbol
                        u = integrand / dv
                        return u, dv
            return None

        return pull_out_dv_rl

    liate_rules = [pull_out_dv(exp), pull_out_dv(sin, cos),
                   pull_out_u(*special_error_functions), pull_out_u(log),
                   pull_out_u(*inverse_trig_functions), pull_out_algebraic,
                   pull_out_u(sin, cos), pull_out_u(sinh, cosh),
                   pull_out_u(exp)]


    dummy = Dummy("temporary")
    # we can integrate log(x), atan(x), and erf by setting dv = 1
    if isinstance(integrand, (log, *inverse_trig_functions,
                              *special_error_functions)):
        integrand = dummy * integrand

    for index, rule in enumerate(liate_rules):
        result = rule(integrand)

        if result:
            u, dv = result

            # Don't pick u to be a constant if possible
            if symbol not in u.free_symbols and not u.has(dummy):
                return None

            u = u.subs(dummy, 1)
            dv = dv.subs(dummy, 1)

            # Don't pick a non-polynomial algebraic to be differentiated
            if rule == pull_out_algebraic and not u.is_polynomial(symbol):
                return None
            # Don't trade one logarithm or special error function for another
            if isinstance(u, (log, *special_error_functions)):
                rec_dv = 1/dv
                if (rec_dv.is_polynomial(symbol) and
                    degree(rec_dv, symbol) == 1):
                    return None

            # Can integrate a polynomial times OrthogonalPolynomial
            if rule == pull_out_algebraic:
                if dv.is_Derivative or dv.has(TrigonometricFunction, HyperbolicFunction) or \
                        isinstance(dv, OrthogonalPolynomial):
                    v_step = yield IntegralInfo(dv, symbol)
                    if v_step.contains_dont_know():
                        return None
                    else:
                        du = u.diff(symbol)
                        v = v_step.eval()
                        return u, dv, v, du, v_step

            # make sure dv is amenable to integration
            accept = False
            cutoff = liate_rules.index(pull_out_algebraic)  # log and inverse trig are usually worth trying
            if index < cutoff:
                accept = True
            elif (rule == pull_out_algebraic and dv.args and
                all(isinstance(a, (sin, cos, exp))
                for a in dv.args)):
                    accept = True
            else:
                for lrule in liate_rules[index + 1:]:
                    r = lrule(integrand)
                    if r and r[0].subs(dummy, 1).equals(dv):
                        accept = True
                        break

            if accept:
                du = u.diff(symbol)
                v_step = yield IntegralInfo(simplify(dv), symbol)
                if not v_step.contains_dont_know():
                    v = v_step.eval()
                    return u, dv, v, du, v_step
    return None


def _parts_rule(integrand, symbol) -> tuple[Expr, Expr, Expr, Expr, Rule] | None:
    """Standalone form of the LIATE search: drive :func:`_parts_rule_gen`
    with a fresh solver."""
    return IntegrationSolver().run_generator(_parts_rule_gen(integrand, symbol))


def parts_rule(integral):
    integrand, symbol = integral
    constant, integrand = integrand.as_coeff_Mul()

    result = yield from _parts_rule_gen(integrand, symbol)

    steps = []
    if result:
        u, dv, v, du, v_step = result
        debug("u : {}, dv : {}, v : {}, du : {}, v_step: {}".format(u, dv, v, du, v_step))
        steps.append(result)

        if isinstance(v, Integral):
            return

        # Set a limit on the number of times u can be used
        if isinstance(u, (sin, cos, exp, sinh, cosh)):
            cachekey = u.xreplace({symbol: _cache_dummy})
            if not (yield PartsUCheck(cachekey)):
                return

        # Try cyclic integration by parts a few times
        for _ in range(4):
            if dv == 1:
                break
            debug("Cyclic integration {} with v: {}, du: {}, integrand: {}".format(_, v, du, integrand))
            coefficient = ((v * du) / integrand).cancel()
            if coefficient == 1:
                break
            if symbol not in coefficient.free_symbols:
                rule = CyclicPartsRule(integrand, symbol,
                    [PartsRule(None, None, u, dv, v_step, None)
                     for (u, dv, v, du, v_step) in steps],
                    (-1) ** len(steps) * coefficient)
                if (constant != 1) and rule:
                    rule = ConstantTimesRule(constant * integrand, symbol, constant, integrand, rule)
                return rule

            # _parts_rule is sensitive to constants, factor it out
            next_constant, next_integrand = (v * du).as_coeff_Mul()
            result = yield from _parts_rule_gen(next_integrand, symbol)

            if result:
                u, dv, v, du, v_step = result
                u *= next_constant
                du *= next_constant
                steps.append((u, dv, v, du, v_step))
            else:
                break

    def make_second_step(steps, integrand):
        if steps:
            u, dv, v, du, v_step = steps[0]
            substep = yield from make_second_step(steps[1:], v * du)
            return PartsRule(integrand, symbol, u, dv, v_step, substep)
        return (yield IntegralInfo(integrand, symbol))

    if steps:
        u, dv, v, du, v_step = steps[0]
        second_step = yield from make_second_step(steps[1:], v * du)
        rule = PartsRule(integrand, symbol, u, dv, v_step, second_step)
        if (constant != 1) and rule:
            rule = ConstantTimesRule(constant * integrand, symbol, constant, integrand, rule)
        return rule


def trig_rule(integral):
    integrand, symbol = integral
    if integrand == sin(symbol):
        return SinRule(integrand, symbol)
    if integrand == cos(symbol):
        return CosRule(integrand, symbol)
    if integrand == sec(symbol)**2:
        return Sec2Rule(integrand, symbol)
    if integrand == csc(symbol)**2:
        return Csc2Rule(integrand, symbol)

    if isinstance(integrand, tan):
        rewritten = sin(*integrand.args) / cos(*integrand.args)
    elif isinstance(integrand, cot):
        rewritten = cos(*integrand.args) / sin(*integrand.args)
    elif isinstance(integrand, sec):
        arg = integrand.args[0]
        rewritten = ((sec(arg)**2 + tan(arg) * sec(arg)) /
                     (sec(arg) + tan(arg)))
    elif isinstance(integrand, csc):
        arg = integrand.args[0]
        rewritten = ((csc(arg)**2 + cot(arg) * csc(arg)) /
                     (csc(arg) + cot(arg)))
    else:
        return

    substep = yield IntegralInfo(rewritten, symbol)
    return RewriteRule(integrand, symbol, rewritten, substep)

def trig_product_rule(integral: IntegralInfo):
    integrand, symbol = integral
    if integrand == sec(symbol) * tan(symbol):
        return SecTanRule(integrand, symbol)
    if integrand == csc(symbol) * cot(symbol):
        return CscCotRule(integrand, symbol)


def trig_product_to_sum_rule(integral: IntegralInfo):
    """
    Rewrite a product of sines and cosines (or hyperbolic sines and
    cosines) of at least two different linear arguments with the
    product-to-sum identities

        sin(A)*sin(B) = cos(A - B)/2 - cos(A + B)/2
        sin(A)*cos(B) = sin(A + B)/2 + sin(A - B)/2
        cos(A)*cos(B) = cos(A - B)/2 + cos(A + B)/2

    (with cosh(A - B) and cosh(A + B) exchanged in the hyperbolic
    versions), e.g. sin(5*x)*cos(3*x) = sin(8*x)/2 + sin(2*x)/2.
    Applying the identity to one pair at a time and recursing reduces
    any such product to a sum of single sines and cosines.
    """
    integrand, x = integral

    coefficient = S.One
    factors = []
    for factor in Mul.make_args(integrand):
        if isinstance(factor, exp):
            # exp(2*x).as_base_exp() would give (E, 2*x)
            base, power = factor, S.One  # type: tuple[Expr, Expr]
        else:
            base, power = factor.as_base_exp()
        if (isinstance(base, (sin, cos, sinh, cosh, exp))
                and power.is_Integer and power > 0
                and not base.args[0].diff(x).has(x)):
            factors.append((type(base), base.args[0], power))
        else:
            # other factors (polynomials, other functions) pass through
            coefficient *= factor

    circular = {f for f, _, _ in factors} & {sin, cos}
    hyperbolic_ = {f for f, _, _ in factors} & {sinh, cosh, exp}
    if circular and hyperbolic_:
        return None
    if not circular and not hyperbolic_ - {exp}:
        # exp alone combines by itself; require a sinh or cosh to pair
        return None
    if sum(power for _, _, power in factors) < 2:
        return None
    if len({argument for _, argument, _ in factors}) < 2:
        # only reduce same-argument powers for a cofactor free of further
        # trigonometric functions: rewriting powers under a trigonometric
        # cofactor can cycle against the expansion rules, and the plain
        # sin(u)**m*cos(u)**n products have dedicated power rules
        if coefficient.has(TrigonometricFunction, HyperbolicFunction):
            return None
        if circular and not coefficient.has(x):
            return None

    def combined(f1, A, f2, B):
        rank = {sin: 0, sinh: 0, cos: 1, cosh: 1, exp: 2}
        if rank[f1] > rank[f2]:
            f1, f2, A, B = f2, f1, B, A
        identities = {
            (sin, sin): (cos(A - B) - cos(A + B))/2,
            (sin, cos): (sin(A + B) + sin(A - B))/2,
            (cos, cos): (cos(A - B) + cos(A + B))/2,
            (sinh, sinh): (cosh(A + B) - cosh(A - B))/2,
            (sinh, cosh): (sinh(A + B) + sinh(A - B))/2,
            (cosh, cosh): (cosh(A + B) + cosh(A - B))/2,
            (sinh, exp): (exp(A + B) - exp(B - A))/2,
            (cosh, exp): (exp(A + B) + exp(B - A))/2,
        }
        return identities[(f1, f2)]

    # combine the first factor with the first one of different argument,
    # or with another same-argument factor in the hyperbolic case
    f1, arg1, _ = factors[0]
    j = next((i for i, (_, argument, _) in enumerate(factors)
              if argument != arg1), None)
    if j is None:
        j = 1 if len(factors) > 1 else 0
    f2, arg2, _ = factors[j]
    rewritten = coefficient*combined(f1, arg1, f2, arg2)
    for i, (f, argument, power) in enumerate(factors):
        used = 2 if i == 0 == j else (1 if i in (0, j) else 0)
        rewritten *= f(argument)**(power - used)
    rewritten = rewritten.expand()
    substep = yield IntegralInfo(rewritten, x)
    return RewriteRule(integrand, x, rewritten, substep)


class TabularPartsRule(AtomicRule):
    """integrate(P(x)*f(a + b*x), x) for a polynomial P and f one of sin,
    cos, sinh, cosh, exp, by tabular integration by parts: alternately
    differentiate P and antidifferentiate f until P vanishes,

        Sum((-1)**k * P^(k)(x) * F_{k+1}(a + b*x) / b**(k+1), (k, 0, deg(P)))

    where F_{k} denotes the k-th antiderivative of f.
    """

    __slots__ = ("poly_part", "func_part")

    poly_part: Expr
    func_part: Expr

    def __init__(
        self, integrand: Expr, variable: Symbol, poly_part: Expr, func_part: Expr
    ) -> None:
        super().__init__(integrand, variable)
        self.poly_part = poly_part
        self.func_part = func_part

    def eval(self) -> Expr:
        if isinstance(self.func_part, Pow) and not isinstance(self.func_part, exp):
            u = cast('Expr', self.func_part.exp)
        else:
            u = cast('Expr', self.func_part.args[0])
        b = u.diff(self.variable)

        def antiderivative(F):
            coefficient, g = F.as_independent(self.variable)
            if isinstance(g, sin):
                return -coefficient*cos(g.args[0])
            if isinstance(g, cos):
                return coefficient*sin(g.args[0])
            if isinstance(g, sinh):
                return coefficient*cosh(g.args[0])
            if isinstance(g, cosh):
                return coefficient*sinh(g.args[0])
            if isinstance(g, Pow) and not isinstance(g, exp):
                # F**u integrates to F**u/log(F) with respect to u
                return coefficient*g/log(g.base)
            return F  # exp

        P = self.poly_part
        F = self.func_part
        total = S.Zero
        sign = S.One
        k = 1
        while P != 0:
            F = antiderivative(F)
            total += sign*P*F/b**k
            P = P.diff(self.variable)
            sign = -sign
            k += 1
        return total.expand()


def tabular_parts_rule(integral: IntegralInfo):
    """Match P(x)*f(a + b*x) for TabularPartsRule, with a Piecewise for
    a possibly vanishing b as in the other linear-argument rules."""
    integrand, x = integral

    func_part = None
    poly_factors = []
    for factor in Mul.make_args(integrand):
        if isinstance(factor, exp):
            # exp(2*x).as_base_exp() would give (E, 2*x)
            base, power = factor, S.One  # type: tuple[Expr, Expr]
        else:
            base, power = factor.as_base_exp()
        if (isinstance(base, (sin, cos, sinh, cosh, exp)) and power == 1
                and base.has(x) and not base.args[0].diff(x).has(x)):
            if func_part is not None:
                return None
            func_part = base
            continue
        if (isinstance(factor, Pow) and not factor.base.has(x)
                and factor.exp.has(x) and not factor.exp.diff(x).has(x)
                and (factor.base.is_positive or factor.base.is_Symbol)):
            # a general-base exponential F**(a + b*x)
            if func_part is not None:
                return None
            func_part = factor
            continue
        poly_factors.append(factor)
    if func_part is None:
        return None
    poly_part = Mul(*poly_factors)
    poly = poly_part.as_poly(x)
    if poly is None or poly.degree() < 1:
        return None

    if isinstance(func_part, Pow) and not isinstance(func_part, exp):
        u = func_part.exp
    else:
        u = func_part.args[0]
    b = u.diff(x)
    generic_step = TabularPartsRule(integrand, x, poly_part, func_part)
    generic_cond = Ne(b, 0)
    if generic_cond is S.true:
        return generic_step
    if isinstance(func_part, Pow) and not isinstance(func_part, exp):
        zero_integrand = poly_part*func_part.base**(u - b*x)
    else:
        zero_integrand = poly_part*func_part.func(u - b*x)
    zero_step = yield IntegralInfo(zero_integrand, x)
    return PiecewiseRule(integrand, x,
                         [(zero_step, Eq(b, 0)), (generic_step, S.true)])


def trig_cmplx_exp_rule(integral: IntegralInfo):
    """
    Strategy that rewrites sin, cos, sinh, and cosh in terms of complex exponentials.
    Useful for integration techniques that handle exponentials better.
    Applies only when the integrand belongs to a class that benefits from exponential rewriting,
    such as combinations involving Gaussian exponentials.

    sin(x)  -> (exp(i*x) - exp(-i*x)) / (2*i)
    cos(x)  -> (exp(i*x) + exp(-i*x)) / 2
    sinh(x) -> (exp(x) - exp(-x)) / 2
    cosh(x) -> (exp(x) + exp(-x)) / 2
    """
    integrand, symbol = integral

    if not integrand.has(exp) and not integrand.has(sin, cos, sinh, cosh):
        return

    a = Wild('a', exclude=[symbol, 0])
    b = Wild('b', exclude=[symbol])
    c = Wild('c', exclude=[symbol])
    # n = Wild('n', exclude=[symbol], properties=[lambda n: n > 0])
    f = WildFunction('f')
    guassian_pattern = exp(a * symbol**2 + b * symbol + c)
    trigexp_over_x_pattern = f*exp(a * symbol)/symbol
    trigexp_over_x_match = integrand.match(trigexp_over_x_pattern)
    if not (any(term.match(guassian_pattern) for term in integrand.atoms(exp))
            or (trigexp_over_x_match and
                trigexp_over_x_match[f].has(sin, cos, sinh, cosh))):
        return

    # Replace trig and hyperbolic functions with their exponential forms
    rewritten = integrand.rewrite([sin, cos, sinh, cosh], exp)

    if rewritten != integrand:
        debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten, symbol))
        steps = yield IntegralInfo(rewritten, symbol)
        return RewriteRule(integrand, symbol, rewritten, steps)


def quadratic_denom_rule(integral):
    integrand, symbol = integral
    if not integrand.is_rational_function(symbol):
        return None
    num, den = integrand.as_numer_denom()
    if den == 1:
        return None
    # Prevents things like c*(c*x + d)**n from hiding the power
    den_const, den_x = den.as_independent(symbol, as_Add=False)
    if den_const != 1:
        num = num / den_const
    den = den_x
    if den.is_Pow:
        q = den.base
        n = den.exp
    else:
        q = den
        n = S.One
    den_poly = Poly(q, symbol)
    num_poly = Poly(num, symbol)
    deg_den = den_poly.degree()
    deg_num = num_poly.degree()
    # TODO: may add n = 1 here and do a general manual rational integration rule
    # instead of letting general substitution rule find the pattern
    if deg_den != 2:
        return None
    if deg_num >= deg_den:
        return None

    def _arctan_match(B, a, c, symbol, degenerate=True):
        # integrates B / a*x**2 + c
        integrand = B / (a*symbol**2 + c)
        pieces = []
        # skips degenerate case if a != 0 or if a = 0 would cause null denominator
        if degenerate and not _if_zero_implies_zero(a, c):
            substituted = integrand.subs(a, 0)
            substep = yield IntegralInfo(substituted, symbol)
            pieces.append((RewriteRule(integrand, symbol, substituted, substep), Eq(a, 0)))
        if degenerate and not _if_zero_implies_zero(c, a):
            substituted = integrand.subs(c, 0)
            substep = yield IntegralInfo(substituted, symbol)
            pieces.append((RewriteRule(integrand, symbol, substituted, substep), Eq(c, 0)))
        if a.is_extended_real and c.is_extended_real:
            positive_cond = c/a > 0
            if positive_cond is not S.true:
                coeff = B/(2*sqrt(-c)*sqrt(a))
                constant = sqrt(-c/a)
                r1 = 1/(symbol-constant)
                r2 = 1/(symbol+constant)
                log_steps = [ReciprocalRule(r1, symbol, symbol-constant),
                            ConstantTimesRule(-r2, symbol, -1, r2, ReciprocalRule(r2, symbol, symbol+constant))]
                rewritten = sub = r1 - r2
                negative_step = AddRule(sub, symbol, log_steps)
                if coeff != 1:
                    rewritten = Mul(coeff, sub, evaluate=False)
                    negative_step = ConstantTimesRule(rewritten, symbol, coeff, sub, negative_step)
                negative_step = RewriteRule(integrand, symbol, rewritten, negative_step)
                if positive_cond is S.false:
                    pieces.append((negative_step, S.true))
                    return PiecewiseRule(integrand, symbol, pieces)
                else:
                    pieces.append((negative_step, c / a < 0))
        general_rule = ArctanRule(integrand, symbol, B, a, c)
        if pieces:
            pieces.append((general_rule, S.true))
            return PiecewiseRule(integrand, symbol, pieces)
        return general_rule


    def _complete_square(B, a, b, c, n, symbol, degenerate_a=True, degenerate_discriminant=True):
        # integrates B / (a*x**2 + b*x + c)**n
        pieces = []
        discriminant = 4*a*c - b**2
        denominator = a*symbol**2 + b*symbol + c
        integrand = B / denominator**n
        # degenerate flags avoid recalculating Piecewise branches recursively
        if degenerate_a and not _if_zero_implies_zero(a, denominator):
            substituted = integrand.subs(a, 0)
            substep = yield IntegralInfo(substituted, symbol)
            pieces.append((RewriteRule(integrand, symbol, substituted, substep), Eq(a, 0)))
        if degenerate_discriminant and not _if_zero_implies_zero(discriminant, denominator):
            u = Dummy("u")
            # we divide by a, Piecewise condition above
            u_func = symbol + b/(2*a)
            rewritten = (B/a**n) * u_func**(-2*n)
            subexpr = (B/a**n) * u**(-2*n)
            substep = yield IntegralInfo(subexpr, u)
            rule = RewriteRule(integrand, symbol, rewritten, URule(rewritten, symbol, u, u_func, substep))
            if discriminant.is_zero:
                if pieces:
                    pieces.append((rule, S.true))
                    return PiecewiseRule(integrand, symbol, pieces)
                return rule
            pieces.append((rule, Eq(discriminant, 0)))
        if n == 1:
            # base case, B / (a*x**2 + b*x + c), solve by substitution with _arctan_match
            u = Dummy("u")
            u_func = symbol + b/(2*a)
            # we put degenerate = False since after substitution, the integrand becomes B/(a*u**2 + discriminant/(4*a)),
            # then the _arctan_match conditions (a != 0 and discriminant !=0) are already computed
            substep = yield from _arctan_match(B, a, discriminant/(4*a), u, degenerate=False)
            general_step = URule(integrand, symbol, u, u_func, substep)
        else:
            # reduction step for B/q**n
            # Differentiate B*T/((n - 1)*discriminant*q**(n - 1)), with T = q'(x), using T**2 = 4*a*q - discriminant,
            # this derivative equals B/q**n - coeff*B/q**(n - 1), so the remaining integral has power n - 1
            T = 2*a*symbol + b
            # we divide by discriminant, Piecewise condition above
            F = B*T / ((n - 1)*discriminant*denominator**(n - 1))
            derivative = Derivative(F, symbol, evaluate=False)
            coeff = 2*a*(2*n - 3) / ((n - 1)*discriminant)
            remainder = B / denominator**(n - 1)
            scaled_remainder = coeff * remainder
            rewritten = derivative + scaled_remainder
            derivative_step = DerivativeRule(derivative, symbol)
            remainder_step = yield from _complete_square(B, a, b, c, n - 1, symbol, degenerate_a=False, degenerate_discriminant=False)
            scaled_step = ConstantTimesRule(scaled_remainder, symbol, coeff, remainder, remainder_step)
            add_step = AddRule(rewritten, symbol, [derivative_step, scaled_step])
            general_step = RewriteRule(integrand, symbol, rewritten, add_step)
        if pieces:
            pieces.append((general_step, S.true))
            return PiecewiseRule(integrand, symbol, pieces)
        return general_step

    def _split_sum(A, B, a, b, c, n, symbol):
        # integrates (A*x + B) / (a*x**2 + b*x + c)**n. Split A*x + B as alpha*q'(x) + beta,
        # then integrate the two terms separately (first by substitution, second with _complete_square)
        pieces = []
        denominator = (a*symbol**2 + b*symbol + c)
        integrand = (A*symbol + B) / denominator**n
        if not _if_zero_implies_zero(a, denominator):
            substituted = integrand.subs(a, 0)
            substep = yield IntegralInfo(substituted, symbol)
            pieces.append((RewriteRule(integrand, symbol, substituted, substep), Eq(a, 0)))
        # we divide by a, Piecewise condition above
        const =  A/(2*a)
        numer1 =  (2*a*symbol + b)
        numer2 = - const*b + B
        qprime_part = numer1 / denominator**n
        u = Dummy('u')
        u_substep = yield IntegralInfo(u**(-n), u)
        step1 = URule(qprime_part, symbol, u, denominator, u_substep)
        if const != 1:
            step1 = ConstantTimesRule(const*qprime_part, symbol, const, qprime_part, step1)
        if numer2.is_zero:
            rewritten = const*qprime_part
            general_step = RewriteRule(integrand, symbol, rewritten, step1)
        else:
            # since degenerate a condition is already computed, degenerate_a = False
            step2 = yield from _complete_square(numer2, a, b, c, n, symbol, degenerate_a=False)
            rewritten = const*qprime_part + numer2/denominator**n
            substeps = AddRule(rewritten, symbol, [step1, step2])
            general_step = RewriteRule(integrand, symbol, rewritten, substeps)
        if pieces:
            pieces.append((general_step, S.true))
            return PiecewiseRule(integrand, symbol, pieces)
        return general_step

    B = num_poly.nth(0)
    a = den_poly.nth(2)
    b = den_poly.nth(1)
    c = den_poly.nth(0)

    normalized_num = num_poly.as_expr()
    normalized_den = den_poly.as_expr()
    normalized_integrand = normalized_num / normalized_den**n

    if b == 0 and deg_num == 0 and n == 1:
        step = yield from _arctan_match(B, a, c, symbol)
    elif deg_num == 1:
        A = num_poly.nth(1)
        step = yield from _split_sum(A, B, a, b, c, n, symbol)
    else:
        step = yield from _complete_square(B, a, b, c, n, symbol)

    if normalized_integrand != integrand:
        step = RewriteRule(integrand, symbol, normalized_integrand, step)

    return step

class RectifyAtanRule(Rule):
    """Wrap the result of a tangent substitution, rewriting each term
    c*atan(w*tan(theta) + v) with numeric v**2 < 4*w as

        c*(theta + atan(((w - 1)*sin(2*theta) + v*(1 + cos(2*theta)))
           /((1 + w) + v*sin(2*theta) + (1 - w)*cos(2*theta))))

    which has the same derivative but no jumps at the poles of the
    tangent, so that the antiderivative is valid on domains of maximum
    extent (D.J. Jeffrey, 1993).  For example 2*atan(3*tan(x/2)) from
    3/(5 - 4*cos(x)) becomes x + 2*atan(sin(x)/(2 - cos(x))).
    """

    __slots__ = ("theta", "substep")

    theta: Expr
    substep: Rule

    def __init__(self, integrand: Expr, variable: Symbol,
                 theta: Expr, substep: Rule) -> None:
        super().__init__(integrand, variable)
        self.theta = theta
        self.substep = substep

    def eval(self) -> Expr:
        result = self.substep.eval()
        theta = self.theta
        if not result.has(atan):
            return result
        T = Dummy("T")
        double = 2*theta

        def rectify(term):
            coefficient, at = term.as_independent(atan)
            if not isinstance(at, atan) or coefficient.has(self.variable):
                return term
            argument = at.args[0]
            arg_poly = argument.subs(tan(theta), T).as_poly(T)
            cotangent = False
            if arg_poly is None or arg_poly.degree() != 1 or arg_poly.has(tan):
                # try an argument linear in cot(theta) = 1/tan(theta)
                arg_poly = argument.subs(tan(theta), 1/T).as_poly(T)
                cotangent = True
                if (arg_poly is None or arg_poly.degree() != 1
                        or arg_poly.has(tan)):
                    return term
            v, w = arg_poly.nth(0), arg_poly.nth(1)
            if v.has(self.variable) or w.has(self.variable):
                return term
            if (v**2 - 4*w).is_negative is not True:
                return term
            if cotangent:
                return coefficient*(-theta + atan(
                    ((w - 1)*sin(double) + v*(1 - cos(double)))
                    /((1 + w) + v*sin(double) - (1 - w)*cos(double))))
            return coefficient*(theta + atan(
                ((w - 1)*sin(double) + v*(1 + cos(double)))
                /((1 + w) + v*sin(double) + (1 - w)*cos(double))))

        return Add(*[rectify(t) for t in Add.make_args(result)])

    def contains_dont_know(self) -> bool:
        return self.substep.contains_dont_know()


def weierstrass_substitution(integral):
    # Look for a rational function in trigonometric functions
    integrand, x = integral

    TRIG = (sin, cos, tan, cot, sec, csc)

    trig_atoms = tuple(ordered(integrand.atoms(*TRIG)))

    if not trig_atoms:
        return None
    if integrand.is_rational_function(*trig_atoms) is not True:
        return None
    masked = integrand.xreplace({atom: Dummy() for atom in trig_atoms})
    # exclude x*sin(x)
    if masked.has(x):
        return None

    variable_atoms = []
    for atom in trig_atoms:
        argument = atom.args[0]
        # exclude sin(cos(x))
        if not argument.is_polynomial(x):
            return None
        argument = Poly(argument, x)
        # exclude sin(2*x**2)
        if argument.degree() > 1:
            return None
        # skip constants
        if argument.degree() == 0:
            continue
        coeff = argument.nth(1)
        phase = argument.nth(0)
        variable_atoms.append((atom, coeff, phase))

    if not variable_atoms:
        return None

    coeff0 = variable_atoms[0][1]
    phase0 = variable_atoms[0][2]
    proportional_phases = True
    ratios = []
    for _, coeff, phase in variable_atoms:
        ratio = (coeff / coeff0).cancel()
        # use this instead of is_rational to avoid taking 3/a
        # if a is declared rational
        if not isinstance(ratio, Rational):
            return None
        ratios.append(ratio)
        if proportional_phases:
            determinant = (coeff0*phase - coeff*phase0).cancel()
            if determinant.is_zero is not True:
                proportional_phases = False

    denominator_lcm = lcm_list([ratio.q for ratio in ratios])

    # Choose the largest common frequency omega such that every
    # coefficient is an integer multiple of omega
    omega = (coeff0 / denominator_lcm).cancel()

    # If all (coefficient, phase) pairs are proportional, absorb the
    # common phase into u. For example, 2*x + 2 and 4*x + 4 become u and 2*u
    if proportional_phases:
        reference_harmonic = (coeff0 / omega).cancel()
        phase_shift = (phase0 / reference_harmonic).cancel()
    else:
        phase_shift = S.Zero
    u_func = omega*x + phase_shift
    u = Dummy("u")
    replacements = {}
    ODD_TRIG = (sin, tan, cot, csc)
    for atom, coeff, phase in variable_atoms:
        harmonic = (coeff / omega).cancel()
        if proportional_phases:
            argument = harmonic*u
        else:
            argument = harmonic*u + phase
        if harmonic.is_negative:
            # Move the minus sign outside to help expand_trig work on the first attempt
            # For example, tan(-2*u + b) would not expand correctly with the minus sign inside
            positive_argument = -argument
            replacement = atom.func(positive_argument)
            if atom.func in ODD_TRIG:
                replacement = -replacement
        else:
            replacement = atom.func(argument)
        replacements[atom] = replacement
    expr_u = integrand.xreplace(replacements)
    # rewrites sin, cos, tan to substitute, for example
    # cos(2*u + 2) as (2*cos(u)**2 - 1)*cos(2) - 2*sin(2)*sin(u)*cos(u)
    expr_u = expand_trig(expr_u)

    zero_replacements = {atom: atom.func(phase) for atom, _, phase in variable_atoms}
    zero_integrand = integrand.xreplace(zero_replacements)

    if omega.is_zero is True:
        zero_substep = yield IntegralInfo(zero_integrand, x)
        return RewriteRule(integrand, x, zero_integrand, zero_substep)
    t = Dummy("t")
    s = Dummy("s")
    c = Dummy("c")

    expr_sc = expr_u.xreplace({
        sin(u): s,
        cos(u): c,
        tan(u): s/c,
        cot(u): c/s,
        sec(u): 1/c,
        csc(u): 1/s,
    })

    def try_sin(w):
        # try t = sin(u_func) when R(s, -c) = -R(s, c), so R/(w*c)
        # is rational in s alone. This avoids poles introduced by tangent substitutions
        transformed = (expr_sc/(w*c)).cancel()
        numerator, denominator = transformed.as_numer_denom()

        numerator = numerator.as_poly(c)
        denominator = denominator.as_poly(c)

        replacements = {s: t}

        # powers of cos must be even so cos**(2*k) can be replaced by (1 - t**2)**k without branch problems
        for polynomial in (numerator, denominator):
            for (power,), coefficient in polynomial.terms():
                if power % 2 and coefficient != 0:
                    return None
                if power:
                    replacements[c**power] = (1 - t**2)**(power // 2)

        numerator = numerator.as_expr().xreplace(replacements)
        denominator = denominator.as_expr().xreplace(replacements)
        transformed = (numerator/denominator).cancel()
        substep = yield IntegralInfo(transformed, t)
        return URule(integrand, x, t, sin(u_func), substep)

    def try_cos(w):
        # try t = cos(u_func) when R(-s, c) = -R(s, c), so -R/(w*s) is rational in c alone,
        # this avoids poles introduced by tangent substitutions
        transformed = (-expr_sc/(w*s)).cancel()
        numerator, denominator = transformed.as_numer_denom()

        numerator = numerator.as_poly(s)
        denominator = denominator.as_poly(s)

        replacements = {c: t}

        # powers of sin must be even so sin**(2*k) can be replaced by (1 - t**2)**k without branch problems
        for polynomial in (numerator, denominator):
            for (power,), coefficient in polynomial.terms():
                if power % 2 and coefficient != 0:
                    return None
                if power:
                    replacements[s**power] = (1 - t**2)**(power // 2)

        numerator = numerator.as_expr().xreplace(replacements)
        denominator = denominator.as_expr().xreplace(replacements)
        transformed = (numerator/denominator).cancel()
        substep = yield IntegralInfo(transformed, t)
        return URule(integrand, x, t, cos(u_func), substep)

    def try_tan(w):
        # try t = tan(u_func) when R(-s, -c) = R(s, c),
        # this produces a lower-degree rational function than the half-angle substitution
        transformed = ((expr_sc*c**2/w).xreplace({s: t*c})).cancel()

        numerator, denominator = transformed.as_numer_denom()

        numerator = numerator.as_poly(c)
        denominator = denominator.as_poly(c)

        replacements = {}

        for polynomial in (numerator, denominator):
            for (power,), _ in polynomial.terms():
                if power % 2:
                    return None
                if power:
                    replacements[c**power] = (1 + t**2)**(-(power // 2))

        numerator = numerator.as_expr().xreplace(replacements)
        denominator = denominator.as_expr().xreplace(replacements)
        transformed = (numerator/denominator).cancel()
        substep = yield IntegralInfo(transformed, t)
        inner = URule(integrand, x, t, tan(u_func), substep)
        return RectifyAtanRule(integrand, x, u_func, inner)

    def try_tan_half(w):
        # fallback to the universal Weierstrass substitution t = tan(u_func/2)
        transformed = expr_u.xreplace({
            sin(u): 2*t/(1 + t**2),
            cos(u): (1 - t**2)/(1 + t**2),
            tan(u): 2*t/(1 - t**2),
            cot(u): (1 - t**2)/(2*t),
            sec(u): (1 + t**2)/(1 - t**2),
            csc(u): (1 + t**2)/(2*t),
        })

        transformed *= (2/(w*(1 + t**2))).cancel()
        substep = yield IntegralInfo(transformed, t)
        inner = URule(integrand, x, t, tan(u_func/2), substep)
        return RectifyAtanRule(integrand, x, u_func/2, inner)

    generic_step = yield from try_sin(omega)
    if generic_step is None:
        generic_step = yield from try_cos(omega)
    if generic_step is None:
        generic_step = yield from try_tan(omega)
    if generic_step is None:
        generic_step = yield from try_tan_half(omega)

    if generic_step is None or generic_step.contains_dont_know():
        return None

    if omega.is_zero is False:
        return generic_step

    zero_substep = yield IntegralInfo(zero_integrand, x)
    zero_step = RewriteRule(integrand, x, zero_integrand, zero_substep)

    return PiecewiseRule(integrand, x, [(zero_step, Eq(omega, 0)), (generic_step, S.true)])

def hyperbolic_rational_substitution(integral):
    """
    Integrate rational functions of hyperbolic functions of a common
    linear argument u = omega*x + phi.

    After rewriting the integrand as a rational function R(sinh(u), cosh(u)),
    the classic sequence of substitutions is tried:

    - R odd in sinh(u):  t = cosh(u), using sinh(u)**2 = cosh(u)**2 - 1
    - R odd in cosh(u):  t = sinh(u), using cosh(u)**2 = sinh(u)**2 + 1
    - R even in both:    t = tanh(u), using cosh(u)**2 = 1/(1 - t**2)
    - otherwise:         t = exp(u), using sinh(u) = (t - 1/t)/2 and
      cosh(u) = (t + 1/t)/2

    Each substitution turns the integral into that of a rational function
    of t.
    """
    integrand, x = integral

    HYPERBOLIC = (sinh, cosh, tanh, coth, sech, csch)
    hyp_atoms = tuple(ordered(integrand.atoms(*HYPERBOLIC)))
    if not hyp_atoms:
        return None
    if integrand.is_rational_function(*hyp_atoms) is not True:
        return None
    masked = integrand.xreplace({atom: Dummy() for atom in hyp_atoms})
    # exclude x*sinh(x)
    if masked.has(x):
        return None

    arguments = {atom.args[0] for atom in hyp_atoms}
    if len(arguments) != 1:
        return None
    u_func, = arguments
    omega = u_func.diff(x)
    # require a linear argument omega*x + phi
    if omega.is_zero or omega.has(x):
        return None

    t = Dummy("t")
    s = Dummy("s")
    c = Dummy("c")

    expr_sc = integrand.xreplace({
        sinh(u_func): s,
        cosh(u_func): c,
        tanh(u_func): s/c,
        coth(u_func): c/s,
        sech(u_func): 1/c,
        csch(u_func): 1/s,
    })

    def _replace_even_powers(transformed, var, square_image):
        # Rewrite transformed, a rational function even in var, replacing
        # var**2 by square_image; return None if odd powers remain
        numerator, denominator = transformed.as_numer_denom()
        numerator = numerator.as_poly(var)
        denominator = denominator.as_poly(var)
        if numerator is None or denominator is None:
            return None
        replacements = {}
        for polynomial in (numerator, denominator):
            for (power,), coefficient in polynomial.terms():
                if power % 2 and coefficient != 0:
                    return None
                if power:
                    replacements[var**power] = square_image**(power // 2)
        numerator = numerator.as_expr().xreplace(replacements)
        denominator = denominator.as_expr().xreplace(replacements)
        return (numerator/denominator).cancel()

    def try_cosh(w):
        # t = cosh(u) when R(-s, c) = -R(s, c), so R/(w*s) is even in s
        # and sinh(u)**2 becomes t**2 - 1; dt = w*sinh(u)*dx
        transformed = (expr_sc/(w*s)).cancel()
        transformed = _replace_even_powers(transformed, s, t**2 - 1)
        if transformed is None:
            return None
        transformed = transformed.xreplace({c: t})
        substep = yield IntegralInfo(transformed, t)
        return URule(integrand, x, t, cosh(u_func), substep)

    def try_sinh(w):
        # t = sinh(u) when R(s, -c) = -R(s, c), so R/(w*c) is even in c
        # and cosh(u)**2 becomes t**2 + 1; dt = w*cosh(u)*dx
        transformed = (expr_sc/(w*c)).cancel()
        transformed = _replace_even_powers(transformed, c, t**2 + 1)
        if transformed is None:
            return None
        transformed = transformed.xreplace({s: t})
        substep = yield IntegralInfo(transformed, t)
        return URule(integrand, x, t, sinh(u_func), substep)

    def try_tanh(w):
        # t = tanh(u) when R(-s, -c) = R(s, c): R(t*c, c)*c**2 is then even
        # in c with cosh(u)**2 = 1/(1 - t**2), and dx = c**2*dt/w
        transformed = ((expr_sc*c**2/w).xreplace({s: t*c})).cancel()
        transformed = _replace_even_powers(transformed, c, 1/(1 - t**2))
        if transformed is None:
            return None
        substep = yield IntegralInfo(transformed, t)
        return URule(integrand, x, t, tanh(u_func), substep)

    def try_exp(w):
        # universal substitution t = exp(u), dx = dt/(w*t)
        transformed = (expr_sc.xreplace({
            s: (t - 1/t)/2,
            c: (t + 1/t)/2,
        })/(w*t)).cancel()
        substep = yield IntegralInfo(transformed, t)
        return URule(integrand, x, t, exp(u_func), substep)

    generic_step = None
    for attempt in (try_cosh, try_sinh, try_tanh, try_exp):
        generic_step = yield from attempt(omega)
        if generic_step is not None and not generic_step.contains_dont_know():
            break
    else:
        return None

    if omega.is_zero is False:
        return generic_step

    phase = u_func - omega*x
    zero_integrand = integrand.xreplace(
        {atom: atom.func(phase) for atom in hyp_atoms})
    zero_substep = yield IntegralInfo(zero_integrand, x)
    zero_step = RewriteRule(integrand, x, zero_integrand, zero_substep)

    return PiecewiseRule(integrand, x, [(zero_step, Eq(omega, 0)), (generic_step, S.true)])

def odd_power_trig_substitution(integral):
    """
    Integrate g(sin(u))*cos(u)**n with n an odd integer by substituting
    t = sin(u), where u is linear in the integration variable and g may
    involve arbitrary (fractional or symbolic) powers; likewise with sin
    and cos exchanged, and for the hyperbolic counterparts.

    For example cos(x)**3/sqrt(sin(x)) is g(sin(x))*cos(x)**3 with
    g(t) = 1/sqrt(t), and the substitution t = sin(x) turns it into
    (1 - t**2)/sqrt(t).

    Also integrate sin(u)**m*cos(u)**n when m + n is an even integer by
    substituting t = tan(u), which gives t**m*(1 + t**2)**(-(m + n)/2 - 1);
    for example sin(x)**(3/2)/cos(x)**(7/2) becomes t**(3/2).  The
    hyperbolic counterpart uses t = tanh(u) and 1 - t**2.

    This complements weierstrass_substitution, which handles integrands
    that are rational functions of the trigonometric functions.
    """
    integrand, x = integral

    trig_atoms = integrand.atoms(sin, cos, tan, cot, sec, csc)
    hyp_atoms = integrand.atoms(sinh, cosh, tanh, coth, sech, csch)
    if trig_atoms and not hyp_atoms:
        atoms = trig_atoms
        sin_f, cos_f, tan_f, cot_f, sec_f, csc_f = sin, cos, tan, cot, sec, csc
        hyperbolic = False
    elif hyp_atoms and not trig_atoms:
        atoms = hyp_atoms
        sin_f, cos_f, tan_f, cot_f, sec_f, csc_f = sinh, cosh, tanh, coth, sech, csch
        hyperbolic = True
    else:
        return None

    arguments = {atom.args[0] for atom in atoms}
    if len(arguments) != 1:
        return None
    u_func, = arguments
    omega = u_func.diff(x)
    # require a linear argument omega*x + phi
    if omega.is_zero or omega.has(x):
        return None

    t = Dummy("t")
    s = Dummy("s")
    c = Dummy("c")
    if hyperbolic:
        # cosh(u)**2 = 1 + t**2 for t = sinh(u); sinh(u)**2 = t**2 - 1 for
        # t = cosh(u); d(cosh(u)) has no minus sign, unlike d(cos(u))
        cos_squared = 1 + t**2
        sin_squared = t**2 - 1
        cos_branch_sign = S.One
    else:
        # cos(u)**2 = 1 - t**2 for t = sin(u) and vice versa
        cos_squared = 1 - t**2
        sin_squared = 1 - t**2
        cos_branch_sign = S.NegativeOne

    # express quotient functions through sin and cos, also inside bases of
    # non-integer powers like (b*sec(u))**(3/2)
    expr_sc = integrand.xreplace({
        tan_f(u_func): sin_f(u_func)/cos_f(u_func),
        cot_f(u_func): cos_f(u_func)/sin_f(u_func),
        sec_f(u_func): 1/cos_f(u_func),
        csc_f(u_func): 1/sin_f(u_func),
    }).xreplace({sin_f(u_func): s, cos_f(u_func): c})
    if expr_sc.has(x):
        return None

    def _attempt(kept, removed, square_image, u_of_t, sign):
        # integrand = g(kept)*removed**power: substitute t = kept(u), using
        # removed**2 = square_image and one power of removed for dt/w
        rest, dependent = expr_sc.as_independent(removed)
        base, power = dependent.as_base_exp()
        if base != removed or power.is_odd is not True:
            return None
        rest = rest.xreplace({kept: t})
        if rest.has(s, c):
            return None
        t_integrand = sign*rest*square_image**((power - 1)/2)/omega
        substep = yield IntegralInfo(t_integrand, t)
        if substep.contains_dont_know():
            return None
        return URule(integrand, x, t, u_of_t, substep)

    def _attempt_tangent():
        # integrand = K*sin(u)**m*cos(u)**n with even integer m + n:
        # substitute t = tan(u), so that sin(u) = t/sqrt(1 + t**2) and
        # cos(u) = 1/sqrt(1 + t**2) give K*t**m*(1 + t**2)**(-(m + n)/2)
        # with dx = dt/(w*(1 + t**2))
        m_exp = n_exp = S.Zero
        coefficient = S.One
        for factor in Mul.make_args(expr_sc):
            base, power = factor.as_base_exp()
            if base == s:
                m_exp += power
            elif base == c:
                n_exp += power
            elif base == s/c:
                m_exp += power
                n_exp -= power
            elif base == c/s:
                n_exp += power
                m_exp -= power
            elif not factor.has(s, c):
                coefficient *= factor
            else:
                return None
        if (m_exp + n_exp).is_even is not True:
            return None
        if hyperbolic:
            t_squared = 1 - t**2
        else:
            t_squared = 1 + t**2
        t_integrand = coefficient*t**m_exp*t_squared**(-(m_exp + n_exp)/2 - 1)/omega
        substep = yield IntegralInfo(t_integrand, t)
        if substep.contains_dont_know():
            return None
        return URule(integrand, x, t, tan_f(u_func), substep)

    generic_step = yield from _attempt(s, c, cos_squared, sin_f(u_func), S.One)
    if generic_step is None:
        generic_step = yield from _attempt(c, s, sin_squared, cos_f(u_func), cos_branch_sign)
    if generic_step is None:
        generic_step = yield from _attempt_tangent()
    if generic_step is None:
        return None

    if omega.is_zero is False:
        return generic_step

    phase = u_func - omega*x
    zero_integrand = integrand.xreplace(
        {atom: atom.func(phase) for atom in atoms})
    zero_substep = yield IntegralInfo(zero_integrand, x)
    zero_step = RewriteRule(integrand, x, zero_integrand, zero_substep)

    return PiecewiseRule(integrand, x, [(zero_step, Eq(omega, 0)), (generic_step, S.true)])

def trig_half_angle_square_rule(integral):
    """
    Rewrite powers of a + b*sin(u), a + b*cos(u) with b = a or b = -a
    (and the hyperbolic counterparts with cosh) using the half-angle
    squares

        1 + cos(u) = 2*cos(u/2)**2       1 - cos(u) = 2*sin(u/2)**2
        1 + sin(u) = 2*cos(u/2 - pi/4)**2
        1 - sin(u) = 2*sin(u/2 - pi/4)**2
        cosh(u) + 1 = 2*cosh(u/2)**2     cosh(u) - 1 = 2*sinh(u/2)**2

    which turns a fractional power p of the sum into a plain power 2*p
    of a single sine or cosine, e.g.

        sqrt(a + a*cos(x)) = sqrt(2*a)*cos(x/2)**(2/2)

    The remaining trigonometric functions of u are rewritten in the
    half-angle argument as well, so that the substitution rules for
    powers of sines and cosines of a common argument apply.
    """
    integrand, x = integral

    for factor in Mul.make_args(integrand):
        base, p = factor.as_base_exp()
        if not (p.is_Rational and not p.is_Integer):
            continue
        trig_atoms = base.atoms(sin, cos, cosh)
        if len(trig_atoms) != 1:
            continue
        g, = trig_atoms
        u = g.args[0]
        if u.diff(x).has(x):
            continue
        g_var = Dummy("g")
        base_poly = base.subs(g, g_var).as_poly(g_var)
        if base_poly is None or base_poly.degree() != 1:
            continue
        c0, c1 = base_poly.nth(0), base_poly.nth(1)
        if c0.has(x) or c1.has(x):
            continue
        if (c0 - c1).cancel() == 0:
            positive_sign = True
        elif (c0 + c1).cancel() == 0:
            positive_sign = False
        else:
            continue
        # base = c0*(1 +- g(u)) = scale*h(v)**2 with the half angle v,
        # choosing the half angle whose cosine or sine is positive around
        # u = 0, as for the 1 +- cos(u) identities
        scale = 2*c0
        if isinstance(g, cos):
            v = u/2
            h = cos if positive_sign else sin
        elif isinstance(g, sin):
            # 1 + sin(u) = 2*cos(u/2 - pi/4)**2
            # 1 - sin(u) = 2*cos(u/2 + pi/4)**2
            v = u/2 - pi/4 if positive_sign else u/2 + pi/4
            h = cos
        else:  # cosh(u) + 1 = 2*cosh(u/2)**2, cosh(u) - 1 = 2*sinh(u/2)**2
            v = u/2
            h = cosh if positive_sign else sinh
            if not positive_sign:
                # base = c1*(cosh(u) - 1), keeping the square positive
                scale = 2*c1
        new_factor = scale**p*h(v)**(2*p)
        # rewrite the rest of the integrand in the half-angle argument
        rest = integrand/factor
        if isinstance(g, sin):
            u_of_v = 2*v + (pi/2 if positive_sign else -pi/2)
        else:
            u_of_v = 2*v
        replacements = {}
        for f in rest.atoms(TrigonometricFunction, HyperbolicFunction):
            if f.args[0] == u:
                replacements[f] = f.func(u_of_v)
        rest = expand_trig(rest.xreplace(replacements))
        rewritten = new_factor*rest
        substep = yield IntegralInfo(rewritten, x)
        if not substep.contains_dont_know():
            return RewriteRule(integrand, x, rewritten, substep)
    return None

def power_substitution_rule(integral):
    """
    Substitute u = x**n when a power x**n appears inside a function
    argument and the integrand divided by the derivative n*x**(n - 1)
    can be written in terms of u alone (assuming x > 0 for the branch
    choice, as usual for such tables):

        Integral(x**3*sinh(a + b*x**2), x)
        -> Integral(u*sinh(a + b*u)/2, u),  u = x**2

    Likewise substitute u = (c + d*x)**(1/q) for a radical of a linear
    expression inside a function argument, inverting x = (u**q - c)/d:

        Integral(sinh(a + b*sqrt(x)), x)
        -> Integral(2*u*sinh(a + b*u), u),  u = sqrt(x)
    """
    integrand, x = integral

    # substituting under a variable binder would capture the bound
    # variable (e.g. the Lambda variable of a RootSum from ratint)
    if integrand.has(RootSum, Lambda):
        return None
    argument_candidates = []
    power_candidates = []
    for f in integrand.atoms(Function):
        if not f.has(x) or isinstance(f, Piecewise):
            continue
        argument = f.args[0]
        # linear arguments are handled by ordinary substitution
        if not argument.has(x) or not argument.diff(x).has(x):
            continue
        pows = [p for p in argument.atoms(Pow) if p.has(x)]
        if len(pows) != 1:
            continue
        p = pows[0]
        base, e = p.base, p.exp
        if e.has(x):
            continue
        if p.base == x and not (e.is_Rational and not e.is_Integer):
            power_candidates.append(e)
        p_var = Dummy("p")
        arg_poly = argument.subs(p, p_var).as_poly(p_var)
        if arg_poly is None or arg_poly.degree() != 1:
            continue
        c0, c1 = arg_poly.nth(0), arg_poly.nth(1)
        if c0.has(x) or c1.has(x):
            continue
        if base == x:
            b0, b1 = S.Zero, S.One
        else:
            base_poly = base.as_poly(x)
            if base_poly is None or base_poly.degree() != 1:
                continue
            b0, b1 = base_poly.nth(0), base_poly.nth(1)
        argument_candidates.append((argument, c0, c1, base, e, b0, b1))

    def _improves(transformed, var):
        # guard against substitution ping-pong: the transformed integrand
        # must not contain radicals of var, inside or outside function
        # arguments, which a later substitution would just undo
        if any(pw.exp.is_Rational and not pw.exp.is_Integer
               for pw in transformed.atoms(Pow) if pw.base.has(var)):
            return False
        return all(argument.is_polynomial(var)
                   for f in transformed.atoms(Function)
                   if f.has(var) and not isinstance(f, Piecewise)
                   for argument in f.args)

    u = Dummy("u")
    # substitute the whole argument u = c0 + c1*(b0 + b1*x)**e, so that
    # the substituted function comes out as a plain f(u)
    for argument, c0, c1, base, e, b0, b1 in ordered(argument_candidates):
        p_of_u = (u - c0)/c1
        x_of_u = (p_of_u**(1/e) - b0)/b1
        du_dx = c1*e*base**(e - 1)*b1
        transformed = (integrand/du_dx).cancel()
        # replace powers of the base whose exponent is a multiple of e
        # exactly, before falling back to the branch-choosing x_of_u
        replacements = {}
        for pw in transformed.atoms(Pow):
            if pw.base == base and not pw.exp.has(x):
                ratio = (pw.exp/e).cancel()
                if ratio.is_Integer:
                    replacements[pw] = p_of_u**ratio
        transformed = powsimp(transformed.xreplace(replacements).subs(x, x_of_u))
        if transformed.has(x) or not _improves(transformed, u):
            continue
        substep = yield IntegralInfo(transformed, u)
        if not substep.contains_dont_know():
            return URule(integrand, x, u, argument, substep)
    # substitute just the power u = x**n (assuming x > 0 for the branch
    # choice, as usual for such tables)
    u_pos = Dummy("u", positive=True)
    for n in ordered(power_candidates):
        transformed = powsimp((integrand/(n*x**(n - 1))).subs(x, u_pos**(1/n)))
        if transformed.has(x) or not _improves(transformed, u_pos):
            continue
        substep = yield IntegralInfo(transformed, u_pos)
        if not substep.contains_dont_know():
            return URule(integrand, x, u_pos, x**n, substep)
    return None

def exp_power_rewrite_rule(integral):
    """
    Rewrite non-integer powers of exponentials as exponentials,
    exp(u)**p = exp(p*u), e.g. sqrt(exp(a + b*x)) = exp(a/2 + b*x/2),
    so that the exponential rules apply; for example
    x*sqrt(exp(a + b*x)) then integrates by parts.
    """
    integrand, x = integral
    replacements = {}
    for pw in integrand.atoms(Pow):
        if (isinstance(pw.base, exp) and pw.base.has(x)
                and not pw.exp.has(x)):
            replacements[pw] = exp((pw.exp*pw.base.args[0]).expand())
    if not replacements:
        return None
    rewritten = integrand.xreplace(replacements)
    substep = yield IntegralInfo(rewritten, x)
    if substep.contains_dont_know():
        return None
    return RewriteRule(integrand, x, rewritten, substep)

class DilogRule(AtomicRule):
    """integrate(log(a + b*x)/(c + d*x), x), non-proportional arguments:

        (log((a*d - b*c)/d)*log(c + d*x)
         - polylog(2, -b*(c + d*x)/(a*d - b*c)))/d

    the classic dilogarithm antiderivative (valid up to the usual branch
    choices of log and polylog).
    """

    __slots__ = ("a", "b", "c", "d")

    a: Expr
    b: Expr
    c: Expr
    d: Expr

    def __init__(self, integrand: Expr, variable: Symbol,
                 a: Expr, b: Expr, c: Expr, d: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def eval(self) -> Expr:
        a, b, c, d = self.a, self.b, self.c, self.d
        x = self.variable
        delta = a*d - b*c
        return (log(delta/d)*log(c + d*x)
                - polylog(2, -b*(c + d*x)/delta))/d


def dilog_rule(integral):
    """Match K*log(a + b*x)/(c + d*x) for DilogRule."""
    integrand, x = integral
    logs = [f for f in integrand.atoms(log) if f.has(x)]
    if len(logs) != 1:
        return None
    L = logs[0]
    w_poly = L.args[0].as_poly(x)
    if w_poly is None or w_poly.degree() != 1:
        return None
    rest = (integrand/L).cancel()
    if rest.has(log):
        return None
    numer, denom = rest.as_numer_denom()
    if numer.has(x):
        return None
    denom_poly = denom.as_poly(x)
    if denom_poly is None or denom_poly.degree() != 1:
        return None
    a, b = w_poly.nth(0), w_poly.nth(1)
    c, d = denom_poly.nth(0), denom_poly.nth(1)
    if (a*d - b*c).cancel() == 0:
        # proportional arguments reduce by an ordinary substitution
        return None
    core = L/(c + d*x)
    step = DilogRule(core, x, a, b, c, d)
    if numer == 1:
        return step
    return ConstantTimesRule(integrand, x, numer, core, step)


def log_times_rational_rule(integral):
    """
    Integrate log(w)*R(x) for a linear w and rational R by splitting R
    into partial fractions: the polynomial and repeated-pole parts
    integrate by parts and the simple-pole parts are dilogarithms.
    """
    integrand, x = integral
    logs = [f for f in integrand.atoms(log) if f.has(x)]
    if len(logs) != 1:
        return None
    L = logs[0]
    R = (integrand/L).cancel()
    if R.has(log) or not R.is_rational_function(x):
        return None
    numer, denom = R.as_numer_denom()
    if denom.as_poly(x) is None or degree(denom, x) < 1:
        return None
    try:
        decomposed = R.apart(x)
    except (PolynomialError, NotImplementedError):
        return None
    if decomposed == R or not isinstance(decomposed, Add):
        return None
    rewritten = Add(*[term*L for term in decomposed.args])
    substep = yield IntegralInfo(rewritten, x)
    if substep.contains_dont_know():
        return None
    return RewriteRule(integrand, x, rewritten, substep)


class ExpOverLinearRule(AtomicRule):
    """integrate(exp(a + b*x)/(c + d*x), x), the shifted exponential
    integral:

        exp(a - b*c/d)*Ei(b*(c + d*x)/d)/d
    """

    __slots__ = ("a", "b", "c", "d")

    a: Expr
    b: Expr
    c: Expr
    d: Expr

    def __init__(self, integrand: Expr, variable: Symbol,
                 a: Expr, b: Expr, c: Expr, d: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def eval(self) -> Expr:
        a, b, c, d = self.a, self.b, self.c, self.d
        x = self.variable
        return exp(a - b*c/d)*Ei(b*(c + d*x)/d)/d


def exp_over_linear_rule(integral):
    """
    Integrate K*exp(a + b*x)/(c + d*x)**n for a positive integer n: the
    n = 1 case is the shifted exponential integral of
    :class:`ExpOverLinearRule` and higher poles reduce with

        exp(w)/q**n = (-exp(w)/(d*(n - 1)*q**(n - 1)))'
                      + b/(d*(n - 1)) * exp(w)/q**(n - 1),  q = c + d*x
    """
    integrand, x = integral

    exps = [f for f in integrand.atoms(exp) if f.has(x)]
    if len(exps) != 1:
        return None
    E = exps[0]
    w_poly = E.args[0].as_poly(x)
    if w_poly is None or w_poly.degree() != 1:
        return None
    rest = powsimp(integrand/E, combine='exp')
    if rest.has(exp):
        return None
    numer, denom = rest.as_numer_denom()
    if numer.has(x):
        return None
    base, n = denom.as_base_exp()
    if base.as_poly(x) is not None and base.as_poly(x).degree() > 1:
        # recover a linear power from an expanded denominator
        constant, factored = factor_terms(denom.factor()).as_independent(x)
        numer, denom = numer/constant, factored
        base, n = denom.as_base_exp()
    base_poly = base.as_poly(x)
    if (base_poly is None or base_poly.degree() != 1
            or not (n.is_Integer and n >= 1)):
        return None
    a0, b0 = w_poly.nth(0), w_poly.nth(1)
    c0, d0 = base_poly.nth(0), base_poly.nth(1)

    core = E/base**n
    if n == 1:
        step: Rule = ExpOverLinearRule(core, x, a0, b0, c0, d0)
    else:
        boundary = -E/(d0*(n - 1)*base**(n - 1))
        coefficient = b0/(d0*(n - 1))
        remainder = E/base**(n - 1)
        substep = yield from exp_over_linear_rule(
            IntegralInfo(remainder, x))
        if substep is None or substep.contains_dont_know():
            return None
        derivative = Derivative(boundary, x)
        rewritten = derivative + coefficient*remainder
        add_step = AddRule(rewritten, x, [
            DerivativeRule(derivative, x),
            ConstantTimesRule(coefficient*remainder, x, coefficient,
                              remainder, substep)])
        step = RewriteRule(core, x, rewritten, add_step)
    if numer != 1:
        step = ConstantTimesRule(integrand, x, numer, core, step)
    return step


def exp_times_rational_rule(integral):
    """
    Integrate exp(w)*R(x) for a linear w and rational R by splitting R
    into partial fractions: the polynomial part integrates by tabular
    parts and the pole terms are (shifted) exponential integrals.
    """
    integrand, x = integral

    exps = [f for f in integrand.atoms(exp) if f.has(x)]
    if len(exps) != 1:
        return None
    E = exps[0]
    w_poly = E.args[0].as_poly(x)
    if w_poly is None or w_poly.degree() != 1:
        return None
    R = powsimp(integrand/E, combine='exp').cancel()
    if R.has(exp) or not R.is_rational_function(x):
        return None
    numer, denom = R.as_numer_denom()
    if denom.as_poly(x) is None or degree(denom, x) < 1:
        return None
    try:
        decomposed = R.apart(x)
    except (PolynomialError, NotImplementedError):
        return None
    if decomposed == R or not isinstance(decomposed, Add):
        return None
    rewritten = Add(*[term*E for term in decomposed.args])
    substep = yield IntegralInfo(rewritten, x)
    if substep.contains_dont_know():
        return None
    return RewriteRule(integrand, x, rewritten, substep)


class LogSinRule(AtomicRule):
    """integrate(log(sin(a + b*x)), x) in terms of the dilogarithm:

        (I*u**2/2 - u*log(2*I) - I*polylog(2, exp(-2*I*u))/2)/b,  u = a + b*x

    (valid up to the usual branch choices).  With a shifted phase this
    also covers log(cos(a + b*x)) = log(sin(a + pi/2 + b*x)).
    """

    __slots__ = ("a", "b")

    a: Expr
    b: Expr

    def __init__(self, integrand: Expr, variable: Symbol,
                 a: Expr, b: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b

    def eval(self) -> Expr:
        u = self.a + self.b*self.variable
        return (S.ImaginaryUnit*u**2/2 - u*log(2*S.ImaginaryUnit)
                - S.ImaginaryUnit*polylog(2, exp(-2*S.ImaginaryUnit*u))/2)/self.b


class LogSinhRule(AtomicRule):
    """integrate(log(sinh(a + b*x)), x) =
    (u**2/2 - u*log(2) + polylog(2, exp(-2*u))/2)/b with u = a + b*x."""

    __slots__ = ("a", "b")

    a: Expr
    b: Expr

    def __init__(self, integrand: Expr, variable: Symbol,
                 a: Expr, b: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b

    def eval(self) -> Expr:
        u = self.a + self.b*self.variable
        return (u**2/2 - u*log(2) + polylog(2, exp(-2*u))/2)/self.b


class LogCoshRule(AtomicRule):
    """integrate(log(cosh(a + b*x)), x) =
    (u**2/2 - u*log(2) + polylog(2, -exp(-2*u))/2)/b with u = a + b*x."""

    __slots__ = ("a", "b")

    a: Expr
    b: Expr

    def __init__(self, integrand: Expr, variable: Symbol,
                 a: Expr, b: Expr) -> None:
        super().__init__(integrand, variable)
        self.a = a
        self.b = b

    def eval(self) -> Expr:
        u = self.a + self.b*self.variable
        return (u**2/2 - u*log(2) + polylog(2, -exp(-2*u))/2)/self.b


def log_of_trig_rule(integral):
    """
    Integrate the logarithm of a product of trigonometric or hyperbolic
    functions of a common linear argument: log(sin(u)), log(cos(u)),
    log(sinh(u)) and log(cosh(u)) close in terms of dilogarithms
    (:class:`LogSinRule` and friends), and a general
    log(k*g1(u)**n1*g2(u)**n2*...) is split into a sum of such
    logarithms first (with the usual branch conventions), e.g.

        log(a*tan(x)**n) = log(a) + n*log(sin(x)) - n*log(cos(x))
    """
    integrand, x = integral

    if not isinstance(integrand, log):
        return None
    coefficient = S.One
    factors = []
    u_common = None
    for factor in Mul.make_args(integrand.args[0]):
        if not factor.has(x):
            coefficient *= factor
            continue
        base, power = factor.as_base_exp()
        if not isinstance(base, (sin, cos, tan, cot, sec, csc,
                                 sinh, cosh, tanh, coth, sech, csch)):
            return None
        if power.has(x):
            return None
        u = base.args[0]
        if u.diff(x).has(x):
            return None
        if u_common is None:
            u_common = u
        elif u != u_common:
            return None
        factors.append((type(base), u, power))
    if not factors:
        return None

    if coefficient == 1 and len(factors) == 1 and factors[0][2] == 1:
        fcls, u, _ = factors[0]
        b = u.diff(x)
        a = u - b*x
        if fcls is sin:
            return LogSinRule(integrand, x, a, b)
        if fcls is cos:
            return LogSinRule(integrand, x, a + pi/2, b)
        if fcls is sinh:
            return LogSinhRule(integrand, x, a, b)
        if fcls is cosh:
            return LogCoshRule(integrand, x, a, b)

    quotient_split = {
        sin: ((sin, 1),), cos: ((cos, 1),),
        tan: ((sin, 1), (cos, -1)), cot: ((cos, 1), (sin, -1)),
        sec: ((cos, -1),), csc: ((sin, -1),),
        sinh: ((sinh, 1),), cosh: ((cosh, 1),),
        tanh: ((sinh, 1), (cosh, -1)), coth: ((cosh, 1), (sinh, -1)),
        sech: ((cosh, -1),), csch: ((sinh, -1),),
    }
    parts = [] if coefficient == 1 else [log(coefficient)]
    for fcls, u, n in factors:
        for g, sign in quotient_split[fcls]:
            parts.append(sign*n*log(g(u)))
    rewritten = Add(*parts)
    if rewritten == integrand:
        return None
    substep = yield IntegralInfo(rewritten, x)
    if substep.contains_dont_know():
        return None
    return RewriteRule(integrand, x, rewritten, substep)


class PiecewiseContinuousRule(Rule):
    """integrate f(x, sign(w1), Abs(w2), floor(w3), ...) for linear wi
    by replacing each piecewise-continuous factor with a temporary
    constant, integrating, and then restoring the constants one at a
    time while removing the jumps at the breakpoints: for a signum

        G_new = G(x, sign(w), ...) - (J/2)*sign(w),
        J = (G(b, 1, ...) - G(b, -1, ...)) at the breakpoint w(b) = 0,

    and for a floor, whose breakpoints are at w(x) = m for integer m,

        G_new = G(x, floor(w), ...) - Sum(J(m), (m, 1, floor(w))),
        J(m) = (G(x_m, m, ...) - G(x_m, m - 1, ...)) at w(x_m) = m,

    with the sum evaluated in closed form when possible, so that the
    result is continuous on the domain of maximum extent (Jeffrey &
    Rich, ISSAC 1998).
    """

    __slots__ = ("sign_replacements", "floor_replacements", "substep")

    # lists of (temporary symbol, argument w) pairs
    sign_replacements: list[tuple[Symbol, Expr]]
    floor_replacements: list[tuple[Symbol, Expr]]
    substep: Rule

    def __init__(self, integrand: Expr, variable: Symbol,
                 sign_replacements: list[tuple[Symbol, Expr]],
                 floor_replacements: list[tuple[Symbol, Expr]],
                 substep: Rule) -> None:
        super().__init__(integrand, variable)
        self.sign_replacements = sign_replacements
        self.floor_replacements = floor_replacements
        self.substep = substep

    def eval(self) -> Expr:
        from sympy.functions.elementary.complexes import sign
        x = self.variable
        G = self.substep.eval()
        for s, w in self.sign_replacements:
            w_poly = w.as_poly(x)
            breakpoint = -w_poly.nth(0)/w_poly.nth(1)
            jump = (G.subs(s, 1) - G.subs(s, -1)).subs(x, breakpoint)
            G = G.subs(s, sign(w)) - jump/2*sign(w)
        for f, w in self.floor_replacements:
            w_poly = w.as_poly(x)
            c0, c1 = w_poly.nth(0), w_poly.nth(1)
            m = Dummy("m", integer=True)
            x_m = (m - c0)/c1
            jump = (G.subs(f, m) - G.subs(f, m - 1)).subs(x, x_m)
            # signum factors proportional to the summation variable have
            # the sign of w throughout the effective summation range
            # (m = 1..floor(w) for w > 0 and, by the usual summation
            # convention, m = floor(w)+1..0 for w < 0), which restores
            # the S(x) factors of Jeffrey & Rich's results
            sign_fixes = {}
            for atom in jump.atoms(sign):
                g_poly = atom.args[0].as_poly(m)
                if (g_poly is not None and g_poly.degree() == 1
                        and g_poly.nth(0).is_zero
                        and g_poly.nth(1).is_positive):
                    sign_fixes[atom] = sign(w)
            if sign_fixes:
                jump = jump.xreplace(sign_fixes)
            total = Sum(jump, (m, 1, floor(w))).doit()
            G = G.subs(f, floor(w)) - total
        return G

    def contains_dont_know(self) -> bool:
        return self.substep.contains_dont_know()


def piecewise_continuous_rule(integral):
    """
    Match integrands containing Abs, sign, Heaviside, floor, ceiling,
    frac or Mod of linear arguments for
    :class:`PiecewiseContinuousRule`, using the preliminary
    transformations of Jeffrey & Rich (ISSAC 1998):

        Abs(w) -> w*sign(w)           Heaviside(w) -> (1 + sign(w))/2
        ceiling(w) -> -floor(-w)      frac(w) -> w - floor(w)
        Mod(w, y) -> w - y*floor(w/y)
    """
    from sympy.functions.elementary.complexes import sign
    integrand, x = integral

    transformed = integrand
    atoms = transformed.atoms(Abs, sign, Heaviside, floor, ceiling,
                              frac, Mod)
    if not any(a.has(x) for a in atoms):
        return None
    # preliminary transformations to sign and floor
    for atom in transformed.atoms(Mod):
        w, y = atom.args
        if w.has(x) and not y.has(x):
            transformed = transformed.xreplace({atom: w - y*floor(w/y)})
    for atom in transformed.atoms(frac):
        transformed = transformed.xreplace(
            {atom: atom.args[0] - floor(atom.args[0])})
    for atom in transformed.atoms(ceiling):
        transformed = transformed.xreplace({atom: -floor(-atom.args[0])})
    for atom in transformed.atoms(Abs, Heaviside):
        w = atom.args[0]
        if not w.has(x):
            continue
        if isinstance(atom, Abs):
            transformed = transformed.xreplace({atom: w*sign(w)})
        else:
            transformed = transformed.xreplace({atom: (1 + sign(w))/2})

    sign_replacements = []
    for atom in transformed.atoms(sign):
        w = atom.args[0]
        if not w.has(x):
            continue
        w_poly = w.as_poly(x)
        if w_poly is None or w_poly.degree() != 1:
            return None
        s = Dummy("s")
        transformed = transformed.xreplace({atom: s})
        sign_replacements.append((s, w))
    floor_replacements = []
    for atom in transformed.atoms(floor):
        w = atom.args[0]
        if not w.has(x):
            continue
        w_poly = w.as_poly(x)
        if w_poly is None or w_poly.degree() != 1:
            return None
        f = Dummy("f")
        transformed = transformed.xreplace({atom: f})
        floor_replacements.append((f, w))
    if (transformed.has(Abs, sign, Heaviside, floor, ceiling, frac, Mod)
            or not (sign_replacements or floor_replacements)):
        return None
    substep = yield IntegralInfo(transformed.expand(), x)
    if substep.contains_dont_know():
        return None
    return PiecewiseContinuousRule(integrand, x, sign_replacements,
                                   floor_replacements, substep)


class LogPolylogRule(AtomicRule):
    """integrate(log(c*x**r)**m*polylog(k, a*x**n)/x, x) by the
    integration-by-parts ladder that raises the polylog order while
    lowering the logarithm power,

        Sum((-1)**j*(m!/(m - j)!)*r**j*log(c*x**r)**(m - j)
            *polylog(k + 1 + j, a*x**n)/n**(j + 1), (j, 0, m))
    """

    __slots__ = ("w", "m", "r", "k", "a", "n")

    w: Expr  # c*x**r, the argument of the logarithm
    m: int
    r: Expr
    k: Expr
    a: Expr
    n: Expr

    def __init__(self, integrand: Expr, variable: Symbol, w: Expr,
                 m: int, r: Expr, k: Expr, a: Expr, n: Expr) -> None:
        super().__init__(integrand, variable)
        self.w = w
        self.m = m
        self.r = r
        self.k = k
        self.a = a
        self.n = n

    def eval(self) -> Expr:
        from sympy import factorial
        x = self.variable
        L = log(self.w)
        return Add(*[(S.NegativeOne**j*factorial(self.m)/factorial(self.m - j)
                      *self.r**j*L**(self.m - j)
                      *polylog(self.k + 1 + j, self.a*x**self.n)
                      /self.n**(j + 1))
                     for j in range(self.m + 1)])


class LogLogBinomialRule(AtomicRule):
    """integrate(log(c*x**r)**m*log(a + b*x**n)/x, x) for a != 0, using
    log(a + b*x**n) = log(a) + log(1 + b*x**n/a) and the polylogarithm
    ladder of :class:`LogPolylogRule` for the second part (with
    log(1 + z) = -polylog(1, -z), so the result contains polylogarithms
    of order 2 and higher only).
    """

    __slots__ = ("w", "m", "r", "a", "b", "n")

    w: Expr
    m: int
    r: Expr
    a: Expr
    b: Expr
    n: Expr

    def __init__(self, integrand: Expr, variable: Symbol, w: Expr,
                 m: int, r: Expr, a: Expr, b: Expr, n: Expr) -> None:
        super().__init__(integrand, variable)
        self.w = w
        self.m = m
        self.r = r
        self.a = a
        self.b = b
        self.n = n

    def eval(self) -> Expr:
        from sympy import factorial
        x = self.variable
        L = log(self.w)
        m, r, n = self.m, self.r, self.n
        constant_part = log(self.a)*L**(m + 1)/((m + 1)*r)
        ladder = Add(*[(S.NegativeOne**j*factorial(m)/factorial(m - j)
                        *r**j*L**(m - j)
                        *polylog(2 + j, -self.b*x**n/self.a)
                        /n**(j + 1))
                       for j in range(m + 1)])
        return constant_part - ladder


def polylog_ladder_rule(integral):
    """
    Integrate products of logarithm powers with polylogarithms or
    logarithms of binomials, closing in higher polylogarithms:

    - K*log(c*x**r)**m*polylog(k, a*x**n)/x  (LogPolylogRule)
    - K*log(c*x**r)**m*log(a + b*x**n)/x with a != 0, splitting
      log(a + b*x**n) = log(a) - polylog(1, -b*x**n/a)
    - K*log(c*x**r)**m/(c2 + d2*x) with m >= 2, by parts, which reduces
      to the previous case
    """
    integrand, x = integral

    coefficient = S.One
    log_power = None      # (w, m, r) with w = c*x**r
    polylog_part = None   # (k, a, n)
    log_binomial = None   # (a0, b0, n)
    denominator = None    # (c2, d2) for 1/(c2 + d2*x)
    inv_x = False
    for factor in Mul.make_args(integrand):
        if not factor.has(x):
            coefficient *= factor
            continue
        if factor == 1/x:
            inv_x = True
            continue
        base, power = factor.as_base_exp()
        if isinstance(base, log) and power.is_Integer and power >= 1:
            argument = base.args[0]
            c_part, x_part = argument.as_independent(x)
            x_base, x_exp = x_part.as_base_exp()
            if (c_part*x_part == argument and x_base == x
                    and not x_exp.has(x) and log_power is None):
                # a logarithm of c*x**r
                log_power = (argument, int(power), x_exp)
                continue
            # log of a binomial a0 + b0*x**n
            terms = Add.make_args(argument)
            if len(terms) == 2 and log_binomial is None and power == 1:
                consts = [t for t in terms if not t.has(x)]
                xterms = [t for t in terms if t.has(x)]
                if len(consts) == 1 and len(xterms) == 1:
                    b0, x_part = xterms[0].as_independent(x)
                    x_base, x_exp = x_part.as_base_exp()
                    if x_base == x and not x_exp.has(x):
                        log_binomial = (consts[0], b0, x_exp)
                        continue
            return None
        if isinstance(base, polylog) and power == 1 and polylog_part is None:
            k, argument = base.args
            a0, x_part = argument.as_independent(x)
            x_base, x_exp = x_part.as_base_exp()
            if x_base == x and not x_exp.has(x):
                polylog_part = (k, a0, x_exp)
                continue
            return None
        if power == -1 and denominator is None:
            base_poly = base.as_poly(x)
            if base_poly is not None and base_poly.degree() == 1:
                denominator = (base_poly.nth(0), base_poly.nth(1))
                continue
        return None

    if log_power is None:
        return None
    w, m, r = log_power

    if polylog_part is not None and inv_x and log_binomial is None \
            and denominator is None:
        k, a0, n = polylog_part
        core = log(w)**m*polylog(k, a0*x**n)/x
        step: Rule = LogPolylogRule(core, x, w, m, r, k, a0, n)
        if coefficient != 1:
            step = ConstantTimesRule(integrand, x, coefficient, core, step)
        return step

    if log_binomial is not None and inv_x and polylog_part is None \
            and denominator is None:
        a0, b0, n = log_binomial
        if a0.is_zero:
            return None
        core = log(w)**m*log(a0 + b0*x**n)/x
        step = LogLogBinomialRule(core, x, w, m, r, a0, b0, n)
        if coefficient != 1:
            step = ConstantTimesRule(integrand, x, coefficient, core, step)
        return step

    if denominator is not None and not inv_x and polylog_part is None \
            and log_binomial is None and m >= 2:
        c2, d2 = denominator
        u = log(w)**m
        dv = 1/(c2 + d2*x)
        v_step = yield IntegralInfo(dv, x)
        if v_step.contains_dont_know():
            return None
        V = v_step.eval()
        remaining = (coefficient*V*u.diff(x)).expand()
        second_step = yield IntegralInfo(remaining, x)
        if second_step.contains_dont_know():
            return None
        return PartsRule(integrand, x, coefficient*u, dv, v_step,
                         second_step)

    return None


def log_inverse_parts_rule(integral):
    """
    Integrate R(x)*L1(x)*L2(x), for a rational function R and two factors
    that are powers of logarithms or of inverse trigonometric or
    hyperbolic functions, by parts with u = L1*L2 and dv = R dx:

        Integral(R*L1*L2, x) = V*L1*L2 - Integral(V*(L1'*L2 + L1*L2'), x)

    Differentiating a factor lowers its transcendence, so the remaining
    integral has a single such factor and falls to the other rules; for
    example (d + e*x**2)*log(c*x**n)*atanh(a*x) reduces to integrals of a
    rational function times atanh(a*x) or times log(c*x**n).
    """
    integrand, x = integral

    L_CLASSES = (log, asin, acos, atan, acot, asec, acsc,
                 asinh, acosh, atanh, acoth, asech, acsch)
    l_factors = []
    rational_factors = []
    for factor in Mul.make_args(integrand):
        base, power = factor.as_base_exp()
        if (isinstance(base, L_CLASSES) and base.has(x)
                and power.is_Integer and power > 0):
            l_factors.append(factor)
        else:
            rational_factors.append(factor)
    if len(l_factors) < 2:
        return None
    R = Mul(*rational_factors)
    if not R.is_rational_function(x):
        return None
    u = Mul(*l_factors)
    v_step = yield IntegralInfo(R, x)
    if v_step.contains_dont_know():
        return None
    V = v_step.eval()
    remaining = (V*u.diff(x)).expand()
    second_step = yield IntegralInfo(remaining, x)
    if second_step.contains_dont_know():
        return None
    return PartsRule(integrand, x, u, R, v_step, second_step)

def linear_power_product_rule(integral):
    """
    Integrate K*(a + b*x)**m*(c + d*x)**n for non-integer rational
    exponents whose sum m + n is an integer, by rewriting the integrand as

        K*((a + b*x)/(c + d*x))**m*(c + d*x)**(m + n)

    a rational function of x and of the fractional linear power
    ((a + b*x)/(c + d*x))**(1/q), which sqrt_fractional_linear_rule
    integrates by substituting t = ((a + b*x)/(c + d*x))**(1/q).

    For example 1/(sqrt(x + 1)*sqrt(x + 2)) becomes
    sqrt((x + 1)/(x + 2))/(x + 1) and t = sqrt((x + 1)/(x + 2)) leads to
    a rational function of t.
    """
    integrand, x = integral

    coefficient = S.One
    linear_powers = []
    for factor in Mul.make_args(integrand):
        if not factor.has(x):
            coefficient *= factor
            continue
        base, power = factor.as_base_exp()
        if not (power.is_Rational and not power.is_Integer):
            return None
        base_poly = base.as_poly(x)
        if base_poly is None or base_poly.degree() != 1:
            return None
        linear_powers.append((base, power))
    if len(linear_powers) != 2:
        return None
    (base1, m), (base2, n) = linear_powers
    if not (m + n).is_Integer:
        return None
    # exclude proportional bases such as (2*x + 2)**m*(x + 1)**n, which
    # reduce to a single power
    if not (base1/base2).cancel().has(x):
        return None

    rewritten = coefficient*(base1/base2)**m*base2**(m + n)
    substep = yield from sqrt_fractional_linear_rule(IntegralInfo(rewritten, x))
    if substep is None or substep.contains_dont_know():
        return None
    return RewriteRule(integrand, x, rewritten, substep)

def inverse_function_substitution_rule(integral):
    """
    Substitute u = ainv(w) for an inverse trigonometric or hyperbolic
    function ainv of a linear argument w, inverting to x and integrating
    in u; for example with u = asin(a*x)

        Integral(x**2*asin(a*x)**4, x)
        -> Integral(sin(u)**2*u**4*cos(u)/a**3, u)

    This succeeds whenever the resulting trigonometric integral is
    doable, e.g. for a polynomial multiplying a positive integer power of
    the inverse function.

    Since u lies in the range of the inverse function, factors such as
    sqrt(cos(u)**2) arising from sqrt(1 - a**2*x**2) simplify to cos(u)
    on that principal branch.
    """
    integrand, x = integral

    INVERSES = {asin: sin, acos: cos, atan: tan, acot: cot, asec: sec,
                acsc: csc, asinh: sinh, acosh: cosh, atanh: tanh,
                acoth: coth, asech: sech, acsch: csch}
    if integrand.has(RootSum, Lambda):
        # substituting x under a variable binder would capture the bound
        # variable
        return None
    inv_atoms = integrand.atoms(*INVERSES)
    if len(inv_atoms) != 1:
        return None
    inv_atom, = inv_atoms
    w = inv_atom.args[0]
    dw = w.diff(x)
    # require a linear argument w = dw*x + phase
    if dw.is_zero or dw.has(x):
        return None
    fn = INVERSES[type(inv_atom)]
    u = Dummy("u")
    masked = integrand.xreplace({inv_atom: u})
    # only attempt when the rest of the integrand is algebraic in x:
    # remaining transcendental factors such as log(x) rarely survive the
    # transformation, and exploring them can be very slow
    if any(f.has(x) for f in masked.atoms(Function)):
        return None
    x_of_u = ((fn(u) - (w - dw*x))/dw).expand()
    replaced = masked.subs(x, x_of_u)
    replaced = (replaced*manual_diff(x_of_u, u)).trigsimp()
    if replaced.has(x):
        return None
    # simplify sqrt(trig(u)**2) -> trig(u): valid on the principal branch
    # (e.g. cos(u) >= 0 for u = asin(w), sinh(u) >= 0 for u = acosh(w))
    substitutions = {}
    for f in (sin, cos, tan, sec, csc, cot, sinh, cosh, tanh, sech, csch, coth):
        substitutions[sqrt(f(u)**2)] = f(u)
        substitutions[sqrt(f(u)**(-2))] = 1/f(u)
    replaced = manual_subs(replaced, substitutions)
    substep = yield IntegralInfo(replaced, u)
    if substep.contains_dont_know():
        return None
    return URule(integrand, x, u, inv_atom, substep)

def sqrt_fractional_linear_rule(integral : IntegralInfo):
    """
    Substitute common ((a*x + b)/(c*x + d))**(1/n)
    """
    integrand, x = integral
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    c = Wild('c', exclude=[x])
    d = Wild('d', exclude=[x])
    base0 = None
    powers, exps, ratios = [], [], []
    constant_bases_subs = {}
    # use ordered() to ensure a selection of the smallest base0 (eg. first sqrt(x), then cbrt(2x), x chosen)
    for pow_ in ordered(integrand.find((Pow))): # collect all ((a*x + b)/(c*x + d))**(p/q)
        base, exp_ = pow_.base, pow_.exp
        if exp_.is_Integer or x not in base.free_symbols: # skip 1/x and sqrt(2)
            continue
        if not exp_.is_Rational: # exclude x**pi
            return None
        num, den = base.as_numer_denom()
        match_num = num.match(a*x + b)
        match_den = den.match(c*x + d)
        if not match_num or not match_den:
            continue
        aa, bb = match_num[a], match_num[b]
        cc, dd = match_den[c], match_den[d]
        if cc.is_zero and dd.is_zero:
            return None
        det = aa*dd - bb*cc
        if det.is_zero: # constant value as sqrt((5*x + 10)/(2*x +  4))
            const_val = (S(aa) / cc) if not cc.is_zero else (S(bb) / dd)
            constant_bases_subs[base] = const_val
            continue
        if base0 is None:
            base0 = base
            a0, b0, c0, d0 = aa, bb, cc, dd
            powers.append(pow_)
            exps.append(exp_)
            ratios.append(S.One)
        else:
            power_ratio = powsimp(pow_ / Pow(base0, exp_), force=False).cancel()
            if power_ratio.has(x):
                return None
            powers.append(pow_)
            exps.append(exp_)
            ratios.append(power_ratio)
    if base0 is None and not constant_bases_subs:
        return None
    if constant_bases_subs:
        integrand = integrand.xreplace(constant_bases_subs)
    if base0 is None:
        substep = yield IntegralInfo(integrand, x)
        if not substep.contains_dont_know():
            debug("Integral: {} is rewritten with {} on symbol: {}".format(integral.integrand, integrand, x))
            return RewriteRule(integral.integrand, x, integrand, substep)
        return None
    q0: Integer = lcm_list([exp_i.q for exp_i in exps])
    u = Dummy("u")
    u_x = base0**(S.One/q0)
    u_pow = u**q0
    x_u = (b0 - d0*u_pow)/(c0*u_pow - a0)
    dx_u = (q0*(a0*d0 - b0*c0)*u**(q0 - 1))/(c0*u_pow - a0)**2
    subs_dict = {pow_i: ratio_i * u**(q0*exp_i) for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
    substituted = integrand.xreplace(subs_dict).xreplace({x: x_u}) * dx_u
    substep = yield IntegralInfo(substituted, u)
    if not substep.contains_dont_know():
        pieces: list[tuple[Rule, Boolean]] = []
        det = a0*d0 - b0*c0
        _, base0_denom = base0.as_numer_denom()
        # skips bases where constant value (degenerate case) is not possible (det != 0 or det = 0 implies den = 0)
        # (eg. (3*x + 2)/(4*x + 3), (3x + b)/(d), (4*x + 3)/(c*x + c))
        if not (_if_zero_implies_zero(det, base0_denom)):
            d0_implies_c0 = _if_zero_implies_zero(d0, c0)
            c0_implies_d0 = _if_zero_implies_zero(c0, d0)
            # skips constant value a/c if d != 0 or d = 0 implies c = 0 (eg. (3*x + 2)/(c*d*x + d), (3*x + 2)/(c*x + 4))
            # takes a/c if they both imply each other (eg. (a*x + b)/3*x + 4)) (taking b/d would be the same)
            if not d0_implies_c0 or (c0_implies_d0 and d0_implies_c0):
                const_val = a0 / c0
                subs_a = {pow_i: ratio_i * Pow(const_val, exp_i) for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
                simplified_a = integrand.xreplace(subs_a)
                degenerate_step_a = yield IntegralInfo(simplified_a, x)
                pieces.append((degenerate_step_a, (And(Eq(det, 0), Ne(c0, 0)))))
            if not c0_implies_d0:
                const_val = b0 / d0
                subs_b = {pow_i: ratio_i * Pow(const_val, exp_i) for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
                simplified_b = integrand.xreplace(subs_b)
                simplified_b = simplified_b.subs({a0: 0, c0: 0}) # if det = 0, c = 0 and d != 0, a must be 0
                degenerate_step_b = yield IntegralInfo(simplified_b, x)
                pieces.append((degenerate_step_b, (And(Eq(det, 0), Eq(c0, 0)))))
        step: Rule = URule(integrand, x, u, u_x, substep)
        if pieces:
            pieces.append((step, S.true))
            step = PiecewiseRule(integrand, x, pieces)
        if constant_bases_subs:
            debug("Integral: {} is rewritten with {} on symbol: {}".format(integral.integrand, integrand, x))
            return RewriteRule(integral.integrand, x, integrand, step)
        return step
    return None

def euler_substitution_rule(integral : IntegralInfo):
    """
    Substitute common sqrt(a + b*x + c*x**2) terms using Euler substitution.
    """
    integrand, x = integral
    base0 = None
    powers, exps, ratios = [], [], []
    # use ordered() to ensure a selection of the smallest base0 (eg. first sqrt(x**2 + 1), then sqrt(2*x**2 + 2), x**2 + 1 chosen)
    for pow_ in ordered(integrand.find(Pow)): # collect all (a + b*x + c*x**2)**(p/2)
        base, exp_ = pow_.base, pow_.exp
        if exp_.is_Integer or x not in base.free_symbols: # skip 1/x and sqrt(2)
            continue
        if not exp_.is_Rational: # exclude (x**2 + 1)**pi
            return None
        if exp_.q != 2:
            return None
        base_poly = base.as_poly(x)
        if base_poly is None or base_poly.degree() != 2: # exclude cube polynomial roots and other radicals
            return None
        aa = base_poly.nth(0)
        bb = base_poly.nth(1)
        cc = base_poly.nth(2)
        R = base_poly.as_expr()
        if base0 is None:
            base0 = R
            a0, b0, c0 = aa, bb, cc
            powers.append(pow_)
            exps.append(exp_)
            ratios.append(S.One)
        else:
            power_ratio = (powsimp(Pow(R, exp_) / Pow(base0, exp_), force=False)).cancel()
            if power_ratio.has(x):
                return None
            powers.append(pow_)
            exps.append(exp_)
            ratios.append(power_ratio)
    if base0 is None:
        return None

    pieces: list[tuple[Rule, Boolean]] = []
    delta = 4*a0*c0 - b0**2
    # substitution not valid for c0 = 0 and delta = 0
    c_zero_cond = Eq(c0, 0)
    delta_zero_cond = Eq(delta, 0)

    def _delta_zero_step():
        shift = x + b0/(2*c0)
        rewritten_base = c0*shift**2
        subs_dict = {pow_i: ratio_i*(rewritten_base)**exp_i for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
        rewritten = integrand.xreplace(subs_dict)
        step = yield IntegralInfo(rewritten, x)
        return RewriteRule(integrand, x, rewritten, step)

    def _c_zero_step():
        degenerate_integrand = integrand.subs(c0, 0)
        if b0.is_zero:
            step = yield IntegralInfo(degenerate_integrand, x)
        else:
            step = yield from sqrt_fractional_linear_rule(IntegralInfo(degenerate_integrand, x))
            if step is None:
                # since calling directly sqrt_fractional_linear_rule could return None we create a DontKnowRule
                step = DontKnowRule(degenerate_integrand, x)
        return step

    def _general_euler_step():
        s = Dummy("s")
        subs_dict = { pow_i: ratio_i * s**(2*exp_i) for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
        rewritten = integrand.xreplace(subs_dict)
        numer, denom = rewritten.as_numer_denom()
        if numer.as_poly(x, s) is None or denom.as_poly(x, s) is None:
            return None
        # Euler's second substitution (u = sqrt(R) + sqrt(c)*x)
        u = Dummy("u")
        sqrt_c0 = sqrt(c0)
        x_u = (u**2 - a0)/(b0 + 2*sqrt_c0*u)
        s_u = u - sqrt_c0*x_u
        dx_u = 2*(b0*u + sqrt_c0*(u**2 + a0))/(b0 + 2*sqrt_c0*u)**2
        substituted = rewritten.xreplace({x: x_u, s: s_u}) * dx_u
        substep = yield IntegralInfo(substituted, u)
        u_func = sqrt(base0) + sqrt_c0*x
        return URule(integrand, x, u, u_func, substep)

    if delta_zero_cond is S.true:
        general_step = yield from _delta_zero_step()
        if general_step.contains_dont_know():
            return None
    else:
        general_step = yield from _general_euler_step()
        if general_step is None or general_step.contains_dont_know():
            return None
    if c0.is_zero is None:
        pieces.append(((yield from _c_zero_step()), c_zero_cond))
    if delta.is_zero is None:
        pieces.append(((yield from _delta_zero_step()), delta_zero_cond))
    if pieces:
        pieces.append((general_step, S.true))
        general_step = PiecewiseRule(integrand, x, pieces)
    return general_step

def sqrt_quadratic_rule(integral: IntegralInfo, degenerate=True):
    # integrate f(x) * (a + b*x + c*x**2)**(n/2),
    # where f(x) is a polynomial and n is an odd integer
    starting_integrand, x = integral

    f = S.One
    root_base = None
    root_exp: Expr = S.Zero

    # collect radicals
    for factor in Mul.make_args(starting_integrand):
        if not factor.has(x):
            f *= factor
            continue
        base, exp = factor.as_base_exp()
        if exp.is_Integer is True:
            f *= factor
            continue
        # exclude x**pi
        if exp.is_Rational is not True:
            return None
        base_poly = base.as_poly(x)
        # exclude sqrt(log(x))
        if base_poly is None or base_poly.degree() != 2:
            return None
        base = base_poly.as_expr()
        if root_base is None:
            root_base = base
            root_exp = exp
            continue
        reference_power = Pow(root_base, exp)
        ratio = powsimp(factor/reference_power, force=False).cancel()
        if ratio.has(x):
            return None
        f *= ratio
        root_exp += exp

    if root_base is None:
        return None
    f_poly = f.as_poly(x)
    if f_poly is None:
        return None
    n = 2*root_exp
    if n.is_Integer is not True or n.is_odd is not True:
        return None
    root_poly = root_base.as_poly(x)
    a = root_poly.nth(0)
    b = root_poly.nth(1)
    c = root_poly.nth(2)
    root_base = a + b*x + c*x**2
    integrand = f*Pow(root_base, n/2)

    generic_cond = Ne(c, 0)
    if not degenerate or generic_cond is S.true:
        degenerate_step = None
    elif b.is_zero:
        degenerate_step = yield IntegralInfo(f*sqrt(a)**n, x)
    else:
        degenerate_integrand = f*sqrt(a + b*x)**n
        degenerate_step = yield from sqrt_fractional_linear_rule(IntegralInfo(degenerate_integrand, x))
        if degenerate_step is None:
            # since  sqrt_fractional_linear_rule does not guarantee a solution
            # we create a DontKnowRule so that _add_degenerate_step adds the degenerate condition
            degenerate_step = DontKnowRule(degenerate_integrand, x)

    def sqrt_quadratic_denom_rule(numer_poly: Poly, integrand: Expr):
        denom = sqrt(a+b*x+c*x**2)
        deg = numer_poly.degree()
        if deg <= 1:
            # integrand == (d+e*x)/sqrt(a+b*x+c*x**2)
            e, d = numer_poly.all_coeffs() if deg == 1 else (S.Zero, numer_poly.as_expr())
            # rewrite numerator to A*(2*c*x+b) + B
            A = e/(2*c)
            B = d-A*b
            pre_substitute = (2*c*x+b)/denom
            constant_step: Rule | None = None
            linear_step: Rule | None = None
            if A != 0:
                u = Dummy("u")
                pow_rule = PowerRule(1/sqrt(u), u, u, -S.Half)
                linear_step = URule(pre_substitute, x, u, a+b*x+c*x**2, pow_rule)
                if A != 1:
                    linear_step = ConstantTimesRule(A*pre_substitute, x, A, pre_substitute, linear_step)
            if B != 0:
                constant_step = yield from inverse_trig_rule(IntegralInfo(1/denom, x), degenerate=False)
                if B != 1:
                    constant_step = ConstantTimesRule(B/denom, x, B, 1/denom, constant_step)  # type: ignore
            if linear_step and constant_step:
                add = Add(A*pre_substitute, B/denom, evaluate=False)
                step: Rule | None = RewriteRule(integrand, x, add, AddRule(add, x, [linear_step, constant_step]))
            else:
                step = linear_step or constant_step
        else:
            coeffs = numer_poly.all_coeffs()
            step = SqrtQuadraticDenomRule(integrand, x, a, b, c, coeffs)
        return step

    def sqrt_quadratic_reduction_rule(integrand: Expr, n: Expr, const: Expr):
        # Implementation of Gradshteyn & Ryzhik 2.263.3
        k = (-n - 1) // 2
        delta = 4*a*c - b**2
        R = c*x**2 + b*x + a

        def double_root_step():
            square_base = sqrt(c)*x + b/(2*sqrt(c))
            nested = Pow(Pow(square_base, 2, evaluate=False), S(n)/2, evaluate=False)
            rewritten = const*nested
            substep = nested_pow_rule(IntegralInfo(rewritten, x))
            return RewriteRule(integrand, x, rewritten, substep)

        if delta.is_zero is True:
            return double_root_step()

        # we divide by delta, then delta  has to be != 0
        term_denom = (2*k - 1) * delta * (R**(S(2*k - 1)/2))
        constant_term = const*2*(2*c*x+b) / term_denom
        coeff = (8*c*(k-1))/((2*k-1) * delta)
        expr = const * R**(S(1)/2 - k)

        rewrite_expr = Derivative(constant_term, x) + coeff * expr
        derive_expr = Derivative(constant_term, x)
        derive_step = yield IntegralInfo(derive_expr, x)

        if coeff == 0:
            substep = derive_step
        else:
            next_step = yield IntegralInfo(expr, x)
            if not next_step:
                next_step = DontKnowRule(expr, x)

            substep = AddRule(
                rewrite_expr,
                x,
                [
                    derive_step,
                    ConstantTimesRule(
                        coeff * expr,
                        x,
                        coeff,
                        expr,
                        next_step
                    )
                ]
            )

        debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewrite_expr, x))
        nondegenerate_step = RewriteRule(integrand, x, rewrite_expr, substep)
        if delta.is_zero is None:
            return _add_degenerate_step(
                Ne(delta, 0),
                nondegenerate_step,
                double_root_step(),
            )
        return nondegenerate_step

    def sqrt_quadratic_polynomial_reduction_rule():
        # reduce non-constant polynomial numerators by writing f = q*R + r,
        # then split the linear remainder into a multiple of R' and a constant.
        terms = []
        steps = []
        root_base = c*x**2 + b*x + a
        root_poly = Poly(root_base, x)
        quotient, rest = f_poly.div(root_poly)
        if not quotient.is_zero:
            # n is increasing by 2 at each step, we will fall in one of the cases above
            quotient_integrand = quotient.as_expr() * sqrt(root_base)**(n + 2)
            quotient_step = yield from sqrt_quadratic_rule(IntegralInfo(quotient_integrand, x), degenerate=False)
            terms.append(quotient_integrand)
            steps.append(quotient_step)
        if not rest.is_zero:
            # split the linear remainder as A*R' + B, where R' = 2*c*x + b.
            e = rest.nth(1)
            d = rest.nth(0)
            A = e/(2*c)
            B = d - A*b
            if A != 0:
                u = Dummy("u")
                # solved by substitution
                base = (2*c*x + b) * sqrt(root_base)**n
                term = A * base
                pow_rule = PowerRule(u**(S(n)/2), u, u, S(n)/2)
                u_step = URule(base, x, u, root_base, pow_rule)
                if A != 1:
                    u_step = ConstantTimesRule(term, x, A, base, u_step)
                terms.append(term)
                steps.append(u_step)
            if B != 0:
                term = B * sqrt(root_base)**n
                # constant case already managed
                const_step = yield from sqrt_quadratic_reduction_rule(term, n, B)
                terms.append(term)
                steps.append(const_step)
        rewritten = Add(*terms, evaluate=False)
        debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten, x))
        if len(steps) == 1:
            substep = steps[0]
        else:
            substep = AddRule(rewritten, x, steps)
        return RewriteRule(integrand, x, rewritten, substep)

    if n > 0:  # rewrite poly * sqrt(s)**(2*k-1) to poly*s**k / sqrt(s)
        numer_poly = f_poly * (a+b*x+c*x**2)**((n+1)/2)
        rewritten = numer_poly.as_expr()/sqrt(a+b*x+c*x**2)
        substep = yield from sqrt_quadratic_denom_rule(numer_poly, rewritten)
        debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten, x))
        generic_step = RewriteRule(integrand, x, rewritten, substep)
    elif n == -1:
        generic_step = yield from sqrt_quadratic_denom_rule(f_poly, integrand)
    elif f_poly.degree() == 0:
        # The numerator must be a const, the formula assumes this
        generic_step = yield from sqrt_quadratic_reduction_rule(integrand, n, f)
    else:
        generic_step = yield from sqrt_quadratic_polynomial_reduction_rule()
    step = _add_degenerate_step(generic_cond, generic_step, degenerate_step)
    if integrand != starting_integrand:
        return RewriteRule(starting_integrand, x, integrand, step)
    return step


def hyperbolic_rule(integral: tuple[Expr, Symbol]):
    integrand, symbol = integral
    if isinstance(integrand, HyperbolicFunction) and integrand.args[0] == symbol:
        if integrand.func == sinh:
            return SinhRule(integrand, symbol)
        if integrand.func == cosh:
            return CoshRule(integrand, symbol)
        u = Dummy('u')
        if integrand.func == tanh:
            rewritten = sinh(symbol)/cosh(symbol)
            debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten, symbol))
            return RewriteRule(integrand, symbol, rewritten,
                   URule(rewritten, symbol, u, cosh(symbol), ReciprocalRule(1/u, u, u)))
        if integrand.func == coth:
            rewritten = cosh(symbol)/sinh(symbol)
            debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten, symbol))
            return RewriteRule(integrand, symbol, rewritten,
                   URule(rewritten, symbol, u, sinh(symbol), ReciprocalRule(1/u, u, u)))
        else:
            rewritten = integrand.rewrite(tanh)
            debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten, symbol))
            if integrand.func == sech:
                return RewriteRule(integrand, symbol, rewritten,
                       URule(rewritten, symbol, u, tanh(symbol/2),
                       ArctanRule(2/(u**2 + 1), u, S(2), S.One, S.One)))
            if integrand.func == csch:
                debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten, symbol))
                return RewriteRule(integrand, symbol, rewritten,
                       URule(rewritten, symbol, u, tanh(symbol/2),
                       ReciprocalRule(1/u, u, u)))

@cacheit
def make_wilds(symbol):
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    m = Wild('m', exclude=[symbol], properties=[lambda n: isinstance(n, Integer)])
    n = Wild('n', exclude=[symbol], properties=[lambda n: isinstance(n, Integer)])

    return a, b, m, n

@cacheit
def sincos_pattern(symbol):
    a, b, m, n = make_wilds(symbol)
    pattern = sin(a*symbol)**m * cos(b*symbol)**n

    return pattern, a, b, m, n

@cacheit
def tansec_pattern(symbol):
    a, b, m, n = make_wilds(symbol)
    pattern = tan(a*symbol)**m * sec(b*symbol)**n

    return pattern, a, b, m, n

@cacheit
def cotcsc_pattern(symbol):
    a, b, m, n = make_wilds(symbol)
    pattern = cot(a*symbol)**m * csc(b*symbol)**n

    return pattern, a, b, m, n

@cacheit
def heaviside_pattern(symbol):
    m = Wild('m', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    g = Wild('g')
    pattern = Heaviside(m*symbol + b) * g

    return pattern, m, b, g

def uncurry(func):
    def uncurry_rl(args):
        debug("Uncurry: {} with args: {}".format(func, args))
        return func(*args)
    return uncurry_rl

def trig_rewriter(rewrite):
    def trig_rewriter_rl(args):
        a, b, m, n, integrand, symbol = args
        rewritten = rewrite(a, b, m, n, integrand, symbol)
        if rewritten != integrand:
            debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewrite, symbol))
            substep = yield IntegralInfo(rewritten, symbol)
            return RewriteRule(integrand, symbol, rewritten, substep)
    return trig_rewriter_rl

sincos_botheven_condition = uncurry(
    lambda a, b, m, n, i, s: m.is_even and n.is_even and
    m.is_nonnegative and n.is_nonnegative)

sincos_botheven = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (((1 - cos(2*a*symbol)) / 2) ** (m / 2)) *
                                    (((1 + cos(2*b*symbol)) / 2) ** (n / 2)) ))

sincos_sinodd_condition = uncurry(lambda a, b, m, n, i, s: m.is_odd and m >= 3)

sincos_sinodd = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (1 - cos(a*symbol)**2)**((m - 1) / 2) *
                                    sin(a*symbol) *
                                    cos(b*symbol) ** n))

sincos_cosodd_condition = uncurry(lambda a, b, m, n, i, s: n.is_odd and n >= 3)

sincos_cosodd = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (1 - sin(b*symbol)**2)**((n - 1) / 2) *
                                    cos(b*symbol) *
                                    sin(a*symbol) ** m))

tansec_seceven_condition = uncurry(lambda a, b, m, n, i, s: n.is_even and n >= 4)
tansec_seceven = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (1 + tan(b*symbol)**2) ** (n/2 - 1) *
                                    sec(b*symbol)**2 *
                                    tan(a*symbol) ** m ))

tansec_tanodd_condition = uncurry(lambda a, b, m, n, i, s: m.is_odd)
tansec_tanodd = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (sec(a*symbol)**2 - 1) ** ((m - 1) / 2) *
                                     tan(a*symbol) *
                                     sec(b*symbol) ** n ))

tan_tansquared_condition = uncurry(lambda a, b, m, n, i, s: m == 2 and n == 0)
tan_tansquared = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( sec(a*symbol)**2 - 1))

cotcsc_csceven_condition = uncurry(lambda a, b, m, n, i, s: n.is_even and n >= 4)
cotcsc_csceven = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (1 + cot(b*symbol)**2) ** (n/2 - 1) *
                                    csc(b*symbol)**2 *
                                    cot(a*symbol) ** m ))

cotcsc_cotodd_condition = uncurry(lambda a, b, m, n, i, s: m.is_odd)
cotcsc_cotodd = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (csc(a*symbol)**2 - 1) ** ((m - 1) / 2) *
                                    cot(a*symbol) *
                                    csc(b*symbol) ** n ))

def trig_sincos_rule(integral):
    integrand, symbol = integral

    if any(integrand.has(f) for f in (sin, cos)):
        pattern, a, b, m, n = sincos_pattern(symbol)
        match = integrand.match(pattern)
        if not match:
            return

        return multiplexer({
            sincos_botheven_condition: sincos_botheven,
            sincos_sinodd_condition: sincos_sinodd,
            sincos_cosodd_condition: sincos_cosodd
        })(tuple(
            [match.get(i, S.Zero) for i in (a, b, m, n)] +
            [integrand, symbol]))

def trig_tansec_rule(integral):
    integrand, symbol = integral

    integrand = integrand.subs({
        1 / cos(symbol): sec(symbol)
    })

    if any(integrand.has(f) for f in (tan, sec)):
        pattern, a, b, m, n = tansec_pattern(symbol)
        match = integrand.match(pattern)
        if not match:
            return

        return multiplexer({
            tansec_tanodd_condition: tansec_tanodd,
            tansec_seceven_condition: tansec_seceven,
            tan_tansquared_condition: tan_tansquared
        })(tuple(
            [match.get(i, S.Zero) for i in (a, b, m, n)] +
            [integrand, symbol]))

def trig_cotcsc_rule(integral):
    integrand, symbol = integral
    integrand = integrand.subs({
        1 / sin(symbol): csc(symbol),
        1 / tan(symbol): cot(symbol),
        cos(symbol) / tan(symbol): cot(symbol)
    })

    if any(integrand.has(f) for f in (cot, csc)):
        pattern, a, b, m, n = cotcsc_pattern(symbol)
        match = integrand.match(pattern)
        if not match:
            return

        return multiplexer({
            cotcsc_cotodd_condition: cotcsc_cotodd,
            cotcsc_csceven_condition: cotcsc_csceven
        })(tuple(
            [match.get(i, S.Zero) for i in (a, b, m, n)] +
            [integrand, symbol]))

def trig_sindouble_rule(integral):
    integrand, symbol = integral
    a = Wild('a', exclude=[sin(2*symbol)])
    match = integrand.match(sin(2*symbol)*a)
    if match:
        sin_double = 2*sin(symbol)*cos(symbol)/sin(2*symbol)
        rewritten = integrand * sin_double
        debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, integrand * sin_double, symbol))
        substeps = yield IntegralInfo(rewritten, symbol)
        return RewriteRule(integrand, symbol, rewritten, substeps)

def trig_powers_products_rule(integral):
    # The subrules may return a generator (e.g. through multiplexer);
    # delegate to it so its requests reach the solver, and fall through to
    # the next subrule if it produced no result.
    for rule in (trig_sincos_rule, trig_tansec_rule,
                 trig_cotcsc_rule, trig_sindouble_rule):
        result = rule(integral)
        if isinstance(result, GeneratorType):
            result = yield from result
        if result is not None:
            return result

def trig_substitution_rule(integral):
    integrand, symbol = integral
    A = Wild('a', exclude=[0, symbol])
    B = Wild('b', exclude=[0, symbol])
    theta = Dummy("theta")
    target_pattern = A + B*symbol**2

    matches = integrand.find(target_pattern)
    for expr in matches:
        match = expr.match(target_pattern)
        a = match.get(A, S.Zero)
        b = match.get(B, S.Zero)

        a_positive = ((a.is_number and a > 0) or a.is_positive)
        b_positive = ((b.is_number and b > 0) or b.is_positive)
        a_negative = ((a.is_number and a < 0) or a.is_negative)
        b_negative = ((b.is_number and b < 0) or b.is_negative)
        x_func = None
        if a_positive and b_positive:
            # a**2 + b*x**2. Assume sec(theta) > 0, -pi/2 < theta < pi/2
            x_func = (sqrt(a)/sqrt(b)) * tan(theta)
            # Do not restrict the domain: tan(theta) takes on any real
            # value on the interval -pi/2 < theta < pi/2 so x takes on
            # any value
            restriction = True
        elif a_positive and b_negative:
            # a**2 - b*x**2. Assume cos(theta) > 0, -pi/2 < theta < pi/2
            constant = sqrt(a)/sqrt(-b)
            x_func = constant * sin(theta)
            restriction = And(symbol > -constant, symbol < constant)
        elif a_negative and b_positive:
            # b*x**2 - a**2. Assume sin(theta) > 0, 0 < theta < pi
            constant = sqrt(-a)/sqrt(b)
            x_func = constant * sec(theta)
            restriction = And(symbol > -constant, symbol < constant)
        if x_func:
            # Manually simplify sqrt(trig(theta)**2) to trig(theta)
            # Valid due to assumed domain restriction
            substitutions = {}
            for f in [sin, cos, tan,
                      sec, csc, cot]:
                substitutions[sqrt(f(theta)**2)] = f(theta)
                substitutions[sqrt(f(theta)**(-2))] = 1/f(theta)

            replaced = integrand.subs(symbol, x_func).trigsimp()
            replaced = manual_subs(replaced, substitutions)
            if not replaced.has(symbol):
                replaced *= manual_diff(x_func, theta)
                replaced = replaced.trigsimp()
                secants = replaced.find(1/cos(theta))
                if secants:
                    replaced = replaced.xreplace({
                        1/cos(theta): sec(theta)
                    })

                substep = yield IntegralInfo(replaced, theta)
                if not substep.contains_dont_know():
                    return TrigSubstitutionRule(integrand, symbol,
                        theta, x_func, replaced, substep, restriction)

def heaviside_rule(integral):
    integrand, symbol = integral
    pattern, m, b, g = heaviside_pattern(symbol)
    match = integrand.match(pattern)
    if match and 0 != match[g]:
        # f = Heaviside(m*x + b)*g
        substep = yield IntegralInfo(match[g], symbol)
        m, b = match[m], match[b]
        return HeavisideRule(integrand, symbol, m*symbol + b, -b/m, substep)


def dirac_delta_rule(integral: IntegralInfo):
    integrand, x = integral
    if len(integrand.args) == 1:
        n = S.Zero
    else:
        n = integrand.args[1] # type: ignore
    if not n.is_Integer or n < 0:
        return
    a, b = Wild('a', exclude=[x]), Wild('b', exclude=[x, 0])
    match = integrand.args[0].match(a+b*x)
    if not match:
        return
    a, b = match[a], match[b]
    generic_cond = Ne(b, 0)
    if generic_cond is S.true:
        degenerate_step = None
    else:
        degenerate_step = ConstantRule(DiracDelta(a, n), x)
    generic_step = DiracDeltaRule(integrand, x, n, a, b)
    return _add_degenerate_step(generic_cond, generic_step, degenerate_step)


def substitution_rule(integral):
    integrand, symbol = integral

    u_var = Dummy("u")
    substitutions = find_substitutions(integrand, symbol, u_var)
    count = 0
    if substitutions:
        # Ask the solver whether all alternative substitutions should be
        # collected or only the first workable one is wanted.
        branch = yield BranchQuery()
        if branch:
            debug("List of Substitution Rules")
        ways = []
        factored_integrand = integrand.factor()
        _, denom_integrand = factored_integrand.as_numer_denom()
        for u_func, c, substituted in substitutions:
            subrule = yield IntegralInfo(substituted, u_var)
            count = count + 1
            if branch:
                debug("Rule {}: {}".format(count, subrule))

            if subrule.contains_dont_know():
                continue

            if simplify(c - 1) != 0:
                _, denom_c = c.as_numer_denom()
                if subrule:
                    subrule = ConstantTimesRule(c * substituted, u_var, c, substituted, subrule)

                if denom_c.free_symbols:
                    pieces = []
                    factors_denom_c = factor_list(denom_c)[1]
                    for pole, _ in factors_denom_c:
                        # only substitute poles introduced by the constant c if they were not already poles of the original integrand
                        if not _if_zero_implies_zero(pole, denom_integrand):
                            rewritten_integral = manual_subs(factored_integrand, pole, 0)
                            debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten_integral, symbol))
                            # additional check not to replace a if it is not valid (for example ln(a*x))
                            if rewritten_integral.has(S.ComplexInfinity, S.Infinity, S.NegativeInfinity, S.NaN):
                                continue
                            substep = yield IntegralInfo(rewritten_integral, symbol)

                            if substep:
                                substep = RewriteRule(integrand, symbol, rewritten_integral, substep)
                                pieces.append((
                                    substep,
                                    Eq(pole, 0)
                                ))
                    if pieces:
                        pieces.append((subrule, True))
                        subrule = PiecewiseRule(substituted, symbol, pieces)

            rule = URule(integrand, symbol, u_var, u_func, subrule)
            if not branch:
                return rule
            ways.append(rule)

        if len(ways) > 1:
            return AlternativeRule(integrand, symbol, ways)
        elif ways:
            return ways[0]


def partial_fractions_rule(integral):
    integrand, symbol = integral
    if not integrand.is_rational_function(symbol):
        return

    rewritten = integrand.apart(symbol)
    if rewritten == integrand:
        # If apart cannot decompose the rational function any further,
        # use ratint as the final fallback for rational integration.
        return RatintRule(integrand, symbol)

    substep = yield IntegralInfo(rewritten, symbol)
    if not isinstance(substep, DontKnowRule):
        return RewriteRule(integrand, symbol, rewritten, substep)

cancel_rule = rewriter(
    # lambda integrand, symbol: integrand.is_algebraic_expr(),
    # lambda integrand, symbol: isinstance(integrand, Mul),
    lambda integrand, symbol: not integrand.is_rational_function(symbol),
    lambda integrand, symbol: integrand.cancel())

distribute_expand_rule = rewriter(
    lambda integrand, symbol: (
        (isinstance(integrand, (Pow, Mul))
         or all(arg.is_Pow or arg.is_polynomial(symbol) for arg in integrand.args))
        and not integrand.is_rational_function(symbol)
    ),
    lambda integrand, symbol: integrand.expand()
)

trig_expand_rule = rewriter(
    # If there are trig functions with different arguments, expand them
    lambda integrand, symbol: (
        len({a.args[0] for a in integrand.atoms(TrigonometricFunction)}) > 1),
    lambda integrand, symbol: integrand.expand(trig=True))

def derivative_rule(integral):
    integrand = integral[0]
    diff_variables = integrand.variables
    undifferentiated_function = integrand.expr
    integrand_variables = undifferentiated_function.free_symbols

    if integral.symbol in integrand_variables:
        if integral.symbol in diff_variables:
            return DerivativeRule(*integral)
        else:
            return DontKnowRule(integrand, integral.symbol)
    else:
        return ConstantRule(*integral)

def rewrites_rule(integral):
    integrand, symbol = integral

    if integrand.match(1/cos(symbol)):
        rewritten = integrand.subs(1/cos(symbol), sec(symbol))
        debug("Integral: {} is rewritten with {} on symbol: {}".format(integrand, rewritten, symbol))
        substep = yield IntegralInfo(rewritten, symbol)
        return RewriteRule(integrand, symbol, rewritten, substep)

def fallback_rule(integral):
    return DontKnowRule(*integral)


# Subproblems are compared with the same dummy variable in place of the
# integration variable for them to match (used for loop detection and for
# counting "u" choices of integration by parts).
_cache_dummy = Dummy("z")

# Canonical replacements for the other Dummy variables of a subproblem,
# keyed by (position of first appearance, assumptions).
_canonical_dummies: dict[tuple, Dummy] = {}


def _integral_cache_key(integrand, symbol):
    """Cache key of a subproblem and the replacements that produced it:
    the integration variable and any other Dummy variables are replaced by
    canonical ones, so that structurally identical subproblems that differ
    only in the identity of their dummy variables (as minted by the
    substitution rules) compare equal."""
    key = integrand.xreplace({symbol: _cache_dummy})
    replacements = {symbol: _cache_dummy}
    index = 0
    for node in preorder_traversal(key):
        if (isinstance(node, Dummy) and node is not _cache_dummy
                and node not in replacements):
            assumptions = tuple(sorted(node.assumptions0.items()))
            canonical = _canonical_dummies.get((index, assumptions))
            if canonical is None:
                canonical = Dummy("z%d" % index, **node.assumptions0)
                _canonical_dummies[(index, assumptions)] = canonical
            replacements[node] = canonical
            index += 1
    if len(replacements) > 1:
        key = key.xreplace(replacements)
    return key, replacements


def _rule_xreplace(rule, mapping):
    """A copy of a Rule tree with the mapping applied to every Expr
    field, to rename the variable and dummies of a cached solution."""
    def map_value(v):
        if isinstance(v, Rule):
            return _rule_xreplace(v, mapping)
        if isinstance(v, Basic):
            return v.xreplace(mapping)
        if isinstance(v, tuple):
            return tuple(map_value(i) for i in v)
        if isinstance(v, list):
            return [map_value(i) for i in v]
        return v
    return type(rule)(*[map_value(getattr(rule, s))
                        for s in rule._get_slots()])


def _integral_key(integral):
    integrand, symbol = integral

    if symbol not in integrand.free_symbols:
        return Number
    for cls in (Symbol, TrigonometricFunction, OrthogonalPolynomial):
        if isinstance(integrand, cls):
            return cls
    return type(integrand)


def _integral_is_subclass(*klasses):
    def _check(integral):
        k = _integral_key(integral)
        return k and issubclass(k, klasses)
    return _check


class IntegrationSolver:
    """Performs the recursive calls on behalf of the integration rules and
    owns all the state of a single integration run.

    A rule function either returns a :class:`Rule` (or ``None``) directly,
    or is a generator. A generator rule does not recurse itself: whenever
    it needs a subproblem solved, it yields an :class:`IntegralInfo` with
    the parameters of the recursive call, and the solver sends the
    resulting :class:`Rule` back into the generator at the yield
    expression. Since a generator may yield any number of times, a rule can
    propose further subproblems after seeing the result of a previous one;
    the assembled :class:`Rule` is the generator's return value. Rules may
    also yield other requests that read or update per-run solver state (see
    :class:`PartsUCheck` and :class:`BranchQuery`).

    Everything that influences how the nested rules are applied (the
    loop-detection set, the integration-by-parts ``u`` counter and the
    options passed to :func:`integral_steps`) is stored on the solver
    instance instead of in module-global variables.
    """

    def __init__(self, max_depth: int | None = None, branch: bool = False,
                 max_calls: int | None = None,
                 soft_max_calls: int | None = None, **other_options):
        # Hard limit on the depth of nested subproblems (None = unlimited):
        # ``max_depth=1`` allows only rules that need no subintegrals.
        self.max_depth = max_depth
        # Total budget of subproblems dispatched by this solver (None =
        # unlimited); once exhausted, further subproblems are reported as
        # unsolvable, bounding the search on hard integrands.
        self.max_calls = max_calls
        # Soft budget: beyond this many subproblems the expensive
        # exploratory rules (substitution search, integration by parts,
        # expansions) are skipped and only the cheap closed-form rules
        # keep running, so hard integrands degrade gracefully instead of
        # exploding.
        self.soft_max_calls = soft_max_calls
        self._calls = 0
        # Whether to collect all applicable rules into AlternativeRule
        # nodes (True) or to keep only the first workable rule (False).
        self.branch = branch
        self.options = other_options
        # Subproblems on the current recursion path, to break cyclic integrals.
        self._active: set[Expr] = set()
        # Subproblems already found unsolvable, so that failing subtrees are
        # not explored again from sibling branches.  Failures that were
        # caused by a cycle break are not recorded, since those subproblems
        # may well be solvable on a different recursion path.
        self._dont_know: set[Expr] = set()
        # Complete solutions by cache key, with the replacements that
        # canonicalized them, for reuse on identical subproblems.
        self._solved: dict[Expr, tuple[Rule, dict]] = {}
        self._loop_breaks = 0
        # Uses of each "u" by integration by parts, to avoid infinite repetition.
        self._parts_u_count: dict[Expr, int] = defaultdict(int)
        self._strategy: Callable[[IntegralInfo], Rule] | None = None

    def solve(self, integrand, symbol) -> Rule:
        """Solve a (sub)integral by dispatching the rules on it.

        This replaces the recursive calls to ``integral_steps`` that the
        rules used to perform themselves.
        """
        # Every ancestor on the recursion path holds exactly one entry in
        # ``_active``, so its size is the current depth.
        if self.max_depth is not None and len(self._active) >= self.max_depth:
            # Recursion budget exhausted: report the subproblem as
            # unsolvable, so that callers fall back to shallower candidates.
            return DontKnowRule(integrand, symbol)
        self._calls += 1
        if self.max_calls is not None and self._calls > self.max_calls:
            return DontKnowRule(integrand, symbol)
        cachekey, canonical_map = _integral_cache_key(integrand, symbol)
        if cachekey in self._dont_know:
            return DontKnowRule(integrand, symbol)
        hit = self._solved.get(cachekey)
        if hit is not None:
            rule, stored_map = hit
            from_canonical = {c: e for e, c in canonical_map.items()}
            mapping = {orig: from_canonical[canonical]
                       for orig, canonical in stored_map.items()
                       if orig != from_canonical[canonical]}
            return _rule_xreplace(rule, mapping) if mapping else rule
        if cachekey in self._active:
            # Stop this attempt, because it leads around in a loop
            self._loop_breaks += 1
            return DontKnowRule(integrand, symbol)
        self._active.add(cachekey)
        loop_breaks = self._loop_breaks
        try:
            if self._strategy is None:
                self._strategy = self._build_strategy()
            result = self._strategy(IntegralInfo(integrand, symbol))
        finally:
            self._active.discard(cachekey)
        if isinstance(result, DontKnowRule):
            if loop_breaks == self._loop_breaks:
                # The failure is unconditional (no cycle break was
                # involved), so this subproblem need not be tried again.
                self._dont_know.add(cachekey)
        elif result is not None and not result.contains_dont_know():
            # Complete solutions do not depend on the recursion path and
            # can be reused for identical subproblems (up to renaming of
            # the variable and dummies).
            self._solved[cachekey] = (result, canonical_map)
        return result

    def run(self, rule, integral):
        """Apply a single rule to ``integral``, driving it if it is a
        generator, and return its result."""
        result = rule(integral)
        if isinstance(result, GeneratorType):
            result = self.run_generator(result)
        return result

    def run_generator(self, gen):
        """Iterate over the requests yielded by a generator rule, perform
        the recursive calls on its behalf, and return the rule's result."""
        response = None
        while True:
            try:
                request = gen.send(response)
            except StopIteration as e:
                return e.value
            if isinstance(request, IntegralInfo):
                response = self.solve(request.integrand, request.symbol)
            elif isinstance(request, PartsUCheck):
                response = self._check_parts_u(request.u_key)
            elif isinstance(request, BranchQuery):
                response = self.branch
            else:
                raise ValueError(
                    "unknown request from integration rule: %s" % (request,))

    def _check_parts_u(self, u_key) -> bool:
        # Set a limit on the number of times u can be used
        if self._parts_u_count[u_key] > 2:
            return False
        self._parts_u_count[u_key] += 1
        return True

    def _build_strategy(self):
        # Each leaf rule is wrapped with ``run`` so that, by the time the
        # combinators from sympy.strategies see a result, any generator has
        # already been driven to completion and only a Rule or None remains.
        def w(rule):
            @wraps(rule)
            def rule_runner(integral):
                return self.run(rule, integral)
            return rule_runner

        def expensive(rule):
            # skipped once the soft call budget is exhausted
            def guarded(integral):
                if (self.soft_max_calls is not None
                        and self._calls > self.soft_max_calls):
                    return None
                return self.run(rule, integral)
            return guarded

        return do_one(
            null_safe(w(special_function_rule)),
            null_safe(switch(_integral_key, {
                Pow: do_one(null_safe(w(power_rule)), null_safe(w(inverse_trig_rule)),
                            null_safe(w(quadratic_denom_rule)),
                            null_safe(w(sqrt_quadratic_rule)),
                            null_safe(w(sqrt_fractional_linear_rule)),
                            null_safe(w(euler_substitution_rule))),
                Symbol: w(power_rule),
                exp: w(exp_rule),
                Add: w(add_rule),
                Mul: do_one(null_safe(w(mul_rule)), null_safe(w(trig_product_rule)),
                            null_safe(w(trig_product_to_sum_rule)),
                            null_safe(w(tabular_parts_rule)),
                            null_safe(w(heaviside_rule)), null_safe(w(quadratic_denom_rule)),
                            null_safe(w(sqrt_quadratic_rule)),
                            null_safe(w(sqrt_fractional_linear_rule)),
                            null_safe(w(euler_substitution_rule)),
                            null_safe(w(trig_cmplx_exp_rule))),
                Derivative: w(derivative_rule),
                TrigonometricFunction: w(trig_rule),
                Heaviside: w(heaviside_rule),
                DiracDelta: w(dirac_delta_rule),
                OrthogonalPolynomial: w(orthogonal_poly_rule),
                Number: w(constant_rule)
            })),
            do_one(
                null_safe(w(trig_rule)),
                null_safe(w(hyperbolic_rule)),
                null_safe(expensive(alternatives(
                    w(rewrites_rule),
                    # cheap specific closed-form rules come before the
                    # general searches, so that table entries win over
                    # (and cost less than) exploratory substitutions
                    w(dilog_rule),
                    w(exp_over_linear_rule),
                    w(log_of_trig_rule),
                    w(polylog_ladder_rule),
                    w(substitution_rule),
                    condition(
                        _integral_is_subclass(Mul, Pow),
                        w(partial_fractions_rule)),
                    condition(
                        _integral_is_subclass(Mul, Pow),
                        w(cancel_rule)),
                    condition(
                        _integral_is_subclass(Mul),
                        w(combine_power_rule)),
                    condition(
                        _integral_is_subclass(Mul, log,
                                              *inverse_trig_functions,
                                              *special_error_functions),
                        w(parts_rule)),
                    condition(
                        _integral_is_subclass(Mul, Pow),
                        w(distribute_expand_rule)),
                    w(trig_powers_products_rule),
                    w(trig_expand_rule),
                    w(power_substitution_rule),
                    w(exp_power_rewrite_rule),
                    w(exp_times_rational_rule),
                    w(piecewise_continuous_rule),
                    branch=self.branch
                ))),
                null_safe(condition(_integral_is_subclass(Mul, Pow), expensive(nested_pow_rule))),
                null_safe(condition(_integral_is_subclass(Mul, Pow), expensive(weierstrass_substitution))),
                null_safe(condition(_integral_is_subclass(Mul, Pow), expensive(hyperbolic_rational_substitution))),
                null_safe(condition(_integral_is_subclass(Mul, Pow), expensive(odd_power_trig_substitution))),
                null_safe(condition(_integral_is_subclass(Mul, Pow), expensive(linear_power_product_rule))),
                null_safe(condition(_integral_is_subclass(Mul, Pow), expensive(inverse_function_substitution_rule))),
                null_safe(condition(_integral_is_subclass(Mul, Pow), expensive(trig_half_angle_square_rule))),
                null_safe(condition(_integral_is_subclass(Mul), expensive(log_times_rational_rule))),
                null_safe(condition(_integral_is_subclass(Mul), expensive(exp_over_linear_rule))),
                null_safe(condition(_integral_is_subclass(Mul), expensive(exp_times_rational_rule))),
                null_safe(condition(_integral_is_subclass(Mul), expensive(log_inverse_parts_rule))),
                null_safe(expensive(power_substitution_rule)),
                null_safe(expensive(trig_substitution_rule))
            ),
            w(fallback_rule))


def integral_steps(integrand, symbol, **options):
    """Returns the steps needed to compute an integral.

    Explanation
    ===========

    This function attempts to mirror what a student would do by hand as
    closely as possible.

    SymPy Gamma uses this to provide a step-by-step explanation of an
    integral. The code it uses to format the results of this function can be
    found at
    https://github.com/sympy/sympy_gamma/blob/master/app/logic/intsteps.py.

    By default, only the first rule that successfully applies at each step
    is kept, even if multiple rules are applicable. Passing ``branch=True``
    causes all applicable rules to be retained and combined into an
    ``AlternativeRule``, so that alternative solution paths are preserved
    rather than discarded.

    Examples
    ========

    >>> from sympy import exp, sin
    >>> from sympy.integrals.manualintegrate import integral_steps
    >>> from sympy.abc import x
    >>> print(repr(integral_steps(exp(x) / (1 + exp(2 * x)), x))) \
    # doctest: +NORMALIZE_WHITESPACE
    URule(integrand=exp(x)/(exp(2*x) + 1), variable=x, u_var=_u, u_func=exp(x),
    substep=ArctanRule(integrand=1/(_u**2 + 1), variable=_u, a=1, b=1, c=1))
    >>> print(repr(integral_steps(sin(x), x))) \
    # doctest: +NORMALIZE_WHITESPACE
    SinRule(integrand=sin(x), variable=x)
    >>> print(repr(integral_steps((x**2 + 3)**2, x))) \
    # doctest: +NORMALIZE_WHITESPACE
    RewriteRule(integrand=(x**2 + 3)**2, variable=x, rewritten=x**4 + 6*x**2 + 9,
    substep=AddRule(integrand=x**4 + 6*x**2 + 9, variable=x,
    substeps=[PowerRule(integrand=x**4, variable=x, base=x, exp=4),
    ConstantTimesRule(integrand=6*x**2, variable=x, constant=6, other=x**2,
    substep=PowerRule(integrand=x**2, variable=x, base=x, exp=2)),
    ConstantRule(integrand=9, variable=x)]))

    Parameters
    ==========

    integrand : Expr
        The expression to integrate.
    symbol : Symbol
        The variable of integration.
    branch : bool, optional
        If True, collect all applicable rules into an ``AlternativeRule``
        instead of returning only the first workable one. Defaults to False.
    max_depth : int, optional
        Hard limit on the depth of nested subproblems; deeper subproblems
        are reported as ``DontKnowRule``. Defaults to None (unlimited).
    max_calls : int, optional
        Total budget of subproblems the solver may dispatch; once
        exhausted, further subproblems are reported as ``DontKnowRule``,
        bounding the search on hard integrands. Defaults to None
        (unlimited).
    soft_max_calls : int, optional
        Soft budget: beyond this many subproblems the expensive
        exploratory rules (substitution search, integration by parts,
        expansions) are skipped while the cheap closed-form rules keep
        running, so hard integrands degrade gracefully. Defaults to
        None (disabled).

    Returns
    =======

    rule : Rule
        The first step; most rules have substeps that must also be
        considered. These substeps can be evaluated using ``manualintegrate``
        to obtain a result.

    """
    return IntegrationSolver(**options).solve(integrand, symbol)


def _combine_logs(result, symbol):
    """Combine sums of logarithms with rationally related coefficients,

        alpha*log(f1) + beta*log(f2) -> (m/n)*log(f1**p*f2**q),

    which for anti-derivatives is always permissible and always
    preferable: the singular points of the collected logarithm are a
    subset of those of the sum (D.J. Jeffrey, 1993).  Applied to the
    top-level sum and inside Piecewise branches.
    """
    if isinstance(result, Piecewise):
        return Piecewise(*[(_combine_logs(e, symbol), c)
                           for e, c in result.args])
    if not result.has(log):
        return result
    log_terms = []
    others = []
    queue = list(Add.make_args(result))
    while queue:
        term = queue.pop()
        coefficient, rest = term.as_independent(log)
        if coefficient.has(symbol):
            others.append(term)
        elif isinstance(rest, log) and rest.args[0].has(symbol):
            log_terms.append((coefficient, rest.args[0]))
        elif (isinstance(rest, Add) and rest.has(log)
                and rest.has(TrigonometricFunction, HyperbolicFunction, exp)
                and all(isinstance(t.as_independent(log)[1], log)
                        and not t.as_independent(log)[0].has(symbol)
                        for t in rest.args)):
            # distribute constants over sums of logarithms of
            # transcendental arguments, e.g. sqrt(3)*(log(A) - log(B))/6
            queue.extend(coefficient*t for t in rest.args)
        else:
            others.append(term)
    if len(log_terms) < 2:
        return result
    # group terms whose coefficient ratios are rational
    groups: list[list[tuple[Expr, Expr]]] = []
    for coefficient, argument in log_terms:
        for group in groups:
            ratio = (coefficient/group[0][0]).cancel()
            if ratio.is_Rational:
                group.append((coefficient, argument))
                break
        else:
            groups.append([(coefficient, argument)])
    combined = []
    for group in groups:
        # only collect logarithms of transcendental arguments, where the
        # collected form removes cancelling singularities; sums of
        # logarithms of rational arguments are left in their usual form
        if (len(group) < 2 or not any(
                argument.has(TrigonometricFunction, HyperbolicFunction, exp)
                for _, argument in group)):
            combined.extend(c*log(argument) for c, argument in group)
            continue
        c0 = group[0][0]
        ratios = [(c/c0).cancel() for c, _ in group]
        n = lcm_list([r.q for r in ratios])
        powers = [r*n for r in ratios]
        if any(abs(power) > 4 for power in powers):
            # avoid unwieldy powers inside the collected logarithm
            combined.extend(c*log(argument) for c, argument in group)
            continue
        collected = Mul(*[argument**int(power)
                          for (_, argument), power in zip(group, powers)]).cancel()
        if collected.has(TrigonometricFunction, HyperbolicFunction):
            simplified = collected.trigsimp()
            if simplified.count_ops() <= collected.count_ops():
                collected = simplified
        combined.append(c0/n*log(collected))
    return Add(*(others + combined))


def manualintegrate(f, var, max_depth=None, max_calls=None,
                    soft_max_calls=None):
    """manualintegrate(f, var)

    Explanation
    ===========

    Compute indefinite integral of a single variable using an algorithm that
    resembles what a student would do by hand.

    Unlike :func:`~.integrate`, var can only be a single symbol.

    If ``max_depth`` is given, subproblems nested more than ``max_depth``
    integrals deep are treated as unsolvable; parts of the result may then
    be returned as unevaluated ``Integral``\\ s.

    Examples
    ========

    >>> from sympy import sin, cos, tan, exp, log, integrate
    >>> from sympy.integrals.manualintegrate import manualintegrate
    >>> from sympy.abc import x
    >>> manualintegrate(1 / x, x)
    log(x)
    >>> integrate(1/x)
    log(x)
    >>> manualintegrate(log(x), x)
    x*log(x) - x
    >>> integrate(log(x))
    x*log(x) - x
    >>> manualintegrate(exp(x) / (1 + exp(2 * x)), x)
    atan(exp(x))
    >>> integrate(exp(x) / (1 + exp(2 * x)))
    RootSum(4*w**2 + 1, Lambda(w, w*log(2*w + exp(x))))
    >>> manualintegrate(cos(x)**4 * sin(x), x)
    -cos(x)**5/5
    >>> integrate(cos(x)**4 * sin(x), x)
    -cos(x)**5/5
    >>> manualintegrate(cos(x)**4 * sin(x)**3, x)
    cos(x)**7/7 - cos(x)**5/5
    >>> integrate(cos(x)**4 * sin(x)**3, x)
    cos(x)**7/7 - cos(x)**5/5
    >>> manualintegrate(tan(x), x)
    -log(cos(x))
    >>> integrate(tan(x), x)
    -log(cos(x))

    See Also
    ========

    sympy.integrals.integrals.integrate
    sympy.integrals.integrals.Integral.doit
    sympy.integrals.integrals.Integral
    """
    result = integral_steps(f, var, max_depth=max_depth,
                            max_calls=max_calls,
                            soft_max_calls=soft_max_calls).eval()
    # If we got Piecewise with two parts, put generic first
    if isinstance(result, Piecewise) and len(result.args) == 2:
        cond = result.args[0][1]
        if isinstance(cond, Eq) and result.args[1][1] == True:
            result = result.func(
                (result.args[1][0], Ne(*cond.args)),
                (result.args[0][0], True))
    # Factor terms like erf(x)*sin(x) that may have been expanded
    def _has_erf_trig_mul(expr):
        for sub in expr.find(Mul):
            if sub.has(erf, erfc, erfi) and sub.has(sin, cos, sinh, cosh):
                return True
        return False
    if _has_erf_trig_mul(f) and _has_erf_trig_mul(result):
        result = factor_terms(result)
    # collect sums of logarithms into single logarithms, which are
    # continuous on wider domains (Jeffrey 1993)
    result = _combine_logs(result, var)
    return result
