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
from sympy.polys.polytools import degree, lcm_list, gcd_list, Poly
from sympy.simplify.radsimp import fraction
from sympy.simplify.simplify import simplify
from sympy.simplify.powsimp import powsimp
from sympy.solvers.solvers import solve
from sympy.strategies.core import switch, do_one, null_safe, condition
from sympy.utilities.iterables import iterable
from sympy.utilities.misc import debug

if TYPE_CHECKING:
    from sympy.core.expr import Expr


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
            step = inverse_trig_rule(IntegralInfo(1/s, x), degenerate=False)
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
        step = sqrt_quadratic_rule(IntegralInfo(self.integrand, self.variable), degenerate=False)
        return step.eval()


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
                substep = integral_steps(rewritten, symbol)
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
                return RewriteRule(integrand, symbol, rewritten, integral_steps(rewritten, symbol))
    return _proxy_rewriter

def multiplexer(conditions):
    """Apply the rule that matches the condition, else None"""
    def multiplexer_rl(expr):
        for key, rule in conditions.items():
            if key(expr):
                return rule(expr)
    return multiplexer_rl

def alternatives(*rules):
    """Strategy that makes an AlternativeRule out of multiple possible results."""
    def _alternatives(integral):
        alts = []
        count = 0
        debug("List of Alternative Rules")
        for rule in rules:
            count = count + 1
            debug("Rule {}: {}".format(count, rule))

            result = rule(integral)
            if (result and not isinstance(result, DontKnowRule) and
                result != integral and result not in alts):
                alts.append(result)
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


def powsimp_rule(integral):
    """
    Strategy that simplifies the exponent of a power.
    exp(a*x**2) * exp(b*x) -> exp((a*x**2 + b*x))
    For example, this is useful for the ErfRule.
    """
    integrand, symbol = integral
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    k = Wild('k', exclude=[symbol])

    match = integrand.match(k**(a*symbol**2) * k**(b*symbol))

    if not match:
        return

    simplified = powsimp(integrand, combine='exp')

    if simplified != integrand:
        steps = integral_steps(simplified, symbol)
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
_wilds = []
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
        degenerate_step = sqrt_fractional_linear_rule(IntegralInfo((a + b * symbol) ** exp, symbol))

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
    results = [integral_steps(g, symbol)
              for g in integrand.as_ordered_terms()]
    return None if None in results else AddRule(integrand, symbol, results)


def mul_rule(integral: IntegralInfo):
    integrand, symbol = integral

    # Constant times function case
    coeff, f = integrand.as_independent(symbol)
    if coeff != 1:
        next_step = integral_steps(f, symbol)
        if next_step is not None:
            return ConstantTimesRule(integrand, symbol, coeff, f, next_step)


special_error_functions = (erf, erfc, erfi, fresnelc, fresnels, Ci, Chi, Si, Shi, Ei, li)


def _parts_rule(integrand, symbol) -> tuple[Expr, Expr, Expr, Expr, Rule] | None:
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

    liate_rules = [pull_out_u(*special_error_functions), pull_out_u(log),
                   pull_out_u(*inverse_trig_functions), pull_out_algebraic,
                   pull_out_u(sin, cos), pull_out_u(sinh, cosh),
                   pull_out_u(exp)]


    dummy = Dummy("temporary")
    # we can integrate log(x) and atan(x) by setting dv = 1
    if isinstance(integrand, (log, *inverse_trig_functions)):
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
            # Don't trade one logarithm for another
            if isinstance(u, log):
                rec_dv = 1/dv
                if (rec_dv.is_polynomial(symbol) and
                    degree(rec_dv, symbol) == 1):
                    return None

            # Can integrate a polynomial times OrthogonalPolynomial
            if rule == pull_out_algebraic:
                if dv.is_Derivative or dv.has(TrigonometricFunction, HyperbolicFunction) or \
                        isinstance(dv, OrthogonalPolynomial):
                    v_step = integral_steps(dv, symbol)
                    if v_step.contains_dont_know():
                        return None
                    else:
                        du = u.diff(symbol)
                        v = v_step.eval()
                        return u, dv, v, du, v_step

            # make sure dv is amenable to integration
            accept = False
            if index < 2:  # log and inverse trig are usually worth trying
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
                v_step = integral_steps(simplify(dv), symbol)
                if not v_step.contains_dont_know():
                    v = v_step.eval()
                    return u, dv, v, du, v_step
    return None


def parts_rule(integral):
    integrand, symbol = integral
    constant, integrand = integrand.as_coeff_Mul()

    result = _parts_rule(integrand, symbol)

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
            if _parts_u_cache[cachekey] > 2:
                return
            _parts_u_cache[cachekey] += 1

        # Try cyclic integration by parts a few times
        for _ in range(4):
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
            result = _parts_rule(next_integrand, symbol)

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
            return PartsRule(integrand, symbol, u, dv, v_step, make_second_step(steps[1:], v * du))
        return integral_steps(integrand, symbol)

    if steps:
        u, dv, v, du, v_step = steps[0]
        rule = PartsRule(integrand, symbol, u, dv, v_step, make_second_step(steps[1:], v * du))
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

    return RewriteRule(integrand, symbol, rewritten, integral_steps(rewritten, symbol))

def trig_product_rule(integral: IntegralInfo):
    integrand, symbol = integral
    if integrand == sec(symbol) * tan(symbol):
        return SecTanRule(integrand, symbol)
    if integrand == csc(symbol) * cot(symbol):
        return CscCotRule(integrand, symbol)


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
        steps = integral_steps(rewritten, symbol)
        return RewriteRule(integrand, symbol, rewritten, steps)


def quadratic_denom_rule(integral):
    integrand, symbol = integral
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    c = Wild('c', exclude=[symbol])

    match = integrand.match(a / (b * symbol ** 2 + c))

    if match:
        a, b, c = match[a], match[b], match[c]
        general_rule = ArctanRule(integrand, symbol, a, b, c)
        if b.is_extended_real and c.is_extended_real:
            positive_cond = c/b > 0
            if positive_cond is S.true:
                return general_rule
            coeff = a/(2*sqrt(-c)*sqrt(b))
            constant = sqrt(-c/b)
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
                return negative_step
            return PiecewiseRule(integrand, symbol, [(general_rule, positive_cond), (negative_step, S.true)])

        power = PowerRule(integrand, symbol, symbol, -2)
        if b != 1:
            power = ConstantTimesRule(integrand, symbol, 1/b, symbol**-2, power)

        return PiecewiseRule(integrand, symbol, [(general_rule, Ne(c, 0)), (power, True)])

    d = Wild('d', exclude=[symbol])
    match2 = integrand.match(a / (b * symbol ** 2 + c * symbol + d))
    if match2:
        b, c =  match2[b], match2[c]
        if b.is_zero:
            return
        u = Dummy('u')
        u_func = symbol + c/(2*b)
        integrand2 = integrand.subs(symbol, u - c / (2*b))
        next_step = integral_steps(integrand2, u)
        if next_step:
            return URule(integrand2, symbol, u, u_func, next_step)
        else:
            return
    e = Wild('e', exclude=[symbol])
    match3 = integrand.match((a* symbol + b) / (c * symbol ** 2 + d * symbol + e))
    if match3:
        a, b, c, d, e = match3[a], match3[b], match3[c], match3[d], match3[e]
        if c.is_zero:
            return
        denominator = c * symbol**2 + d * symbol + e
        const =  a/(2*c)
        numer1 =  (2*c*symbol+d)
        numer2 = - const*d + b
        u = Dummy('u')
        step1 = URule(integrand, symbol,
                      u, denominator, integral_steps(u**(-1), u))
        if const != 1:
            step1 = ConstantTimesRule(const*numer1/denominator, symbol,
                                      const, numer1/denominator, step1)
        if numer2.is_zero:
            return step1
        step2 = integral_steps(numer2/denominator, symbol)
        substeps = AddRule(integrand, symbol, [step1, step2])
        rewriten = const*numer1/denominator+numer2/denominator
        return RewriteRule(integrand, symbol, rewriten, substeps)

    return


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
    bases, qs, ratios = [], [], []
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
            bases.append(base)
            ratios.append(S.One)
            qs.append(exp_.q)
        else:
            K = (base / base0).cancel()
            if K.has(x): # cannot substitute both sqrt(x) and sqrt(x + 1)
                return None
            bases.append(base)
            ratios.append(K)
            qs.append(exp_.q)
    if base0 is None and not constant_bases_subs:
        return None
    if constant_bases_subs:
        integrand = integrand.subs(constant_bases_subs)
    if base0 is None:
        substep = integral_steps(integrand, x)
        if not substep.contains_dont_know():
            return RewriteRule(integral.integrand, x, integrand, substep)
        return None
    q0: Integer = lcm_list(qs)
    u = Dummy("u")
    u_x = base0**(S.One/q0)
    u_pow = u**q0
    x_u = (b0 - d0*u_pow)/(c0*u_pow - a0)
    dx_u = (q0*(a0*d0 - b0*c0)*u**(q0 - 1))/(c0*u_pow - a0)**2
    subs_dict = {}
    for base_i, ratio_i, q_i in zip(bases, ratios, qs):
        subs_dict[base_i**(S.One/q_i)] = (ratio_i)**(S.One/q_i) * u**(q0/q_i)
    substituted = integrand.subs(subs_dict).subs(x, x_u) * dx_u
    substep = integral_steps(substituted, u)
    if not substep.contains_dont_know():
        step: Rule = URule(integrand, x, u, u_x, substep)
        # in these cases, deteminant would be 0 only if both c0 and d0 were 0 (null denom), no need of Piecewise
        if (c0.is_zero and (a0.is_zero is False)) or (d0.is_zero and (b0.is_zero is False)):
            if constant_bases_subs:
                return RewriteRule(integral.integrand, x, integrand, step)
            else:
                return step
        generic_cond = Ne(a0*d0 - b0*c0, 0)
        if generic_cond is not S.true:
            pieces: list[tuple[Rule, Boolean]] = [(step, generic_cond)]
            cond_c0 = Ne(c0, 0)
            if cond_c0 is not S.false:
                const_val = a0 / c0
                subs_a = {base_i: ratio_i * const_val for base_i, ratio_i in zip(bases, ratios)}
                simplified_a = integrand.subs(subs_a)
                degenerate_step_a = integral_steps(simplified_a, x)
                pieces.append((degenerate_step_a, cond_c0))
            if cond_c0 is not S.true:
                const_val = b0 / d0
                subs_b = {base_i: ratio_i * const_val for base_i, ratio_i in zip(bases, ratios)}
                simplified_b = integrand.subs(subs_b)
                simplified_b.subs({a0: 0, c0: 0}) # if det=0 and c0=0 and d0!=0, a0=0
                degenerate_step_b = integral_steps(simplified_b, x)
                pieces.append((degenerate_step_b, S.true))
            step = PiecewiseRule(integrand, x, pieces)
        if constant_bases_subs:
            return RewriteRule(integral.integrand, x, integrand, step)
        else:
            return step
    return None


def sqrt_quadratic_rule(integral: IntegralInfo, degenerate=True):
    integrand, x = integral
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    c = Wild('c', exclude=[x, 0])
    f = Wild('f')
    n = Wild('n', properties=[lambda n: n.is_Integer and n.is_odd])
    match = integrand.match(f*sqrt(a+b*x+c*x**2)**n)
    if not match:
        return
    a, b, c, f, n = match[a], match[b], match[c], match[f], match[n]
    f_poly = f.as_poly(x)
    if f_poly is None:
        return

    generic_cond = Ne(c, 0)
    if not degenerate or generic_cond is S.true:
        degenerate_step = None
    elif b.is_zero:
        degenerate_step = integral_steps(f*sqrt(a)**n, x)
    else:
        degenerate_step = sqrt_fractional_linear_rule(IntegralInfo(f*sqrt(a+b*x)**n, x))

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
                constant_step = inverse_trig_rule(IntegralInfo(1/denom, x), degenerate=False)
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

    def sqrt_quadratic_reduction_rule(integrand: Expr, n: int):
        # Implementation of Gradshteyn & Ryzhik 2.263.3
        k = (-n - 1) // 2
        delta = 4*a*c - b**2
        R = c*x**2 + b*x + a

        term_denom = (2*k - 1) * delta * (R**(S(2*k - 1)/2))
        constant_term = f*2*(2*c*x+b) / term_denom
        coeff = (8*c*(k-1))/((2*k-1) * delta)
        expr = f * R**(S(1)/2 - k)

        rewrite_expr = Derivative(constant_term, x) + coeff * expr
        derive_expr = Derivative(constant_term, x)
        derive_step = integral_steps(derive_expr, x)

        if coeff == 0:
            substep = derive_step
        else:
            next_step = integral_steps(expr, x)
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

        return RewriteRule(integrand, x, rewrite_expr, substep)

    if n > 0:  # rewrite poly * sqrt(s)**(2*k-1) to poly*s**k / sqrt(s)
        numer_poly = f_poly * (a+b*x+c*x**2)**((n+1)/2)
        rewritten = numer_poly.as_expr()/sqrt(a+b*x+c*x**2)
        substep = sqrt_quadratic_denom_rule(numer_poly, rewritten)
        generic_step = RewriteRule(integrand, x, rewritten, substep)
    elif n == -1:
        generic_step = sqrt_quadratic_denom_rule(f_poly, integrand)
    elif f_poly.degree() == 0:
        # The numerator must be a const, the formula assumes this
        generic_step = sqrt_quadratic_reduction_rule(integrand, n)
    else:
        # Handle non-constant numerators (eg. x / R**(-3/2))
        # This requires splitting the integral as A*(2ax+b) + B form
        return None
    return _add_degenerate_step(generic_cond, generic_step, degenerate_step)


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
            return RewriteRule(integrand, symbol, rewritten,
                   URule(rewritten, symbol, u, cosh(symbol), ReciprocalRule(1/u, u, u)))
        if integrand.func == coth:
            rewritten = cosh(symbol)/sinh(symbol)
            return RewriteRule(integrand, symbol, rewritten,
                   URule(rewritten, symbol, u, sinh(symbol), ReciprocalRule(1/u, u, u)))
        else:
            rewritten = integrand.rewrite(tanh)
            if integrand.func == sech:
                return RewriteRule(integrand, symbol, rewritten,
                       URule(rewritten, symbol, u, tanh(symbol/2),
                       ArctanRule(2/(u**2 + 1), u, S(2), S.One, S.One)))
            if integrand.func == csch:
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
        return func(*args)
    return uncurry_rl

def trig_rewriter(rewrite):
    def trig_rewriter_rl(args):
        a, b, m, n, integrand, symbol = args
        rewritten = rewrite(a, b, m, n, integrand, symbol)
        if rewritten != integrand:
            return RewriteRule(integrand, symbol, rewritten, integral_steps(rewritten, symbol))
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
        return integral_steps(integrand * sin_double, symbol)

def trig_powers_products_rule(integral):
    return do_one(null_safe(trig_sincos_rule),
                  null_safe(trig_tansec_rule),
                  null_safe(trig_cotcsc_rule),
                  null_safe(trig_sindouble_rule))(integral)

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

                substep = integral_steps(replaced, theta)
                if not substep.contains_dont_know():
                    return TrigSubstitutionRule(integrand, symbol,
                        theta, x_func, replaced, substep, restriction)

def heaviside_rule(integral):
    integrand, symbol = integral
    pattern, m, b, g = heaviside_pattern(symbol)
    match = integrand.match(pattern)
    if match and 0 != match[g]:
        # f = Heaviside(m*x + b)*g
        substep = integral_steps(match[g], symbol)
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
        debug("List of Substitution Rules")
        ways = []
        for u_func, c, substituted in substitutions:
            subrule = integral_steps(substituted, u_var)
            count = count + 1
            debug("Rule {}: {}".format(count, subrule))

            if subrule.contains_dont_know():
                continue

            if simplify(c - 1) != 0:
                _, denom = c.as_numer_denom()
                if subrule:
                    subrule = ConstantTimesRule(c * substituted, u_var, c, substituted, subrule)

                if denom.free_symbols:
                    piecewise = []
                    could_be_zero = []

                    if isinstance(denom, Mul):
                        could_be_zero = denom.args
                    else:
                        could_be_zero.append(denom)

                    for expr in could_be_zero:
                        if not fuzzy_not(expr.is_zero):
                            substep = integral_steps(manual_subs(integrand, expr, 0), symbol)

                            if substep:
                                piecewise.append((
                                    substep,
                                    Eq(expr, 0)
                                ))
                    piecewise.append((subrule, True))
                    subrule = PiecewiseRule(substituted, symbol, piecewise)

            ways.append(URule(integrand, symbol, u_var, u_func, subrule))

        if len(ways) > 1:
            return AlternativeRule(integrand, symbol, ways)
        elif ways:
            return ways[0]


partial_fractions_rule = rewriter(
    lambda integrand, symbol: integrand.is_rational_function(),
    lambda integrand, symbol: integrand.apart(symbol))

cancel_rule = rewriter(
    # lambda integrand, symbol: integrand.is_algebraic_expr(),
    # lambda integrand, symbol: isinstance(integrand, Mul),
    lambda integrand, symbol: True,
    lambda integrand, symbol: integrand.cancel())

distribute_expand_rule = rewriter(
    lambda integrand, symbol: (
        isinstance(integrand, (Pow, Mul)) or all(arg.is_Pow or arg.is_polynomial(symbol) for arg in integrand.args)),
    lambda integrand, symbol: integrand.expand())

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
        return RewriteRule(integrand, symbol, rewritten, integral_steps(rewritten, symbol))

def fallback_rule(integral):
    return DontKnowRule(*integral)


# Cache is used to break cyclic integrals.
# Need to use the same dummy variable in cached expressions for them to match.
# Also record "u" of integration by parts, to avoid infinite repetition.
_integral_cache: dict[Expr, Expr | None] = {}
_parts_u_cache: dict[Expr, int] = defaultdict(int)
_cache_dummy = Dummy("z")

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

    Returns
    =======

    rule : Rule
        The first step; most rules have substeps that must also be
        considered. These substeps can be evaluated using ``manualintegrate``
        to obtain a result.

    """
    cachekey = integrand.xreplace({symbol: _cache_dummy})
    if cachekey in _integral_cache:
        if _integral_cache[cachekey] is None:
            # Stop this attempt, because it leads around in a loop
            return DontKnowRule(integrand, symbol)
        else:
            # TODO: This is for future development, as currently
            # _integral_cache gets no values other than None
            return (_integral_cache[cachekey].xreplace(_cache_dummy, symbol),
                symbol)
    else:
        _integral_cache[cachekey] = None

    integral = IntegralInfo(integrand, symbol)

    def key(integral):
        integrand = integral.integrand

        if symbol not in integrand.free_symbols:
            return Number
        for cls in (Symbol, TrigonometricFunction, OrthogonalPolynomial):
            if isinstance(integrand, cls):
                return cls
        return type(integrand)

    def integral_is_subclass(*klasses):
        def _integral_is_subclass(integral):
            k = key(integral)
            return k and issubclass(k, klasses)
        return _integral_is_subclass

    result = do_one(
        null_safe(special_function_rule),
        null_safe(switch(key, {
            Pow: do_one(null_safe(power_rule), null_safe(inverse_trig_rule),
                        null_safe(quadratic_denom_rule),
                        null_safe(sqrt_quadratic_rule),
                        null_safe(sqrt_fractional_linear_rule)),
            Symbol: power_rule,
            exp: exp_rule,
            Add: add_rule,
            Mul: do_one(null_safe(mul_rule), null_safe(trig_product_rule),
                        null_safe(heaviside_rule), null_safe(quadratic_denom_rule),
                        null_safe(sqrt_quadratic_rule),
                        null_safe(sqrt_fractional_linear_rule),
                        null_safe(powsimp_rule),
                        null_safe(trig_cmplx_exp_rule)),
            Derivative: derivative_rule,
            TrigonometricFunction: trig_rule,
            Heaviside: heaviside_rule,
            DiracDelta: dirac_delta_rule,
            OrthogonalPolynomial: orthogonal_poly_rule,
            Number: constant_rule
        })),
        do_one(
            null_safe(trig_rule),
            null_safe(hyperbolic_rule),
            null_safe(alternatives(
                rewrites_rule,
                substitution_rule,
                condition(
                    integral_is_subclass(Mul, Pow),
                    partial_fractions_rule),
                condition(
                    integral_is_subclass(Mul, Pow),
                    cancel_rule),
                condition(
                    integral_is_subclass(Mul, log,
                    *inverse_trig_functions),
                    parts_rule),
                condition(
                    integral_is_subclass(Mul, Pow),
                    distribute_expand_rule),
                trig_powers_products_rule,
                trig_expand_rule
            )),
            null_safe(condition(integral_is_subclass(Mul, Pow), nested_pow_rule)),
            null_safe(trig_substitution_rule)
        ),
        fallback_rule)(integral)
    del _integral_cache[cachekey]
    return result


def manualintegrate(f, var):
    """manualintegrate(f, var)

    Explanation
    ===========

    Compute indefinite integral of a single variable using an algorithm that
    resembles what a student would do by hand.

    Unlike :func:`~.integrate`, var can only be a single symbol.

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
    result = integral_steps(f, var).eval()
    # Clear the cache of u-parts
    _parts_u_cache.clear()
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
    return result
