"""Integration method that emulates by-hand techniques.

This module also provides functionality to get the steps used to evaluate a
particular integral, in the ``integral_steps`` function. This will return
a ``Rule`` representing the first step of computing the antiderivative;
its ``subgoals`` (solved via :class:`ManualSolver`) represent further work.

Architecture
============

Integration is modeled as an AND-OR search over a directed acyclic
hypergraph, built lazily while solving:

- Nodes are goals: ``IntegralInfo(integrand, symbol)``.
- Hyperedges are ``Rule`` instances: each has one target ``goal``, a
  rewritten ``result`` expression (which may embed unresolved
  ``Integral(...)`` atoms), and zero or more ``subgoals`` that must be
  solved before ``result`` is the final answer.
- AND is a single ``Rule`` with multiple subgoals. OR is implicit: several
  ``Rule`` candidates proposed for the same goal; solving any one solves
  the goal.

A *proposer* is a function ``(goal, solver) -> Rule | list[Rule] | None``
that proposes candidate rules for a goal without recursing back into the
solver to test subgoal solvability -- that's the search's job. The one
exception is ``PartsRule``: the antiderivative ``v`` of ``dv`` is a real
solver goal, but by-parts's *second* integral (``du*v``) cannot be stated
until ``v``'s value is known, so ``Rule.expand()`` lets a rule derive one
more goal once an earlier one has actually resolved.

``ManualSolver`` drives a depth-first AND-OR search with a solver-side
cycle guard (``open_path``): revisiting a goal already being expanded on
the current search path fails that specific candidate, not the whole goal.
Genuinely cyclic integration patterns (by-parts loops like
``exp(x)*sin(x)``) are handled by a dedicated proposer
(``exp_trig_cyclic_rule``) that computes the closed form directly via
algebra, never by the solver "solving an equation for I".

The integrator can be extended with new heuristics and evaluation
techniques: add a ``Rule`` subclass, then write a proposer function that
accepts an ``IntegralInfo`` (and the solver) and returns a ``Rule``
instance, a list of them, or ``None``; register it in ``PROPOSER_TABLE``
(keyed by the dispatch type from ``key()``) or ``FALLBACK_PROPOSERS``. To
enable simple substitutions, add the match to ``find_substitutions``.

"""

from __future__ import annotations
from typing import NamedTuple, Callable, Sequence, TYPE_CHECKING
from abc import ABC
from collections import defaultdict
from collections.abc import Mapping
from enum import Enum, auto

from sympy.core.add import Add
from sympy.core.cache import cacheit
from sympy.core.containers import Dict, Tuple as STuple
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
from .rationaltools import ratint
from sympy.logic.boolalg import And, Boolean
from sympy.ntheory.factor_ import primefactors
from sympy.polys.polytools import degree, factor_list, lcm_list, gcd_list, Poly
from sympy.simplify.radsimp import fraction
from sympy.simplify.simplify import simplify
from sympy.simplify.powsimp import powsimp
from sympy.solvers.solvers import solve
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


class IntegralInfo(NamedTuple):
    integrand: Expr
    symbol: Symbol


# ---------------------------------------------------------------------------
# Node identity: canonicalization
# ---------------------------------------------------------------------------
#
# Rules introduce fresh Dummy symbols (u-substitution's u, TrigSubstitution's
# theta, ...). Goals identical up to dummy-renaming must collapse to one
# node, or the solver's goal-keyed dicts (and the open_path cycle guard)
# never see the sharing/cycles they need to.
#
# canon() must be a *pure, idempotent* function of an expression's shape:
# the same logical goal reached through two different Python objects has to
# canonicalize to two EQUAL results. SymPy Dummy equality is by identity
# (dummy_index), not by name, so minting a *fresh* Dummy("_d0") on every
# call would make canon() silently non-idempotent -- two calls on the same
# input would produce two unequal outputs, corrupting every goal-keyed dict
# lookup. Reuse a fixed, pre-created pool instead.

_CANON_POOL = [Dummy(f"_d{i}") for i in range(256)]


def canon(expr: Expr) -> Expr:
    """Deterministically rename every free Dummy symbol in expr to a
    canonical sequence (_d0, _d1, ...), assigned in a fixed (preorder)
    traversal order. Non-Dummy symbols are untouched."""
    seen: dict[Dummy, Dummy] = {}

    def walk(e):
        if isinstance(e, Dummy):
            if e not in seen:
                if len(seen) >= len(_CANON_POOL):
                    raise RuntimeError("canon(): dummy pool exhausted")
                seen[e] = _CANON_POOL[len(seen)]
        for a in getattr(e, "args", ()):
            walk(a)

    walk(expr)
    return expr.xreplace(seen) if seen else expr


def canonical_goal(info: IntegralInfo) -> IntegralInfo:
    """Canonicalizes integrand and symbol together (the bound variable may
    itself be a Dummy from a prior substitution), so a shared dummy renames
    consistently across both."""
    combined = canon(STuple(info.integrand, info.symbol))
    return IntegralInfo(combined[0], combined[1])


# ---------------------------------------------------------------------------
# Rule
# ---------------------------------------------------------------------------

class Rule(ABC):
    """A hyperedge: ``goal`` is the IntegralInfo this rule targets,
    ``result`` is the rewritten expression (may embed ``Integral(...)``
    atoms referencing ``subgoals``), and ``subgoals`` are the goals that
    must be solved before ``result`` is the final answer.

    ``subgoals`` is authoritative data supplied by the proposer that built
    the rule -- it is never derived by scanning ``result``.
    """

    __slots__ = ('goal', 'result', 'subgoals')

    def __init__(
        self,
        goal: IntegralInfo,
        result: Expr,
        subgoals: tuple[IntegralInfo, ...] = (),
    ) -> None:
        self.goal = goal
        self.result = result
        self.subgoals = subgoals

    def eval(self, resolved: Mapping[IntegralInfo, Expr]) -> None:
        """Mutates self.result/self.subgoals in place. No return value.

        Default: generic substitution. For each subgoal g, replace the
        literal atom Integral(g.integrand, g.symbol) inside self.result
        with resolved[g], then clear self.subgoals to mark this rule fully
        resolved. Override only for rules whose composition isn't plain
        substitution (USubstitutionRule, TrigSubstitutionRule,
        HeavisideRule, PartsRule)."""
        for g in self.subgoals:
            atom = Integral(g.integrand, g.symbol)
            self.result = self.result.xreplace({atom: resolved[g]})
        self.subgoals = ()

    def expand(self, resolved: Mapping[IntegralInfo, Expr]) -> IntegralInfo | None:
        """Return one more goal to solve before this rule is AND-satisfied,
        derived from whatever has resolved so far, or None once there's
        nothing more needed. Default: no dynamic expansion -- every
        subgoal was already declared upfront. Only PartsRule overrides
        this: its second integral (du*v) needs v's concrete value before
        it can even be stated."""
        return None

    def contains_dont_know(self) -> bool:
        """True if this (already-resolved) rule's value still contains an
        unevaluated Integral -- i.e. some subgoal, however deep, could not
        be solved.

        This has to be a structural check on `result`, not a walk of a
        Rule tree: once `_satisfy()` resolves a rule's subgoals, their Rule
        objects aren't kept around (only the flattened `.result` Expr,
        possibly with embedded `Integral(...)` atoms from deeper
        failures). Only call this after the rule has actually been through
        ManualSolver._satisfy (i.e. on a `winning_rule`/`solver.solve()`
        return value) -- `.result` is a placeholder before that.
        """
        return self.result.has(Integral)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self.goal == other.goal and self.result == other.result
                and self.subgoals == other.subgoals)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(goal={self.goal!r})"


class DontKnowRule(Rule):
    """Placeholder for a goal the solver could not resolve."""

    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, Integral(goal.integrand, goal.symbol))


class GoalStatus(Enum):
    OPEN = auto()
    SOLVED = auto()
    FAILED = auto()


# ---------------------------------------------------------------------------
# Rule subclasses: zero-subgoal closed forms
# ---------------------------------------------------------------------------

class ConstantRule(Rule):
    """integrate(a, x)  ->  a*x"""

    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, goal.integrand * goal.symbol)


class PowerRule(Rule):
    """integrate(x**a, x)"""

    __slots__ = ("base", "exp")

    def __init__(self, goal: IntegralInfo, base: Expr, exp: Expr) -> None:
        result = Piecewise(
            ((base**(exp + 1))/(exp + 1), Ne(exp, -1)),
            (log(base), True),
        )
        super().__init__(goal, result)
        self.base = base
        self.exp = exp


class NestedPowRule(Rule):
    """integrate((x**a)**b, x)"""

    __slots__ = ("base", "exp")

    def __init__(self, goal: IntegralInfo, base: Expr, exp: Expr) -> None:
        m = base * goal.integrand
        result = Piecewise((m / (exp + 1), Ne(exp, -1)),
                            (m * log(base), True))
        super().__init__(goal, result)
        self.base = base
        self.exp = exp


class TrigRule(Rule, ABC):
    __slots__ = ()


class SinRule(TrigRule):
    """integrate(sin(x), x) -> -cos(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, -cos(goal.symbol))


class CosRule(TrigRule):
    """integrate(cos(x), x) -> sin(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, sin(goal.symbol))


class SecTanRule(TrigRule):
    """integrate(sec(x)*tan(x), x) -> sec(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, sec(goal.symbol))


class CscCotRule(TrigRule):
    """integrate(csc(x)*cot(x), x) -> -csc(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, -csc(goal.symbol))


class Sec2Rule(TrigRule):
    """integrate(sec(x)**2, x) -> tan(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, tan(goal.symbol))


class Csc2Rule(TrigRule):
    """integrate(csc(x)**2, x) -> -cot(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, -cot(goal.symbol))


class HyperbolicRule(Rule, ABC):
    __slots__ = ()


class SinhRule(HyperbolicRule):
    """integrate(sinh(x), x) -> cosh(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, cosh(goal.symbol))


class CoshRule(HyperbolicRule):
    """integrate(cosh(x), x) -> sinh(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, sinh(goal.symbol))


class ExpRule(Rule):
    """integrate(a**x, x) -> a**x/ln(a)"""

    __slots__ = ("base", "exp")

    def __init__(self, goal: IntegralInfo, base: Expr, exp: Expr) -> None:
        super().__init__(goal, goal.integrand / log(base))
        self.base = base
        self.exp = exp


class ReciprocalRule(Rule):
    """integrate(1/x, x) -> ln(x)"""

    __slots__ = ("base",)

    def __init__(self, goal: IntegralInfo, base: Expr) -> None:
        super().__init__(goal, log(base))
        self.base = base


class ArcsinRule(Rule):
    """integrate(1/sqrt(1-x**2), x) -> asin(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, asin(goal.symbol))


class ArcsinhRule(Rule):
    """integrate(1/sqrt(1+x**2), x) -> asinh(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, asinh(goal.symbol))


class ReciprocalSqrtQuadraticRule(Rule):
    """integrate(1/sqrt(a+b*x+c*x**2), x) -> log(2*sqrt(c)*sqrt(a+b*x+c*x**2)+b+2*c*x)/sqrt(c)"""

    __slots__ = ("a", "b", "c")

    def __init__(self, goal: IntegralInfo, a: Expr, b: Expr, c: Expr) -> None:
        x = goal.symbol
        result = log(2*sqrt(c)*sqrt(a + b*x + c*x**2) + b + 2*c*x)/sqrt(c)
        super().__init__(goal, result)
        self.a, self.b, self.c = a, b, c


class ArctanRule(Rule):
    """integrate(a/(b*x**2+c), x) -> a/b / sqrt(c/b) * atan(x/sqrt(c/b))"""

    __slots__ = ("a", "b", "c")

    def __init__(self, goal: IntegralInfo, a: Expr, b: Expr, c: Expr) -> None:
        x = goal.symbol
        result = a/b / sqrt(c/b) * atan(x/sqrt(c/b))
        super().__init__(goal, result)
        self.a, self.b, self.c = a, b, c


class RatintRule(Rule):
    """Integrate a rational function using ``ratint`` as a fallback."""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        super().__init__(goal, ratint(goal.integrand, goal.symbol))


class OrthogonalPolyRule(Rule, ABC):
    __slots__ = ("n",)


class JacobiRule(OrthogonalPolyRule):
    __slots__ = ("a", "b")

    def __init__(self, goal: IntegralInfo, n: Expr, a: Expr, b: Expr) -> None:
        x = goal.symbol
        result = Piecewise(
            (2*jacobi(n + 1, a - 1, b - 1, x)/(n + a + b), Ne(n + a + b, 0)),
            (x, Eq(n, 0)),
            ((a + b + 2)*x**2/4 + (a - b)*x/2, Eq(n, 1)))
        super().__init__(goal, result)
        self.n, self.a, self.b = n, a, b


class GegenbauerRule(OrthogonalPolyRule):
    __slots__ = ("a",)

    def __init__(self, goal: IntegralInfo, n: Expr, a: Expr) -> None:
        x = goal.symbol
        result = Piecewise(
            (gegenbauer(n + 1, a - 1, x)/(2*(a - 1)), Ne(a, 1)),
            (chebyshevt(n + 1, x)/(n + 1), Ne(n, -1)),
            (S.Zero, True))
        super().__init__(goal, result)
        self.n, self.a = n, a


class ChebyshevTRule(OrthogonalPolyRule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, n: Expr) -> None:
        x = goal.symbol
        result = Piecewise(
            ((chebyshevt(n + 1, x)/(n + 1) -
              chebyshevt(n - 1, x)/(n - 1))/2, Ne(Abs(n), 1)),
            (x**2/2, True))
        super().__init__(goal, result)
        self.n = n


class ChebyshevURule(OrthogonalPolyRule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, n: Expr) -> None:
        x = goal.symbol
        result = Piecewise(
            (chebyshevt(n + 1, x)/(n + 1), Ne(n, -1)),
            (S.Zero, True))
        super().__init__(goal, result)
        self.n = n


class LegendreRule(OrthogonalPolyRule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, n: Expr) -> None:
        x = goal.symbol
        result = (legendre(n + 1, x) - legendre(n - 1, x))/(2*n + 1)
        super().__init__(goal, result)
        self.n = n


class HermiteRule(OrthogonalPolyRule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, n: Expr) -> None:
        x = goal.symbol
        result = hermite(n + 1, x)/(2*(n + 1))
        super().__init__(goal, result)
        self.n = n


class LaguerreRule(OrthogonalPolyRule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, n: Expr) -> None:
        x = goal.symbol
        result = laguerre(n, x) - laguerre(n + 1, x)
        super().__init__(goal, result)
        self.n = n


class AssocLaguerreRule(OrthogonalPolyRule):
    __slots__ = ("a",)

    def __init__(self, goal: IntegralInfo, n: Expr, a: Expr) -> None:
        result = -assoc_laguerre(n + 1, a - 1, goal.symbol)
        super().__init__(goal, result)
        self.n, self.a = n, a


class IRule(Rule, ABC):
    __slots__ = ("a", "b")

    def __init__(self, goal: IntegralInfo, a: Expr, b: Expr) -> None:
        super().__init__(goal, self._closed_form(goal, a, b))
        self.a, self.b = a, b

    @staticmethod
    def _closed_form(goal: IntegralInfo, a: Expr, b: Expr) -> Expr:
        raise NotImplementedError


class CiRule(IRule):
    __slots__ = ()

    @staticmethod
    def _closed_form(goal, a, b):
        x = goal.symbol
        return cos(b)*Ci(a*x) - sin(b)*Si(a*x)


class ChiRule(IRule):
    __slots__ = ()

    @staticmethod
    def _closed_form(goal, a, b):
        x = goal.symbol
        return cosh(b)*Chi(a*x) + sinh(b)*Shi(a*x)


class EiRule(IRule):
    __slots__ = ()

    @staticmethod
    def _closed_form(goal, a, b):
        x = goal.symbol
        return exp(b)*Ei(a*x)


class SiRule(IRule):
    __slots__ = ()

    @staticmethod
    def _closed_form(goal, a, b):
        x = goal.symbol
        return sin(b)*Ci(a*x) + cos(b)*Si(a*x)


class ShiRule(IRule):
    __slots__ = ()

    @staticmethod
    def _closed_form(goal, a, b):
        x = goal.symbol
        return sinh(b)*Chi(a*x) + cosh(b)*Shi(a*x)


class LiRule(IRule):
    __slots__ = ()

    @staticmethod
    def _closed_form(goal, a, b):
        return li(a*goal.symbol + b)/a


class ErfRule(Rule):
    __slots__ = ("a", "b", "c")

    def __init__(self, goal: IntegralInfo, a: Expr, b: Expr, c: Expr) -> None:
        x = goal.symbol
        if a.is_extended_real:
            result = Piecewise(
                (sqrt(S.Pi)/sqrt(-a)/2 * exp(c - b**2/(4*a)) *
                    erf((-2*a*x - b)/(2*sqrt(-a))), a < 0),
                (sqrt(S.Pi)/sqrt(a)/2 * exp(c - b**2/(4*a)) *
                    erfi((2*a*x + b)/(2*sqrt(a))), True))
        else:
            result = sqrt(S.Pi)/sqrt(a)/2 * exp(c - b**2/(4*a)) * \
                    erfi((2*a*x + b)/(2*sqrt(a)))
        super().__init__(goal, result)
        self.a, self.b, self.c = a, b, c


class FresnelCRule(Rule):
    __slots__ = ("a", "b", "c")

    def __init__(self, goal: IntegralInfo, a: Expr, b: Expr, c: Expr) -> None:
        x = goal.symbol
        result = sqrt(S.Pi)/sqrt(2*a) * (
            cos(b**2/(4*a) - c)*fresnelc((2*a*x + b)/sqrt(2*a*S.Pi)) +
            sin(b**2/(4*a) - c)*fresnels((2*a*x + b)/sqrt(2*a*S.Pi)))
        super().__init__(goal, result)
        self.a, self.b, self.c = a, b, c


class FresnelSRule(Rule):
    __slots__ = ("a", "b", "c")

    def __init__(self, goal: IntegralInfo, a: Expr, b: Expr, c: Expr) -> None:
        x = goal.symbol
        result = sqrt(S.Pi)/sqrt(2*a) * (
            cos(b**2/(4*a) - c)*fresnels((2*a*x + b)/sqrt(2*a*S.Pi)) -
            sin(b**2/(4*a) - c)*fresnelc((2*a*x + b)/sqrt(2*a*S.Pi)))
        super().__init__(goal, result)
        self.a, self.b, self.c = a, b, c


class PolylogRule(Rule):
    __slots__ = ("a", "b")

    def __init__(self, goal: IntegralInfo, a: Expr, b: Expr) -> None:
        super().__init__(goal, polylog(b + 1, a * goal.symbol))
        self.a, self.b = a, b


class UpperGammaRule(Rule):
    __slots__ = ("a", "e")

    def __init__(self, goal: IntegralInfo, a: Expr, e: Expr) -> None:
        x = goal.symbol
        result = x**e * (-a*x)**(-e) * uppergamma(e + 1, -a*x)/a
        super().__init__(goal, result)
        self.a, self.e = a, e


class EllipticFRule(Rule):
    __slots__ = ("a", "d")

    def __init__(self, goal: IntegralInfo, a: Expr, d: Expr) -> None:
        result = elliptic_f(goal.symbol, d/a)/sqrt(a)
        super().__init__(goal, result)
        self.a, self.d = a, d


class EllipticERule(Rule):
    __slots__ = ("a", "d")

    def __init__(self, goal: IntegralInfo, a: Expr, d: Expr) -> None:
        result = elliptic_e(goal.symbol, d/a)*sqrt(a)
        super().__init__(goal, result)
        self.a, self.d = a, d


class DiracDeltaRule(Rule):
    __slots__ = ("n", "a", "b")

    def __init__(self, goal: IntegralInfo, n: Expr, a: Expr, b: Expr) -> None:
        x = goal.symbol
        if n == 0:
            result = Heaviside(a + b*x)/b
        else:
            result = DiracDelta(a + b*x, n - 1)/b
        super().__init__(goal, result)
        self.n, self.a, self.b = n, a, b


class DerivativeRule(Rule):
    """integrate(f'(x), x) -> f(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo) -> None:
        integrand = goal.integrand
        assert isinstance(integrand, Derivative)
        variable_count = list(integrand.variable_count)
        for i, (var, count) in enumerate(variable_count):
            if var == goal.symbol:
                variable_count[i] = (var, count - 1)
                break
        result = Derivative(integrand.expr, *variable_count)
        super().__init__(goal, result)


class ExpTrigCyclicRule(Rule):
    """integrate(exp(a*x)*{sin,cos,sinh,cosh}(b*x), x)

    Replaces CyclicPartsRule: the closed form of applying integration by
    parts twice and solving the resulting linear equation for the original
    integral is computed directly via algebra at proposal time (see
    exp_trig_cyclic_rule) -- zero subgoals, no recursion."""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, closed_form: Expr) -> None:
        super().__init__(goal, closed_form)


# ---------------------------------------------------------------------------
# Rule subclasses: subgoals declared upfront
# ---------------------------------------------------------------------------

class ConstantTimesRule(Rule):
    """integrate(a*f(x), x)  ->  a*integrate(f(x), x)"""

    __slots__ = ("constant", "other")

    def __init__(self, goal: IntegralInfo, constant: Expr, other: Expr) -> None:
        subgoal = IntegralInfo(other, goal.symbol)
        result = constant * Integral(other, goal.symbol)
        super().__init__(goal, result, (subgoal,))
        self.constant = constant
        self.other = other


class AddRule(Rule):
    """integrate(f(x) + g(x), x) -> integrate(f(x), x) + integrate(g(x), x)"""

    __slots__ = ("terms",)

    def __init__(self, goal: IntegralInfo, terms: Sequence[Expr]) -> None:
        subgoals = tuple(IntegralInfo(t, goal.symbol) for t in terms)
        result = Add(*(Integral(t, goal.symbol) for t in terms))
        super().__init__(goal, result, subgoals)
        self.terms = terms


class PiecewiseRule(Rule):
    """integrate(Piecewise((f, cond1), (g, cond2), ...), x)

    Each branch whose expression still needs solving is embedded as
    Integral(...) and listed in subgoals; branches that are already a
    concrete Expr (e.g. a degenerate case computed directly, with no
    search needed) are left as-is and never touched by eval()."""

    __slots__ = ("conditions", "branches")

    def __init__(
        self,
        goal: IntegralInfo,
        branches: Sequence[tuple[Expr, bool | Boolean]],
        subgoals: tuple[IntegralInfo, ...] = (),
    ) -> None:
        # subgoals is explicit, authoritative data from the proposer: each
        # branch expression is either already-concrete (degenerate cases
        # computed eagerly, no search needed) or an Integral(...) atom
        # matching one entry of `subgoals`. `branches` is kept around (not
        # just folded into `result`) so degenerate-case combinators can
        # merge/rebuild Piecewise structure without re-parsing result.args.
        result = Piecewise(*branches)
        super().__init__(goal, result, subgoals)
        self.conditions = tuple(cond for _, cond in branches)
        self.branches = tuple(branches)


class RewriteRule(Rule):
    """Rewrite integrand to another (single-subgoal) form that is easier to
    handle."""

    __slots__ = ("rewritten",)

    def __init__(self, goal: IntegralInfo, rewritten: Expr) -> None:
        subgoal = IntegralInfo(rewritten, goal.symbol)
        result = Integral(rewritten, goal.symbol)
        super().__init__(goal, result, (subgoal,))
        self.rewritten = rewritten


class CompleteSquareRule(RewriteRule):
    """Rewrite a+b*x+c*x**2 to a-b**2/(4*c) + c*(x+b/(2*c))**2"""
    __slots__ = ()


class HeavisideRule(Rule):
    """integrate(Heaviside(m*x+b)*g(x), x)
    -> Heaviside(harg) * (result - result.subs(x, ibnd))

    Needs a custom eval(): after the subgoal (integrating g) resolves,
    the antiderivative must be re-substituted at x=ibnd to enforce
    continuity, which is not expressible as plain atom substitution."""

    __slots__ = ("harg", "ibnd", "_subgoal")

    def __init__(self, goal: IntegralInfo, harg: Expr, ibnd: Expr, g: Expr) -> None:
        subgoal = IntegralInfo(g, goal.symbol)
        result = Heaviside(harg) * Integral(g, goal.symbol)
        super().__init__(goal, result, (subgoal,))
        self.harg = harg
        self.ibnd = ibnd
        self._subgoal = subgoal

    def eval(self, resolved: Mapping[IntegralInfo, Expr]) -> None:
        # If we are integrating over x and the integrand has the form
        #       Heaviside(m*x+b)*g(x) == Heaviside(harg)*g(symbol)
        # then there needs to be continuity at -b/m == ibnd,
        # so we subtract the appropriate term.
        base_result = resolved[self._subgoal]
        variable = self.goal.symbol
        self.result = Heaviside(self.harg) * (
            base_result - base_result.subs(variable, self.ibnd))
        self.subgoals = ()


class USubstitutionRule(Rule):
    """integrate(f(g(x))*g'(x), x) -> integrate(f(u), u), u = g(x)"""

    __slots__ = ("u_var", "u_func")

    def __init__(self, goal: IntegralInfo, u_var: Symbol, u_func: Expr, substituted: Expr) -> None:
        subgoal = IntegralInfo(substituted, u_var)
        result = Integral(substituted, u_var)
        super().__init__(goal, result, (subgoal,))
        self.u_var = u_var
        self.u_func = u_func

    def eval(self, resolved: Mapping[IntegralInfo, Expr]) -> None:
        super().eval(resolved)  # generic substitution first, clears subgoals
        if self.u_func.is_Pow:
            base, exp_ = self.u_func.as_base_exp()
            if exp_ == -1:
                # avoid needless -log(1/x) from substitution
                self.result = self.result.subs(log(self.u_var), -log(base))
        self.result = self.result.subs(self.u_var, self.u_func)


class TrigSubstitutionRule(Rule):
    """integrate(f(x), x), x = theta-substituted trig function

    Needs a custom eval(): inverts the substitution via solve(), applies
    trigsimp(), and wraps in a Piecewise gated on `restriction` (which is
    not always a catch-all True branch, unlike most other Piecewise
    results) -- none of that is plain atom substitution."""

    __slots__ = ("theta", "func", "rewritten", "restriction", "_subgoal")

    def __init__(
        self,
        goal: IntegralInfo,
        theta: Expr,
        func: Expr,
        rewritten: Expr,
        restriction: bool | Boolean,
        substep_goal: IntegralInfo,
    ) -> None:
        result = Integral(rewritten, theta)
        super().__init__(goal, result, (substep_goal,))
        self.theta = theta
        self.func = func
        self.rewritten = rewritten
        self.restriction = restriction
        self._subgoal = substep_goal

    def eval(self, resolved: Mapping[IntegralInfo, Expr]) -> None:
        theta, func, x = self.theta, self.func, self.goal.symbol
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
        substep_result = resolved[self._subgoal]
        self.result = Piecewise(
            (substep_result.subs(substitution).trigsimp(), self.restriction)  # type: ignore
        )
        self.subgoals = ()


class PartsRule(Rule):
    """integrate(u(x)*v'(x), x) -> u(x)*v(x) - integrate(u'(x)*v(x), x)

    `v` is declared as a real subgoal (`v_goal`) instead of being computed
    eagerly via a recursive call before this Rule even exists. The second
    subgoal (integrate du*v) genuinely cannot be stated until v has a
    concrete value, so it is not declared upfront: `expand()` is called by
    the solver once v_goal resolves, and only then derives and hands back
    the second goal. `eval()` combines u*v - second from the resolved
    mapping -- no rule needs to reach back into the solver.
    """

    __slots__ = ("u", "dv", "du", "v_goal", "second_goal")

    def __init__(self, goal: IntegralInfo, u: Expr, dv: Expr, du: Expr) -> None:
        v_goal = IntegralInfo(dv, goal.symbol)
        placeholder = Integral(goal.integrand, goal.symbol)
        super().__init__(goal, placeholder, (v_goal,))
        self.u = u
        self.dv = dv
        self.du = du
        self.v_goal = v_goal
        self.second_goal: IntegralInfo | None = None

    def expand(self, resolved: Mapping[IntegralInfo, Expr]) -> IntegralInfo | None:
        if self.second_goal is None:
            v = resolved[self.v_goal]
            # NOT canonicalized here: canonicalization is purely internal
            # bookkeeping for the solver's own memoization/cycle-detection
            # dicts (done inside ManualSolver.solve()); the goal handed to
            # proposers -- and therefore the variable names appearing in
            # the returned value -- must stay exactly as constructed here.
            self.second_goal = IntegralInfo(self.du * v, self.goal.symbol)
            self.subgoals = (self.second_goal,)
            return self.second_goal
        return None

    def eval(self, resolved: Mapping[IntegralInfo, Expr]) -> None:
        v = resolved[self.v_goal]
        second = resolved[self.second_goal]
        self.result = self.u * v - second
        self.subgoals = ()


class SqrtQuadraticDenomRule(Rule):
    """integrate(poly(x)/sqrt(a+b*x+c*x**2), x)"""

    __slots__ = ("a", "b", "c", "coeffs")

    def __init__(
        self,
        goal: IntegralInfo,
        a: Expr,
        b: Expr,
        c: Expr,
        coeffs: list[Expr],
        result: Expr,
    ) -> None:
        super().__init__(goal, result)
        self.a, self.b, self.c, self.coeffs = a, b, c, coeffs


class SqrtQuadraticRule(Rule):
    """integrate(sqrt(a+b*x+c*x**2), x) -- thin identity wrapper kept for
    parity with the original class catalog; sqrt_quadratic_rule proposes
    the real (possibly multi-subgoal) rule directly."""

    __slots__ = ("a", "b", "c")

    def __init__(self, goal: IntegralInfo, a: Expr, b: Expr, c: Expr, result: Expr) -> None:
        super().__init__(goal, result)
        self.a, self.b, self.c = a, b, c


# ---------------------------------------------------------------------------
# Substitution-finding helpers (unchanged: already pure / non-recursive)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# ManualSolver
# ---------------------------------------------------------------------------
#
# Design note (deliberate departure from the design doc, per instruction to
# follow the existing implementation where they conflict): the old
# do_one/switch dispatch commits to the FIRST proposer that returns a Rule
# at all, without checking whether that rule's subgoals ultimately resolve
# -- an unsolvable piece surfaces as a literal embedded Integral(...) atom
# in the final answer (see e.g. manualintegrate(atan(x)*log(x), x), which
# is expected to return a partially-evaluated expression with exactly such
# a term). A handful of proposers (substitution_rule, parts_rule,
# sqrt_fractional_linear_rule, euler_substitution_rule, the rewriter-based
# rules) DO recursively pre-verify a candidate before deciding whether to
# propose it at all, using `solver.solve(...).contains_dont_know()`.
#
# To reproduce this exactly, `solve()` always returns *some* Rule
# (DontKnowRule at worst -- goals never "fail" outright), and dispatch is
# simple first-match: PROPOSER_TABLE[key(goal)] entries, then
# FALLBACK_PROPOSERS, in order; the first proposer to return non-None wins,
# no backtracking across proposers based on downstream success. Proposers
# that need old-style pre-verification call solver.solve() themselves.


def _key(goal: IntegralInfo):
    """Dispatch key for a goal, matching the pre-existing key() logic."""
    integrand = goal.integrand
    if goal.symbol not in integrand.free_symbols:
        return Number
    for cls in (Symbol, TrigonometricFunction, OrthogonalPolynomial):
        if isinstance(integrand, cls):
            return cls
    return type(integrand)


class ManualSolver:
    """Drives the AND-OR search. `proposer_table`/`fallback_proposers` are
    supplied by the caller (see PROPOSER_TABLE/FALLBACK_PROPOSERS below) so
    the dispatch data is fully external to this class."""

    def __init__(self, proposer_table: dict, fallback_proposers=()):
        self.proposer_table = proposer_table
        self.fallback_proposers = tuple(fallback_proposers)

        self.goals: dict[IntegralInfo, GoalStatus] = {}
        self.winning_rule: dict[IntegralInfo, Rule] = {}
        self.resolved: dict[IntegralInfo, Expr] = {}
        self.open_path: list[IntegralInfo] = []

        # was module-level _parts_u_cache: limits how many times a given
        # "u" (in u(x)*v'(x) by-parts) may be reused along a search path.
        self.parts_u_usage: dict[Expr, int] = defaultdict(int)
        self.cache_dummy = Dummy("z")

    def solve(self, goal: IntegralInfo) -> Rule:
        """Always returns a Rule for `goal` (DontKnowRule at worst).

        Memoization (`winning_rule`/`resolved`) is keyed by the RAW goal,
        not the canonical one -- a solved Rule's `.result` is a concrete
        Expr built out of whatever Dummy objects the caller's `goal`
        actually contained (proposers always see the raw goal, never the
        canonicalized one). Two goals that are only *canonically* equal
        (e.g. two distinct `Dummy("u")` instances from unrelated
        u-substitutions, structurally identical after renaming) do NOT
        share a cache entry: reusing goal A's cached result for goal B
        would hand back an Expr mentioning goal A's dummy, which nothing
        downstream would know to substitute -- exactly the class of bug
        this comment is here to prevent reintroducing. Canonical goals are
        used ONLY for the open_path cycle guard below, which is a pure
        membership check with no value reuse, so it's safe there.
        """
        if goal in self.winning_rule:
            return self.winning_rule[goal]
        cgoal = canonical_goal(goal)
        if cgoal in self.open_path:
            # solver-side search cycle: stop this path, don't loop forever.
            return DontKnowRule(goal)

        self.goals[cgoal] = GoalStatus.OPEN
        self.open_path.append(cgoal)
        try:
            rule = self._propose(goal)
            if rule is None:
                rule = DontKnowRule(goal)
            self._satisfy(rule)
            self.goals[cgoal] = (GoalStatus.FAILED if rule.contains_dont_know()
                                  else GoalStatus.SOLVED)
            self.winning_rule[goal] = rule
            self.resolved[goal] = rule.result
            return rule
        finally:
            self.open_path.pop()

    def extract(self, goal: IntegralInfo) -> Expr:
        """Convenience: solve `goal` and return its resolved Expr."""
        return self.solve(goal).result

    def _propose(self, goal: IntegralInfo) -> Rule | None:
        # special_function_rule is tried unconditionally, before anything
        # else, exactly like the original do_one(special_function_rule, ...)
        # (defined later in this module; resolved by name at call time).
        produced = special_function_rule(goal, self)
        if produced is not None:
            return produced
        proposers = self.proposer_table.get(_key(goal), ()) + self.fallback_proposers
        for proposer in proposers:
            debug(f"Trying proposer {proposer.__name__} on {goal.integrand}")
            produced = proposer(goal, self)
            if produced is not None:
                return produced
        return None

    def satisfy(self, rule: Rule) -> Expr:
        """Resolve rule's subgoals (looping expand() for dynamically-growing
        rules like PartsRule), call eval(), and return the resulting value.

        Proposers that build a composite/custom-eval Rule (USubstitutionRule,
        PartsRule, HeavisideRule, TrigSubstitutionRule, RewriteRule,
        ConstantTimesRule, AddRule, a PiecewiseRule with real subgoals, ...)
        purely as a deterministic sub-computation -- i.e. they need the
        Rule's *value* right now, not to hand the Rule back to the solver --
        MUST go through this rather than reading `.result` off a freshly
        constructed Rule directly: `.result` right after construction is
        often just an unevaluated placeholder (e.g. USubstitutionRule's
        `Integral(substituted, u_var)`, before the u_var -> u_func
        back-substitution and log(1/x) fixup that only happen inside its
        eval()).
        """
        self._satisfy(rule)
        return rule.result

    def _satisfy(self, rule: Rule) -> None:
        """Resolve rule's subgoals (looping expand() for dynamically-growing
        rules like PartsRule), then call eval()."""
        resolved_map: dict[IntegralInfo, Expr] = {}
        frontier = rule.subgoals
        while True:
            for g in frontier:
                resolved_map[g] = self.extract(g)
            rule.subgoals = frontier
            next_goal = rule.expand(resolved_map)
            if next_goal is None:
                break
            frontier = (next_goal,)
        rule.eval(resolved_map)


# ---------------------------------------------------------------------------
# Proposers
# ---------------------------------------------------------------------------
#
# Proposer = Callable[[IntegralInfo, ManualSolver], Rule | None]
#
# A proposer takes one goal (and the solver, for the handful of proposers
# that need to pre-verify a candidate the way the original code did) and
# returns a Rule for that goal, or None.

Proposer = Callable[[IntegralInfo, ManualSolver], "Rule | None"]


def constant_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    return ConstantRule(goal)


def power_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    base, expt = integrand.as_base_exp()

    if symbol not in expt.free_symbols and isinstance(base, Symbol):
        if simplify(expt + 1) == 0:
            return ReciprocalRule(goal, base)
        return PowerRule(goal, base, expt)
    elif symbol not in base.free_symbols and isinstance(expt, Symbol):
        rule = ExpRule(goal, base, expt)

        if fuzzy_not(log(base).is_zero):
            return rule
        elif log(base).is_zero:
            return ConstantRule(IntegralInfo(S.One, symbol))

        return PiecewiseRule(goal, [
            (rule.result, Ne(log(base), 0)),
            (ConstantRule(IntegralInfo(S.One, symbol)).result, True),
        ])
    return None


def exp_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if isinstance(integrand.args[0], Symbol):
        return ExpRule(goal, E, integrand.args[0])
    return None


def combine_power_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    """
    Strategy that simplifies the exponent of a power.
    exp(a*x**2) * exp(b*x) -> exp((a*x**2 + b*x))
    For example, this is useful for the ErfRule.
    """
    integrand, symbol = goal
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    k = Wild('k', exclude=[symbol])
    rest = Wild('rest')

    match = integrand.match(rest * k**(a*symbol**2) * k**(b*symbol))

    if not match:
        return None

    simplified = powsimp(integrand, combine='exp')

    if simplified != integrand:
        return RewriteRule(goal, simplified)
    return None


def orthogonal_poly_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    orthogonal_poly_classes = {
        jacobi: JacobiRule,
        gegenbauer: GegenbauerRule,
        chebyshevt: ChebyshevTRule,
        chebyshevu: ChebyshevURule,
        legendre: LegendreRule,
        hermite: HermiteRule,
        laguerre: LaguerreRule,
        assoc_laguerre: AssocLaguerreRule,
    }
    orthogonal_poly_var_index = {
        jacobi: 3,
        gegenbauer: 2,
        assoc_laguerre: 2,
    }
    integrand, symbol = goal
    for klass in orthogonal_poly_classes:
        if isinstance(integrand, klass):
            var_index = orthogonal_poly_var_index.get(klass, 1)
            if (integrand.args[var_index] is symbol and not
                    any(v.has(symbol) for v in integrand.args[:var_index])):
                return orthogonal_poly_classes[klass](goal, *integrand.args[:var_index])
    return None


_special_function_patterns: list[tuple[type, Expr, Callable | None, tuple]] = []
_wilds: list[Wild] = []
_sf_symbol = Dummy('x')


def special_function_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if not _special_function_patterns:
        a = Wild('a', exclude=[_sf_symbol], properties=[lambda x: not x.is_zero])
        b = Wild('b', exclude=[_sf_symbol])
        c = Wild('c', exclude=[_sf_symbol])
        d = Wild('d', exclude=[_sf_symbol], properties=[lambda x: not x.is_zero])
        e = Wild('e', exclude=[_sf_symbol], properties=[
            lambda x: not (x.is_nonnegative and x.is_integer)])
        _wilds.extend((a, b, c, d, e))
        # patterns consist of a SymPy class, a wildcard expr, an optional
        # condition coded as a lambda (when Wild properties are not enough),
        # followed by an applicable rule
        linear_pattern = a*_sf_symbol + b
        quadratic_pattern = a*_sf_symbol**2 + b*_sf_symbol + c
        _special_function_patterns.extend((
            (Mul, exp(linear_pattern, evaluate=False)/_sf_symbol, None, EiRule),
            (Mul, cos(linear_pattern, evaluate=False)/_sf_symbol, None, CiRule),
            (Mul, cosh(linear_pattern, evaluate=False)/_sf_symbol, None, ChiRule),
            (Mul, sin(linear_pattern, evaluate=False)/_sf_symbol, None, SiRule),
            (Mul, sinh(linear_pattern, evaluate=False)/_sf_symbol, None, ShiRule),
            (Pow, 1/log(linear_pattern, evaluate=False), None, LiRule),
            (exp, exp(quadratic_pattern, evaluate=False), None, ErfRule),
            (sin, sin(quadratic_pattern, evaluate=False), None, FresnelSRule),
            (cos, cos(quadratic_pattern, evaluate=False), None, FresnelCRule),
            (Mul, _sf_symbol**e*exp(a*_sf_symbol, evaluate=False), None, UpperGammaRule),
            (Mul, polylog(b, a*_sf_symbol, evaluate=False)/_sf_symbol, None, PolylogRule),
            (Pow, 1/sqrt(a - d*sin(_sf_symbol, evaluate=False)**2),
                lambda a, d: a != d, EllipticFRule),
            (Pow, sqrt(a - d*sin(_sf_symbol, evaluate=False)**2),
                lambda a, d: a != d, EllipticERule),
        ))
    _integrand = integrand.subs(symbol, _sf_symbol)
    for type_, pattern, constraint, rule in _special_function_patterns:
        if isinstance(_integrand, type_):
            match = _integrand.match(pattern)
            if match:
                wild_vals = tuple(match.get(w) for w in _wilds
                                  if match.get(w) is not None)
                if constraint is None or constraint(*wild_vals):
                    return rule(goal, *wild_vals)
    return None


def _add_degenerate_step(
    goal: IntegralInfo,
    generic_cond: bool | Boolean,
    generic_step: Rule,
    degenerate_step: Rule | None,
) -> Rule:
    if degenerate_step is None:
        return generic_step
    if isinstance(generic_step, PiecewiseRule):
        branches = [(expr, (cond & generic_cond).simplify())
                    for expr, cond in generic_step.branches]
    else:
        branches = [(generic_step.result, generic_cond)]
    subgoals = list(generic_step.subgoals)
    if isinstance(degenerate_step, PiecewiseRule):
        branches += list(degenerate_step.branches)
    else:
        branches.append((degenerate_step.result, S.true))
    subgoals += list(degenerate_step.subgoals)
    return PiecewiseRule(goal, branches, tuple(subgoals))


def nested_pow_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    # nested (c*(a+b*x)**d)**e
    integrand, x = goal

    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    pattern = a_ + b_*x
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
        return None
    if generic_cond is S.true:
        degenerate_step = None
    else:
        # equivalent with subs(b, 0) but no need to find b
        degenerate_step = ConstantRule(IntegralInfo(integrand.subs(x, 0), x))
    generic_step = NestedPowRule(goal, base, exp_)
    return _add_degenerate_step(goal, generic_cond, generic_step, degenerate_step)


def inverse_trig_rule(goal: IntegralInfo, solver: ManualSolver, degenerate: bool = True) -> Rule | None:
    """
    Set degenerate=False on recursive call where coefficient of quadratic term
    is assumed non-zero.
    """
    integrand, symbol = goal
    base, exp = integrand.as_base_exp()
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    c = Wild('c', exclude=[symbol, 0])
    match = base.match(a + b*symbol + c*symbol**2)

    if not match:
        return None

    def make_inverse_trig(RuleClass, a, sign_a, c, sign_c, h) -> Rule:
        # RuleClass is always ArcsinRule or ArcsinhRule here: both are
        # zero-subgoal closed forms depending only on the (substituted)
        # variable, so the whole composition (substitution + constant
        # factor + "complete the square" relabeling) is deterministic --
        # no search needed, just build the closed-form value directly.
        # (The old URule/ConstantTimesRule/CompleteSquareRule wrapper chain
        # existed only to shape the now-unsupported step-by-step display
        # tree; it never changed the computed value.)
        quadratic_base = sqrt(c/a)*(symbol - h)
        constant = 1/sqrt(c)
        value = RuleClass(IntegralInfo(S.Zero, quadratic_base)).result
        if constant != 1:
            value = constant * value
        return Rule(goal, value)

    a, b, c = [match.get(i, S.Zero) for i in (a, b, c)]
    generic_cond = Ne(c, 0)
    if not degenerate or generic_cond is S.true:
        degenerate_step = None
    elif b.is_zero:
        degenerate_step = ConstantRule(IntegralInfo(a ** exp, symbol))
    else:
        degenerate_step = sqrt_fractional_linear_rule(
            IntegralInfo((a + b * symbol) ** exp, symbol), solver)

    if simplify(2*exp + 1) == 0:
        h, k = -b/(2*c), a - b**2/(4*c)  # rewrite base to k + c*(symbol-h)**2
        non_square_cond = Ne(k, 0)
        square_step = None
        if non_square_cond is not S.true:
            square_step = NestedPowRule(
                IntegralInfo(1/sqrt(c*(symbol-h)**2), symbol), symbol-h, S.NegativeOne)
        if non_square_cond is S.false:
            return square_step
        generic_step: Rule = ReciprocalSqrtQuadraticRule(goal, a, b, c)
        step = _add_degenerate_step(goal, non_square_cond, generic_step, square_step)
        if k.is_real and c.is_real:
            # list of (rule, condition)
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
                branches = [(r.result, c) for r, c in rules]
                subgoals = tuple(g for r, _ in rules for g in r.subgoals)
                step = PiecewiseRule(goal, branches, subgoals)
            else:
                step = generic_step
        return _add_degenerate_step(goal, generic_cond, step, degenerate_step)
    if exp == S.Half:
        step = sqrt_quadratic_rule(goal, solver, degenerate=False)
        return _add_degenerate_step(goal, generic_cond, step, degenerate_step)
    return None


def add_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if not isinstance(integrand, Add):
        return None
    return AddRule(goal, integrand.as_ordered_terms())


def mul_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    # Constant times function case
    coeff, f = integrand.as_independent(symbol)
    if coeff != 1:
        return ConstantTimesRule(goal, coeff, f)
    return None


def exp_trig_cyclic_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    """integrate(exp(a*x)*{sin,cos,sinh,cosh}(b*x), x)

    Replaces CyclicPartsRule (design doc section 5.6): pattern-matches and
    writes the closed form obtained (mathematically) by applying
    integration by parts twice and solving the resulting linear equation
    for the original integral -- directly via algebra, at proposal time.
    Zero subgoals, no recursion: this is exactly the kind of "productive
    cycle" the solver's cycle guard is not expected to unwind on its own.
    """
    integrand, symbol = goal
    a = Wild('a', exclude=[symbol, 0])
    b = Wild('b', exclude=[symbol, 0])
    for trig in (sin, cos, sinh, cosh):
        match = integrand.match(exp(a*symbol) * trig(b*symbol))
        if not match:
            continue
        aa, bb = match[a], match[b]
        if trig in (sinh, cosh):
            denom = aa**2 - bb**2
            if denom == 0:
                continue
            if trig is sinh:
                closed = exp(aa*symbol)*(aa*sinh(bb*symbol) - bb*cosh(bb*symbol))/denom
            else:
                closed = exp(aa*symbol)*(aa*cosh(bb*symbol) - bb*sinh(bb*symbol))/denom
        else:
            denom = aa**2 + bb**2
            if trig is sin:
                closed = exp(aa*symbol)*(aa*sin(bb*symbol) - bb*cos(bb*symbol))/denom
            else:
                closed = exp(aa*symbol)*(aa*cos(bb*symbol) + bb*sin(bb*symbol))/denom
        return ExpTrigCyclicRule(goal, closed)
    return None


special_error_functions = (erf, erfc, erfi, fresnelc, fresnels, Ci, Chi, Si, Shi, Ei, li)


def _parts_rule(
    integrand: Expr, symbol: Symbol, solver: ManualSolver | None = None,
) -> tuple[Expr, Expr, Expr] | None:
    """Picks a LIATE-ordered (u, dv) by-parts split, returning (u, dv, du)
    or None if none is structurally viable.

    Ported from the original recursive-verify algorithm: still calls
    solver.solve() to check that a candidate `dv` is actually integrable
    before committing to it (exactly mirroring the old
    integral_steps()/contains_dont_know() pre-check -- this is one of the
    few proposers, per the user, that keeps that pattern rather than the
    design doc's "never recurse" ideal, since removing it would change
    which split gets picked and risk breaking exact-equality tests). It no
    longer eagerly extracts `v`'s value, though: PartsRule's v_goal/
    expand() now own that (see the design review's finding #1 on
    PartsRule's self-contradiction) -- the solver.solve() call here is
    purely a pre-check, memoized, so PartsRule's later re-solve of the same
    goal is free.
    """
    if solver is None:
        solver = ManualSolver(PROPOSER_TABLE, FALLBACK_PROPOSERS)

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
                    u = Mul(*args)  # type: ignore
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
                    v_rule = solver.solve(IntegralInfo(dv, symbol))
                    if v_rule.contains_dont_know():
                        return None
                    du = u.diff(symbol)
                    return u, dv, du

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
                simplified_dv = simplify(dv)
                v_rule = solver.solve(IntegralInfo(simplified_dv, symbol))
                if not v_rule.contains_dont_know():
                    return u, simplified_dv, du
    return None


def parts_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    constant, stripped = integrand.as_coeff_Mul()
    if constant != 1:
        # let the solver re-dispatch on the stripped integrand (which will
        # have constant==1 there); this function then finds the (u, dv)
        # split as normal.
        return ConstantTimesRule(goal, constant, stripped)

    result = _parts_rule(integrand, symbol, solver)
    if result is None:
        return None
    u, dv, du = result

    # Set a limit on the number of times u can be used (was module-level
    # _parts_u_cache, now solver-instance-scoped).
    if isinstance(u, (sin, cos, exp, sinh, cosh)):
        cachekey = u.xreplace({symbol: solver.cache_dummy})
        if solver.parts_u_usage[cachekey] > 2:
            return None
        solver.parts_u_usage[cachekey] += 1

    return PartsRule(goal, u, dv, du)


def trig_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if integrand == sin(symbol):
        return SinRule(goal)
    if integrand == cos(symbol):
        return CosRule(goal)
    if integrand == sec(symbol)**2:
        return Sec2Rule(goal)
    if integrand == csc(symbol)**2:
        return Csc2Rule(goal)

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
        return None

    return RewriteRule(goal, rewritten)


def trig_product_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if integrand == sec(symbol) * tan(symbol):
        return SecTanRule(goal)
    if integrand == csc(symbol) * cot(symbol):
        return CscCotRule(goal)
    return None


def trig_cmplx_exp_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
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
    integrand, symbol = goal

    if not integrand.has(exp) and not integrand.has(sin, cos, sinh, cosh):
        return None

    a = Wild('a', exclude=[symbol, 0])
    b = Wild('b', exclude=[symbol])
    c = Wild('c', exclude=[symbol])
    f = WildFunction('f')
    guassian_pattern = exp(a * symbol**2 + b * symbol + c)
    trigexp_over_x_pattern = f*exp(a * symbol)/symbol
    trigexp_over_x_match = integrand.match(trigexp_over_x_pattern)
    if not (any(term.match(guassian_pattern) for term in integrand.atoms(exp))
            or (trigexp_over_x_match and
                trigexp_over_x_match[f].has(sin, cos, sinh, cosh))):
        return None

    # Replace trig and hyperbolic functions with their exponential forms
    rewritten = integrand.rewrite([sin, cos, sinh, cosh], exp)

    if rewritten != integrand:
        return RewriteRule(goal, rewritten)
    return None


def quadratic_denom_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
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

    def _arctan_match(B, a, c, symbol, degenerate=True) -> Rule:
        # integrates B / a*x**2 + c
        integrand = B / (a*symbol**2 + c)
        inner_goal = IntegralInfo(integrand, symbol)
        pieces: list[tuple[Expr, bool | Boolean]] = []
        subgoals: list[IntegralInfo] = []
        # skips degenerate case if a != 0 or if a = 0 would cause null denominator
        if degenerate and not _if_zero_implies_zero(a, c):
            substituted = integrand.subs(a, 0)
            pieces.append((Integral(substituted, symbol), Eq(a, 0)))
            subgoals.append(IntegralInfo(substituted, symbol))
        if degenerate and not _if_zero_implies_zero(c, a):
            substituted = integrand.subs(c, 0)
            pieces.append((Integral(substituted, symbol), Eq(c, 0)))
            subgoals.append(IntegralInfo(substituted, symbol))
        if a.is_extended_real and c.is_extended_real:
            positive_cond = c/a > 0
            if positive_cond is not S.true:
                coeff = B/(2*sqrt(-c)*sqrt(a))
                constant = sqrt(-c/a)
                # log(x - constant) - log(x + constant), computed directly:
                # reciprocal-of-linear isn't reachable through the general
                # Pow dispatch (power_rule only handles base**expr with a
                # bare Symbol base), so this is built as a closed form here,
                # exactly like the original ReciprocalRule construction.
                value = coeff * (log(symbol - constant) - log(symbol + constant))
                if positive_cond is S.false:
                    pieces.append((value, S.true))
                    return PiecewiseRule(inner_goal, pieces, tuple(subgoals))
                else:
                    pieces.append((value, c / a < 0))
        general_value = ArctanRule(inner_goal, B, a, c).result
        if pieces:
            pieces.append((general_value, S.true))
            return PiecewiseRule(inner_goal, pieces, tuple(subgoals))
        return Rule(inner_goal, general_value)

    def _complete_square(B, a, b, c, n, symbol, degenerate_a=True, degenerate_discriminant=True) -> Rule:
        # integrates B / (a*x**2 + b*x + c)**n
        discriminant = 4*a*c - b**2
        denominator = a*symbol**2 + b*symbol + c
        integrand = B / denominator**n
        inner_goal = IntegralInfo(integrand, symbol)
        pieces: list[tuple[Expr, bool | Boolean]] = []
        subgoals: list[IntegralInfo] = []
        # degenerate flags avoid recalculating Piecewise branches recursively
        if degenerate_a and not _if_zero_implies_zero(a, denominator):
            substituted = integrand.subs(a, 0)
            pieces.append((Integral(substituted, symbol), Eq(a, 0)))
            subgoals.append(IntegralInfo(substituted, symbol))
        if degenerate_discriminant and not _if_zero_implies_zero(discriminant, denominator):
            u = Dummy("u")
            # we divide by a, Piecewise condition above
            u_func = symbol + b/(2*a)
            # (B/a**n) * u**(-2*n) w.r.t. u: constant times a bare-Dummy
            # power, exactly what mul_rule+power_rule would deterministically
            # produce via full dispatch -- computed directly instead.
            pow_value = PowerRule(IntegralInfo(u**(-2*n), u), u, -2*n).result
            usub_value = ((B/a**n) * pow_value).subs(u, u_func)
            if discriminant.is_zero:
                if pieces:
                    pieces.append((usub_value, S.true))
                    return PiecewiseRule(inner_goal, pieces, tuple(subgoals))
                return Rule(inner_goal, usub_value)
            pieces.append((usub_value, Eq(discriminant, 0)))
        if n == 1:
            # base case, B / (a*x**2 + b*x + c), solve by substitution with _arctan_match
            u = Dummy("u")
            u_func = symbol + b/(2*a)
            # degenerate=False: after substitution the integrand becomes
            # B/(a*u**2 + discriminant/(4*a)); the _arctan_match conditions
            # (a != 0 and discriminant != 0) are already established by the
            # Piecewise branches above.
            substep = _arctan_match(B, a, discriminant/(4*a), u, degenerate=False)
            value = substep.result.subs(u, u_func)
            general_step: Rule = Rule(inner_goal, value, substep.subgoals)
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
            derivative_value = DerivativeRule(IntegralInfo(derivative, symbol)).result
            remainder_step = _complete_square(B, a, b, c, n - 1, symbol, degenerate_a=False, degenerate_discriminant=False)
            value = derivative_value + coeff * remainder_step.result
            general_step = Rule(inner_goal, value, remainder_step.subgoals)
        if pieces:
            pieces.append((general_step.result, S.true))
            subgoals += list(general_step.subgoals)
            return PiecewiseRule(inner_goal, pieces, tuple(subgoals))
        return general_step

    def _split_sum(A, B, a, b, c, n, symbol) -> Rule:
        # integrates (A*x + B) / (a*x**2 + b*x + c)**n. Split A*x + B as alpha*q'(x) + beta,
        # then integrate the two terms separately (first by substitution, second with _complete_square)
        denominator = (a*symbol**2 + b*symbol + c)
        integrand = (A*symbol + B) / denominator**n
        inner_goal = IntegralInfo(integrand, symbol)
        pieces: list[tuple[Expr, bool | Boolean]] = []
        subgoals: list[IntegralInfo] = []
        if not _if_zero_implies_zero(a, denominator):
            substituted = integrand.subs(a, 0)
            pieces.append((Integral(substituted, symbol), Eq(a, 0)))
            subgoals.append(IntegralInfo(substituted, symbol))
        # we divide by a, Piecewise condition above
        const = A/(2*a)
        numer1 = (2*a*symbol + b)
        numer2 = -const*b + B
        qprime_part = numer1 / denominator**n
        u = Dummy('u')
        # integrate u**(-n) w.r.t. u directly: a bare Dummy base, so this
        # is exactly what power_rule would deterministically produce via
        # full dispatch -- no need to actually go through the solver.
        pow_value = PowerRule(IntegralInfo(u**(-n), u), u, -n).result
        step1_value = pow_value.subs(u, denominator)
        if const != 1:
            step1_value = const * step1_value
        if numer2.is_zero:
            value = step1_value
            general_step: Rule = Rule(inner_goal, value, ())
        else:
            # since degenerate a condition is already computed, degenerate_a = False
            step2 = _complete_square(numer2, a, b, c, n, symbol, degenerate_a=False)
            value = step1_value + step2.result
            general_step = Rule(inner_goal, value, step2.subgoals)
        if pieces:
            pieces.append((general_step.result, S.true))
            subgoals += list(general_step.subgoals)
            return PiecewiseRule(inner_goal, pieces, tuple(subgoals))
        return general_step

    B = num_poly.nth(0)
    a = den_poly.nth(2)
    b = den_poly.nth(1)
    c = den_poly.nth(0)

    normalized_num = num_poly.as_expr()
    normalized_den = den_poly.as_expr()
    normalized_integrand = normalized_num / normalized_den**n

    if b == 0 and deg_num == 0 and n == 1:
        step = _arctan_match(B, a, c, symbol)
    elif deg_num == 1:
        A = num_poly.nth(1)
        step = _split_sum(A, B, a, b, c, n, symbol)
    else:
        step = _complete_square(B, a, b, c, n, symbol)

    # `step`'s value is already the antiderivative of normalized_integrand
    # (== integrand, just possibly written differently after clearing a
    # common denominator factor) -- reuse it directly instead of discarding
    # it and re-dispatching from scratch.
    return Rule(goal, step.result, step.subgoals)


def sqrt_fractional_linear_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    """
    Substitute common ((a*x + b)/(c*x + d))**(1/n)
    """
    integrand, x = goal
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    c = Wild('c', exclude=[x])
    d = Wild('d', exclude=[x])
    base0 = None
    powers, exps, ratios = [], [], []
    constant_bases_subs = {}
    # use ordered() to ensure a selection of the smallest base0 (eg. first sqrt(x), then cbrt(2x), x chosen)
    for pow_ in ordered(integrand.find((Pow))):  # collect all ((a*x + b)/(c*x + d))**(p/q)
        base, exp_ = pow_.base, pow_.exp
        if exp_.is_Integer or x not in base.free_symbols:  # skip 1/x and sqrt(2)
            continue
        if not exp_.is_Rational:  # exclude x**pi
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
        if det.is_zero:  # constant value as sqrt((5*x + 10)/(2*x +  4))
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
        sub_goal = IntegralInfo(integrand, x)
        substep = solver.solve(sub_goal)
        if not substep.contains_dont_know():
            return RewriteRule(goal, integrand)
        return None
    q0: Integer = lcm_list([exp_i.q for exp_i in exps])
    u = Dummy("u")
    u_x = base0**(S.One/q0)
    u_pow = u**q0
    x_u = (b0 - d0*u_pow)/(c0*u_pow - a0)
    dx_u = (q0*(a0*d0 - b0*c0)*u**(q0 - 1))/(c0*u_pow - a0)**2
    subs_dict = {pow_i: ratio_i * u**(q0*exp_i) for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
    substituted = integrand.xreplace(subs_dict).xreplace({x: x_u}) * dx_u
    u_goal = IntegralInfo(substituted, u)
    substep = solver.solve(u_goal)
    if not substep.contains_dont_know():
        pieces: list[tuple[Expr, bool | Boolean]] = []
        subgoals: list[IntegralInfo] = []
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
                pieces.append((Integral(simplified_a, x), (And(Eq(det, 0), Ne(c0, 0)))))
                subgoals.append(IntegralInfo(simplified_a, x))
            if not c0_implies_d0:
                const_val = b0 / d0
                subs_b = {pow_i: ratio_i * Pow(const_val, exp_i) for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
                simplified_b = integrand.xreplace(subs_b)
                simplified_b = simplified_b.subs({a0: 0, c0: 0})  # if det = 0, c = 0 and d != 0, a must be 0
                pieces.append((Integral(simplified_b, x), (And(Eq(det, 0), Eq(c0, 0)))))
                subgoals.append(IntegralInfo(simplified_b, x))
        usub_value = solver.satisfy(USubstitutionRule(IntegralInfo(integrand, x), u, u_x, substituted))
        if pieces:
            pieces.append((usub_value, S.true))
            return PiecewiseRule(goal, pieces, tuple(subgoals))
        return Rule(goal, usub_value)
    return None


def euler_substitution_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    """
    Substitute common sqrt(a + b*x + c*x**2) terms using Euler substitution.
    """
    integrand, x = goal
    base0 = None
    powers, exps, ratios = [], [], []
    # use ordered() to ensure a selection of the smallest base0 (eg. first sqrt(x**2 + 1), then sqrt(2*x**2 + 2), x**2 + 1 chosen)
    for pow_ in ordered(integrand.find(Pow)):  # collect all (a + b*x + c*x**2)**(p/2)
        base, exp_ = pow_.base, pow_.exp
        if exp_.is_Integer or x not in base.free_symbols:  # skip 1/x and sqrt(2)
            continue
        if not exp_.is_Rational:  # exclude (x**2 + 1)**pi
            return None
        if exp_.q != 2:
            return None
        base_poly = base.as_poly(x)
        if base_poly is None or base_poly.degree() != 2:  # exclude cube polynomial roots and other radicals
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

    pieces: list[tuple[Expr, bool | Boolean]] = []
    subgoals: list[IntegralInfo] = []
    delta = 4*a0*c0 - b0**2
    # substitution not valid for c0 = 0 and delta = 0
    c_zero_cond = Eq(c0, 0)
    delta_zero_cond = Eq(delta, 0)

    def _delta_zero_step():
        shift = x + b0/(2*c0)
        rewritten_base = c0*shift**2
        subs_dict = {pow_i: ratio_i*(rewritten_base)**exp_i for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
        rewritten = integrand.xreplace(subs_dict)
        return rewritten

    def _c_zero_step():
        degenerate_integrand = integrand.subs(c0, 0)
        return degenerate_integrand

    def _general_euler_step():
        s = Dummy("s")
        subs_dict = {pow_i: ratio_i * s**(2*exp_i) for pow_i, exp_i, ratio_i in zip(powers, exps, ratios)}
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
        u_func = sqrt(base0) + sqrt_c0*x
        return u, u_func, substituted

    if delta_zero_cond is S.true:
        rewritten = _delta_zero_step()
        general_rule = solver.solve(IntegralInfo(rewritten, x))
        if general_rule.contains_dont_know():
            return None
        general_value = solver.satisfy(RewriteRule(goal, rewritten))
    else:
        general_data = _general_euler_step()
        if general_data is None:
            return None
        u, u_func, substituted = general_data
        substep = solver.solve(IntegralInfo(substituted, u))
        if substep.contains_dont_know():
            return None
        general_value = solver.satisfy(USubstitutionRule(goal, u, u_func, substituted))

    if c0.is_zero is None:
        c_zero_integrand = _c_zero_step()
        pieces.append((Integral(c_zero_integrand, x), c_zero_cond))
        subgoals.append(IntegralInfo(c_zero_integrand, x))
    if delta.is_zero is None:
        delta_zero_integrand = _delta_zero_step()
        pieces.append((Integral(delta_zero_integrand, x), delta_zero_cond))
        subgoals.append(IntegralInfo(delta_zero_integrand, x))
    if pieces:
        pieces.append((general_value, S.true))
        return PiecewiseRule(goal, pieces, tuple(subgoals))
    return Rule(goal, general_value)


def _sqrt_quadratic_denom_value(a, b, c, coeffs, x, solver: ManualSolver) -> Expr:
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
        I0 = S.Zero
    else:
        step = inverse_trig_rule(IntegralInfo(1/s, x), solver, degenerate=False)
        I0 = constant*step.result
    return Add(*(result_coeffs[i]*x**(len(coeffs)-2-i)
                 for i in range(len(result_coeffs))), e/c)*s + I0


def sqrt_quadratic_rule(goal: IntegralInfo, solver: ManualSolver, degenerate: bool = True) -> Rule | None:
    # integrate f(x) * (a + b*x + c*x**2)**(n/2),
    # where f(x) is a polynomial and n is an odd integer
    starting_integrand, x = goal

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
    inner_goal = IntegralInfo(integrand, x)

    generic_cond = Ne(c, 0)
    if not degenerate or generic_cond is S.true:
        degenerate_step: Rule | None = None
    elif b.is_zero:
        degenerate_step = solver.solve(IntegralInfo(f*sqrt(a)**n, x))
    else:
        degenerate_integrand = f*sqrt(a + b*x)**n
        degenerate_step = sqrt_fractional_linear_rule(IntegralInfo(degenerate_integrand, x), solver)
        if degenerate_step is None:
            degenerate_step = DontKnowRule(IntegralInfo(degenerate_integrand, x))

    def sqrt_quadratic_denom_rule(numer_poly: Poly, integrand: Expr) -> Rule:
        denom = sqrt(a+b*x+c*x**2)
        deg = numer_poly.degree()
        if deg <= 1:
            # integrand == (d+e*x)/sqrt(a+b*x+c*x**2)
            e, d = numer_poly.all_coeffs() if deg == 1 else (S.Zero, numer_poly.as_expr())
            # rewrite numerator to A*(2*c*x+b) + B
            A = e/(2*c)
            B = d-A*b
            linear_value: Expr | None = None
            constant_value: Expr | None = None
            if A != 0:
                # d/dx sqrt(a+b*x+c*x**2) = (2*c*x+b) / (2*sqrt(a+b*x+c*x**2))
                # so integral of (2*c*x+b)/denom is 2*sqrt(a+b*x+c*x**2)
                linear_value = 2*sqrt(a+b*x+c*x**2)
                if A != 1:
                    linear_value = A*linear_value
            if B != 0:
                constant_step = inverse_trig_rule(IntegralInfo(1/denom, x), solver, degenerate=False)
                constant_value = constant_step.result
                if B != 1:
                    constant_value = B*constant_value
            if linear_value is not None and constant_value is not None:
                value = linear_value + constant_value
            else:
                value = linear_value if linear_value is not None else constant_value
            return Rule(IntegralInfo(integrand, x), value)
        else:
            coeffs = numer_poly.all_coeffs()
            value = _sqrt_quadratic_denom_value(a, b, c, coeffs, x, solver)
            return SqrtQuadraticDenomRule(IntegralInfo(integrand, x), a, b, c, coeffs, value)

    def sqrt_quadratic_reduction_rule(integrand: Expr, n: Expr, const: Expr) -> Rule:
        # Implementation of Gradshteyn & Ryzhik 2.263.3
        k = (-n - 1) // 2
        delta = 4*a*c - b**2
        R = c*x**2 + b*x + a
        red_goal = IntegralInfo(integrand, x)

        def double_root_step() -> Rule:
            square_base = sqrt(c)*x + b/(2*sqrt(c))
            nested = Pow(Pow(square_base, 2, evaluate=False), S(n)/2, evaluate=False)
            rewritten = const*nested
            substep = nested_pow_rule(IntegralInfo(rewritten, x), solver)
            value = substep.result if substep is not None else Integral(rewritten, x)
            return Rule(red_goal, value)

        if delta.is_zero is True:
            return double_root_step()

        # we divide by delta, then delta  has to be != 0
        term_denom = (2*k - 1) * delta * (R**(S(2*k - 1)/2))
        constant_term = const*2*(2*c*x+b) / term_denom
        coeff = (8*c*(k-1))/((2*k-1) * delta)
        expr = const * R**(S(1)/2 - k)

        derive_expr = Derivative(constant_term, x)
        derive_value = DerivativeRule(IntegralInfo(derive_expr, x)).result

        if coeff == 0:
            rewrite_value = derive_value
        else:
            next_value = solver.extract(IntegralInfo(expr, x))
            rewrite_value = derive_value + coeff * next_value

        nondegenerate_step = Rule(red_goal, rewrite_value)
        if delta.is_zero is None:
            return _add_degenerate_step(
                red_goal,
                Ne(delta, 0),
                nondegenerate_step,
                double_root_step(),
            )
        return nondegenerate_step

    def sqrt_quadratic_polynomial_reduction_rule() -> Rule:
        # reduce non-constant polynomial numerators by writing f = q*R + r,
        # then split the linear remainder into a multiple of R' and a constant.
        terms = []
        values = []
        root_base_ = c*x**2 + b*x + a
        root_poly = Poly(root_base_, x)
        quotient, rest = f_poly.div(root_poly)
        if not quotient.is_zero:
            # n is increasing by 2 at each step, we will fall in one of the cases above
            quotient_integrand = quotient.as_expr() * sqrt(root_base_)**(n + 2)
            quotient_step = sqrt_quadratic_rule(IntegralInfo(quotient_integrand, x), solver, degenerate=False)
            terms.append(quotient_integrand)
            values.append(quotient_step.result if quotient_step is not None else Integral(quotient_integrand, x))
        if not rest.is_zero:
            # split the linear remainder as A*R' + B, where R' = 2*c*x + b.
            e = rest.nth(1)
            d = rest.nth(0)
            A = e/(2*c)
            B = d - A*b
            if A != 0:
                # solved by substitution u = root_base_:
                # integral of (2*c*x+b)*sqrt(root_base_)**n dx
                #   = integral of u**(n/2) du = 2*sqrt(root_base_)**(n+2)/(n+2)
                base = (2*c*x + b) * sqrt(root_base_)**n
                term = A * base
                u_value = 2*sqrt(root_base_)**(n+2) / (n+2)
                values.append(A*u_value if A != 1 else u_value)
                terms.append(term)
            if B != 0:
                term = B * sqrt(root_base_)**n
                const_step = sqrt_quadratic_reduction_rule(term, n, B)
                terms.append(term)
                values.append(const_step.result)
        rewritten = Add(*terms, evaluate=False)
        value = Add(*values)
        return Rule(IntegralInfo(integrand, x), value)

    if n > 0:  # rewrite poly * sqrt(s)**(2*k-1) to poly*s**k / sqrt(s)
        numer_poly = f_poly * (a+b*x+c*x**2)**((n+1)/2)
        rewritten = numer_poly.as_expr()/sqrt(a+b*x+c*x**2)
        substep = sqrt_quadratic_denom_rule(numer_poly, rewritten)
        generic_step: Rule = Rule(inner_goal, substep.result)
    elif n == -1:
        generic_step = sqrt_quadratic_denom_rule(f_poly, integrand)
    elif f_poly.degree() == 0:
        generic_step = sqrt_quadratic_reduction_rule(integrand, n, f)
    else:
        generic_step = sqrt_quadratic_polynomial_reduction_rule()
    step = _add_degenerate_step(inner_goal, generic_cond, generic_step, degenerate_step)
    return Rule(goal, step.result, step.subgoals)


def hyperbolic_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if isinstance(integrand, HyperbolicFunction) and integrand.args[0] == symbol:
        if integrand.func == sinh:
            return SinhRule(goal)
        if integrand.func == cosh:
            return CoshRule(goal)
        if integrand.func == tanh:
            return RewriteRule(goal, sinh(symbol)/cosh(symbol))
        if integrand.func == coth:
            return RewriteRule(goal, cosh(symbol)/sinh(symbol))
        # sech, csch: rewriting via tanh(x/2) produces a rational function
        # of tanh(x/2) that substitution_rule/find_substitutions rediscovers
        # (verified: matches the u = tanh(x/2) substitution the original
        # code built by hand -- e.g. sech(x).rewrite(tanh) ==
        # (1-tanh(x/2)**2)/(tanh(x/2)**2+1), and find_substitutions finds
        # exactly (tanh(x/2), 2, 1/(u**2+1)) for it).
        rewritten = integrand.rewrite(tanh)
        return RewriteRule(goal, rewritten)
    return None


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


def _trig_rewrite_rule(goal: IntegralInfo, rewritten: Expr) -> Rule | None:
    if rewritten != goal.integrand:
        return RewriteRule(goal, rewritten)
    return None


def trig_sincos_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if any(integrand.has(f) for f in (sin, cos)):
        pattern, a, b, m, n = sincos_pattern(symbol)
        match = integrand.match(pattern)
        if not match:
            return None
        aa, bb, mm, nn = (match.get(i, S.Zero) for i in (a, b, m, n))
        if mm.is_even and nn.is_even and mm.is_nonnegative and nn.is_nonnegative:
            rewritten = (((1 - cos(2*aa*symbol)) / 2) ** (mm / 2) *
                         ((1 + cos(2*bb*symbol)) / 2) ** (nn / 2))
        elif mm.is_odd and mm >= 3:
            rewritten = ((1 - cos(aa*symbol)**2)**((mm - 1) / 2) *
                         sin(aa*symbol) * cos(bb*symbol) ** nn)
        elif nn.is_odd and nn >= 3:
            rewritten = ((1 - sin(bb*symbol)**2)**((nn - 1) / 2) *
                         cos(bb*symbol) * sin(aa*symbol) ** mm)
        else:
            return None
        return _trig_rewrite_rule(goal, rewritten)
    return None


def trig_tansec_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    integrand = integrand.subs({
        1 / cos(symbol): sec(symbol)
    })

    if any(integrand.has(f) for f in (tan, sec)):
        pattern, a, b, m, n = tansec_pattern(symbol)
        match = integrand.match(pattern)
        if not match:
            return None
        aa, bb, mm, nn = (match.get(i, S.Zero) for i in (a, b, m, n))
        if mm.is_odd:
            rewritten = ((sec(aa*symbol)**2 - 1) ** ((mm - 1) / 2) *
                         tan(aa*symbol) * sec(bb*symbol) ** nn)
        elif nn.is_even and nn >= 4:
            rewritten = ((1 + tan(bb*symbol)**2) ** (nn/2 - 1) *
                         sec(bb*symbol)**2 * tan(aa*symbol) ** mm)
        elif mm == 2 and nn == 0:
            rewritten = sec(aa*symbol)**2 - 1
        else:
            return None
        return _trig_rewrite_rule(goal, rewritten)
    return None


def trig_cotcsc_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    integrand = integrand.subs({
        1 / sin(symbol): csc(symbol),
        1 / tan(symbol): cot(symbol),
        cos(symbol) / tan(symbol): cot(symbol)
    })

    if any(integrand.has(f) for f in (cot, csc)):
        pattern, a, b, m, n = cotcsc_pattern(symbol)
        match = integrand.match(pattern)
        if not match:
            return None
        aa, bb, mm, nn = (match.get(i, S.Zero) for i in (a, b, m, n))
        if mm.is_odd:
            rewritten = ((csc(aa*symbol)**2 - 1) ** ((mm - 1) / 2) *
                         cot(aa*symbol) * csc(bb*symbol) ** nn)
        elif nn.is_even and nn >= 4:
            rewritten = ((1 + cot(bb*symbol)**2) ** (nn/2 - 1) *
                         csc(bb*symbol)**2 * cot(aa*symbol) ** mm)
        else:
            return None
        return _trig_rewrite_rule(goal, rewritten)
    return None


def trig_sindouble_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    a = Wild('a', exclude=[sin(2*symbol)])
    match = integrand.match(sin(2*symbol)*a)
    if match:
        sin_double = 2*sin(symbol)*cos(symbol)/sin(2*symbol)
        rewritten = integrand * sin_double
        return RewriteRule(goal, rewritten)
    return None


def trig_powers_products_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    for proposer in (trig_sincos_rule, trig_tansec_rule, trig_cotcsc_rule, trig_sindouble_rule):
        result = proposer(goal, solver)
        if result is not None:
            return result
    return None


def trig_substitution_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
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
        restriction: bool | Boolean = True
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

                substep_goal = IntegralInfo(replaced, theta)
                substep = solver.solve(substep_goal)
                if not substep.contains_dont_know():
                    return TrigSubstitutionRule(goal, theta, x_func, replaced,
                                                 restriction, substep_goal)
    return None


def heaviside_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    pattern, m, b, g = heaviside_pattern(symbol)
    match = integrand.match(pattern)
    if match and 0 != match[g]:
        # f = Heaviside(m*x+b)*g
        mm, bb = match[m], match[b]
        return HeavisideRule(goal, mm*symbol + bb, -bb/mm, match[g])
    return None


def dirac_delta_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, x = goal
    if len(integrand.args) == 1:
        n = S.Zero
    else:
        n = integrand.args[1]  # type: ignore
    if not n.is_Integer or n < 0:
        return None
    a, b = Wild('a', exclude=[x]), Wild('b', exclude=[x, 0])
    match = integrand.args[0].match(a+b*x)
    if not match:
        return None
    a, b = match[a], match[b]
    generic_cond = Ne(b, 0)
    if generic_cond is S.true:
        degenerate_step = None
    else:
        degenerate_step = ConstantRule(IntegralInfo(DiracDelta(a, n), x))
    generic_step = DiracDeltaRule(goal, n, a, b)
    return _add_degenerate_step(goal, generic_cond, generic_step, degenerate_step)


def substitution_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal

    u_var = Dummy("u")
    substitutions = find_substitutions(integrand, symbol, u_var)
    if not substitutions:
        return None
    factored_integrand = integrand.factor()
    _, denom_integrand = factored_integrand.as_numer_denom()
    for u_func, c, substituted in substitutions:
        subrule = solver.solve(IntegralInfo(substituted, u_var))
        if subrule.contains_dont_know():
            continue

        sub_value = subrule.result
        if simplify(c - 1) != 0:
            _, denom_c = c.as_numer_denom()
            sub_value = c * sub_value

            if denom_c.free_symbols:
                pieces: list[tuple[Expr, bool | Boolean]] = []
                subgoals: list[IntegralInfo] = []
                factors_denom_c = factor_list(denom_c)[1]
                for pole, _ in factors_denom_c:
                    # only substitute poles introduced by the constant c if they were not already poles of the original integrand
                    if not _if_zero_implies_zero(pole, denom_integrand):
                        rewritten_integral = manual_subs(factored_integrand, pole, 0)
                        # additional check not to replace a if it is not valid (for example ln(a*x))
                        if rewritten_integral.has(S.ComplexInfinity, S.Infinity,
                                                   S.NegativeInfinity, S.NaN):
                            continue
                        pieces.append((Integral(rewritten_integral, symbol), Eq(pole, 0)))
                        subgoals.append(IntegralInfo(rewritten_integral, symbol))
                if pieces:
                    pieces.append((sub_value, True))
                    sub_value = solver.satisfy(
                        PiecewiseRule(IntegralInfo(substituted, symbol), pieces, tuple(subgoals)))

        # same composition USubstitutionRule.eval() performs: avoid needless
        # -log(1/x) from the substitution, then undo it.
        if u_func.is_Pow:
            base, exp_ = u_func.as_base_exp()
            if exp_ == -1:
                sub_value = sub_value.subs(log(u_var), -log(base))
        return Rule(goal, sub_value.subs(u_var, u_func))
    return None


def _rewrite_if_solvable(goal: IntegralInfo, rewritten: Expr, solver: ManualSolver) -> Rule | None:
    """Shared by cancel_rule/partial_fractions_rule/distribute_expand_rule/
    trig_expand_rule: these rewrites are not always productive -- e.g.
    cancel_rule's .cancel() can undo combine_power_rule's exponent merge
    and hand back an ancestor goal, which the open_path cycle guard then
    fails, embedding an unresolved Integral. The original rewriter()
    wrapper guarded against exactly this by pre-checking the rewritten
    form before committing to it -- but with a *shallow* check
    (`isinstance(substep, DontKnowRule)`), not the deep/recursive
    `contains_dont_know()` that substitution_rule/_parts_rule/etc use.
    That means a rewrite whose top-level technique succeeds (e.g. add_rule
    splitting a sum) is accepted even if some of ITS children fail and
    leave their own embedded Integral -- e.g. rewriting (f(x)-f(-x))/x to
    f(x)/x - f(-x)/x is still "useful progress" worth keeping even though
    neither piece, integrating an opaque Function, ever resolves. Only a
    rewrite that fails outright (the cycle guard's direct DontKnowRule, or
    no applicable technique at all) is rejected here, matching that."""
    substep = solver.solve(IntegralInfo(rewritten, goal.symbol))
    if isinstance(substep, DontKnowRule):
        return None
    return Rule(goal, substep.result)


def partial_fractions_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if not integrand.is_rational_function(symbol):
        return None

    rewritten = integrand.apart(symbol)
    if rewritten == integrand:
        # If apart cannot decompose the rational function any further,
        # use ratint as the final fallback for rational integration.
        return RatintRule(goal)
    return _rewrite_if_solvable(goal, rewritten, solver)


def cancel_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if integrand.is_rational_function(symbol):
        return None
    rewritten = integrand.cancel()
    if rewritten != integrand:
        return _rewrite_if_solvable(goal, rewritten, solver)
    return None


def distribute_expand_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    applicable = (
        (isinstance(integrand, (Pow, Mul))
         or all(arg.is_Pow or arg.is_polynomial(symbol) for arg in integrand.args))
        and not integrand.is_rational_function(symbol)
    )
    if not applicable:
        return None
    rewritten = integrand.expand()
    if rewritten != integrand:
        return _rewrite_if_solvable(goal, rewritten, solver)
    return None


def trig_expand_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    # If there are trig functions with different arguments, expand them
    integrand, symbol = goal
    if len({a.args[0] for a in integrand.atoms(TrigonometricFunction)}) > 1:
        rewritten = integrand.expand(trig=True)
        if rewritten != integrand:
            return _rewrite_if_solvable(goal, rewritten, solver)
    return None


def derivative_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand = goal.integrand
    if not isinstance(integrand, Derivative):
        return None
    diff_variables = integrand.variables
    undifferentiated_function = integrand.expr
    integrand_variables = undifferentiated_function.free_symbols

    if goal.symbol in integrand_variables:
        if goal.symbol in diff_variables:
            return DerivativeRule(goal)
        return None
    return ConstantRule(goal)


def rewrites_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    integrand, symbol = goal
    if integrand.match(1/cos(symbol)):
        rewritten = integrand.subs(1/cos(symbol), sec(symbol))
        return RewriteRule(goal, rewritten)
    return None


# ---------------------------------------------------------------------------
# Dispatch: PROPOSER_TABLE / FALLBACK_PROPOSERS
# ---------------------------------------------------------------------------
#
# key(goal) reuses the pre-existing dispatch-key logic (see _key() above).
# This replaces the do_one/switch/condition/multiplexer/alternatives/
# rewriter combinator layer entirely with plain data plus ManualSolver's
# simple first-match dispatch. A few proposers (partial_fractions_rule,
# cancel_rule, combine_power_rule, parts_rule, distribute_expand_rule,
# nested_pow_rule) were only ever tried for a subset of goal types in the
# original code (the old `condition(integral_is_subclass(...), rule)`
# wrappers); rather than hardcoding them into PROPOSER_TABLE (which would
# reorder them relative to rewrites_rule/substitution_rule for e.g. `log`-
# or `atan`-keyed goals), they're kept in FALLBACK_PROPOSERS behind a
# small type gate, preserving the exact original try-order.

def _gate(*keys):
    keys = frozenset(keys)

    def _gate_decorator(proposer: Proposer) -> Proposer:
        def gated(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
            if _key(goal) not in keys:
                return None
            return proposer(goal, solver)
        gated.__name__ = proposer.__name__
        return gated
    return _gate_decorator


_gated_partial_fractions_rule = _gate(Mul, Pow)(partial_fractions_rule)
_gated_cancel_rule = _gate(Mul, Pow)(cancel_rule)
_gated_combine_power_rule = _gate(Mul)(combine_power_rule)
_gated_parts_rule = _gate(Mul, log, *inverse_trig_functions, *special_error_functions)(parts_rule)
_gated_distribute_expand_rule = _gate(Mul, Pow)(distribute_expand_rule)
_gated_nested_pow_rule = _gate(Mul, Pow)(nested_pow_rule)


PROPOSER_TABLE: dict[type, tuple[Proposer, ...]] = {
    Pow: (power_rule, inverse_trig_rule, quadratic_denom_rule, sqrt_quadratic_rule,
          sqrt_fractional_linear_rule, euler_substitution_rule),
    Symbol: (power_rule,),
    exp: (exp_rule,),
    Add: (add_rule,),
    Mul: (mul_rule, trig_product_rule, heaviside_rule, quadratic_denom_rule,
          sqrt_quadratic_rule, sqrt_fractional_linear_rule, euler_substitution_rule,
          trig_cmplx_exp_rule),
    Derivative: (derivative_rule,),
    TrigonometricFunction: (trig_rule,),
    Heaviside: (heaviside_rule,),
    DiracDelta: (dirac_delta_rule,),
    OrthogonalPolynomial: (orthogonal_poly_rule,),
    Number: (constant_rule,),
}

_gated_exp_trig_cyclic_rule = _gate(Mul)(exp_trig_cyclic_rule)

FALLBACK_PROPOSERS: tuple[Proposer, ...] = (
    trig_rule,
    hyperbolic_rule,
    rewrites_rule,
    substitution_rule,
    _gated_partial_fractions_rule,
    _gated_cancel_rule,
    _gated_combine_power_rule,
    # tried immediately before parts_rule: exp(a*x)*{sin,cos,sinh,cosh}(b*x)
    # is exactly the pattern that makes naive by-parts cycle back to a
    # constant multiple of the original goal; this dedicated proposer
    # (design doc section 5.6) computes the closed form directly via
    # algebra instead of relying on the open_path cycle guard, which would
    # otherwise just fail that candidate and leave an unevaluated Integral.
    _gated_exp_trig_cyclic_rule,
    _gated_parts_rule,
    _gated_distribute_expand_rule,
    trig_powers_products_rule,
    trig_expand_rule,
    _gated_nested_pow_rule,
    trig_substitution_rule,
)


def fallback_rule(goal: IntegralInfo, solver: ManualSolver) -> Rule | None:
    return DontKnowRule(goal)


def integral_steps(integrand, symbol, **options):
    """Returns the first step needed to compute an integral.

    Explanation
    ===========

    This function attempts to mirror what a student would do by hand as
    closely as possible.

    SymPy Gamma used to use this to provide a step-by-step explanation of
    an integral; the step-by-step ``Rule`` tree this returns is no longer
    guaranteed to have any particular shape (see ``ManualSolver`` above),
    only that ``.contains_dont_know()`` accurately reports whether the
    integral fully resolved.

    Examples
    ========

    >>> from sympy import exp, sin
    >>> from sympy.integrals.manualintegrate import integral_steps
    >>> from sympy.abc import x
    >>> integral_steps(sin(x), x)
    SinRule(goal=IntegralInfo(integrand=sin(x), symbol=x))

    Returns
    =======

    rule : Rule
        The first step; it may have subgoals that must also be solved
        (already resolved into ``rule.result`` by the time this returns).

    """
    solver = ManualSolver(PROPOSER_TABLE, FALLBACK_PROPOSERS)
    return solver.solve(IntegralInfo(integrand, symbol))


def _reorder_piecewise(result: Expr) -> Expr:
    # If we got Piecewise with two parts, put generic first
    if isinstance(result, Piecewise) and len(result.args) == 2:
        cond = result.args[0][1]
        if isinstance(cond, Eq) and result.args[1][1] == True:  # noqa: E712
            result = result.func(
                (result.args[1][0], Ne(*cond.args)),
                (result.args[0][0], True))
    return result


def _has_erf_trig_mul(expr: Expr) -> bool:
    for sub in expr.find(Mul):
        if sub.has(erf, erfc, erfi) and sub.has(sin, cos, sinh, cosh):
            return True
    return False


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
    solver = ManualSolver(PROPOSER_TABLE, FALLBACK_PROPOSERS)
    result = solver.solve(IntegralInfo(f, var)).result
    result = _reorder_piecewise(result)
    # Factor terms like erf(x)*sin(x) that may have been expanded
    if _has_erf_trig_mul(f) and _has_erf_trig_mul(result):
        result = factor_terms(result)
    return result
