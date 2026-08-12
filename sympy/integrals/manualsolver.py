"""Proof-of-concept AND-OR / hypergraph search engine for manual integration.

This is a prototype, not a replacement for :mod:`sympy.integrals.manualintegrate`
(which is left untouched). It implements just enough of a goal-directed
rewrite system -- ``Rule``, ``Proposer``, ``ManualSolver`` -- to demonstrate
two fixes identified while reviewing a proposed redesign of
``manualintegrate``:

1. **Integration by parts without eager recursion.** In the naive version of
   this design, ``PartsRule``'s antiderivative ``v`` (of ``dv``) would need
   to be known *before* the rule can even state its second subgoal
   (``integrate(du*v, x)``) -- a self-contradiction, since subgoals are
   supposed to be declared before they're solved. Here, ``v`` is a genuine
   solver goal (``PartsRule.v_goal``); the second subgoal is only derived,
   via the new ``Rule.expand()`` hook, once the solver has actually
   resolved ``v``. No rule ever reaches back into the solver.

2. **Custom (non-generic) ``eval()`` overrides.** Beyond the always-cited
   ``USubstitutionRule`` (log(1/x) fixup), ``HeavisideRule`` also needs more
   than plain atom substitution -- it re-substitutes into an already-
   resolved subgoal's value -- and is implemented here too.

Only a small, deliberately non-exhaustive rule/proposer set is implemented:
enough to run end to end on a handful of representative integrals and
cross-check the results against the real ``manualintegrate``. Run this
module directly (``python -m sympy.integrals.manualsolver``) for a demo.
"""

from __future__ import annotations

from abc import ABC
from collections import defaultdict
from enum import Enum, auto
from typing import Callable, Mapping, TYPE_CHECKING

from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.symbol import Dummy, Symbol, Wild
from sympy.core.containers import Tuple as STuple
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.trigonometric import sin, cos, atan
from sympy.functions.special.delta_functions import Heaviside
from .integrals import Integral
from .manualintegrate import (
    IntegralInfo, find_substitutions, manualintegrate as reference_manualintegrate,
)

if TYPE_CHECKING:
    from sympy.core.expr import Expr


# ---------------------------------------------------------------------------
# 1. Node identity: canonicalization (design doc section 2.1)
# ---------------------------------------------------------------------------

# canon() must be a pure, idempotent function of an expression's shape: the
# SAME logical goal reached via two different Python objects has to
# canonicalize to two EQUAL results, or every goal-keyed dict lookup in the
# solver breaks. SymPy Dummy equality is by identity (dummy_index), not by
# name, so minting a *fresh* Dummy("_d0") on every call would silently make
# canon() non-idempotent -- two calls on the same input would produce two
# unequal outputs. Reuse a fixed, pre-created pool instead.
_CANON_POOL = [Dummy(f"_d{i}") for i in range(64)]


def _rename_dummies(expr):
    seen: dict[Dummy, Dummy] = {}

    def walk(e):
        if isinstance(e, Dummy):
            if e not in seen:
                if len(seen) >= len(_CANON_POOL):
                    raise RuntimeError("canon() dummy pool exhausted")
                seen[e] = _CANON_POOL[len(seen)]
        for a in getattr(e, "args", ()):
            walk(a)

    walk(expr)
    return expr.xreplace(seen) if seen else expr


def canonical_goal(info: IntegralInfo) -> IntegralInfo:
    combined = _rename_dummies(STuple(info.integrand, info.symbol))
    return IntegralInfo(combined[0], combined[1])


# ---------------------------------------------------------------------------
# 2. Rule base class
# ---------------------------------------------------------------------------

class Rule(ABC):
    """goal/result/subgoals, exactly as in the design doc, plus one addition:
    `expand()`, used only by rules (currently just PartsRule) whose later
    subgoals cannot be known until an earlier subgoal is actually resolved."""

    __slots__ = ("goal", "result", "subgoals")

    def __init__(self, goal: IntegralInfo, result: Expr, subgoals=()):
        self.goal = goal
        self.result = result
        self.subgoals = tuple(subgoals)

    def eval(self, resolved: Mapping[IntegralInfo, Expr]) -> None:
        """Default: generic atom substitution."""
        for g in self.subgoals:
            atom = Integral(g.integrand, g.symbol)
            self.result = self.result.xreplace({atom: resolved[g]})
        self.subgoals = ()

    def expand(self, resolved: Mapping[IntegralInfo, Expr]):
        """Return one more IntegralInfo goal to solve before this rule is
        considered AND-satisfied, or None if there's nothing more needed.
        `resolved` accumulates every subgoal resolved so far (not just the
        latest batch), so a rule can compute its next goal from an earlier
        one's concrete value. Default: no dynamic expansion."""
        return None

    def __repr__(self):
        return f"{type(self).__name__}(goal={self.goal.integrand!r})"


# --- zero-subgoal atomic rules ---------------------------------------------

class ConstantRule(Rule):
    """integrate(a, x) -> a*x"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo):
        super().__init__(goal, goal.integrand * goal.symbol)


class PowerRule(Rule):
    """integrate(x**n, x) -> x**(n+1)/(n+1), n != -1"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, base, exp_):
        super().__init__(goal, base ** (exp_ + 1) / (exp_ + 1))


class ReciprocalRule(Rule):
    """integrate(1/x, x) -> log(x)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, base):
        super().__init__(goal, log(base))


class ExpRule(Rule):
    """integrate(a**x, x) -> a**x/log(a)"""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, base, exp_):
        super().__init__(goal, goal.integrand / log(base))


class SinRule(Rule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo):
        super().__init__(goal, -cos(goal.symbol))


class CosRule(Rule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo):
        super().__init__(goal, sin(goal.symbol))


class ExpTrigCyclicRule(Rule):
    """Design doc section 5.6: closed form computed directly via algebra at
    proposal time, no recursion, zero subgoals -- CyclicPartsRule's
    replacement."""
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, closed_form):
        super().__init__(goal, closed_form)


# --- composite rules (subgoals declared upfront) ----------------------------

class AddRule(Rule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, terms):
        subgoals = tuple(IntegralInfo(t, goal.symbol) for t in terms)
        result = Add(*(Integral(t, goal.symbol) for t in terms))
        super().__init__(goal, result, subgoals)


class ConstantTimesRule(Rule):
    __slots__ = ()

    def __init__(self, goal: IntegralInfo, constant, other):
        subgoal = IntegralInfo(other, goal.symbol)
        result = constant * Integral(other, goal.symbol)
        super().__init__(goal, result, (subgoal,))


# --- USubstitutionRule: confirmed custom eval() (design doc section 3.2) ---

class USubstitutionRule(Rule):
    __slots__ = ("u_var", "u_func", "constant")

    def __init__(self, goal: IntegralInfo, u_var, u_func, substituted, constant=S.One):
        subgoal = IntegralInfo(substituted, u_var)
        result = constant * Integral(substituted, u_var)
        super().__init__(goal, result, (subgoal,))
        self.u_var = u_var
        self.u_func = u_func
        self.constant = constant

    def eval(self, resolved):
        super().eval(resolved)  # generic substitution first, clears subgoals
        if self.u_func.is_Pow:
            base, exp_ = self.u_func.as_base_exp()
            if exp_ == -1:
                # avoid needless -log(1/x) from substitution
                self.result = self.result.subs(log(self.u_var), -log(base))
        self.result = self.result.subs(self.u_var, self.u_func)


# --- HeavisideRule: another confirmed custom eval(), missing from the doc --

class HeavisideRule(Rule):
    """integrate(Heaviside(m*x+b)*g(x), x)
    -> Heaviside(harg) * (result - result.subs(x, ibnd))

    `eval()` here does a SECOND substitution on the already-resolved
    subgoal's value (result.subs(variable, ibnd)) -- this is not expressible
    as plain atom-for-atom substitution, so like USubstitutionRule it must
    override eval() rather than relying on Rule's default."""
    __slots__ = ("harg", "ibnd", "_subgoal")

    def __init__(self, goal: IntegralInfo, harg, ibnd, g_expr):
        subgoal = IntegralInfo(g_expr, goal.symbol)
        result = Heaviside(harg) * Integral(g_expr, goal.symbol)
        super().__init__(goal, result, (subgoal,))
        self.harg = harg
        self.ibnd = ibnd
        self._subgoal = subgoal

    def eval(self, resolved):
        base_result = resolved[self._subgoal]
        variable = self.goal.symbol
        self.result = Heaviside(self.harg) * (
            base_result - base_result.subs(variable, self.ibnd)
        )
        self.subgoals = ()


# --- PartsRule: the headline fix --------------------------------------------

class PartsRule(Rule):
    """integrate(u*dv, x) -> u*v - integrate(du*v, x)

    THE FIX (design doc review, finding #1): `v` is declared as a real
    subgoal (`v_goal`) instead of being computed eagerly via a recursive
    integral_steps() call before the Rule even exists. The *second* subgoal
    (integrate du*v) genuinely cannot be stated until v has a concrete
    value -- so it isn't declared upfront at all. Instead `expand()` is
    called by the solver once v_goal resolves; only then does it derive
    and hand back the second goal. `eval()` still only needs a resolved
    mapping (v and the second integral's value) to combine into the final
    u*v - second answer -- no rule needs to reach back into the solver.
    """
    __slots__ = ("u", "dv", "du", "v_goal", "second_goal")

    def __init__(self, goal: IntegralInfo, u, dv, du):
        v_goal = IntegralInfo(dv, goal.symbol)
        placeholder = Integral(goal.integrand, goal.symbol)
        super().__init__(goal, placeholder, (v_goal,))
        self.u = u
        self.dv = dv
        self.du = du
        self.v_goal = v_goal
        self.second_goal: IntegralInfo | None = None

    def expand(self, resolved):
        if self.second_goal is None:
            v = resolved[self.v_goal]
            self.second_goal = canonical_goal(IntegralInfo(self.du * v, self.goal.symbol))
            self.subgoals = (self.second_goal,)
            return self.second_goal
        return None

    def eval(self, resolved):
        v = resolved[self.v_goal]
        second = resolved[self.second_goal]
        self.result = self.u * v - second
        self.subgoals = ()


# ---------------------------------------------------------------------------
# 3. Proposers -- deliberately small, not the full LIATE/dispatch catalog
# ---------------------------------------------------------------------------

Proposer = Callable[[IntegralInfo, "ManualSolver"], "list[Rule] | Rule | None"]


def constant_rule(goal, solver):
    if goal.symbol not in goal.integrand.free_symbols:
        return ConstantRule(goal)
    return None


def power_rule(goal, solver):
    integrand, symbol = goal
    if not (integrand.is_Pow or isinstance(integrand, Symbol) or type(integrand) is exp):
        return None
    base, expt = integrand.as_base_exp()
    if symbol not in expt.free_symbols and isinstance(base, Symbol):
        if (expt + 1) == 0:
            return ReciprocalRule(goal, base)
        return PowerRule(goal, base, expt)
    if symbol not in base.free_symbols and isinstance(expt, Symbol):
        return ExpRule(goal, base, expt)
    return None


def trig_rule(goal, solver):
    integrand, symbol = goal
    if integrand == sin(symbol):
        return SinRule(goal)
    if integrand == cos(symbol):
        return CosRule(goal)
    return None


def add_rule(goal, solver):
    integrand, symbol = goal
    if not isinstance(integrand, Add):
        return None
    return AddRule(goal, integrand.as_ordered_terms())


def mul_constant_rule(goal, solver):
    integrand, symbol = goal
    if not isinstance(integrand, Mul):
        return None
    coeff, f = integrand.as_independent(symbol)
    if coeff == 1:
        return None
    return ConstantTimesRule(goal, coeff, f)


def substitution_rule(goal, solver):
    """One OR-candidate Rule per structural substitution match -- no
    recursive pre-filtering. find_substitutions() itself never calls back
    into the solver (it only differentiates/cancels), so it's already
    compliant with the "no recursion" proposer rule and is reused as-is."""
    integrand, symbol = goal
    u_var = Dummy("u")
    subs = find_substitutions(integrand, symbol, u_var)
    if not subs:
        return None
    return [
        USubstitutionRule(goal, u_var, u_func, substituted, constant=c)
        for u_func, c, substituted in subs
    ]


def exp_trig_cyclic_rule(goal, solver):
    """design doc section 5.6 -- retires CyclicPartsRule. Pattern-matches
    exp(a*x)*{sin,cos}(b*x) and writes the closed form directly."""
    integrand, symbol = goal
    a = Wild('a', exclude=[symbol, 0])
    b = Wild('b', exclude=[symbol, 0])
    for trig in (sin, cos):
        match = integrand.match(exp(a * symbol) * trig(b * symbol))
        if not match:
            continue
        aa, bb = match[a], match[b]
        if trig is sin:
            closed = (exp(aa * symbol) * (aa * sin(bb * symbol) - bb * cos(bb * symbol))
                      / (aa ** 2 + bb ** 2))
        else:
            closed = (exp(aa * symbol) * (aa * cos(bb * symbol) + bb * sin(bb * symbol))
                      / (aa ** 2 + bb ** 2))
        return ExpTrigCyclicRule(goal, closed)
    return None


def parts_rule(goal, solver):
    """A deliberately tiny stand-in for the real LIATE-based _parts_rule:
    handles integrand = log(x)/atan(x) (dv=1 trick) and integrand =
    algebraic * {sin, cos, exp} (u=algebraic, dv=rest), including allowing
    u = exp(...) when nothing else qualifies, purely so the demo can show
    a plain (non-cyclic-aware) PartsRule hitting the solver's open_path
    cycle guard on exp(x)*sin(x) and failing -- which is exactly why
    exp_trig_cyclic_rule exists as a dedicated non-recursive proposer.
    No recursive solvability probing anywhere in this function.
    """
    integrand, symbol = goal
    if isinstance(integrand, (log, atan)):
        u, dv = integrand, S.One
    elif isinstance(integrand, Mul):
        alg = [a for a in integrand.args
               if not a.has(symbol) or a.is_polynomial(symbol)]
        rest = [a for a in integrand.args if a not in alg]
        if alg and rest:
            u, dv = Mul(*alg), Mul(*rest)
        else:
            exp_factors = [a for a in integrand.args if a.func is exp]
            if len(exp_factors) != 1:
                return None
            u = exp_factors[0]
            rest = [a for a in integrand.args if a is not u]
            dv = Mul(*rest)
        if symbol not in u.free_symbols:
            return None  # never pick a constant u
    else:
        return None

    du = u.diff(symbol)
    cache_key = u.xreplace({symbol: solver.cache_dummy})
    if solver.parts_u_usage[cache_key] > 2:
        return None
    solver.parts_u_usage[cache_key] += 1
    return PartsRule(goal, u, dv, du)


def heaviside_rule(goal, solver):
    integrand, symbol = goal
    if not integrand.has(Heaviside):
        return None
    m = Wild('m', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    g = Wild('g')
    match = integrand.match(Heaviside(m * symbol + b) * g)
    if not match or match[g] == 0:
        return None
    mm, bb, gg = match[m], match[b], match[g]
    return HeavisideRule(goal, mm * symbol + bb, -bb / mm, gg)


# ---------------------------------------------------------------------------
# 4. ManualSolver: AND-OR search with a solver-side cycle guard
# ---------------------------------------------------------------------------

class GoalStatus(Enum):
    OPEN = auto()
    SOLVED = auto()
    FAILED = auto()


def _key(goal: IntegralInfo):
    integrand = goal.integrand
    if goal.symbol not in integrand.free_symbols:
        return Number
    if integrand.is_Pow or isinstance(integrand, Symbol):
        return "pow"
    return type(integrand)


class ManualSolver:
    def __init__(self, proposer_table: dict, fallback_proposers=()):
        self.proposer_table = proposer_table
        self.fallback_proposers = list(fallback_proposers)

        self.goals: dict[IntegralInfo, GoalStatus] = {}
        self.winning_rule: dict[IntegralInfo, Rule] = {}
        self.candidates: dict[IntegralInfo, list[Rule]] = {}
        self.resolved: dict[IntegralInfo, Expr] = {}
        self.open_path: list[IntegralInfo] = []
        self.parts_u_usage = defaultdict(int)
        self.cache_dummy = Dummy("z")

    def _generate_candidates(self, goal: IntegralInfo) -> list[Rule]:
        proposers = list(self.proposer_table.get(_key(goal), ())) + self.fallback_proposers
        out: list[Rule] = []
        for proposer in proposers:
            produced = proposer(goal, self)
            if produced is None:
                continue
            out.extend(produced if isinstance(produced, list) else [produced])
        return out

    def solve(self, goal: IntegralInfo) -> bool:
        cgoal = canonical_goal(goal)
        if cgoal in self.goals:
            return self.goals[cgoal] is GoalStatus.SOLVED
        if cgoal in self.open_path:
            # solver-side search cycle (design doc section 6.3): this
            # specific expansion path fails, not the whole goal.
            return False

        self.goals[cgoal] = GoalStatus.OPEN
        self.open_path.append(cgoal)
        try:
            candidates = self.candidates.setdefault(cgoal, self._generate_candidates(goal))
            for rule in candidates:
                value = self._try_rule(rule)
                if value is not None:
                    self.goals[cgoal] = GoalStatus.SOLVED
                    self.winning_rule[cgoal] = rule
                    self.resolved[cgoal] = value
                    return True
            self.goals[cgoal] = GoalStatus.FAILED
            return False
        finally:
            self.open_path.pop()

    def _try_rule(self, rule: Rule):
        """AND across (possibly dynamically-growing) subgoals. POC
        simplification: subgoals are resolved eagerly here (not in a
        separate later extract() pass) because rules with expand() -- i.e.
        PartsRule -- need a concrete value to derive their next goal. Rules
        that never override expand() behave exactly like the doc's
        two-phase solve-then-extract design."""
        accumulated: dict[IntegralInfo, Expr] = {}
        frontier = rule.subgoals
        while True:
            for g in frontier:
                if not self.solve(g):
                    return None
                accumulated[g] = self.resolved[canonical_goal(g)]
            rule.subgoals = frontier
            next_goal = rule.expand(accumulated)
            if next_goal is None:
                break
            frontier = (next_goal,)
        rule.eval(accumulated)
        return rule.result


def solve_integral(f, var, proposer_table, fallback_proposers=()):
    solver = ManualSolver(proposer_table, fallback_proposers)
    root = canonical_goal(IntegralInfo(f, var))
    if solver.solve(root):
        return solver.resolved[root], solver
    return Integral(f, var), solver


# ---------------------------------------------------------------------------
# 5. Demo
# ---------------------------------------------------------------------------

FULL_TABLE = {
    Number: [constant_rule],
    "pow": [power_rule],
    Add: [add_rule],
    Mul: [mul_constant_rule, heaviside_rule, substitution_rule, exp_trig_cyclic_rule, parts_rule],
    exp: [power_rule],
    log: [parts_rule],
    atan: [parts_rule],
    Heaviside: [heaviside_rule],
}
FALLBACK = [trig_rule]


def _demo(label, f, var):
    got, _ = solve_integral(f, var, FULL_TABLE, FALLBACK)
    ref = reference_manualintegrate(f, var)
    ok = (got - ref).simplify() == 0
    status = "OK " if ok else "MISMATCH"
    print(f"[{status}] {label}")
    print(f"    poc:  {got}")
    print(f"    real: {ref}")


if __name__ == "__main__":
    from sympy.abc import x

    print("=== basic atomic / composite rules ===")
    _demo("sin(x)", sin(x), x)
    _demo("x**2 + sin(x)", x**2 + sin(x), x)
    _demo("3*sin(x)", 3*sin(x), x)

    print("\n=== USubstitutionRule (custom eval): x*exp(-x**2) ===")
    _demo("x*exp(-x**2)", x*exp(-x**2), x)

    print("\n=== PartsRule with v as a real solver goal: log(x) ===")
    _demo("log(x)", log(x), x)

    print("\n=== PartsRule dynamic second-subgoal, algebraic*trig: x*sin(x) ===")
    _demo("x*sin(x)", x*sin(x), x)

    print("\n=== HeavisideRule (custom eval): Heaviside(x - 1)*x ===")
    _demo("Heaviside(x-1)*x", Heaviside(x - 1)*x, x)

    print("\n=== cyclic-by-parts: exp(x)*sin(x) ===")
    print("-- plain PartsRule only (no dedicated cyclic rule) --")
    table_no_cyclic = dict(FULL_TABLE)
    table_no_cyclic[Mul] = [mul_constant_rule, substitution_rule, parts_rule]
    got_fail, solver_fail = solve_integral(exp(x)*sin(x), x, table_no_cyclic, FALLBACK)
    print(f"    result: {got_fail}   (expected to stay an unevaluated Integral --")
    print("     every PartsRule expansion eventually re-hits the same goal,")
    print("     which the open_path cycle guard fails rather than solving)")
    print("-- with exp_trig_cyclic_rule added --")
    _demo("exp(x)*sin(x)", exp(x)*sin(x), x)
