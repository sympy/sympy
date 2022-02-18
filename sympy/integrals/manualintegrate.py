"""Integration method that emulates by-hand techniques.

This module also provides functionality to get the steps used to evaluate a
particular integral, in the ``integral_steps`` function. This will return
nested namedtuples representing the integration rules used. The
``manualintegrate`` function computes the integral using those steps given
an integrand; given the steps, ``_manualintegrate`` will evaluate them.

The integrator can be extended with new heuristics and evaluation
techniques. To do so, write a function that accepts an ``IntegralInfo``
object and returns either a namedtuple representing a rule or
``None``. Then, write another function that accepts the namedtuple's fields
and returns the antiderivative, and decorate it with
``@evaluates(namedtuple_type)``.  If the new technique requires a new
match, add the key and call to the antiderivative function to integral_steps.
To enable simple substitutions, add the match to find_substitutions.

"""

from typing import Dict as tDict, Optional
from collections import namedtuple, defaultdict
from collections.abc import Mapping
from functools import reduce

from sympy.core.add import Add
from sympy.core.cache import cacheit
from sympy.core.containers import Dict
from sympy.core.expr import Expr
from sympy.core.function import Derivative
from sympy.core.logic import fuzzy_not
from sympy.core.mul import Mul
from sympy.core.numbers import Integer, Number, E
from sympy.core.power import Pow
from sympy.core.relational import Eq, Ne, Gt, Lt
from sympy.core.singleton import S
from sympy.core.symbol import Dummy, Symbol, Wild
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.hyperbolic import (cosh, sinh, acosh, asinh,
                                                   acoth, atanh)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import (TrigonometricFunction,
    cos, sin, tan, cot, csc, sec, acos, asin, atan, acot, acsc, asec)
from sympy.functions.special.delta_functions import Heaviside
from sympy.functions.special.error_functions import (erf, erfi, fresnelc,
    fresnels, Ci, Chi, Si, Shi, Ei, li)
from sympy.functions.special.gamma_functions import uppergamma
from sympy.functions.special.elliptic_integrals import elliptic_e, elliptic_f
from sympy.functions.special.polynomials import (chebyshevt, chebyshevu,
    legendre, hermite, laguerre, assoc_laguerre, gegenbauer, jacobi,
    OrthogonalPolynomial)
from sympy.functions.special.zeta_functions import polylog
from .integrals import Integral
from sympy.logic.boolalg import And
from sympy.ntheory.factor_ import divisors
from sympy.polys.polytools import degree
from sympy.simplify.radsimp import fraction
from sympy.simplify.simplify import simplify
from sympy.solvers.solvers import solve
from sympy.strategies.core import switch, do_one, null_safe, condition
from sympy.utilities.iterables import iterable
from sympy.utilities.misc import debug


def Rule(name, props=""):
    # GOTCHA: namedtuple class name not considered!
    def __eq__(self, other):
        return self.__class__ == other.__class__ and tuple.__eq__(self, other)
    __neq__ = lambda self, other: not __eq__(self, other)
    cls = namedtuple(name, props + " context symbol")
    cls.__eq__ = __eq__
    cls.__ne__ = __neq__
    return cls

ConstantRule = Rule("ConstantRule", "constant")
ConstantTimesRule = Rule("ConstantTimesRule", "constant other substep")
PowerRule = Rule("PowerRule", "base exp")
AddRule = Rule("AddRule", "substeps")
URule = Rule("URule", "u_var u_func constant substep")
PartsRule = Rule("PartsRule", "u dv v_step second_step")
CyclicPartsRule = Rule("CyclicPartsRule", "parts_rules coefficient")
TrigRule = Rule("TrigRule", "func arg")
ExpRule = Rule("ExpRule", "base exp")
ReciprocalRule = Rule("ReciprocalRule", "func")
ArcsinRule = Rule("ArcsinRule")
InverseHyperbolicRule = Rule("InverseHyperbolicRule", "func")
AlternativeRule = Rule("AlternativeRule", "alternatives")
DontKnowRule = Rule("DontKnowRule")
DerivativeRule = Rule("DerivativeRule")
RewriteRule = Rule("RewriteRule", "rewritten substep")
PiecewiseRule = Rule("PiecewiseRule", "subfunctions")
HeavisideRule = Rule("HeavisideRule", "harg ibnd substep")
TrigSubstitutionRule = Rule("TrigSubstitutionRule",
                            "theta func rewritten substep restriction")
ArctanRule = Rule("ArctanRule", "a b c")
ArccothRule = Rule("ArccothRule", "a b c")
ArctanhRule = Rule("ArctanhRule", "a b c")
JacobiRule = Rule("JacobiRule", "n a b")
GegenbauerRule = Rule("GegenbauerRule", "n a")
ChebyshevTRule = Rule("ChebyshevTRule", "n")
ChebyshevURule = Rule("ChebyshevURule", "n")
LegendreRule = Rule("LegendreRule", "n")
HermiteRule = Rule("HermiteRule", "n")
LaguerreRule = Rule("LaguerreRule", "n")
AssocLaguerreRule = Rule("AssocLaguerreRule", "n a")
CiRule = Rule("CiRule", "a b")
ChiRule = Rule("ChiRule", "a b")
EiRule = Rule("EiRule", "a b")
SiRule = Rule("SiRule", "a b")
ShiRule = Rule("ShiRule", "a b")
ErfRule = Rule("ErfRule", "a b c")
FresnelCRule = Rule("FresnelCRule", "a b c")
FresnelSRule = Rule("FresnelSRule", "a b c")
LiRule = Rule("LiRule", "a b")
PolylogRule = Rule("PolylogRule", "a b")
UpperGammaRule = Rule("UpperGammaRule", "a e")
EllipticFRule = Rule("EllipticFRule", "a d")
EllipticERule = Rule("EllipticERule", "a d")

IntegralInfo = namedtuple('IntegralInfo', 'integrand symbol')

evaluators = {}
def evaluates(rule):
    def _evaluates(func):
        func.rule = rule
        evaluators[rule] = func
        return func
    return _evaluates

def contains_dont_know(rule):
    if isinstance(rule, DontKnowRule):
        return True
    else:
        for val in rule:
            if isinstance(val, tuple):
                if contains_dont_know(val):
                    return True
            elif isinstance(val, list):
                if any(contains_dont_know(i) for i in val):
                    return True
    return False

def manual_diff(f, symbol):
    """Derivative of f in form expected by find_substitutions

    SymPy's derivatives for some trig functions (like cot) aren't in a form
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
            return sum([manual_diff(arg, symbol) for arg in f.args])
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
        if symbol not in substituted.free_symbols:
            # replaced everything already
            return False
        debug("substituted: {}, u: {}, u_var: {}".format(substituted, u, u_var))
        substituted = manual_subs(substituted, u, u_var).cancel()

        if symbol not in substituted.free_symbols:
            # avoid increasing the degree of a rational function
            if integrand.is_rational_function(symbol) and substituted.is_rational_function(u_var):
                deg_before = max([degree(t, symbol) for t in integrand.as_numer_denom()])
                deg_after = max([degree(t, u_var) for t in substituted.as_numer_denom()])
                if deg_after > deg_before:
                    return False
            return substituted.as_independent(u_var, as_Add=False)

        # special treatment for substitutions u = (a*x+b)**(1/n)
        if (isinstance(u, Pow) and (1/u.exp).is_Integer and
            Abs(u.exp) < 1):
                a = Wild('a', exclude=[symbol])
                b = Wild('b', exclude=[symbol])
                match = u.base.match(a*symbol + b)
                if match:
                    a, b = [match.get(i, S.Zero) for i in (a, b)]
                    if a != 0 and b != 0:
                        substituted = substituted.subs(symbol,
                            (u_var**(1/u.exp) - b)/a)
                        return substituted.as_independent(u_var, as_Add=False)

        return False

    def possible_subterms(term):
        if isinstance(term, (TrigonometricFunction,
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
            r = []
            if term.args[1].is_constant(symbol):
                r.append(term.args[0])
            elif term.args[0].is_constant(symbol):
                r.append(term.args[1])
            if term.args[1].is_Integer:
                r.extend([term.args[0]**d for d in divisors(term.args[1])
                    if 1 < d < abs(term.args[1])])
                if term.args[0].is_Add:
                    r.extend([t for t in possible_subterms(term.args[0])
                        if t.is_Pow])
            return r
        elif isinstance(term, Add):
            r = []
            for arg in term.args:
                r.append(arg)
                r.extend(possible_subterms(arg))
            return r
        return []

    for u in possible_subterms(integrand):
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
                    return RewriteRule(
                        rewritten,
                        substep,
                        integrand, symbol)
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
                return RewriteRule(
                    rewritten,
                    integral_steps(rewritten, symbol),
                    integrand, symbol)
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
            doable = [rule for rule in alts if not contains_dont_know(rule)]
            if doable:
                return AlternativeRule(doable, *integral)
            else:
                return AlternativeRule(alts, *integral)
    return _alternatives

def constant_rule(integral):
    return ConstantRule(integral.integrand, *integral)

def power_rule(integral):
    integrand, symbol = integral
    base, expt = integrand.as_base_exp()

    if symbol not in expt.free_symbols and isinstance(base, Symbol):
        if simplify(expt + 1) == 0:
            return ReciprocalRule(base, integrand, symbol)
        return PowerRule(base, expt, integrand, symbol)
    elif symbol not in base.free_symbols and isinstance(expt, Symbol):
        rule = ExpRule(base, expt, integrand, symbol)

        if fuzzy_not(log(base).is_zero):
            return rule
        elif log(base).is_zero:
            return ConstantRule(1, 1, symbol)

        return PiecewiseRule([
            (rule, Ne(log(base), 0)),
            (ConstantRule(1, 1, symbol), True)
        ], integrand, symbol)

def exp_rule(integral):
    integrand, symbol = integral
    if isinstance(integrand.args[0], Symbol):
        return ExpRule(E, integrand.args[0], integrand, symbol)


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
                    args = integrand.args[:var_index] + (integrand, symbol)
                    return orthogonal_poly_classes[klass](*args)


def special_function_rule(integral):
    integrand, symbol = integral
    a = Wild('a', exclude=[symbol], properties=[lambda x: not x.is_zero])
    b = Wild('b', exclude=[symbol])
    c = Wild('c', exclude=[symbol])
    d = Wild('d', exclude=[symbol], properties=[lambda x: not x.is_zero])
    e = Wild('e', exclude=[symbol], properties=[
        lambda x: not (x.is_nonnegative and x.is_integer)])
    wilds = (a, b, c, d, e)
    # patterns consist of a SymPy class, a wildcard expr, an optional
    # condition coded as a lambda (when Wild properties are not enough),
    # followed by an applicable rule
    patterns = (
        (Mul, exp(a*symbol + b)/symbol, None, EiRule),
        (Mul, cos(a*symbol + b)/symbol, None, CiRule),
        (Mul, cosh(a*symbol + b)/symbol, None, ChiRule),
        (Mul, sin(a*symbol + b)/symbol, None, SiRule),
        (Mul, sinh(a*symbol + b)/symbol, None, ShiRule),
        (Pow, 1/log(a*symbol + b), None, LiRule),
        (exp, exp(a*symbol**2 + b*symbol + c), None, ErfRule),
        (sin, sin(a*symbol**2 + b*symbol + c), None, FresnelSRule),
        (cos, cos(a*symbol**2 + b*symbol + c), None, FresnelCRule),
        (Mul, symbol**e*exp(a*symbol), None, UpperGammaRule),
        (Mul, polylog(b, a*symbol)/symbol, None, PolylogRule),
        (Pow, 1/sqrt(a - d*sin(symbol)**2),
            lambda a, d: a != d, EllipticFRule),
        (Pow, sqrt(a - d*sin(symbol)**2),
            lambda a, d: a != d, EllipticERule),
    )
    for p in patterns:
        if isinstance(integrand, p[0]):
            match = integrand.match(p[1])
            if match:
                wild_vals = tuple(match.get(w) for w in wilds
                                  if match.get(w) is not None)
                if p[2] is None or p[2](*wild_vals):
                    args = wild_vals + (integrand, symbol)
                    return p[3](*args)


def inverse_trig_rule(integral):
    integrand, symbol = integral
    base, exp = integrand.as_base_exp()
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    match = base.match(a + b*symbol**2)

    if not match:
        return

    def negative(x):
        return x.is_negative or x.could_extract_minus_sign()

    def ArcsinhRule(integrand, symbol):
        return InverseHyperbolicRule(asinh, integrand, symbol)

    def ArccoshRule(integrand, symbol):
        return InverseHyperbolicRule(acosh, integrand, symbol)

    def make_inverse_trig(RuleClass, base_exp, a, sign_a, b, sign_b):
        u_var = Dummy("u")
        current_base = base
        current_symbol = symbol
        constant = u_func = u_constant = substep = None
        factored = integrand
        if a != 1:
            constant = a**base_exp
            current_base = sign_a + sign_b * (b/a) * current_symbol**2
            factored = current_base ** base_exp
        if (b/a) != 1:
            u_func = sqrt(b/a) * symbol
            u_constant = sqrt(a/b)
            current_symbol = u_var
            current_base = sign_a + sign_b * current_symbol**2

        substep = RuleClass(current_base ** base_exp, current_symbol)
        if u_func is not None:
            if u_constant != 1 and substep is not None:
                substep = ConstantTimesRule(
                    u_constant, current_base ** base_exp, substep,
                    u_constant * current_base ** base_exp, symbol)
            substep = URule(u_var, u_func, u_constant, substep, factored, symbol)
        if constant is not None and substep is not None:
            substep = ConstantTimesRule(constant, factored, substep, integrand, symbol)
        return substep

    a, b = [match.get(i, S.Zero) for i in (a, b)]
    # list of (rule, base_exp, a, sign_a, b, sign_b, condition)
    possibilities = []

    if simplify(2*exp + 1) == 0:
        possibilities.append((ArcsinRule, exp, a, 1, -b, -1, And(a > 0, b < 0)))
        possibilities.append((ArcsinhRule, exp, a, 1, b, 1, And(a > 0, b > 0)))
        possibilities.append((ArccoshRule, exp, -a, -1, b, 1, And(a < 0, b > 0)))

    possibilities = [p for p in possibilities if p[-1] is not S.false]
    if a.is_number and b.is_number:
        possibility = [p for p in possibilities if p[-1] is S.true]
        if len(possibility) == 1:
            return make_inverse_trig(*possibility[0][:-1])
    elif possibilities:
        return PiecewiseRule(
            [(make_inverse_trig(*p[:-1]), p[-1]) for p in possibilities],
            integrand, symbol)

def add_rule(integral):
    integrand, symbol = integral
    results = [integral_steps(g, symbol)
              for g in integrand.as_ordered_terms()]
    return None if None in results else AddRule(results, integrand, symbol)

def mul_rule(integral):
    integrand, symbol = integral

    # Constant times function case
    coeff, f = integrand.as_independent(symbol)
    next_step = integral_steps(f, symbol)

    if coeff != 1 and next_step is not None:
        return ConstantTimesRule(
            coeff, f,
            next_step,
            integrand, symbol)

def _parts_rule(integrand, symbol):
    # LIATE rule:
    # log, inverse trig, algebraic, trigonometric, exponential
    def pull_out_algebraic(integrand):
        integrand = integrand.cancel().together()
        # iterating over Piecewise args would not work here
        algebraic = ([] if isinstance(integrand, Piecewise)
            else [arg for arg in integrand.args if arg.is_algebraic_expr(symbol)])
        if algebraic:
            u = Mul(*algebraic)
            dv = (integrand / u).cancel()
            return u, dv

    def pull_out_u(*functions):
        def pull_out_u_rl(integrand):
            if any(integrand.has(f) for f in functions):
                args = [arg for arg in integrand.args
                        if any(isinstance(arg, cls) for cls in functions)]
                if args:
                    u = reduce(lambda a,b: a*b, args)
                    dv = integrand / u
                    return u, dv

        return pull_out_u_rl

    liate_rules = [pull_out_u(log), pull_out_u(*inverse_trig_functions),
                   pull_out_algebraic, pull_out_u(sin, cos),
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
                return

            u = u.subs(dummy, 1)
            dv = dv.subs(dummy, 1)

            # Don't pick a non-polynomial algebraic to be differentiated
            if rule == pull_out_algebraic and not u.is_polynomial(symbol):
                return
            # Don't trade one logarithm for another
            if isinstance(u, log):
                rec_dv = 1/dv
                if (rec_dv.is_polynomial(symbol) and
                    degree(rec_dv, symbol) == 1):
                        return

            # Can integrate a polynomial times OrthogonalPolynomial
            if rule == pull_out_algebraic and isinstance(dv, OrthogonalPolynomial):
                    v_step = integral_steps(dv, symbol)
                    if contains_dont_know(v_step):
                        return
                    else:
                        du = u.diff(symbol)
                        v = _manualintegrate(v_step)
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
                if not contains_dont_know(v_step):
                    v = _manualintegrate(v_step)
                    return u, dv, v, du, v_step


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
                rule = CyclicPartsRule(
                    [PartsRule(u, dv, v_step, None, None, None)
                     for (u, dv, v, du, v_step) in steps],
                    (-1) ** len(steps) * coefficient,
                    integrand, symbol
                )
                if (constant != 1) and rule:
                    rule = ConstantTimesRule(constant, integrand, rule,
                                             constant * integrand, symbol)
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
            return PartsRule(u, dv, v_step,
                             make_second_step(steps[1:], v * du),
                             integrand, symbol)
        else:
            steps = integral_steps(integrand, symbol)
            if steps:
                return steps
            else:
                return DontKnowRule(integrand, symbol)

    if steps:
        u, dv, v, du, v_step = steps[0]
        rule = PartsRule(u, dv, v_step,
                         make_second_step(steps[1:], v * du),
                         integrand, symbol)
        if (constant != 1) and rule:
            rule = ConstantTimesRule(constant, integrand, rule,
                                     constant * integrand, symbol)
        return rule


def trig_rule(integral):
    integrand, symbol = integral
    if isinstance(integrand, (sin, cos)):
        arg = integrand.args[0]

        if not isinstance(arg, Symbol):
            return  # perhaps a substitution can deal with it

        if isinstance(integrand, sin):
            func = 'sin'
        else:
            func = 'cos'

        return TrigRule(func, arg, integrand, symbol)

    if integrand == sec(symbol)**2:
        return TrigRule('sec**2', symbol, integrand, symbol)
    elif integrand == csc(symbol)**2:
        return TrigRule('csc**2', symbol, integrand, symbol)

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

    return RewriteRule(
        rewritten,
        integral_steps(rewritten, symbol),
        integrand, symbol
    )

def trig_product_rule(integral):
    integrand, symbol = integral

    sectan = sec(symbol) * tan(symbol)
    q = integrand / sectan

    if symbol not in q.free_symbols:
        rule = TrigRule('sec*tan', symbol, sectan, symbol)
        if q != 1 and rule:
            rule = ConstantTimesRule(q, sectan, rule, integrand, symbol)

        return rule

    csccot = -csc(symbol) * cot(symbol)
    q = integrand / csccot

    if symbol not in q.free_symbols:
        rule = TrigRule('csc*cot', symbol, csccot, symbol)
        if q != 1 and rule:
            rule = ConstantTimesRule(q, csccot, rule, integrand, symbol)

        return rule

def quadratic_denom_rule(integral):
    integrand, symbol = integral
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    c = Wild('c', exclude=[symbol])

    match = integrand.match(a / (b * symbol ** 2 + c))

    if match:
        a, b, c = match[a], match[b], match[c]
        if b.is_extended_real and c.is_extended_real:
            return PiecewiseRule([(ArctanRule(a, b, c, integrand, symbol), Gt(c / b, 0)),
                                (ArccothRule(a, b, c, integrand, symbol), And(Gt(symbol ** 2, -c / b), Lt(c / b, 0))),
                                (ArctanhRule(a, b, c, integrand, symbol), And(Lt(symbol ** 2, -c / b), Lt(c / b, 0))),
            ], integrand, symbol)
        else:
            return ArctanRule(a, b, c, integrand, symbol)

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
            return URule(u, u_func, None, next_step, integrand2, symbol)
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
        step1 = URule(u,
                      denominator,
                      const,
                      integral_steps(u**(-1), u),
                      integrand,
                      symbol)
        if const != 1:
            step1 = ConstantTimesRule(const,
                                      numer1/denominator,
                                      step1,
                                      const*numer1/denominator,
                                      symbol)
        if numer2.is_zero:
            return step1
        step2 = integral_steps(numer2/denominator, symbol)
        substeps = AddRule([step1, step2], integrand, symbol)
        rewriten = const*numer1/denominator+numer2/denominator
        return RewriteRule(rewriten, substeps, integrand, symbol)

    return

def root_mul_rule(integral):
    integrand, symbol = integral
    a = Wild('a', exclude=[symbol])
    b = Wild('b', exclude=[symbol])
    c = Wild('c')
    match = integrand.match(sqrt(a * symbol + b) * c)

    if not match:
        return

    a, b, c = match[a], match[b], match[c]
    d = Wild('d', exclude=[symbol])
    e = Wild('e', exclude=[symbol])
    f = Wild('f')
    recursion_test = c.match(sqrt(d * symbol + e) * f)
    if recursion_test:
        return

    u = Dummy('u')
    u_func = sqrt(a * symbol + b)
    integrand = integrand.subs(u_func, u)
    integrand = integrand.subs(symbol, (u**2 - b) / a)
    integrand = integrand * 2 * u / a
    next_step = integral_steps(integrand, u)
    if next_step:
        return URule(u, u_func, None, next_step, integrand, symbol)

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
            return RewriteRule(
                rewritten,
                integral_steps(rewritten, symbol),
                integrand, symbol)
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
                if not contains_dont_know(substep):
                    return TrigSubstitutionRule(
                        theta, x_func, replaced, substep, restriction,
                        integrand, symbol)

def heaviside_rule(integral):
    integrand, symbol = integral
    pattern, m, b, g = heaviside_pattern(symbol)
    match = integrand.match(pattern)
    if match and 0 != match[g]:
        # f = Heaviside(m*x + b)*g
        v_step = integral_steps(match[g], symbol)
        result = _manualintegrate(v_step)
        m, b = match[m], match[b]
        return HeavisideRule(m*symbol + b, -b/m, result, integrand, symbol)

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

            if contains_dont_know(subrule):
                continue

            if simplify(c - 1) != 0:
                _, denom = c.as_numer_denom()
                if subrule:
                    subrule = ConstantTimesRule(c, substituted, subrule, substituted, u_var)

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
                    subrule = PiecewiseRule(piecewise, substituted, symbol)

            ways.append(URule(u_var, u_func, c,
                              subrule,
                              integrand, symbol))

        if len(ways) > 1:
            return AlternativeRule(ways, integrand, symbol)
        elif ways:
            return ways[0]

    elif integrand.has(exp):
        u_func = exp(symbol)
        c = 1
        substituted = integrand / u_func.diff(symbol)
        substituted = substituted.subs(u_func, u_var)

        if symbol not in substituted.free_symbols:
            return URule(u_var, u_func, c,
                         integral_steps(substituted, u_var),
                         integrand, symbol)

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
        all(arg.is_Pow or arg.is_polynomial(symbol) for arg in integrand.args)
        or isinstance(integrand, Pow)
        or isinstance(integrand, Mul)),
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
        return ConstantRule(integral.integrand, *integral)

def rewrites_rule(integral):
    integrand, symbol = integral

    if integrand.match(1/cos(symbol)):
        rewritten = integrand.subs(1/cos(symbol), sec(symbol))
        return RewriteRule(rewritten, integral_steps(rewritten, symbol), integrand, symbol)

def fallback_rule(integral):
    return DontKnowRule(*integral)

# Cache is used to break cyclic integrals.
# Need to use the same dummy variable in cached expressions for them to match.
# Also record "u" of integration by parts, to avoid infinite repetition.
_integral_cache = {}  # type: tDict[Expr, Optional[Expr]]
_parts_u_cache = defaultdict(int)  # type: tDict[Expr, int]
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
    URule(u_var=_u, u_func=exp(x), constant=1,
    substep=PiecewiseRule(subfunctions=[(ArctanRule(a=1, b=1, c=1, context=1/(_u**2 + 1), symbol=_u), True),
        (ArccothRule(a=1, b=1, c=1, context=1/(_u**2 + 1), symbol=_u), False),
        (ArctanhRule(a=1, b=1, c=1, context=1/(_u**2 + 1), symbol=_u), False)],
    context=1/(_u**2 + 1), symbol=_u), context=exp(x)/(exp(2*x) + 1), symbol=x)
    >>> print(repr(integral_steps(sin(x), x))) \
    # doctest: +NORMALIZE_WHITESPACE
    TrigRule(func='sin', arg=x, context=sin(x), symbol=x)
    >>> print(repr(integral_steps((x**2 + 3)**2, x))) \
    # doctest: +NORMALIZE_WHITESPACE
    RewriteRule(rewritten=x**4 + 6*x**2 + 9,
    substep=AddRule(substeps=[PowerRule(base=x, exp=4, context=x**4, symbol=x),
        ConstantTimesRule(constant=6, other=x**2,
            substep=PowerRule(base=x, exp=2, context=x**2, symbol=x),
                context=6*x**2, symbol=x),
        ConstantRule(constant=9, context=9, symbol=x)],
    context=x**4 + 6*x**2 + 9, symbol=x), context=(x**2 + 3)**2, symbol=x)


    Returns
    =======

    rule : namedtuple
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
        elif isinstance(integrand, TrigonometricFunction):
            return TrigonometricFunction
        elif isinstance(integrand, Derivative):
            return Derivative
        else:
            for cls in (Pow, Symbol, exp, log,
                        Add, Mul, *inverse_trig_functions,
                        Heaviside, OrthogonalPolynomial):
                if isinstance(integrand, cls):
                    return cls


    def integral_is_subclass(*klasses):
        def _integral_is_subclass(integral):
            k = key(integral)
            return k and issubclass(k, klasses)
        return _integral_is_subclass

    result = do_one(
        null_safe(special_function_rule),
        null_safe(switch(key, {
            Pow: do_one(null_safe(power_rule), null_safe(inverse_trig_rule), \
                              null_safe(quadratic_denom_rule)),
            Symbol: power_rule,
            exp: exp_rule,
            Add: add_rule,
            Mul: do_one(null_safe(mul_rule), null_safe(trig_product_rule), \
                              null_safe(heaviside_rule), null_safe(quadratic_denom_rule), \
                              null_safe(root_mul_rule)),
            Derivative: derivative_rule,
            TrigonometricFunction: trig_rule,
            Heaviside: heaviside_rule,
            OrthogonalPolynomial: orthogonal_poly_rule,
            Number: constant_rule
        })),
        do_one(
            null_safe(trig_rule),
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
            null_safe(trig_substitution_rule)
        ),
        fallback_rule)(integral)
    del _integral_cache[cachekey]
    return result

@evaluates(ConstantRule)
def eval_constant(constant, integrand, symbol):
    return constant * symbol

@evaluates(ConstantTimesRule)
def eval_constanttimes(constant, other, substep, integrand, symbol):
    return constant * _manualintegrate(substep)

@evaluates(PowerRule)
def eval_power(base, exp, integrand, symbol):
    return Piecewise(
        ((base**(exp + 1))/(exp + 1), Ne(exp, -1)),
        (log(base), True),
        )

@evaluates(ExpRule)
def eval_exp(base, exp, integrand, symbol):
    return integrand / log(base)

@evaluates(AddRule)
def eval_add(substeps, integrand, symbol):
    return sum(map(_manualintegrate, substeps))

@evaluates(URule)
def eval_u(u_var, u_func, constant, substep, integrand, symbol):
    result = _manualintegrate(substep)
    if u_func.is_Pow and u_func.exp == -1:
        # avoid needless -log(1/x) from substitution
        result = result.subs(log(u_var), -log(u_func.base))
    return result.subs(u_var, u_func)

@evaluates(PartsRule)
def eval_parts(u, dv, v_step, second_step, integrand, symbol):
    v = _manualintegrate(v_step)

    return u * v - _manualintegrate(second_step)

@evaluates(CyclicPartsRule)
def eval_cyclicparts(parts_rules, coefficient, integrand, symbol):
    coefficient = 1 - coefficient
    result = []

    sign = 1
    for rule in parts_rules:
        result.append(sign * rule.u * _manualintegrate(rule.v_step))
        sign *= -1

    return Add(*result) / coefficient

@evaluates(TrigRule)
def eval_trig(func, arg, integrand, symbol):
    if func == 'sin':
        return -cos(arg)
    elif func == 'cos':
        return sin(arg)
    elif func == 'sec*tan':
        return sec(arg)
    elif func == 'csc*cot':
        return csc(arg)
    elif func == 'sec**2':
        return tan(arg)
    elif func == 'csc**2':
        return -cot(arg)

@evaluates(ArctanRule)
def eval_arctan(a, b, c, integrand, symbol):
    return a / b * 1 / sqrt(c / b) * atan(symbol / sqrt(c / b))

@evaluates(ArccothRule)
def eval_arccoth(a, b, c, integrand, symbol):
    return - a / b * 1 / sqrt(-c / b) * acoth(symbol / sqrt(-c / b))

@evaluates(ArctanhRule)
def eval_arctanh(a, b, c, integrand, symbol):
    return - a / b * 1 / sqrt(-c / b) * atanh(symbol / sqrt(-c / b))

@evaluates(ReciprocalRule)
def eval_reciprocal(func, integrand, symbol):
    return log(func)

@evaluates(ArcsinRule)
def eval_arcsin(integrand, symbol):
    return asin(symbol)

@evaluates(InverseHyperbolicRule)
def eval_inversehyperbolic(func, integrand, symbol):
    return func(symbol)

@evaluates(AlternativeRule)
def eval_alternative(alternatives, integrand, symbol):
    return _manualintegrate(alternatives[0])

@evaluates(RewriteRule)
def eval_rewrite(rewritten, substep, integrand, symbol):
    return _manualintegrate(substep)

@evaluates(PiecewiseRule)
def eval_piecewise(substeps, integrand, symbol):
    return Piecewise(*[(_manualintegrate(substep), cond)
                             for substep, cond in substeps])

@evaluates(TrigSubstitutionRule)
def eval_trigsubstitution(theta, func, rewritten, substep, restriction, integrand, symbol):
    func = func.subs(sec(theta), 1/cos(theta))
    func = func.subs(csc(theta), 1/sin(theta))
    func = func.subs(cot(theta), 1/tan(theta))

    trig_function = list(func.find(TrigonometricFunction))
    assert len(trig_function) == 1
    trig_function = trig_function[0]
    relation = solve(symbol - func, trig_function)
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
    elif isinstance(trig_function, tan):
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
        (_manualintegrate(substep).subs(substitution).trigsimp(), restriction)
    )

@evaluates(DerivativeRule)
def eval_derivativerule(integrand, symbol):
    # isinstance(integrand, Derivative) should be True
    variable_count = list(integrand.variable_count)
    for i, (var, count) in enumerate(variable_count):
        if var == symbol:
            variable_count[i] = (var, count-1)
            break
    return Derivative(integrand.expr, *variable_count)

@evaluates(HeavisideRule)
def eval_heaviside(harg, ibnd, substep, integrand, symbol):
    # If we are integrating over x and the integrand has the form
    #       Heaviside(m*x+b)*g(x) == Heaviside(harg)*g(symbol)
    # then there needs to be continuity at -b/m == ibnd,
    # so we subtract the appropriate term.
    return Heaviside(harg)*(substep - substep.subs(symbol, ibnd))

@evaluates(JacobiRule)
def eval_jacobi(n, a, b, integrand, symbol):
    return Piecewise(
        (2*jacobi(n + 1, a - 1, b - 1, symbol)/(n + a + b), Ne(n + a + b, 0)),
        (symbol, Eq(n, 0)),
        ((a + b + 2)*symbol**2/4 + (a - b)*symbol/2, Eq(n, 1)))

@evaluates(GegenbauerRule)
def eval_gegenbauer(n, a, integrand, symbol):
    return Piecewise(
        (gegenbauer(n + 1, a - 1, symbol)/(2*(a - 1)), Ne(a, 1)),
        (chebyshevt(n + 1, symbol)/(n + 1), Ne(n, -1)),
        (S.Zero, True))

@evaluates(ChebyshevTRule)
def eval_chebyshevt(n, integrand, symbol):
    return Piecewise(((chebyshevt(n + 1, symbol)/(n + 1) -
        chebyshevt(n - 1, symbol)/(n - 1))/2, Ne(Abs(n), 1)),
        (symbol**2/2, True))

@evaluates(ChebyshevURule)
def eval_chebyshevu(n, integrand, symbol):
    return Piecewise(
        (chebyshevt(n + 1, symbol)/(n + 1), Ne(n, -1)),
        (S.Zero, True))

@evaluates(LegendreRule)
def eval_legendre(n, integrand, symbol):
    return (legendre(n + 1, symbol) - legendre(n - 1, symbol))/(2*n + 1)

@evaluates(HermiteRule)
def eval_hermite(n, integrand, symbol):
    return hermite(n + 1, symbol)/(2*(n + 1))

@evaluates(LaguerreRule)
def eval_laguerre(n, integrand, symbol):
    return laguerre(n, symbol) - laguerre(n + 1, symbol)

@evaluates(AssocLaguerreRule)
def eval_assoclaguerre(n, a, integrand, symbol):
    return -assoc_laguerre(n + 1, a - 1, symbol)

@evaluates(CiRule)
def eval_ci(a, b, integrand, symbol):
    return cos(b)*Ci(a*symbol) - sin(b)*Si(a*symbol)

@evaluates(ChiRule)
def eval_chi(a, b, integrand, symbol):
    return cosh(b)*Chi(a*symbol) + sinh(b)*Shi(a*symbol)

@evaluates(EiRule)
def eval_ei(a, b, integrand, symbol):
    return exp(b)*Ei(a*symbol)

@evaluates(SiRule)
def eval_si(a, b, integrand, symbol):
    return sin(b)*Ci(a*symbol) + cos(b)*Si(a*symbol)

@evaluates(ShiRule)
def eval_shi(a, b, integrand, symbol):
    return sinh(b)*Chi(a*symbol) + cosh(b)*Shi(a*symbol)

@evaluates(ErfRule)
def eval_erf(a, b, c, integrand, symbol):
    if a.is_extended_real:
        return Piecewise(
            (sqrt(S.Pi/(-a))/2 * exp(c - b**2/(4*a)) *
                erf((-2*a*symbol - b)/(2*sqrt(-a))), a < 0),
            (sqrt(S.Pi/a)/2 * exp(c - b**2/(4*a)) *
                erfi((2*a*symbol + b)/(2*sqrt(a))), True))
    else:
        return sqrt(S.Pi/a)/2 * exp(c - b**2/(4*a)) * \
                erfi((2*a*symbol + b)/(2*sqrt(a)))

@evaluates(FresnelCRule)
def eval_fresnelc(a, b, c, integrand, symbol):
    return sqrt(S.Pi/(2*a)) * (
        cos(b**2/(4*a) - c)*fresnelc((2*a*symbol + b)/sqrt(2*a*S.Pi)) +
        sin(b**2/(4*a) - c)*fresnels((2*a*symbol + b)/sqrt(2*a*S.Pi)))

@evaluates(FresnelSRule)
def eval_fresnels(a, b, c, integrand, symbol):
    return sqrt(S.Pi/(2*a)) * (
        cos(b**2/(4*a) - c)*fresnels((2*a*symbol + b)/sqrt(2*a*S.Pi)) -
        sin(b**2/(4*a) - c)*fresnelc((2*a*symbol + b)/sqrt(2*a*S.Pi)))

@evaluates(LiRule)
def eval_li(a, b, integrand, symbol):
    return li(a*symbol + b)/a

@evaluates(PolylogRule)
def eval_polylog(a, b, integrand, symbol):
    return polylog(b + 1, a*symbol)

@evaluates(UpperGammaRule)
def eval_uppergamma(a, e, integrand, symbol):
    return symbol**e * (-a*symbol)**(-e) * uppergamma(e + 1, -a*symbol)/a

@evaluates(EllipticFRule)
def eval_elliptic_f(a, d, integrand, symbol):
    return elliptic_f(symbol, d/a)/sqrt(a)

@evaluates(EllipticERule)
def eval_elliptic_e(a, d, integrand, symbol):
    return elliptic_e(symbol, d/a)*sqrt(a)

@evaluates(DontKnowRule)
def eval_dontknowrule(integrand, symbol):
    return Integral(integrand, symbol)

def _manualintegrate(rule):
    evaluator = evaluators.get(rule.__class__)
    if not evaluator:

        raise ValueError("Cannot evaluate rule %s" % repr(rule))
    return evaluator(*rule)

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
    RootSum(4*_z**2 + 1, Lambda(_i, _i*log(2*_i + exp(x))))
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
    result = _manualintegrate(integral_steps(f, var))
    # Clear the cache of u-parts
    _parts_u_cache.clear()
    # If we got Piecewise with two parts, put generic first
    if isinstance(result, Piecewise) and len(result.args) == 2:
        cond = result.args[0][1]
        if isinstance(cond, Eq) and result.args[1][1] == True:
            result = result.func(
                (result.args[1][0], Ne(*cond.args)),
                (result.args[0][0], True))
    return result
