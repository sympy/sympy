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
``@evaluates(namedtuple_type)``.

"""
from __future__ import print_function, division

from collections import namedtuple

import sympy

from sympy.core.compatibility import reduce
from sympy.functions.elementary.trigonometric import TrigonometricFunction
from sympy.simplify import fraction
from sympy.strategies.core import (switch, identity, do_one, null_safe,
                                   condition, tryit)

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
LogRule = Rule("LogRule", "func")
ArctanRule = Rule("ArctanRule")
AlternativeRule = Rule("AlternativeRule", "alternatives")
DontKnowRule = Rule("DontKnowRule")
DerivativeRule = Rule("DerivativeRule")
RewriteRule = Rule("RewriteRule", "rewritten substep")
PiecewiseRule = Rule("PiecewiseRule", "subfunctions")

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
        if isinstance(f, sympy.tan):
            return arg.diff(symbol) * sympy.sec(arg)**2
        elif isinstance(f, sympy.cot):
            return -arg.diff(symbol) * sympy.csc(arg)**2
        elif isinstance(f, sympy.sec):
            return arg.diff(symbol) * sympy.sec(arg) * sympy.tan(arg)
        elif isinstance(f, sympy.csc):
            return -arg.diff(symbol) * sympy.csc(arg) * sympy.cot(arg)
        elif isinstance(f, sympy.Add):
            return sum([manual_diff(arg, symbol) for arg in f.args])
    return f.diff(symbol)

# Method based on that on SIN, described in "Symbolic Integration: The
# Stormy Decade"

def find_substitutions(integrand, symbol, u_var):
    results = []

    def test_subterm(u, u_diff):
        substituted = integrand / u_diff
        if symbol not in substituted.free_symbols:
            # replaced everything already
            return False

        substituted = substituted.subs(u, u_var).cancel()
        if symbol not in substituted.free_symbols:
            _, denom = integrand.as_numer_denom()
            _, denom2 = substituted.as_numer_denom()
            denom2 = denom2.subs(u_var, u)
            constant = denom/denom2
            return (constant, substituted / constant)
        return False

    def possible_subterms(term):
        if any(isinstance(term, cls)
               for cls in (sympy.sin, sympy.cos, sympy.tan,
                           sympy.asin, sympy.acos, sympy.atan,
                           sympy.exp, sympy.log)):
            return [term.args[0]]
        elif isinstance(term, sympy.Mul):
            r = []
            for u in term.args:
                numer, denom = fraction(u)
                if numer == 1:
                    r.append(denom)
                    r.extend(possible_subterms(denom))
                else:
                    r.append(u)
                    r.extend(possible_subterms(u))
            return r
        elif isinstance(term, sympy.Pow):
            if term.args[1].is_constant(symbol):
                return [term.args[0]]
            elif term.args[0].is_constant(symbol):
                return [term.args[1]]
        elif isinstance(term, sympy.Add):
            return term.args
        return []

    for u in possible_subterms(integrand):
        if u == symbol:
            continue
        u_diff = manual_diff(u, symbol)
        new_integrand = test_subterm(u, u_diff)
        if new_integrand is not False:
            constant, new_integrand = new_integrand
            substitution = (u, constant, new_integrand)
            if substitution not in results:
                results.append(substitution)

    return results

def rewriter(condition, rewrite):
    """Strategy that rewrites an integrand."""
    def _rewriter(integral):
        integrand, symbol = integral
        if condition(*integral):
            rewritten = rewrite(*integral)
            if rewritten != integrand:
                substep = integral_steps(rewritten, symbol)
                if not isinstance(substep, DontKnowRule):
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
        for rule in rules:
            result = rule(integral)
            if (result and not isinstance(result, DontKnowRule) and
                result != integral and result not in alts):
                alts.append(result)
        if len(alts) == 1:
            return alts[0]
        elif len(alts) > 1:
            return AlternativeRule(alts, *integral)
    return _alternatives

def constant_rule(integral):
    integrand, symbol = integral
    return ConstantRule(integral.integrand, *integral)

def power_rule(integral):
    integrand, symbol = integral
    base, exp = integrand.as_base_exp()

    if symbol not in exp.free_symbols and isinstance(base, sympy.Symbol):
        if sympy.simplify(exp + 1) == 0:
            return LogRule(base, integrand, symbol)
        return PowerRule(base, exp, integrand, symbol)
    elif symbol not in base.free_symbols and isinstance(exp, sympy.Symbol):
        rule = ExpRule(base, exp, integrand, symbol)

        if sympy.ask(~sympy.Q.zero(sympy.log(base))):
            return rule
        elif sympy.ask(sympy.Q.zero(sympy.log(base))):
            return ConstantRule(1, 1, symbol)

        return PiecewiseRule([
            (ConstantRule(1, 1, symbol), sympy.Eq(sympy.log(base), 0)),
            (rule, True)
        ], integrand, symbol)

def exp_rule(integral):
    integrand, symbol = integral
    if isinstance(integrand.args[0], sympy.Symbol):
        return ExpRule(sympy.E, integrand.args[0], integrand, symbol)

def arctan_rule(integral):
    integrand, symbol = integral
    base, exp = integrand.as_base_exp()

    if sympy.simplify(exp + 1) == 0:
        a = sympy.Wild('a', exclude=[symbol])
        b = sympy.Wild('b', exclude=[symbol])
        match = base.match(a + b*symbol**2)
        if match:
            a, b = match[a], match[b]

            if ((isinstance(a, sympy.Number) and a < 0) or
                (isinstance(b, sympy.Number) and b < 0)):
                return
            if (sympy.ask(sympy.Q.negative(a) | sympy.Q.negative(b) |
                          sympy.Q.is_true(a <= 0) | sympy.Q.is_true(b <= 0))):
                return

            if a != 1 or b != 1:
                b_condition = b >= 0
                u_var = sympy.Dummy("u")
                rewritten = (sympy.Integer(1) / a) * (base / a) ** (-1)
                u_func = sympy.sqrt(sympy.sympify(b) / a) * symbol
                constant = 1 / sympy.sqrt(sympy.sympify(b) / a)
                substituted = rewritten.subs(u_func, u_var)

                if a == b:
                    substep = ArctanRule(integrand, symbol)
                else:
                    subrule = ArctanRule(substituted, u_var)
                    if constant != 1:
                        b_condition = b > 0
                        subrule = ConstantTimesRule(
                            constant, substituted, subrule,
                            substituted, symbol)

                    substep = URule(u_var, u_func, constant,
                                    subrule,
                                    integrand, symbol)

                if a != 1:
                    other = (base / a) ** (-1)
                    substep = ConstantTimesRule(
                        sympy.Integer(1) / a, other, substep,
                        integrand, symbol)

                return PiecewiseRule([
                    (substep, sympy.And(a > 0, b_condition))
                ], integrand, symbol)

            return ArctanRule(integrand, symbol)

def add_rule(integral):
    integrand, symbol = integral
    return AddRule(
        [integral_steps(g, symbol)
         for g in integrand.as_ordered_terms()],
        integrand, symbol)

def mul_rule(integral):
    integrand, symbol = integral
    args = integrand.args

    # Constant times function case
    coeff, f = integrand.as_independent(symbol)

    if coeff != 1:
        return ConstantTimesRule(
            coeff, f,
            integral_steps(f, symbol),
            integrand, symbol)

def _parts_rule(integrand, symbol):
    # LIATE rule:
    # log, inverse trig, algebraic (polynomial), trigonometric, exponential
    def pull_out_polys(integrand):
        integrand = integrand.together()
        polys = [arg for arg in integrand.args if arg.is_polynomial(symbol)]
        if polys:
            u = sympy.Mul(*polys)
            dv = integrand / u
            return u, dv

    def pull_out_u(*functions):
        def pull_out_u_rl(integrand):
            if any([integrand.has(f) for f in functions]):
                args = [arg for arg in integrand.args
                        if any(isinstance(arg, cls) for cls in functions)]
                if args:
                    u = reduce(lambda a,b: a*b, args)
                    dv = integrand / u
                    return u, dv

        return pull_out_u_rl

    liate_rules = [pull_out_u(sympy.log), pull_out_u(sympy.atan),
                   pull_out_polys, pull_out_u(sympy.sin, sympy.cos),
                   pull_out_u(sympy.exp)]


    dummy = sympy.Dummy("temporary")
    # we can integrate log(x) and atan(x) by setting dv = 1
    if isinstance(integrand, sympy.log) or isinstance(integrand, sympy.atan):
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

            for rule in liate_rules[index + 1:]:
                r = rule(integrand)
                # make sure dv is amenable to integration
                if r and r[0].subs(dummy, 1) == dv:
                    du = u.diff(symbol)
                    v_step = integral_steps(dv, symbol)
                    v = _manualintegrate(v_step)

                    return u, dv, v, du, v_step

def parts_rule(integral):
    integrand, symbol = integral
    constant, integrand = integrand.as_coeff_Mul()

    result = _parts_rule(integrand, symbol)

    steps = []
    if result:
        u, dv, v, du, v_step = result
        steps.append(result)

        if isinstance(v, sympy.Integral):
            return

        while True:
            if symbol not in (integrand / (v * du)).cancel().free_symbols:
                coefficient = ((v * du) / integrand).cancel()
                rule = CyclicPartsRule(
                    [PartsRule(u, dv, v_step, None, None, None)
                     for (u, dv, v, du, v_step) in steps],
                    (-1) ** len(steps) * coefficient,
                    integrand, symbol
                )
                if constant != 1:
                    rule = ConstantTimesRule(constant, integrand, rule,
                                             constant * integrand, symbol)
                return rule

            result = _parts_rule(v * du, symbol)

            if result:
                u, dv, v, du, v_step = result
                steps.append(result)
            else:
                break

    def make_second_step(steps, integrand):
        if steps:
            u, dv, v, du, v_step = steps[0]
            return PartsRule(u, dv, v_step,
                             make_second_step(steps[1:], v * du),
                             integrand, symbol)
        else:
            return integral_steps(integrand, symbol)

    if steps:
        u, dv, v, du, v_step = steps[0]
        rule = PartsRule(u, dv, v_step,
                         make_second_step(steps[1:], v * du),
                         integrand, symbol)
        if constant != 1:
            rule = ConstantTimesRule(constant, integrand, rule,
                                     constant * integrand, symbol)
        return rule


def trig_rule(integral):
    integrand, symbol = integral
    if isinstance(integrand, sympy.sin) or isinstance(integrand, sympy.cos):
        arg = integrand.args[0]

        if not isinstance(arg, sympy.Symbol):
            return  # perhaps a substitution can deal with it

        if isinstance(integrand, sympy.sin):
            func = 'sin'
        else:
            func = 'cos'

        return TrigRule(func, arg, integrand, symbol)

    if isinstance(integrand, sympy.tan):
        rewritten = sympy.sin(*integrand.args) / sympy.cos(*integrand.args)
    elif isinstance(integrand, sympy.cot):
        rewritten = sympy.cos(*integrand.args) / sympy.sin(*integrand.args)
    elif isinstance(integrand, sympy.sec):
        arg = integrand.args[0]
        rewritten = ((sympy.sec(arg)**2 + sympy.tan(arg) * sympy.sec(arg)) /
                     (sympy.sec(arg) + sympy.tan(arg)))
    elif isinstance(integrand, sympy.csc):
        arg = integrand.args[0]
        rewritten = ((sympy.csc(arg)**2 + sympy.cot(arg) * sympy.csc(arg)) /
                     (sympy.csc(arg) + sympy.cot(arg)))
    return RewriteRule(
        rewritten,
        integral_steps(rewritten, symbol),
        integrand, symbol
    )

def trig_product_rule(integral):
    integrand, symbol = integral

    sectan = sympy.sec(symbol) * sympy.tan(symbol)
    q = integrand / sectan

    if symbol not in q.free_symbols:
        rule = TrigRule('sec*tan', symbol, sectan, symbol)
        if q != 1:
            rule = ConstantTimesRule(q, sectan, rule, integrand, symbol)

        return rule

    csccot = -sympy.csc(symbol) * sympy.cot(symbol)
    q = integrand / csccot

    if symbol not in q.free_symbols:
        rule = TrigRule('csc*cot', symbol, csccot, symbol)
        if q != 1:
            rule = ConstantTimesRule(q, csccot, rule, integrand, symbol)

        return rule

@sympy.cacheit
def make_wilds(symbol):
    a = sympy.Wild('a', exclude=[symbol])
    b = sympy.Wild('b', exclude=[symbol])
    m = sympy.Wild('m', exclude=[symbol], properties=[lambda n: isinstance(n, sympy.Integer)])
    n = sympy.Wild('n', exclude=[symbol], properties=[lambda n: isinstance(n, sympy.Integer)])

    return a, b, m, n

@sympy.cacheit
def sincos_pattern(symbol):
    a, b, m, n = make_wilds(symbol)
    pattern = sympy.sin(a*symbol)**m * sympy.cos(b*symbol)**n

    return pattern, a, b, m, n

@sympy.cacheit
def tansec_pattern(symbol):
    a, b, m, n = make_wilds(symbol)
    pattern = sympy.tan(a*symbol)**m * sympy.sec(b*symbol)**n

    return pattern, a, b, m, n

@sympy.cacheit
def cotcsc_pattern(symbol):
    a, b, m, n = make_wilds(symbol)
    pattern = sympy.cot(a*symbol)**m * sympy.csc(b*symbol)**n

    return pattern, a, b, m, n

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

sincos_botheven_condition = uncurry(lambda a, b, m, n, i, s: m.is_even and n.is_even)

sincos_botheven = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (((1 - sympy.cos(2*a*symbol)) / 2) ** (m / 2)) *
                                    (((1 + sympy.cos(2*b*symbol)) / 2) ** (n / 2)) ))

sincos_sinodd_condition = uncurry(lambda a, b, m, n, i, s: m.is_odd and m >= 3)

sincos_sinodd = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (1 - sympy.cos(a*symbol)**2)**((m - 1) / 2) *
                                    sympy.sin(a*symbol) *
                                    sympy.cos(b*symbol) ** n))

sincos_cosodd_condition = uncurry(lambda a, b, m, n, i, s: n.is_odd and n >= 3)

sincos_cosodd = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (1 - sympy.sin(b*symbol)**2)**((n - 1) / 2) *
                                    sympy.cos(b*symbol) *
                                    sympy.sin(a*symbol) ** m))

tansec_seceven_condition = uncurry(lambda a, b, m, n, i, s: n.is_even and n >= 4)
tansec_seceven = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (1 + sympy.tan(b*symbol)**2) ** (n/2 - 1) *
                                    sympy.sec(b*symbol)**2 *
                                    sympy.tan(a*symbol) ** m ))

tansec_tanodd_condition = uncurry(lambda a, b, m, n, i, s: m.is_odd)
tansec_tanodd = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (sympy.sec(a*symbol)**2 - 1) ** ((m - 1) / 2) *
                                     sympy.tan(a*symbol) *
                                     sympy.sec(b*symbol) ** n ))

cotcsc_csceven_condition = uncurry(lambda a, b, m, n, i, s: n.is_even and n >= 4)
cotcsc_csceven = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (1 + sympy.cot(b*symbol)**2) ** (n/2 - 1) *
                                    sympy.csc(b*symbol)**2 *
                                    sympy.cot(a*symbol) ** m ))

cotcsc_cotodd_condition = uncurry(lambda a, b, m, n, i, s: m.is_odd)
cotcsc_cotodd = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( (sympy.csc(a*symbol)**2 - 1) ** ((m - 1) / 2) *
                                    sympy.cot(a*symbol) *
                                    sympy.csc(b*symbol) ** n ))

def trig_powers_products_rule(integral):
    integrand, symbol = integral

    if any(integrand.has(f) for f in (sympy.sin, sympy.cos)):
        pattern, a, b, m, n = sincos_pattern(symbol)
        match = integrand.match(pattern)

        if match:
            a, b, m, n = match.get(a, 0),match.get(b, 0), match.get(m, 0), match.get(n, 0)
            return multiplexer({
                sincos_botheven_condition: sincos_botheven,
                sincos_sinodd_condition: sincos_sinodd,
                sincos_cosodd_condition: sincos_cosodd
            })((a, b, m, n, integrand, symbol))

    integrand = integrand.subs({
        1 / sympy.cos(symbol): sympy.sec(symbol)
    })

    if any(integrand.has(f) for f in (sympy.tan, sympy.sec)):
        pattern, a, b, m, n = tansec_pattern(symbol)
        match = integrand.match(pattern)

        if match:
            a, b, m, n = match.get(a, 0),match.get(b, 0), match.get(m, 0), match.get(n, 0)
            return multiplexer({
                tansec_tanodd_condition: tansec_tanodd,
                tansec_seceven_condition: tansec_seceven
            })((a, b, m, n, integrand, symbol))

    integrand = integrand.subs({
        1 / sympy.sin(symbol): sympy.csc(symbol),
        1 / sympy.tan(symbol): sympy.cot(symbol),
        sympy.cos(symbol) / sympy.tan(symbol): sympy.cot(symbol)
    })

    if any(integrand.has(f) for f in (sympy.cot, sympy.csc)):
        pattern, a, b, m, n = cotcsc_pattern(symbol)
        match = integrand.match(pattern)

        if match:
            a, b, m, n = match.get(a, 0),match.get(b, 0), match.get(m, 0), match.get(n, 0)
            return multiplexer({
                cotcsc_cotodd_condition: cotcsc_cotodd,
                cotcsc_csceven_condition: cotcsc_csceven
            })((a, b, m, n, integrand, symbol))

def substitution_rule(integral):
    integrand, symbol = integral

    u_var = sympy.Dummy("u")
    substitutions = find_substitutions(integrand, symbol, u_var)
    if substitutions:
        ways = []
        for u_func, c, substituted in substitutions:
            subrule = integral_steps(substituted, u_var)
            if contains_dont_know(subrule):
                continue

            if sympy.simplify(c - 1) != 0:
                _, denom = c.as_numer_denom()
                subrule = ConstantTimesRule(c, substituted, subrule, substituted, symbol)

                if denom.free_symbols:
                    piecewise = []
                    could_be_zero = []

                    if isinstance(denom, sympy.Mul):
                        could_be_zero = denom.args
                    else:
                        could_be_zero.append(denom)

                    for expr in could_be_zero:
                        if not sympy.ask(~sympy.Q.zero(expr)):
                            substep = integral_steps(integrand.subs(expr, 0), symbol)

                            if substep:
                                piecewise.append((
                                    substep,
                                    sympy.Eq(expr, 0)
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

    elif integrand.has(sympy.exp):
        u_func = sympy.exp(symbol)
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

distribute_expand_rule = rewriter(
    lambda integrand, symbol: (
        all(arg.is_Pow or arg.is_polynomial(symbol) for arg in integrand.args)
        or isinstance(integrand, sympy.Pow)
        or isinstance(integrand, sympy.Mul)),
    lambda integrand, symbol: integrand.expand())

def derivative_rule(integral):
    variables = integral[0].args[1:]

    if variables[-1] == integral.symbol:
        return DerivativeRule(*integral)
    else:
        return ConstantRule(integral.integrand, *integral)

def fallback_rule(integral):
    return DontKnowRule(*integral)

# Cache is used to break cyclic integrals
_integral_cache = {}
def integral_steps(integrand, symbol, **options):
    """Returns the steps needed to compute an integral.

    This function attempts to mirror what a student would do by hand as
    closely as possible.

    SymPy Gamma uses this to provide a step-by-step explanation of an
    integral. The code it uses to format the results of this function can be
    found at
    https://github.com/sympy/sympy_gamma/blob/master/app/logic/intsteps.py.

    Examples
    ========

    >>> from sympy import exp, sin, cos
    >>> from sympy.integrals.manualintegrate import integral_steps
    >>> from sympy.abc import x
    >>> print(repr(integral_steps(exp(x) / (1 + exp(2 * x)), x))) \
    # doctest: +NORMALIZE_WHITESPACE
    URule(u_var=_u, u_func=exp(x), constant=1,
        substep=ArctanRule(context=1/(_u**2 + 1), symbol=_u),
        context=exp(x)/(exp(2*x) + 1), symbol=x)
    >>> print(repr(integral_steps(sin(x), x))) \
    # doctest: +NORMALIZE_WHITESPACE
    TrigRule(func='sin', arg=x, context=sin(x), symbol=x)
    >>> print(repr(integral_steps((x**2 + 3)**2 , x))) \
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
    cachekey = (integrand, symbol)
    if cachekey in _integral_cache:
        if _integral_cache[cachekey] is None:
            # cyclic integral! null_safe will eliminate that path
            return None
        else:
            return _integral_cache[cachekey]
    else:
        _integral_cache[cachekey] = None

    integral = IntegralInfo(integrand, symbol)

    def key(integral):
        integrand = integral.integrand

        if isinstance(integrand, TrigonometricFunction):
            return TrigonometricFunction
        elif isinstance(integrand, sympy.Derivative):
            return sympy.Derivative
        elif symbol not in integrand.free_symbols:
            return sympy.Number
        else:
            for cls in (sympy.Pow, sympy.Symbol, sympy.exp, sympy.log,
                        sympy.Add, sympy.Mul, sympy.atan):
                if isinstance(integrand, cls):
                    return cls

    def integral_is_subclass(*klasses):
        def _integral_is_subclass(integral):
            k = key(integral)
            return k and issubclass(k, klasses)
        return _integral_is_subclass

    result = do_one(
        null_safe(switch(key, {
            sympy.Pow: do_one(null_safe(power_rule), null_safe(arctan_rule)),
            sympy.Symbol: power_rule,
            sympy.exp: exp_rule,
            sympy.Add: add_rule,
            sympy.Mul: do_one(null_safe(mul_rule), null_safe(trig_product_rule)),
            sympy.Derivative: derivative_rule,
            TrigonometricFunction: trig_rule,
            sympy.Number: constant_rule
        })),
        null_safe(
            alternatives(
                substitution_rule,
                condition(
                    integral_is_subclass(sympy.Mul, sympy.log, sympy.atan),
                    parts_rule),
                condition(
                    integral_is_subclass(sympy.Mul, sympy.Pow),
                    partial_fractions_rule),
                condition(
                    integral_is_subclass(sympy.Mul, sympy.Pow),
                    distribute_expand_rule),
                trig_powers_products_rule
            )
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
    return (base ** (exp + 1)) / (exp + 1)

@evaluates(ExpRule)
def eval_exp(base, exp, integrand, symbol):
    return integrand / sympy.ln(base)

@evaluates(AddRule)
def eval_add(substeps, integrand, symbol):
    return sum(map(_manualintegrate, substeps))

@evaluates(URule)
def eval_u(u_var, u_func, constant, substep, integrand, symbol):
    result = _manualintegrate(substep)
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

    return sympy.Add(*result) / coefficient

@evaluates(TrigRule)
def eval_trig(func, arg, integrand, symbol):
    if func == 'sin':
        return -sympy.cos(arg)
    elif func == 'cos':
        return sympy.sin(arg)
    elif func == 'sec*tan':
        return sympy.sec(arg)
    elif func == 'csc*cot':
        return sympy.csc(arg)

@evaluates(LogRule)
def eval_log(func, integrand, symbol):
    return sympy.ln(func)

@evaluates(ArctanRule)
def eval_arctan(integrand, symbol):
    return sympy.atan(symbol)

@evaluates(AlternativeRule)
def eval_alternative(alternatives, integrand, symbol):
    return _manualintegrate(alternatives[0])

@evaluates(RewriteRule)
def eval_rewrite(rewritten, substep, integrand, symbol):
    return _manualintegrate(substep)

@evaluates(PiecewiseRule)
def eval_piecewise(substeps, integrand, symbol):
    return sympy.Piecewise(*[(_manualintegrate(substep), cond)
                             for substep, cond in substeps])

@evaluates(DerivativeRule)
def eval_derivativerule(integrand, symbol):
    # isinstance(integrand, Derivative) should be True
    if len(integrand.args) == 2:
        return integrand.args[0]
    else:
        return sympy.Derivative(integrand.args[0], *integrand.args[1:-1])

@evaluates(DontKnowRule)
def eval_dontknowrule(integrand, symbol):
    return sympy.Integral(integrand, symbol)

def _manualintegrate(rule):
    evaluator = evaluators.get(rule.__class__)
    if not evaluator:
        raise ValueError("Cannot evaluate rule %s" % rule)
    return evaluator(*rule)

def manualintegrate(f, var):
    """manualintegrate(f, var)

    Compute indefinite integral of a single variable using an algorithm that
    resembles what a student would do by hand.

    Unlike ``integrate``, var can only be a single symbol.

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
    -log(sin(x)**2 - 1)/2

    See Also
    ========

    sympy.integrals.integrals.integrate
    sympy.integrals.integrals.Integral.doit
    sympy.integrals.integrals.Integral
    """
    return _manualintegrate(integral_steps(f, var))
