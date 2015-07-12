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

from __future__ import division, print_function

from collections import namedtuple

import sympy

from sympy.core.compatibility import reduce
from sympy.functions.elementary.trigonometric import TrigonometricFunction
from sympy.strategies.core import switch, do_one, null_safe, condition

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
ArctanRule = Rule("ArctanRule")
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
        elif isinstance(f, sympy.Mul):
            if len(f.args) == 2 and isinstance(f.args[0], sympy.Number):
                return f.args[0] * manual_diff(f.args[1], symbol)
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
            return substituted.as_independent(u_var, as_Add=False)

        return False

    def possible_subterms(term):
        if isinstance(term, (TrigonometricFunction,
                             sympy.asin, sympy.acos, sympy.atan,
                             sympy.exp, sympy.log, sympy.Heaviside)):
            return [term.args[0]]
        elif isinstance(term, sympy.Mul):
            r = []
            for u in term.args:
                r.append(u)
                r.extend(possible_subterms(u))
            return r
        elif isinstance(term, sympy.Pow):
            if term.args[1].is_constant(symbol):
                return [term.args[0]]
            elif term.args[0].is_constant(symbol):
                return [term.args[1]]
        elif isinstance(term, sympy.Add):
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
        elif alts:
            doable = [rule for rule in alts if not contains_dont_know(rule)]
            if doable:
                return AlternativeRule(doable, *integral)
            else:
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
            return ReciprocalRule(base, integrand, symbol)
        return PowerRule(base, exp, integrand, symbol)
    elif symbol not in base.free_symbols and isinstance(exp, sympy.Symbol):
        rule = ExpRule(base, exp, integrand, symbol)

        if sympy.log(base).is_nonzero:
            return rule
        elif sympy.log(base).is_zero:
            return ConstantRule(1, 1, symbol)

        return PiecewiseRule([
            (ConstantRule(1, 1, symbol), sympy.Eq(sympy.log(base), 0)),
            (rule, True)
        ], integrand, symbol)

def exp_rule(integral):
    integrand, symbol = integral
    if isinstance(integrand.args[0], sympy.Symbol):
        return ExpRule(sympy.E, integrand.args[0], integrand, symbol)

def inverse_trig_rule(integral):
    integrand, symbol = integral
    base, exp = integrand.as_base_exp()
    a = sympy.Wild('a', exclude=[symbol])
    b = sympy.Wild('b', exclude=[symbol])
    match = base.match(a + b*symbol**2)

    if not match:
        return

    def negative(x):
        return x.is_negative or x.could_extract_minus_sign()

    def ArcsinhRule(integrand, symbol):
        return InverseHyperbolicRule(sympy.asinh, integrand, symbol)

    def ArccoshRule(integrand, symbol):
        return InverseHyperbolicRule(sympy.acosh, integrand, symbol)

    def make_inverse_trig(RuleClass, base_exp, a, sign_a, b, sign_b):
        u_var = sympy.Dummy("u")
        current_base = base
        current_symbol = symbol
        constant = u_func = u_constant = substep = None
        factored = integrand
        if a != 1:
            constant = a**base_exp
            current_base = sign_a + sign_b * (b/a) * current_symbol**2
            factored = current_base ** base_exp
        if (b/a) != 1:
            u_func = sympy.sqrt(b/a) * symbol
            u_constant = sympy.sqrt(a/b)
            current_symbol = u_var
            current_base = sign_a + sign_b * current_symbol**2

        substep = RuleClass(current_base ** base_exp, current_symbol)
        if u_func is not None:
            if u_constant != 1:
                substep = ConstantTimesRule(
                    u_constant, current_base ** base_exp, substep,
                    u_constant * current_base ** base_exp, symbol)
            substep = URule(u_var, u_func, u_constant, substep, factored, symbol)
        if constant is not None:
            substep = ConstantTimesRule(constant, factored, substep, integrand, symbol)
        return substep

    a, b = match[a], match[b]

    # list of (rule, base_exp, a, sign_a, b, sign_b, condition)
    possibilities = []

    if sympy.simplify(exp + 1) == 0 and not (negative(a) or negative(b)):
        possibilities.append((ArctanRule, exp, a, 1, b, 1, sympy.And(a > 0, b > 0)))
    elif sympy.simplify(2*exp + 1) == 0:
        possibilities.append((ArcsinRule, exp, a, 1, -b, -1, sympy.And(a > 0, b < 0)))
        possibilities.append((ArcsinhRule, exp, a, 1, b, 1, sympy.And(a > 0, b > 0)))
        possibilities.append((ArccoshRule, exp, -a, -1, b, 1, sympy.And(a < 0, b > 0)))

    possibilities = [p for p in possibilities if p[-1] is not sympy.false]
    if a.is_number and b.is_number:
        possibility = [p for p in possibilities if p[-1] is sympy.true]
        if len(possibility) == 1:
            return make_inverse_trig(*possibility[0][:-1])
    elif possibilities:
        return PiecewiseRule(
            [(make_inverse_trig(*p[:-1]), p[-1]) for p in possibilities],
            integrand, symbol)

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

    liate_rules = [pull_out_u(sympy.log), pull_out_u(sympy.atan, sympy.asin, sympy.acos),
                   pull_out_polys, pull_out_u(sympy.sin, sympy.cos),
                   pull_out_u(sympy.exp)]


    dummy = sympy.Dummy("temporary")
    # we can integrate log(x) and atan(x) by setting dv = 1
    if isinstance(integrand, (sympy.log, sympy.atan, sympy.asin, sympy.acos)):
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

    if integrand == sympy.sec(symbol)**2:
        return TrigRule('sec**2', symbol, integrand, symbol)
    elif integrand == sympy.csc(symbol)**2:
        return TrigRule('csc**2', symbol, integrand, symbol)

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
    else:
        return

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

@sympy.cacheit
def heaviside_pattern(symbol):
    m = sympy.Wild('m', exclude=[symbol])
    b = sympy.Wild('b', exclude=[symbol])
    g = sympy.Wild('g')
    pattern = sympy.Heaviside(m*symbol + b) * g

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

tan_tansquared_condition = uncurry(lambda a, b, m, n, i, s: m == 2 and n == 0)
tan_tansquared = trig_rewriter(
    lambda a, b, m, n, i, symbol: ( sympy.sec(a*symbol)**2 - 1))

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

def trig_sincos_rule(integral):
    integrand, symbol = integral

    if any(integrand.has(f) for f in (sympy.sin, sympy.cos)):
        pattern, a, b, m, n = sincos_pattern(symbol)
        match = integrand.match(pattern)

        if match:
            a, b, m, n = match.get(a, 0), match.get(b, 0), match.get(m, 0), match.get(n, 0)
            return multiplexer({
                sincos_botheven_condition: sincos_botheven,
                sincos_sinodd_condition: sincos_sinodd,
                sincos_cosodd_condition: sincos_cosodd
            })((a, b, m, n, integrand, symbol))

def trig_tansec_rule(integral):
    integrand, symbol = integral

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
                tansec_seceven_condition: tansec_seceven,
                tan_tansquared_condition: tan_tansquared
            })((a, b, m, n, integrand, symbol))

def trig_cotcsc_rule(integral):
    integrand, symbol = integral
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

def trig_powers_products_rule(integral):
    return do_one(null_safe(trig_sincos_rule),
                  null_safe(trig_tansec_rule),
                  null_safe(trig_cotcsc_rule))(integral)

def trig_substitution_rule(integral):
    integrand, symbol = integral
    a = sympy.Wild('a', exclude=[0, symbol])
    b = sympy.Wild('b', exclude=[0, symbol])
    theta = sympy.Dummy("theta")

    matches = integrand.find(a + b*symbol**2)
    if matches:
        for expr in matches:
            match = expr.match(a + b*symbol**2)
            a = match[a]
            b = match[b]

            a_positive = ((a.is_number and a > 0) or a.is_positive)
            b_positive = ((b.is_number and b > 0) or b.is_positive)
            x_func = None
            if a_positive and b_positive:
                # a**2 + b*x**2. Assume sec(theta) > 0, -pi/2 < theta < pi/2
                x_func = (sympy.sqrt(a)/sympy.sqrt(b)) * sympy.tan(theta)
                # Do not restrict the domain: tan(theta) takes on any real
                # value on the interval -pi/2 < theta < pi/2 so x takes on
                # any value
                restriction = True
            elif a_positive and not b_positive:
                # a**2 - b*x**2. Assume cos(theta) > 0, -pi/2 < theta < pi/2
                constant = sympy.sqrt(a)/sympy.sqrt(-b)
                x_func = constant * sympy.sin(theta)
                restriction = sympy.And(symbol > -constant, symbol < constant)
            elif not a_positive and b_positive:
                # b*x**2 - a**2. Assume sin(theta) > 0, 0 < theta < pi
                constant = sympy.sqrt(-a)/sympy.sqrt(b)
                x_func = constant * sympy.sec(theta)
                restriction = sympy.And(symbol > -constant, symbol < constant)
            if x_func:
                # Manually simplify sqrt(trig(theta)**2) to trig(theta)
                # Valid due to assumed domain restriction
                substitutions = {}
                for f in [sympy.sin, sympy.cos, sympy.tan,
                          sympy.sec, sympy.csc, sympy.cot]:
                    substitutions[sympy.sqrt(f(theta)**2)] = f(theta)
                    substitutions[sympy.sqrt(f(theta)**(-2))] = 1/f(theta)

                replaced = integrand.subs(symbol, x_func).trigsimp()
                replaced = replaced.subs(substitutions)
                if not replaced.has(symbol):
                    replaced *= manual_diff(x_func, theta)
                    replaced = replaced.trigsimp()
                    secants = replaced.find(1/sympy.cos(theta))
                    if secants:
                        replaced = replaced.xreplace({
                            1/sympy.cos(theta): sympy.sec(theta)
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
                subrule = ConstantTimesRule(c, substituted, subrule, substituted, u_var)

                if denom.free_symbols:
                    piecewise = []
                    could_be_zero = []

                    if isinstance(denom, sympy.Mul):
                        could_be_zero = denom.args
                    else:
                        could_be_zero.append(denom)

                    for expr in could_be_zero:
                        if not expr.is_nonzero:
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

def rewrites_rule(integral):
    integrand, symbol = integral

    if integrand.match(1/sympy.cos(symbol)):
        rewritten = integrand.subs(1/sympy.cos(symbol), sympy.sec(symbol))
        return RewriteRule(rewritten, integral_steps(rewritten, symbol), integrand, symbol)

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
                        sympy.Add, sympy.Mul, sympy.atan, sympy.asin, sympy.acos, sympy.Heaviside):
                if isinstance(integrand, cls):
                    return cls

    def integral_is_subclass(*klasses):
        def _integral_is_subclass(integral):
            k = key(integral)
            return k and issubclass(k, klasses)
        return _integral_is_subclass

    result = do_one(
        null_safe(switch(key, {
            sympy.Pow: do_one(null_safe(power_rule), null_safe(inverse_trig_rule)),
            sympy.Symbol: power_rule,
            sympy.exp: exp_rule,
            sympy.Add: add_rule,
            sympy.Mul: do_one(null_safe(mul_rule), null_safe(trig_product_rule), \
                null_safe(heaviside_rule)),
            sympy.Derivative: derivative_rule,
            TrigonometricFunction: trig_rule,
            sympy.Heaviside: heaviside_rule,
            sympy.Number: constant_rule
        })),
        do_one(
            null_safe(trig_rule),
            null_safe(alternatives(
                rewrites_rule,
                substitution_rule,
                condition(
                    integral_is_subclass(sympy.Mul, sympy.Pow),
                    partial_fractions_rule),
                condition(
                    integral_is_subclass(sympy.Mul, sympy.log, sympy.atan, sympy.asin, sympy.acos),
                    parts_rule),
                condition(
                    integral_is_subclass(sympy.Mul, sympy.Pow),
                    distribute_expand_rule),
                trig_powers_products_rule
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
    elif func == 'sec**2':
        return sympy.tan(arg)
    elif func == 'csc**2':
        return -sympy.cot(arg)

@evaluates(ReciprocalRule)
def eval_reciprocal(func, integrand, symbol):
    return sympy.ln(func)

@evaluates(ArctanRule)
def eval_arctan(integrand, symbol):
    return sympy.atan(symbol)

@evaluates(ArcsinRule)
def eval_arcsin(integrand, symbol):
    return sympy.asin(symbol)

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
    return sympy.Piecewise(*[(_manualintegrate(substep), cond)
                             for substep, cond in substeps])

@evaluates(TrigSubstitutionRule)
def eval_trigsubstitution(theta, func, rewritten, substep, restriction, integrand, symbol):
    func = func.subs(sympy.sec(theta), 1/sympy.cos(theta))

    trig_function = list(func.find(TrigonometricFunction))
    assert len(trig_function) == 1
    trig_function = trig_function[0]
    relation = sympy.solve(symbol - func, trig_function)
    assert len(relation) == 1
    numer, denom = sympy.fraction(relation[0])

    if isinstance(trig_function, sympy.sin):
        opposite = numer
        hypotenuse = denom
        adjacent = sympy.sqrt(denom**2 - numer**2)
        inverse = sympy.asin(relation[0])
    elif isinstance(trig_function, sympy.cos):
        adjacent = numer
        hypotenuse = denom
        opposite = sympy.sqrt(denom**2 - numer**2)
        inverse = sympy.acos(relation[0])
    elif isinstance(trig_function, sympy.tan):
        opposite = numer
        adjacent = denom
        hypotenuse = sympy.sqrt(denom**2 + numer**2)
        inverse = sympy.atan(relation[0])

    substitution = [
        (sympy.sin(theta), opposite/hypotenuse),
        (sympy.cos(theta), adjacent/hypotenuse),
        (sympy.tan(theta), opposite/adjacent),
        (theta, inverse)
    ]
    return sympy.Piecewise(
        (_manualintegrate(substep).subs(substitution).trigsimp(), restriction)
    )

@evaluates(DerivativeRule)
def eval_derivativerule(integrand, symbol):
    # isinstance(integrand, Derivative) should be True
    if len(integrand.args) == 2:
        return integrand.args[0]
    else:
        return sympy.Derivative(integrand.args[0], *integrand.args[1:-1])

@evaluates(HeavisideRule)
def eval_heaviside(harg, ibnd, substep, integrand, symbol):
    # If we are integrating over x and the integrand has the form
    #       Heaviside(m*x+b)*g(x) == Heaviside(harg)*g(symbol)
    # then there needs to be continuity at -b/m == ibnd,
    # so we subtract the appropriate term.
    return sympy.Heaviside(harg)*(substep - substep.subs(symbol, ibnd))

@evaluates(DontKnowRule)
def eval_dontknowrule(integrand, symbol):
    return sympy.Integral(integrand, symbol)

def _manualintegrate(rule):
    evaluator = evaluators.get(rule.__class__)
    if not evaluator:
        raise ValueError("Cannot evaluate rule %s" % repr(rule))
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
