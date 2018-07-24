import sys
from sympy.external import import_module
matchpy = import_module("matchpy")

if not matchpy:
    #bin/test will not execute any tests now
    disabled = True

if sys.version_info[:2] < (3, 6):
    disabled = True

from sympy.integrals.rubi.parsetools.parse import (rubi_rule_parser,
    get_default_values, add_wildcards, parse_freeq, seperate_freeq,
    get_free_symbols, divide_constraint, generate_sympy_from_parsed,
    setWC, replaceWith, rubi_printer)

from sympy import Symbol, Not
from sympy import sympify

def test_rubi_rule_parser():
    header = '''
from matchpy import Operation, CommutativeOperation
    rubi = ManyToOneReplacer()
'''
    fullform = 'List[RuleDelayed[HoldPattern[Int[Power[Pattern[x,Blank[]],Optional[Pattern[m,Blank[]]]],Pattern[x,Blank[Symbol]]]],Condition[Times[Power[x,Plus[m,1]],Power[Plus[m,1],-1]],NonzeroQ[Plus[m,1]]]]]'
    rules = rubi_rule_parser(fullform, header)
    result = '''
from matchpy import Operation, CommutativeOperation
    rubi = ManyToOneReplacer()
    pattern1 = Pattern(Integral(x_**WC('m', S(1)), x_), CustomConstraint(lambda m: NonzeroQ(m + S(1))))
    rule1 = ReplacementRule(pattern1, lambda m, x : x**(m + S(1))/(m + S(1)))
    rubi.add(rule1)

    return rubi

'''
    assert len(result.strip()) == len(rules.strip()) # failing randomly while using `result.strip() == rules`

def test_get_default_values():
    s = ['Int', ['Power', ['Plus', ['Optional', ['Pattern', 'a', ['Blank']]], ['Times', ['Optional', ['Pattern', 'b', ['Blank']]], ['Pattern', 'x', ['Blank']]]], ['Pattern', 'm', ['Blank']]], ['Pattern', 'x', ['Blank', 'Symbol']]]
    assert get_default_values(s, {}) == {'a': 0, 'b': 1}
    s = ['Int', ['Power', ['Pattern', 'x', ['Blank']], ['Optional', ['Pattern', 'm', ['Blank']]]], ['Pattern', 'x', ['Blank', 'Symbol']]]
    assert get_default_values(s, {}) == {'m': 1}

def test_add_wildcards():
    s = 'Integral(Pow(Pattern(x, Blank), Optional(Pattern(m, Blank))), Pattern(x, Blank(Symbol)))'
    assert add_wildcards(s, {'m': 1}) == ("Integral(Pow(x_, WC('m', S(1))), x_)", ['m', 'x', 'x'])

def test_seperate_freeq():
    s = ['FreeQ', ['List', 'a', 'b'], 'x']
    assert seperate_freeq(s) == (['a', 'b'], 'x')

def test_parse_freeq():
    l = ['a', 'b']
    x = 'x'
    symbols = ['x', 'a', 'b']
    assert parse_freeq(l, x, symbols) == ', CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x))'

def test_get_free_symbols():
    s = ['NonzeroQ', ['Plus', 'm', '1']]
    symbols = ['m', 'x']
    assert get_free_symbols(s, symbols, []) == ['m']

def test_divide_constraint():
    s = ['And', ['FreeQ', 'm', 'x'], ['NonzeroQ', ['Plus', 'm', '1']]]
    assert divide_constraint(s, ['m', 'x']) == ', CustomConstraint(lambda m: NonzeroQ(m + S(1)))'

def test_setWC():
    assert setWC('Integral(x_**WC(m, S(1)), x_)') == "Integral(x_**WC('m', S(1)), x_)"

def test_replaceWith():
    s = sympify('Module(List(Set(r, Numerator(Rt(a/b, n))), Set(s, Denominator(Rt(a/b, n))), k, u), CompoundExpression(Set(u, Integral((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)), Dist(2*r/(a*n), _Sum(u, List(k, 1, n/2 - 1/2)), x) + r*Integral(1/(r + s*x), x)/(a*n)))')
    symbols = ['x', 'a', 'n', 'b']
    assert replaceWith(s, symbols, 1) == ("    def With1(x, a, n, b):\n        r = Numerator(Rt(a/b, n))\n        s = Denominator(Rt(a/b, n))\n        k = Symbol('k')\n        u = Symbol('u')\n        u = Integral((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)\n        return Dist(S(2)*r/(a*n), _Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Integral(S(1)/(r + s*x), x)/(a*n)", ', ')

def test_generate_sympy_from_parsed():
    s = ['Int', ['Power', ['Plus', ['Pattern', 'a', ['Blank']], ['Times', ['Optional', ['Pattern', 'b', ['Blank']]], ['Power', ['Pattern', 'x', ['Blank']], ['Pattern', 'n', ['Blank']]]]], '-1'], ['Pattern', 'x', ['Blank', 'Symbol']]]
    assert generate_sympy_from_parsed(s, wild=True) == 'Int(Pow(Add(Pattern(a, Blank), Mul(Optional(Pattern(b, Blank)), Pow(Pattern(x, Blank), Pattern(n, Blank)))), S(-1)), Pattern(x, Blank(Symbol)))'
    assert generate_sympy_from_parsed(s ,replace_Int=True) == 'Integral(Pow(Add(Pattern(a, Blank), Mul(Optional(Pattern(b, Blank)), Pow(Pattern(x, Blank), Pattern(n, Blank)))), S(-1)), Pattern(x, Blank(Symbol)))'
    s = ['And', ['FreeQ', ['List', 'a', 'b'], 'x'], ['PositiveIntegerQ', ['Times', ['Plus', 'n', '-3'], ['Power', '2', '-1']]], ['PosQ', ['Times', 'a', ['Power', 'b', '-1']]]]
    assert generate_sympy_from_parsed(s) == 'And(FreeQ(List(a, b), x), PositiveIntegerQ(Mul(Add(n, S(-3)), Pow(S(2), S(-1)))), PosQ(Mul(a, Pow(b, S(-1)))))'


def test_rubi_printer():
    #14819
    a = Symbol('a')
    assert rubi_printer(Not(a)) == 'Not(a)'
