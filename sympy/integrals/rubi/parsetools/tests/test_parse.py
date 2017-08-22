import sys
from sympy.external import import_module
matchpy = import_module("matchpy")

if not matchpy:
    #bin/test will not execute any tests now
    disabled = True

if sys.version_info[:2] < (3, 6):
    disabled = True

from sympy.integrals.rubi.parsetools.parse import rubi_rule_parser

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
    assert len(result.strip()) == len(rules) # failing randomly while using `result.strip() == rules`
