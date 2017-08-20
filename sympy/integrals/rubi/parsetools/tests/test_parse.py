from sympy.integrals.rubi.parsetools.parse import rubi_rule_parser

def test_rubi_rule_parser():
    header = '''
from matchpy import Operation, CommutativeOperation
    rubi = ManyToOneReplacer()
    '''
    fullform = 'List[RuleDelayed[HoldPattern[Int[Power[Pattern[x,Blank[]],Optional[Pattern[m,Blank[]]]],Pattern[x,Blank[Symbol]]]],Condition[Times[Power[x,Plus[m,1]],Power[Plus[m,1],-1]],And[FreeQ[m,x],NonzeroQ[Plus[m,1]]]]]]'
    rules = rubi_rule_parser(fullform, header)
    result = '''
from matchpy import Operation, CommutativeOperation
    rubi = ManyToOneReplacer()
        pattern1 = Pattern(Integral(x_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + 1)))
    rule1 = ReplacementRule(pattern1, lambda x, m : x**(m + 1)/(m + 1))
    rubi.add(rule1)

    return rubi
'''
    assert result.strip() == rules
