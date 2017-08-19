'''
Parser for FullForm[Downvalues[]] of Mathematica rules.

This parser is customised to parse the output in MatchPy rules format. Multiple
`Constraints` are divided into individual `Constraints` because it helps the
MatchPy's `ManyToOneReplacer` to backtrack earlier and improve the speed.

Parsed output is formatted into readable format by using `sympify` and print the
expression using `sstr`. This replaces `And`, `Mul`, 'Pow' by their respective
symbols.

References
==========
[1] http://reference.wolfram.com/language/ref/FullForm.html
[2] http://reference.wolfram.com/language/ref/DownValues.html
[3] https://gist.github.com/Upabjojr/bc07c49262944f9c1eb0
'''

import re
from sympy import sympify
from sympy.printing import sstr

replacements = dict(
        Times="Mul",
        Plus="Add",
        Power="Pow",
        Sum="_Sum" # Temporarily replace `Sum` because it can raise errors while sympifying
)

parsed = '''
from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    Pattern, ReplacementRule, ManyToOneReplacer = matchpy.Pattern, matchpy.ReplacementRule, matchpy.ManyToOneReplacer
    from matchpy import Operation, CommutativeOperation, AssociativeOperation, OneIdentityOperation, CustomConstraint
    from matchpy.expressions.functions import register_operation_iterator, register_operation_factory

    from sympy import Pow, Add, Integral, Basic, Mul, S, Or, And, symbols
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf)
    from sympy.functions.elementary.hyperbolic import acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch
    from sympy.functions.elementary.trigonometric import atan, acsc, asin, acot, acos, asec
    from sympy.integrals.rubi.symbol import WC
    from sympy.integrals.rubi.utility_function import (Int, Set, With, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcCsch, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, FixSimplify, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, SimpFixFactor, Simp, SimpHelp, SplitProduct, SplitSum)

    Operation.register(Integral)
    register_operation_iterator(Integral, lambda a: (a._args[0],) + a._args[1], lambda a: len((a._args[0],) + a._args[1]))

    Operation.register(Pow)
    OneIdentityOperation.register(Pow)
    register_operation_iterator(Pow, lambda a: a._args, lambda a: len(a._args))

    Operation.register(Add)
    OneIdentityOperation.register(Add)
    CommutativeOperation.register(Add)
    AssociativeOperation.register(Add)
    register_operation_iterator(Add, lambda a: a._args, lambda a: len(a._args))

    Operation.register(Mul)
    OneIdentityOperation.register(Mul)
    CommutativeOperation.register(Mul)
    AssociativeOperation.register(Mul)
    register_operation_iterator(Mul, lambda a: a._args, lambda a: len(a._args))

    Operation.register(sin)
    register_operation_iterator(sin, lambda a: a._args, lambda a: len(a._args))

    Operation.register(cos)
    register_operation_iterator(cos, lambda a: a._args, lambda a: len(a._args))

    Operation.register(tan)
    register_operation_iterator(tan, lambda a: a._args, lambda a: len(a._args))

    Operation.register(cot)
    register_operation_iterator(cot, lambda a: a._args, lambda a: len(a._args))

    Operation.register(csc)
    register_operation_iterator(csc, lambda a: a._args, lambda a: len(a._args))

    Operation.register(sec)
    register_operation_iterator(sec, lambda a: a._args, lambda a: len(a._args))

    Operation.register(sinh)
    register_operation_iterator(sinh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(cosh)
    register_operation_iterator(cosh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(tanh)
    register_operation_iterator(tanh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(coth)
    register_operation_iterator(coth, lambda a: a._args, lambda a: len(a._args))

    Operation.register(csch)
    register_operation_iterator(csch, lambda a: a._args, lambda a: len(a._args))

    Operation.register(sech)
    register_operation_iterator(sech, lambda a: a._args, lambda a: len(a._args))

    Operation.register(asin)
    register_operation_iterator(asin, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acos)
    register_operation_iterator(acos, lambda a: a._args, lambda a: len(a._args))

    Operation.register(atan)
    register_operation_iterator(atan, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acot)
    register_operation_iterator(acot, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acsc)
    register_operation_iterator(acsc, lambda a: a._args, lambda a: len(a._args))

    Operation.register(asec)
    register_operation_iterator(asec, lambda a: a._args, lambda a: len(a._args))

    Operation.register(asinh)
    register_operation_iterator(asinh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acosh)
    register_operation_iterator(acosh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(atanh)
    register_operation_iterator(atanh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acoth)
    register_operation_iterator(acoth, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acsch)
    register_operation_iterator(acsch, lambda a: a._args, lambda a: len(a._args))

    Operation.register(asech)
    register_operation_iterator(asech, lambda a: a._args, lambda a: len(a._args))

    def sympy_op_factory(old_operation, new_operands, variable_name):
         return type(old_operation)(*new_operands)

    register_operation_factory(Basic, sympy_op_factory)

    A_, B_, C_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, z_ = [WC(i) for i in 'ABCabcdefghijklmnpqrtuvswxz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, jn_, mn_, non2_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'n1', 'n2', 'n3',' Pq', 'Pm', ' Px', 'Qm', 'Qr', 'jn', 'mn', 'non2']]
    p, q = symbols('p q')

def rubi_object():
    rubi = ManyToOneReplacer()

'''

def parse_full_form(wmexpr):
    '''
    Parses FullForm[Downvalues[]] generated by Mathematica
    '''
    out = []
    stack = [out]
    generator = re.finditer(r'[\[\],]', wmexpr)
    last_pos = 0
    for match in generator:
        if match is None:
            break
        position = match.start()
        last_expr = wmexpr[last_pos:position].replace(',', '').replace(']', '').replace('[', '').strip()

        if match.group() == ',':
            if last_expr != '':
                stack[-1].append(last_expr)
        elif match.group() == ']':
            if last_expr != '':
                stack[-1].append(last_expr)
            stack.pop()
            current_pos = stack[-1]
        elif match.group() == '[':
            stack[-1].append([last_expr])
            stack.append(stack[-1][-1])
        last_pos = match.end()
    return out[0]

def get_default_values(parsed, default_values={}):
    '''
    Returns Optional variables and their values in the pattern.
    '''

    if not isinstance(parsed, list):
        return default_values

    if parsed[0] == "Times": # find Default arguments for "Times"
        for i in parsed[1:]:
            if i[0] == "Optional":
                default_values[(i[1][1])] = 1

    if parsed[0] == "Plus": # find Default arguments for "Plus"
        for i in parsed[1:]:
            if i[0] == "Optional":
                default_values[(i[1][1])] = 0

    if parsed[0] == "Power": # find Default arguments for "Power"
        for i in parsed[1:]:
            if i[0] == "Optional":
                default_values[(i[1][1])] = 1

    if len(parsed) == 1:
        return default_values

    for i in parsed:
        default_values = get_default_values(i, default_values)

    return default_values

def add_wildcards(string, optional={}):
    '''
    Replaces `Pattern(variable)` by `variable` in `string`.
    Returns the free symbols present in the string.
    '''

    symbols = [] # stores symbols present in the expression

    p = r'(Optional\(Pattern\((\w+), Blank\)\))'
    matches = re.findall(p, string)
    for i in matches:
        string = string.replace(i[0], "WC('{}', S({}))".format(i[1], optional[i[1]]))
        symbols.append(i[1])

    p = r'(Pattern\((\w+), Blank\))'
    matches = re.findall(p, string)
    for i in matches:
        string = string.replace(i[0], i[1] + '_')
        symbols.append(i[1])

    p = r'(Pattern\((\w+), Blank\(Symbol\)\))'
    matches = re.findall(p, string)
    for i in matches:
        string = string.replace(i[0], i[1] + '_')
        symbols.append(i[1])

    return string, symbols

def seperate_freeq(s, variables=[], x=None):
    '''
    Helper function for parse_freeq().
    '''
    if s[0] == 'FreeQ':
        if len(s[1]) == 1:
            variables = [s[1]]
        else:
            variables = s[1][1:]
        x = s[2]
    else:
        for i in s[1:]:
            variables, x = seperate_freeq(i, variables, x)
    return variables, x

def parse_freeq(l, x):
    '''
    Divies FreeQ() constraint into smaller constraints.
    '''
    res = []
    for i in l:
        res.append(('CustomConstraint(lambda {}, {}: FreeQ({}, {}))').format(i, x, i, x))
    if res != []:
        return ', ' + ', '.join(res)
    return ''

def generate_sympy_from_parsed(parsed, wild=False, symbols=[], replace_Int=False):
    '''
    Parses list into Python syntax.

    Parameters
    ==========
    wild : When set to True, the symbols are replaced as wild symbols.
    symbols : Symbols already present in the pattern.
    replace_Int : When set to True, `Int` is replaced by `Integral`(used to parse pattern in the rule).
    '''
    out = ""

    if not isinstance(parsed, list):
        try: #return S(number) if parsed is Number
            float(parsed)
            return "S({})".format(parsed)
        except:
            pass
        if parsed in symbols:
            if wild:
                return parsed + '_'
        return parsed

    if parsed[0] == 'Rational':
        return 'S({})/S({})'.format(generate_sympy_from_parsed(parsed[1], wild=wild, symbols=symbols, replace_Int=replace_Int), generate_sympy_from_parsed(parsed[2], wild=wild, symbols=symbols, replace_Int=replace_Int))
    elif parsed[0] != 'FreeQ': # FreeQ is parsed seperately
        if parsed[0] in replacements:
            out += replacements[parsed[0]]
        else:
            if parsed[0] == 'Int' and replace_Int:
                out += 'Integral'
            else:
                out += parsed[0]

        if len(parsed) == 1:
            return out

        result = [generate_sympy_from_parsed(i, wild=wild, symbols=symbols, replace_Int=replace_Int) for i in parsed[1:]]
        if '' in result:
            result.remove('')

        out += "("
        out += ", ".join(result)
        out += ")"

    return out

def get_free_symbols(s, symbols, free_symbols=[]):
    '''
    Returns free symbols present in `s`.
    '''
    if not isinstance(s, list):
        if s in symbols:
            free_symbols.append(s)
        return free_symbols

    for i in s:
        free_symbols = get_free_symbols(i, symbols, free_symbols)

    return free_symbols

def _divide_constriant(s, symbols):
    # Creates a CustomConstraint of the form `CustomConstraint(lambda a, x: FreeQ(a, x))`
    if s[0] == 'FreeQ':
        return ''
    lambda_symbols = list(set(get_free_symbols(s, symbols, [])))
    return 'CustomConstraint(lambda {}: {})'.format(', '.join(lambda_symbols), sstr(sympify(generate_sympy_from_parsed(s)), sympy_integers=True))

def divide_constraint(s, symbols):
    '''
    Divides multiple constraints into smaller constraints.
    '''
    if s[0] == 'And':
        result = [_divide_constriant(i, symbols) for i in s[1:]]
    else:
        result = [_divide_constriant(s, symbols)]

    r = ['']
    for i in result:
        if i != '':
            r.append(i)

    return ', '.join(r)

def setWC(string):
    '''
    Replaces `WC(a, b)` by `WC('a', S(b))`
    '''
    p = r'(WC\((\w+), ([-+]?\d)\))'
    matches = re.findall(p, string)
    for i in matches:
        string = string.replace(i[0], "WC('{}', S({}))".format(i[1], i[2]))

    return string

def downvalues_rules(r, parsed):
    '''
    Function which generates parsed rules by substituting all possible
    combinations of default values.
    '''
    res = []

    index = 0
    for i in r:
        #print('parsing rule {}'.format(r.index(i) + 1)) # uncomment to display number of rules parsed

        # Parse Pattern
        if i[1][1][0] == 'Condition':
            p = i[1][1][1].copy()
        else:
            p = i[1][1].copy()

        optional = get_default_values(p, {})
        pattern = generate_sympy_from_parsed(p.copy(), replace_Int=True)
        pattern, free_symbols = add_wildcards(pattern, optional=optional)
        free_symbols = list(set(free_symbols)) #remove common symbols

        # Parse Transformed Expression and Constraints
        if i[2][0] == 'Condition': # parse rules without constraints seperately
            constriant = divide_constraint(i[2][2], free_symbols) # seperate And constraints into indivdual constraints
            FreeQ_vars, FreeQ_x = seperate_freeq(i[2][2].copy()) # seperate FreeQ into indivdual constraints
            transformed = generate_sympy_from_parsed(i[2][1].copy(), symbols=free_symbols)
        else:
            constriant = ''
            FreeQ_vars, FreeQ_x = [], []
            transformed = generate_sympy_from_parsed(i[2].copy(), symbols=free_symbols)

        FreeQ_constraint = parse_freeq(FreeQ_vars, FreeQ_x)
        pattern = sympify(pattern)
        pattern = sstr(pattern, sympy_integers=True)
        pattern = setWC(pattern)
        transformed = sympify(transformed)

        index += 1
        parsed = parsed + '    pattern' + str(index) +' = Pattern(' + pattern + '' + FreeQ_constraint + '' + constriant + ')'
        parsed = parsed + '\n    ' + 'rule' + str(index) +' = ReplacementRule(' + 'pattern' + sstr(index, sympy_integers=True) + ', lambda ' + ', '.join(free_symbols) + ' : ' + str(transformed) + ')\n    '
        parsed = parsed + 'rubi.add(rule'+ str(index) +')\n\n'

    parsed = parsed + '    return rubi\n'

    # Replace modified functions such as `_Sum`
    parsed = parsed.replace('_Sum', 'Sum')

    return parsed

def rubi_rule_parser(fullform, header=None):
    '''
    Parses rules in MatchPy format.

    Parameters
    ==========
    fullform : FullForm of the rule as string.
    header : Header imports for the file. Uses default imports if None.

    Examples
    ========
    >>> from sympy.integrals.rubi.parsetools.parse import rubi_rule_parser
    >>> header = \'''
    ... from matchpy import Operation, CommutativeOperation
    ... def rubi_object():
    ...     rubi = ManyToOneReplacer()
    ...
    ... \'''
    >>> fullform = 'List[RuleDelayed[HoldPattern[Int[Power[Pattern[x,Blank[]],Optional[Pattern[m,Blank[]]]],Pattern[x,Blank[Symbol]]]],Condition[Times[Power[x,Plus[m,1]],Power[Plus[m,1],-1]],And[FreeQ[m,x],NonzeroQ[Plus[m,1]]]]]]'
    >>> rules = rubi_rule_parser(fullform, header)
    >>> print(rules)
    from matchpy import Operation, CommutativeOperation
    def rubi_object():
        rubi = ManyToOneReplacer()
    <BLANKLINE>
        pattern1 = Pattern(Integral(x_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda m: NonzeroQ(m + 1)))
        rule1 = ReplacementRule(pattern1, lambda m, x : x**(m + 1)/(m + 1))
        rubi.add(rule1)
    <BLANKLINE>
        return rubi
    '''
    if not header: # use default header values
        header = parsed

    res = parse_full_form(fullform)
    rules = []

    for i in res: # separate all rules
        if i[0] == 'RuleDelayed':
            rules.append(i)

    return downvalues_rules(rules, header).strip()
