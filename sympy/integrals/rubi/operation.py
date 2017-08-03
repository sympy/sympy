from sympy.external import import_module
matchpy = import_module("matchpy")

if matchpy:
    Operation, Arity = matchpy.Operation, matchpy.Arity
else:
    Operation = object
    class Arity(object):
        nullary = (0, True)
        unary = (1, True)
        binary = (2, True)
        ternary = (3, True)
        polyadic = (2, False)
        variadic = (0, False)

class Int(Operation):
    name = "Int"
    arity = Arity.binary

class Mul(Operation):
    name = "Mul"
    arity = Arity.variadic
    associative = True
    commutative = True
    one_identity = True

class Add(Operation):
    name = "Add"
    arity = Arity.variadic
    associative = True
    commutative = True
    one_identity = True

class Pow(Operation):
    name = "Pow"
    arity = Arity.variadic
    one_identity = True

class And(Operation):
    name = "And"
    arity = Arity.variadic
    one_identity = True

class Or(Operation):
    name = "Or"
    arity = Arity.variadic
    one_identity = True

class ZeroQ(Operation):
    name = "ZeroQ"
    arity = Arity.unary

class NonzeroQ(Operation):
    name = "NonzeroQ"
    arity = Arity.unary

class List(Operation):
    name = "List"
    arity = Arity.variadic

class Log(Operation):
    name = "Log"
    arity = Arity.unary

class RemoveContent(Operation):
    name = "RemoveContent"
    arity = Arity.binary

class PositiveIntegerQ(Operation):
    name = "PositiveIntegerQ"
    arity = Arity.variadic

class NegativeIntegerQ(Operation):
    name = "NegativeIntegerQ"
    arity = Arity.variadic

class PositiveQ(Operation):
    name = "PositiveQ"
    arity = Arity.unary

class IntegerQ(Operation):
    name = "IntegerQ"
    arity = Arity.unary

class IntegersQ(Operation):
    name = "IntegersQ"
    arity = Arity.variadic

class PosQ(Operation):
    name = "PosQ"
    arity = Arity.unary

class NegQ(Operation):
    name = "NegQ"
    arity = Arity.unary

class FracPart(Operation):
    name = "FracPart"
    arity = Arity.unary

class IntPart(Operation):
    name = "IntPart"
    arity = Arity.unary

class RationalQ(Operation):
    name = "RationalQ"
    arity = Arity.variadic

class Subst(Operation):
    name = "Subst"
    arity = Arity.variadic

class LinearQ(Operation):
    name = "LinearQ"
    arity = Arity.binary

class Sqrt(Operation):
    name = "Sqrt"
    arity = Arity.unary

class NegativeQ(Operation):
    name = "NegativeQ"
    arity = Arity.unary

class ArcCosh(Operation):
    name = "ArcCosh"
    arity = Arity.unary

class Rational(Operation):
    name = "Rational"
    arity = Arity.binary

class Less(Operation):
    name = "Less"
    arity = Arity.variadic

class Not(Operation):
    name = "Not"
    arity = Arity.unary

class Simplify(Operation):
    name = "Simplify"
    arity = Arity.unary

class Denominator(Operation):
    name = "Denominator"
    arity = Arity.unary

class Coefficient(Operation):
    name = "Coefficient"
    arity = Arity.ternary

class SumSimplerQ(Operation):
    name = "SumSimplerQ"
    arity = Arity.binary

class Equal(Operation):
    name = "Equal"
    arity = Arity.binary

class Unequal(Operation):
    name = "Unequal"
    arity = Arity.binary

class SimplerQ(Operation):
    name = "SimplerQ"
    arity = Arity.binary

class LessEqual(Operation):
    name = "LessEqual"
    arity = Arity.variadic

class IntLinearcQ(Operation):
    name = "IntLinearcQ"
    arity = Arity.variadic

class Greater(Operation):
    name = "Greater"
    arity = Arity.variadic

class GreaterEqual(Operation):
    name = "GreaterEqual"
    arity = Arity.variadic

class FractionQ(Operation):
    name = "FractionQ"
    arity = Arity.variadic

class ExpandIntegrand(Operation):
    name = "ExpandIntegrand"
    arity = Arity.binary

class With(Operation):
    name = "With"
    arity = Arity.binary

class Set(Operation):
    name = "Set"
    arity = Arity.binary

class Hypergeometric2F1(Operation):
    name = "Hypergeometric2F1"
    arity = Arity.variadic

class TogetherSimplify(Operation):
    name = "TogetherSimplify"
    arity = Arity.unary

class Inequality(Operation):
    name = "Inequality"
    arity = Arity.variadic
