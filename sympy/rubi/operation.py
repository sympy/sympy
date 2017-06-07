from matchpy import Operation, Arity

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
    arity = Arity.binary
    associative = False
    commutative = False

class And(Operation):
    name = "And"
    arity = Arity.variadic
    one_identity = True

class Or(Operation):
    name = "Or"
    arity = Arity.variadic
    one_identity = True

class FreeQ(Operation):
    name = "FreeQ"
    arity = Arity.binary

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
    name = "log"
    arity = Arity.unary

class RemoveContent(Operation):
    name = "RemoveContent"
    arity = Arity.binary

class PositiveIntegerQ(Operation):
    name = "PositiveIntegerQ"
    arity = Arity.unary

class NegativeIntegerQ(Operation):
    name = "NegativeIntegerQ"
    arity = Arity.unary

class PositiveQ(Operation):
    name = "PositiveQ"
    arity = Arity.unary

class IntegerQ(Operation):
    name = "IntegerQ"
    arity = Arity.unary

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
    arity = Arity.unary

class Subst(Operation):
    name = "Subst"
    arity = Arity.variadic

class LinearQ(Operation):
    name = "LinearQ"
    arity = Arity.binary

class Sqrt(Operation):
    name = "Sqrt"
    arity = Arity.unary

class ArcCosh(Operation):
    name = "ArcCosh"
    arity = Arity.unary









class ConstantSymbol_parse(Operation):
    name = "ConstantSymbol"
    arity = Arity.unary
