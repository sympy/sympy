from matchpy import Operation, Arity

class Int(Operation):
    name = "Int"
    arity = Arity.binary
    associative = False
    def __str__(self):
        return 'Int({}, {})'.format(self.operands[0], self.operands[1])

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

class Log(Operation):
    name = "log"
    arity = Arity.unary
    def __str__(self):
        return 'log({})'.format(self.operands[0])

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

class NonzeroQ(Operation):
    name = "NonzeroQ"
    arity = Arity.unary

class List(Operation):
    name = "list"
    arity = Arity.variadic
    def __str__(self):
        try:
            res = str([s.name for s in self.operands]).replace("'", "")
        except:
            res = str([s.variable_name for s in self.operands]).replace("'", "")
        return 'list('+ res +')'

class RemoveContent(Operation):
    name = "RemoveContent"
    arity = Arity.binary


class Log_parse(Operation):
    name = 'Log'
    arity = Arity.unary

class List_parse(Operation):
    name = "List"
    arity = Arity.variadic

class ConstantSymbol_parse(Operation):
    name = "ConstantSymbol"
    arity = Arity.unary
