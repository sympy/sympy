from sympy.core.expr import Expr
from sympy.core.decorators import call_highest_priority
from sympy.core.sympify import _sympify
from sympy.physics.quantum import Operator, KetBase, BraBase, OuterProduct, InnerProduct
from sympy.core.numbers import Number
from sympy import Symbol

def _Qrules_(class1, class2):
    if issubclass(class1, (Number, Symbol)):
        return class2
    elif issubclass(class2, (Number, Symbol)):
        return class1
    elif issubclass(class1, (Operator, OuterProduct)):
        if issubclass(class2, (Operator, OuterProduct)):
            return class2
        elif issubclass(class2, KetBase):
            return class2
        elif issubclass(class2, BraBase):
            raise Exception("Operator*<Bra| is not allowed")
    elif issubclass(class2, (Operator, OuterProduct)):
        if issubclass(class1, (Operator, OuterProduct)):
            return class2
        elif issubclass(class1, BraBase):
            return class1
        elif issubclass(class1, KetBase):
            raise Exception("|Ket>*Operator is not allowed")
    elif issubclass(class1, KetBase):
        if issubclass(class2, BraBase):
            return Number
        else:
            raise Exception("%s*%s is not allowed" % (class1.__name__, class2.__name__))
    elif issubclass(class1, InnerProduct):
        return class2
    elif issubclass(class2, InnerProduct):
        return class1
    else:
        raise Exception("%s*%s is not allowed" % (class1.__name__, class2.__name__))

class QAssocOp(Expr):
    _op_priority = 100.0
    name = 'QAssocOp'

    #Mul and add need Expand Methods as well as Identity methods. I need to figure out how to set what something evaluates to
    def __new__(cls, *args, **assumptions):
	if len(args) == 0:
            return cls.identity #Every bin-op must define identity (something...something...group theory)
        if len(args) == 1:
            return _sympify(args[0])
        parts = cls.flatten(map(_sympify, args))
        return Expr.__new__(cls, *parts, **assumptions)

    @classmethod
    def flatten(cls, seq):
        print 'flatten', seq
        #determine if this will work
        validationMeth = '_validate_%s' % cls.name
        print validationMeth
        for i in range(len(seq)-1):
            if hasattr(seq[i], validationMeth):
                getattr(seq[i], validationMeth)(seq[i+1])
            elif hasattr(seq[i+1], validationMeth):
                getattr(seq[i+1], validationMeth)(seq[i]) #do rvalidate FIXME
            else:
                _Qrules_(seq[i].__class__, seq[i+1].__class__)
        # apply associativity, no commutativity property is used
        new_seq = []
        while seq:
            o = seq.pop(0)
            if o.__class__ is cls: # classes must match exactly
                seq = list(o[:]) + seq
                continue
            new_seq.append(o)
        # c_part, nc_part, order_symbols
        return new_seq
        
    def __getitem__(self, number):
        return self.args[number]   

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        from sympy.physics.qmul import QMul
        return QMul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        from sympy.physics.qmul import QMul
        return QMul(other, self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return QAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return QAdd(other, self)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return QPow(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        return QPow(other, self)

