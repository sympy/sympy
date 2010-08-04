from sympy.core.expr import Expr
from sympy.core.decorators import call_highest_priority
from sympy.core.sympify import _sympify
from sympy.physics.quantum import Operator, KetBase, BraBase, OuterProduct, InnerProduct
from sympy.core.numbers import Number
from sympy import Symbol

class QRules(object):
#also need to keep track of hilbert space
    def __init__(self, qclass, hilbert_space):
        self.qclass = qclass
        self.hilbert_space = hilbert_space

    def __getitem__(self, number):
        if number == 1:
            return self.hilbert_space
        elif number == 0:
            return self.qclass
        else:
            raise Exception("Index out of bounds")

    @staticmethod
    def _rules_QMul(Object1, Object2):
        if isinstance(Object1, QRules):
            class1 = Object1 
        elif hasattr(Object1, 'eval_to'): 
            class1 = Object1.eval_to
        else:
            class1 = QRules(Object1.__class__, Object1.hilbert_space)
    
        if isinstance(Object2, QRules):
            class2 = Object2 
        elif hasattr(Object2, 'eval_to'): 
            class2 = Object2.eval_to
        else:
            class2 = QRules(Object2.__class__, Object2.hilbert_space)

        if class1.hilbert_space != class2.hilbert_space:
            raise Exception("Hilbert Spaces do not match")

        if (not issubclass(class1[0], (Operator, OuterProduct, KetBase, BraBase))) and (not issubclass(class2[0], (Operator, OuterProduct,    KetBase, BraBase))):
            return QRules(Number, None)    
        elif issubclass(class1[0], (Number, Symbol)):
            return class2
        elif issubclass(class2[0], (Number, Symbol)):
            return class1
        elif issubclass(class1[0], (Operator, OuterProduct)):
            if issubclass(class2[0], (Operator, OuterProduct)):
                return class2
            elif issubclass(class2[0], KetBase):
                return class2
        elif issubclass(class2[0], (Operator, OuterProduct)):
            if issubclass(class1[0], (Operator, OuterProduct)):
                return class2
            elif issubclass(class1[0], BraBase):
                return class1
        elif issubclass(class1[0], KetBase) and issubclass(class2[0], BraBase):
            return QRules(OuterProduct, class1.hilbert_space)
        elif issubclass(class1[0], BraBase) and issubclass(class2[0], KetBase):
            return QRules(InnerProduct, class1.hilbert_space)
        elif issubclass(class1[0], InnerProduct):
            return class2
        elif issubclass(class2[0], InnerProduct):
            return class1
        raise Exception("%s*%s is not allowed" % (class1[0].__name__, class2[0].__name__))

    @staticmethod    
    def _rules_QAdd(Object1, Object2):
        if isinstance(Object1, QRules):
            class1 = Object1 
        elif hasattr(Object1, 'eval_to'): 
            class1 = Object1.eval_to
        else:
            class1 = QRules(Object1.__class__, Object1.hilbert_space) 
    
        if isinstance(Object2, type):
            class2 = Object2
        elif hasattr(Object2, 'eval_to'):
            class2 = Object2.eval_to
        else:
            class2 = QRules(Object2.__class__, Object2.hilbert_space)
            
        if (not issubclass(class1, (Operator, OuterProduct, KetBase, BraBase))) and (not issubclass(class2, (Operator, OuterProduct,    KetBase, BraBase))):
            return QRules(Number, None)
        elif class1 == class2:
            return class1
        else:
            raise Exception("Can't add (%s + %s)" % (class1[0].__name__, class2[0].__name__))

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
        #determine if this will work
        validationMeth = '_validate_%s' % cls.name
        rules = getattr(QRules, '_rules_%s' % cls.name)
        result = 0
        for i in range(len(seq)):
            if result:
                result = rules(result, seq[i])
            else:
                if hasattr(seq[i], 'eval_to'):
                    result = seq[i].eval_to
                else:
                    result = QRules(seq[i].__class__, seq[i].hilbert_space)
                
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

    @property
    def hilbert_space(self):
        return self.args[0].hilbert_space 
        
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
        from sympy.physics.qadd import QAdd    
        return QAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        from sympy.physics.qadd import QAdd       
        return QAdd(other, self)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return QPow(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        return QPow(other, self)

