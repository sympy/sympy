from sympy.core.expr import Expr
from sympy.core.decorators import call_highest_priority
from sympy.core.sympify import sympify
from sympy.physics.quantum import Operator, KetBase, BraBase, OuterProduct, InnerProduct, StateBase
from sympy.core.numbers import Number
from sympy import Symbol
from sympy.core.mul import Mul
from sympy.printing.str import sstr

class QAssocOp(Expr):
    _op_priority = 100.0
    __slots__ = ['evaluates', 'hilbert_space']

    #Mul and add need Expand Methods as well as Identity methods. I need to figure out how to set what something evaluates
    def __new__(cls, *args, **assumptions):
        if len(args) == 1:
            return sympify(args[0])
        return cls.instantiate(map(sympify, args))

    @classmethod
    def instantiate(cls, seq):
        #determine if this will work flatten needs to flattening (pull out non quantum parts as they belong to an abelian group
        rules = getattr(cls, '_rules_%s' % cls.__name__)
        for i in range(len(seq)-1):
            result = rules(seq[i], seq[i+1])
        return result

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

    def _sympystr(self, printer, *args):
        string = ''
        length = len(self.args)
        for i in range(length):
            string = string + sstr(self.args[i])
            if i != length-1:
                string = string + self.__class__.binop
        return string 
    

    def _new_rawargs(self, *args):
        """create new instance of own class with args exactly as provided by caller

           This is handy when we want to optimize things, e.g.

           >>> from sympy import Mul, symbols
           >>> from sympy.abc import x, y
           >>> e = Mul(3,x,y)
           >>> e.args
           (3, x, y)
           >>> Mul(*e.args[1:])
           x*y
           >>> e._new_rawargs(*e.args[1:])  # the same as above, but faster
           x*y

        """
        obj = Expr.__new__(type(self), *args)  # NB no assumptions for Add/Mul

        return obj

