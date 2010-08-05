from sympy.physics.qoperations import QAssocOp
from sympy.physics.quantum import *
from sympy.core.expr import Expr
from sympy.core.cache import cacheit
from sympy.core.mul import Mul

class QMul(QAssocOp):
    binop = '*'

    @classmethod
    def _rules_QMul(cls, Object1, Object2):
        #built in flatten
        from sympy.physics.quantum import StateBase, Operator, OuterProduct, KetBase, BraBase, OuterProduct, InnerProduct
        if (not isinstance(Object1, (Operator, OuterProduct, StateBase, QAssocOp))) and (not isinstance(Object2, (Operator, OuterProduct, StateBase, QAssocOp))):
            return Mul(Object1, Object2)
        elif not isinstance(Object1, (Operator, OuterProduct, StateBase, QAssocOp)):
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object2.hilbert_space
            retVal.evaluates = Object2.evaluates 
            return retVal
        elif not isinstance(Object2, (Operator, OuterProduct, StateBase, QAssocOp)):
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object1.hilbert_space
            retVal.evaluates = Object1.evaluates 
            return retVal
            
        if Object1.hilbert_space != Object2.hilbert_space:
            raise Exception("Hilbert Spaces do not match") 
         
        if issubclass(Object1.evaluates, (Operator, OuterProduct)):
            if issubclass(Object2.evaluates, (Operator, OuterProduct)):
                retVal = cls.QMulflatten(Object1, Object2)
                retVal.hilbert_space = Object2.hilbert_space
                retVal.evaluates = Object2.evaluates
                return retVal
            elif issubclass(Object2.evaluates, KetBase):
                retVal = cls.QMulflatten(Object1, Object2)            
                retVal.hilbert_space = Object2.hilbert_space
                retVal.evaluates = Object2.evaluates
                return retVal
        elif issubclass(Object2.evaluates, (Operator, OuterProduct)):
            if issubclass(Object1.evaluates, (Operator, OuterProduct)):
                retVal = cls.QMulflatten(Object1, Object2)
                retVal.hilbert_space = Object2.hilbert_space
                retVal.evaluates = Object2.evaluates
                return retVal
            elif issubclass(Object1.evaluates, BraBase):
                retVal = cls.QMulflatten(Object1, Object2)
                retVal.hilbert_space = Object2.hilbert_space
                retVal.evaluates = Object1.evaluates
                return Object1
        #figure out inner and outer products
        elif issubclass(Object1.evaluates, KetBase) and issubclass(Object2.evaluates, BraBase):
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object2.hilbert_space
            retVal.evaluates = OuterProduct           
            return retVal
        elif issubclass(Object1.evaluates, BraBase) and issubclass(Object2.evaluates, KetBase):
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object2.hilbert_space
            retVal.evaluates = InnerProduct
            return retVal
        raise Exception("%s*%s is not allowed" % (Object1.evaluates.__name__, Object2.evaluates.__name__))  

    @classmethod
    def QMulflatten(cls, Object1, Object2):
        # get Non-QUantum stuff out
        from sympy.physics.quantum import StateBase, Operator, OuterProduct, KetBase, BraBase, OuterProduct, InnerProduct
        from sympy.core.mul import Mul
        Qseq = []
        Eseq = []
        seq = [Object1, Object2]
        while seq:
            o = seq.pop(0)                
            if o.__class__ is cls: # classes must match exactly
                seq = list(o[:]) + seq
                continue
            if not isinstance(o, (Operator, OuterProduct, StateBase, QAssocOp)):
                Eseq.append(o)
            else:
                Qseq.append(o)
        # c_part, nc_part, order_symbols
        return Expr.__new__(cls, Mul(*Eseq), *Qseq)

    @cacheit
    def as_two_terms(self):
        args = self.args

        if len(args) == 1:
            return S.One, self
        elif len(args) == 2:
            return args

        else:
            return args[0], self._new_rawargs(*args[1:])

    @property
    def identity(self):
        return S.One
                
