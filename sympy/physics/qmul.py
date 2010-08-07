from sympy.physics.qoperations import QAssocOp
from sympy.physics.quantum import *
from sympy.core.expr import Expr
from sympy.core.cache import cacheit
from sympy.core.mul import Mul
from sympy.physics.quantumbasic import QuantumError, QuantumBasic

class QMul(QAssocOp):
    binop = '*'

    @classmethod
    def _rules_QMul(cls, Object1, Object2):
        #built in flatten
        from sympy.physics.quantum import StateBase, Operator, OuterProduct, KetBase, BraBase, OuterProduct, InnerProduct
        if (not isinstance(Object1, QuantumBasic)) and (not isinstance(Object2, QuantumBasic)):
            return Mul(Object1, Object2)
        elif not isinstance(Object1, QuantumBasic):
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object2.hilbert_space
            retVal.evaluates = Object2.evaluates 
            return retVal
        elif not isinstance(Object2, QuantumBasic):
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
                return retVal
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

    def _eval_dagger(self):
        newargs = []
        for item in reversed(self.args):
            newargs.append(Dagger(item))
        return QMul(*newargs)

    @staticmethod
    def _expandsums(sums):
        """
        Helper function for _eval_expand_mul.

        sums must be a list of instances of Basic.
        """
        from sympy.utilities.iterables import make_list
        from sympy.physics.qadd import QAdd
        L = len(sums)
        if L == 1:
            return sums[0].args
        terms = []
        left = QMul._expandsums(sums[:L//2])
        right = QMul._expandsums(sums[L//2:])

        terms = [QMul(a, b) for a in left for b in right]
        added = QAdd(*terms)
        return make_list(added, QAdd) #it may have collapsed down to one term

    def _eval_expand_mul(self, deep=True, **hints):
        from sympy.physics.qadd import QAdd
        plain, sums, rewrite = [], [], False
        for factor in self.args:
            if deep:
                term = factor.expand(deep=deep, **hints)
                if term != factor:
                    factor = term
                    rewrite = True

            if factor.is_Add or isinstance(factor, QAdd):
                sums.append(factor)
                rewrite = True
            else:
                if factor.is_commutative:
                    plain.append(factor)
                else:
                    Wrapper = QuantumBasic
                    sums.append(Wrapper(factor))

        if not rewrite:
            return self
        else:
            if sums:
                terms = QMul._expandsums(sums)
                plain = QMul(*plain)
                return QAdd(*[QMul(plain, term) for term in terms])
            else:
                return QMul(*plain)     
