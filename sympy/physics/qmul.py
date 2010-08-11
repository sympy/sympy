from sympy.physics.qoperations import QAssocOp
from sympy.physics.quantum import StateBase, Operator, OuterProduct, KetBase, BraBase, OuterProduct, InnerProduct
from sympy.core.expr import Expr
from sympy.core.cache import cacheit
from sympy.core.mul import Mul
from sympy.physics.quantumbasic import QuantumError, QuantumBasic
from sympy import (
    Expr, Basic, sympify, Add, Mul, Pow, 
    I, Function, Integer, S, sympify, Matrix, oo
)
from sympy.printing.pretty.stringpict import prettyForm

class QMul(QAssocOp):
    binopPretty = prettyForm(u'\u00B7')
    binop = '*'

    @classmethod
    def _rules_QMul(cls, Object1, Object2):
        """
            This method is called by new to instantiate a QMul, InnerProduct, or OuterProduct class
            Applies rules of what can and can't be added together by checking types and hilbert spaces
            Returns a Mul, QMul, InnerProduct or OuterProduct object on success
            Raises and exception if input violates quantum shape rules
        """
        if Object1 == 1:
            return Object2
        if Object2 == 1:
            return Object1
        elif Object2 == 0 or Object1 == 0:
            return 0
        
        if (not isinstance(Object1, QuantumBasic)) and (not isinstance(Object2, QuantumBasic)):
            return Mul(Object1, Object2)
        elif not isinstance(Object1, QuantumBasic):
            if Object1 == S.One:
                return Object2
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object2.hilbert_space
            retVal.evaluates = Object2.evaluates 
            return retVal
        elif not isinstance(Object2, QuantumBasic):
            if Object2 == S.One:
                return Object1
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object1.hilbert_space
            retVal.evaluates = Object1.evaluates 
            return retVal

        if Object1.hilbert_space != Object2.hilbert_space:
            raise Exception("Hilbert Spaces do not match")

        if issubclass(Object1.evaluates, InnerProduct):
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object2.hilbert_space
            retVal.evaluates = Object2.evaluates
            return retVal
        elif issubclass(Object2.evaluates, InnerProduct):
            retVal = cls.QMulflatten(Object1, Object2)
            retVal.hilbert_space = Object1.hilbert_space
            retVal.evaluates = Object1.evaluates
            return retVal         
        elif issubclass(Object1.evaluates, Operator):
            if issubclass(Object2.evaluates, (Operator, InnerProduct)):
                retVal = cls.QMulflatten(Object1, Object2)
                retVal.hilbert_space = Object1.hilbert_space
                retVal.evaluates = Object1.evaluates
                return retVal
            elif issubclass(Object2.evaluates, KetBase):
                retVal = cls.QMulflatten(Object1, Object2)            
                retVal.hilbert_space = Object2.hilbert_space
                retVal.evaluates = Object2.evaluates
                return retVal
        elif issubclass(Object2.evaluates, Operator):
            if issubclass(Object1.evaluates, (Operator, InnerProduct)):
                retVal = cls.QMulflatten(Object1, Object2)
                retVal.hilbert_space = Object2.hilbert_space
                retVal.evaluates = Object2.evaluates
                return retVal
            elif issubclass(Object1.evaluates, BraBase):
                retVal = cls.QMulflatten(Object1, Object2)
                retVal.hilbert_space = Object1.hilbert_space
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
        """
            Flattens out QMul objects.
            Places Non-Quantum objects at front of Qmul in Mul object and Quantum objects behind it
        """
        Qseq = []
        Eseq = []
        seq = [Object1, Object2]
        while seq:
            o = seq.pop(0)                
            if o.__class__ is cls: # classes must match exactly
                seq = list(o[:]) + seq
                continue
            if not isinstance(o, (QuantumBasic)):
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
            return args[0], self._new_rawargs(self.evaluates, self.hilbert_space, *args[1:]) 

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
        """
            A facsimile of _eval_expand_mul in regular Mul
        """
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
                    sums.append(factor)

        if not rewrite:
            return self
        else:
            if sums:
                terms = QMul._expandsums(sums)
                plain = QMul(*plain)
                return QAdd(*[QMul(plain, term) for term in terms])
            else:
                return QMul(*plain)     
