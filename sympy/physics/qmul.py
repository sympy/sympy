from sympy.physics.qoperations import QAssocOp
from sympy.physics.quantum import Operator, OuterProduct, KetBase, BraBase,\
InnerProduct, Dagger
from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.cache import cacheit
from sympy.core.mul import Mul
from sympy.physics.qexpr import QuantumError, QExpr
from sympy import  Mul, S
from sympy.printing.pretty.stringpict import prettyForm
from sympy.physics.hilbert import HilbertSpaceException
from sympy.physics.qpow import QPow

class QMul(QAssocOp):
    binopPretty = prettyForm(u'\u00B7')
    binop = '*'

    @classmethod
    def _rules_QMul(cls, object1, object2):
        """
            This method is called by new to instantiate a QMul, InnerProduct, or
            OuterProduct class Applies rules of what can and can't be added
            together by checking types and hilbert spaces Returns a Mul, QMul,
            InnerProduct or OuterProduct object on success Raises and exception
            if input violates quantum shape rules
        """
        if object1 is S.One:
            return object2
        if object2 is S.One:
            return object1
        elif object2 is S.Zero or object1 is S.Zero:
            return 0

        if (not isinstance(object1, QExpr)) and (not isinstance(object2,\
        QExpr)):
            return Mul(object1, object2)
        elif not isinstance(object1, QExpr):
            if object1 == S.One:
                return object2
            retVal = cls.QMulflatten(object1, object2)
            retVal.hilbert_space = object2.hilbert_space
            retVal.acts_like = object2.acts_like
            return retVal
        elif not isinstance(object2, QExpr):
            if object2 == S.One:
                return object1
            retVal = cls.QMulflatten(object1, object2)
            retVal.hilbert_space = object1.hilbert_space
            retVal.acts_like = object1.acts_like
            return retVal

        if object1.hilbert_space != object2.hilbert_space:
            raise HilbertSpaceException("Hilbert Spaces do not match")

        if issubclass(object1.acts_like, InnerProduct):
            retVal = cls.QMulflatten(object1, object2)
            retVal.hilbert_space = object2.hilbert_space
            retVal.acts_like = object2.acts_like
            return retVal
        elif issubclass(object2.acts_like, InnerProduct):
            retVal = cls.QMulflatten(object1, object2)
            retVal.hilbert_space = object1.hilbert_space
            retVal.acts_like = object1.acts_like
            return retVal
        elif issubclass(object1.acts_like, Operator):
            if issubclass(object2.acts_like, (Operator, InnerProduct)):
                retVal = cls.QMulflatten(object1, object2)
                retVal.hilbert_space = object1.hilbert_space
                retVal.acts_like = object1.acts_like
                return retVal
            elif issubclass(object2.acts_like, KetBase):
                retVal = cls.QMulflatten(object1, object2)
                retVal.hilbert_space = object2.hilbert_space
                retVal.acts_like = object2.acts_like
                return retVal
        elif issubclass(object2.acts_like, Operator):
            if issubclass(object1.acts_like, (Operator, InnerProduct)):
                retVal = cls.QMulflatten(object1, object2)
                retVal.hilbert_space = object2.hilbert_space
                retVal.acts_like = object2.acts_like
                return retVal
            elif issubclass(object1.acts_like, BraBase):
                retVal = cls.QMulflatten(object1, object2)
                retVal.hilbert_space = object1.hilbert_space
                retVal.acts_like = object1.acts_like
                return retVal
        #figure out inner and outer products
        elif issubclass(object1.acts_like, KetBase) and\
        issubclass(object2.acts_like, BraBase):
            retVal = cls.QMulflatten(object1, object2)
            retVal.hilbert_space = object2.hilbert_space
            retVal.acts_like = OuterProduct
            return retVal
        elif issubclass(object1.acts_like, BraBase) and\
        issubclass(object2.acts_like, KetBase):
            retVal = cls.QMulflatten(object1, object2)
            retVal.hilbert_space = object2.hilbert_space
            retVal.acts_like = InnerProduct
            return retVal
        raise QuantumError("%s*%s is not allowed" % (object1.acts_like.__name__,\
        object2.acts_like.__name__))

    @classmethod
    def QMulflatten(cls, object1, object2):
        """
            Flattens out QMul objects.
            Places Non-Quantum objects at front of Qmul in Mul object and
            Quantum objects behind it
        """
        Qseq = []
        Eseq = []
        seq = [object1, object2]
        last_argument = None
        while seq:
            o = seq.pop(0)
            if o.__class__ is cls: # classes must match exactly
                seq = list(o[:]) + seq
                continue
            if not isinstance(o, (QExpr)):
                Eseq.append(o)
            else:
                #now, figure out if we should combine terms inside a Pow
                if not last_argument is None:
                    if isinstance(last_argument, QPow):
                        basel = last_argument.base
                        powerl = last_argument.exp
                    else:
                        basel = last_argument
                        powerl = 1

                    if isinstance(o, QPow):
                        baseo = o.base
                        powero = o.exp
                    else:
                        baseo = o
                        powero = 1

                    #for now we won't combine anything questionable into a QPow
                    #(e.g. complicated expressions whose exponent's act_like ==
                    #InnerProduct)
                    if baseo == basel and not isinstance(powero, QExpr)\
                    and not isinstance(powero, QExpr):
                        Qseq.pop(-1)
                        o = baseo**(powero+powerl)

                Qseq.append(o)
                last_argument = o
        nqpart = Mul(*Eseq)
        if nqpart == 1:
            if len(Qseq) ==1:
                return Qseq[0]
            return Expr.__new__(cls, *Qseq)
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
            return args[0], self.__class__._new_rawargs(self.acts_like,\
            self.hilbert_space, *args[1:])

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
                    Wrapper = Basic
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
