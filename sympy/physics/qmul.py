from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.cache import cacheit
from sympy.core.mul import Mul
from sympy.core.basic import S
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.qoperations import QAssocOp
from sympy.physics.qexpr import QuantumError, QExpr
from sympy.physics.hilbert import HilbertSpaceException
from sympy.physics.qpow import QPow


class QMul(QAssocOp):

    binop = '*'
    binop_pretty = prettyForm(u'\u00B7')

    @classmethod
    def _apply_rules(cls, obj1, obj2):
        """Apply rules to simplify a new QMul instance.

        This method is called by ``eval`` to instantiate a QMul class and
        applies rules of what can and can't be added together by checking
        types and hilbert spaces.
        """

        from sympy.physics.quantum import (
            Operator, OuterProduct, KetBase, BraBase,
            InnerProduct
        )
        
        if obj1 is S.One:
            return obj2
        if obj2 is S.One:
            return obj1
        elif obj2 is S.Zero or obj1 is S.Zero:
            return 0

        obj1_is_qexpr = isinstance(obj1, QExpr)
        obj2_is_qexpr = isinstance(obj2, QExpr)

        if (not obj1_is_qexpr) and (not obj2_is_qexpr):
            return Mul(obj1, obj2)
        elif not obj2_is_qexpr:
            if obj1 == S.One:
                return obj2
            result = cls.flatten(obj1, obj2)
            result.hilbert_space = obj2.hilbert_space
            result.acts_like = obj2.acts_like
            return result
        elif not obj2_is_qexpr:
            if obj2 == S.One:
                return obj1
            result = cls.flatten(obj1, obj2)
            result.hilbert_space = obj1.hilbert_space
            result.acts_like = obj1.acts_like
            return result

        if obj1.hilbert_space != obj2.hilbert_space:
            raise HilbertSpaceException("Hilbert Spaces do not match")

        if issubclass(obj1.acts_like, InnerProduct):
            result = cls.flatten(obj1, obj2)
            result.hilbert_space = obj2.hilbert_space
            result.acts_like = obj2.acts_like
            return result
        elif issubclass(obj2.acts_like, InnerProduct):
            result = cls.flatten(obj1, obj2)
            result.hilbert_space = obj1.hilbert_space
            result.acts_like = obj1.acts_like
            return result
        elif issubclass(obj1.acts_like, Operator):
            if issubclass(obj2.acts_like, (Operator, InnerProduct)):
                result = cls.flatten(obj1, obj2)
                result.hilbert_space = obj1.hilbert_space
                result.acts_like = obj1.acts_like
                return result
            elif issubclass(obj2.acts_like, KetBase):
                result = cls.flatten(obj1, obj2)
                result.hilbert_space = obj2.hilbert_space
                result.acts_like = obj2.acts_like
                return result
        elif issubclass(obj2.acts_like, Operator):
            if issubclass(obj1.acts_like, (Operator, InnerProduct)):
                result = cls.flatten(obj1, obj2)
                result.hilbert_space = obj2.hilbert_space
                result.acts_like = obj2.acts_like
                return result
            elif issubclass(obj1.acts_like, BraBase):
                result = cls.flatten(obj1, obj2)
                result.hilbert_space = obj1.hilbert_space
                result.acts_like = obj1.acts_like
                return result

        # Figure out inner and outer products.
        elif issubclass(obj1.acts_like, KetBase) and\
        issubclass(obj2.acts_like, BraBase):
            result = cls.flatten(obj1, obj2)
            result.hilbert_space = obj2.hilbert_space
            result.acts_like = OuterProduct
            return result
        elif issubclass(obj1.acts_like, BraBase) and\
        issubclass(obj2.acts_like, KetBase):
            result = cls.flatten(obj1, obj2)
            result.hilbert_space = obj2.hilbert_space
            result.acts_like = InnerProduct
            return result

        raise QuantumError("%s*%s is not allowed" % (
            obj1.acts_like.__name__, obj2.acts_like.__name__
        ))

    @classmethod
    def flatten(cls, obj1, obj2):
        """
            Flattens out QMul objects.
            Places Non-Quantum objects at front of Qmul in Mul object and
            Quantum objects behind it
        """
        Qseq = []
        Eseq = []
        seq = [obj1, obj2]
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
        from sympy.physics.quantum import Dagger
        newargs = [Dagger(item) for item in reversed(self.args)]
        return QMul(*newargs)

    @staticmethod
    def _expandsums(sums):
        """Helper function for _eval_expand_mul.

        The argument ``sums`` must be a list of instances of Basic.
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
        """A facsimile of _eval_expand_mul in regular Mul."""
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
