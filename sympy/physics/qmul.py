from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.cache import cacheit
from sympy.core.mul import Mul
from sympy.core.basic import S
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.qoperations import QAssocOp
from sympy.physics.qexpr import QuantumError, QExpr
from sympy.physics.hilbert import HilbertSpaceError
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

        e_part, q_part = cls._flatten_to_parts(obj1, obj2)

        # Try to create InnerProduct or OuterProduct, but we first need
        # to make sure both obj1 and obj2 are QExprs.
        if len(q_part)==2 and obj1_is_qexpr and obj2_is_qexpr:
            if issubclass(obj1.acts_like, KetBase) and \
                 issubclass(obj2.acts_like, BraBase):
                try:
                    result = OuterProduct(q_part[0], q_part[1])
                except (TypeError, QuantumError):
                    pass
                else:
                    q_part = [result]
            elif issubclass(obj1.acts_like, BraBase) and \
                 issubclass(obj2.acts_like, KetBase):
                try:
                    result = InnerProduct(q_part[0], q_part[1])
                except (TypeError, QuantumError):
                    pass
                else:
                    e_part.append(result)
                    q_part = []

        result = cls._new_from_parts(e_part, q_part)
        if not isinstance(result, QExpr):
            return result

        if not obj1_is_qexpr:
            result.hilbert_space = obj2.hilbert_space
            result.acts_like = obj2.acts_like
            return result
        elif not obj2_is_qexpr:
            result.hilbert_space = obj1.hilbert_space
            result.acts_like = obj1.acts_like
            return result

        # Both obj1 and obj2 are QExprs from here on out.
        obj1_acts_like = obj1.acts_like
        obj2_acts_like = obj2.acts_like

        if obj1.hilbert_space != obj2.hilbert_space:
            raise HilbertSpaceError(
                "Incompatible hilbert spaces: %r and %r" % (obj1, obj2)
            )

        result.hilbert_space = obj1.hilbert_space
        
        if issubclass(obj1_acts_like, InnerProduct):
            result.acts_like = obj2_acts_like
            return result
        elif issubclass(obj2_acts_like, InnerProduct):
            result.acts_like = obj1_acts_like
            return result

        if issubclass(obj1_acts_like, Operator):
            if issubclass(obj2_acts_like, Operator):
                result.acts_like = obj1_acts_like
                return result
            elif issubclass(obj2_acts_like, KetBase):
                result.acts_like = obj2_acts_like
                return result
        elif issubclass(obj2_acts_like, Operator):
            if issubclass(obj1_acts_like, Operator):
                result.acts_like = obj2_acts_like
                return result
            elif issubclass(obj1_acts_like, BraBase):
                result.acts_like = obj1_acts_like
                return result

        if issubclass(obj1_acts_like, KetBase) and \
           issubclass(obj2_acts_like, BraBase):
            result.acts_like = OuterProduct
            return result
        elif issubclass(obj1_acts_like, BraBase) and \
             issubclass(obj2_acts_like, KetBase):
            result.acts_like = InnerProduct
            return result

        raise TypeError("Can't multiply: %s and %s" % (
            obj1.acts_like.__name__, obj2.acts_like.__name__
        ))

    @classmethod
    def flatten(cls, obj1, obj2):
        e_part, q_part = cls._flatten_to_parts(obj1, obj2)
        return cls._new_from_parts(e_part, q_part)

    @classmethod
    def _flatten_to_parts(cls, obj1, obj2):
        """Flattens out QMul args.

        This method splits the args into two parts: objects that are not QExpr
        instances and object that are non-commutative QExpr instances.
        """
        e_part = []  # Expr parts
        q_part = []  # Non-commuting QExpr parts (Operator, State, etc.)
        seq = [obj1, obj2]
        last_argument = None
        while seq:
            o = seq.pop(0)

            # Here we flatten args that are also QMuls - the classes must
            # match exactly.
            if o.__class__ is cls:
                seq = list(o[:]) + seq
                continue

            if not isinstance(o, QExpr):
                e_part.append(o)
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
                        q_part.pop(-1)
                        o = baseo**(powero+powerl)

                q_part.append(o)
                last_argument = o
        return e_part, q_part

    @classmethod
    def _new_from_parts(cls, e_part, q_part):
        e_part = Mul(*e_part)
        # Handle the case where we only have a commuting Expr part.
        if len(q_part) == 0:
            return e_part
        # Now the case where we only have a QExpr part.
        if e_part == 1:
            if len(q_part) == 1:
                return q_part[0]
            return Expr.__new__(cls, *q_part, **{'commutative': False})
        # Both Expr and QExpr parts
        return Expr.__new__(cls, e_part, *q_part, **{'commutative': False})

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
