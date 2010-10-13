from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.basic import S
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.qoperations import QAssocOp
from sympy.physics.qexpr import QuantumError, QExpr
from sympy.physics.hilbert import HilbertSpaceError


class QAdd(QAssocOp):
    """
        Quantum Add operation
    """
    binop = ' + '
    binop_pretty =  prettyForm(u' \u002B ')

    @classmethod
    def _is_qadd_allowed(cls, e):
        """Is the expression an instance that can be QAdd'd."""
        from sympy.physics.qpow import QPow
        from sympy.physics.quantum import (
            StateBase, Operator, Dagger, Commutator
        )
        _allowed = (
            StateBase, Operator, QAssocOp, QPow, Dagger, 
            Commutator
        )
        return isinstance(e, _allowed)

    @classmethod
    def _acts_the_same(cls, e1, e2):
        """Do the two expressions act the same."""
        return issubclass(e1.acts_like, e2.acts_like) \
        or issubclass(e2.acts_like, e1.acts_like)

    @classmethod
    def _apply_rules(cls, obj1, obj2):
        """Apply rules to simplify a new QAdd instance.

        This method is called by ``eval`` to instantiate a QAdd class and
        applies rules of what can and can't be added together by checking
        types and hilbert spaces.
        """
        if obj2 is S.Zero:
            return obj1
        elif obj1 is S.Zero:
            return obj2

        if not isinstance(obj1, QExpr) and not isinstance(obj2, QExpr):
            return Add(obj1, obj2)

        if not cls._is_qadd_allowed(obj1) or \
           not cls._is_qadd_allowed(obj2):
            raise TypeError("Can't add: %s and %s" % (
                obj1.__class__.__name__, obj2.__class__.__name__
            ))

        if obj1.hilbert_space != obj2.hilbert_space:
            raise HilbertSpaceError(
                "Incompatible hilbert spaces: %r and %r" % (obj1, obj2)
            )

        if cls._acts_the_same(obj1, obj2):
            result = cls.flatten([obj1, obj2])
            result.hilbert_space = obj1.hilbert_space
            result.acts_like = obj1.acts_like
            return result
        else:
            raise TypeError("Can't add two objects that act like: %s and %s"\
            % (obj1.acts_like.__name__, obj2.acts_like.__name__))

    @classmethod
    def flatten(cls, seq):
        """Flatten additive expressions.

        This simplifies sums by combining like terms. This assumes
        associativity for multiplication and assoc. and comm. for addition.
        Each element of the seq must be a QExpr instance.
        """
        from sympy.physics.qmul import QMul
        terms = {}      # term -> coeff
                        # e.g. x**2 -> 5   for ... + 5*x**2 + ...

        for o in seq:
            if isinstance(o, QAdd):
                seq.extend(o.args)
                continue
            elif isinstance(o, QMul):
                c = o.args[0]
                if not isinstance(c, QExpr):
                    if c is S.One:
                        s = o
                    else:
                        s = o.as_two_terms()[1]
                else:
                    c = S.One
                    s = o
            else:
                c = S.One
                s = o


            # now we have:
            # o = c*s, where
            #
            # c is a Number
            # s is a QExpr

            # let's collect terms with the same s, so e.g.
            # 2*x**2 + 3*x**2  ->  5*x**2
            if s in terms:
                terms[s] += c
            else:
                terms[s] = c


        # now let's construct new args:
        # [2*x**2, x**3, 7*x**4, pi, ...]
        newseq = []
        noncommutative = False
        for s,c in terms.items():
            # 0*s
            if c is S.Zero:
                continue
            # 1*s
            elif c is S.One:
                newseq.append(s)
            # c*s
            else:
                if isinstance(s, QMul):
                    # Mul, already keeps its arguments in perfect order.
                    # so we can simply put c in slot0 and go the fast way.
                    cs = QMul._new_rawargs(s.acts_like, s.hilbert_space, *((c,)\
                    + s.args)) #figure out rawargs F
                    newseq.append(cs)
                else:
                    # alternatively we have to call all Mul's machinery (slow)
                    newseq.append(QMul(c,s))

            noncommutative = noncommutative or not s.is_commutative

        # order args canonically
        # Currently we sort things using hashes, as it is quite fast. A better
        # solution is not to sort things at all - but this needs some more
        # fixing.
        newseq.sort(key=hash)

        if len(newseq)==1:
            return newseq[0]
        else:
            if noncommutative:
                return Expr.__new__(cls, *newseq, **{'commutative': False})
            else:
                return Expr.__new__(cls, *newseq)

    @property
    def identity(self):
        return S.Zero

    def _eval_dagger(self):
        from sympy.physics.quantum import Dagger
        newargs = [Dagger(item) for item in self.args]
        return QAdd(*newargs)

    def _eval_expand_mul(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_mul'):
                newterm = term._eval_expand_mul(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

