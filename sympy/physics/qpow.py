from sympy import Expr, Pow, S, sympify
from sympy.physics.quantum import InnerProduct, OuterProduct, Operator,\
Dagger
from sympy.core.decorators import call_highest_priority
from sympy.physics.qexpr import QuantumError, QExpr
from sympy.printing.str import sstr

class QPow(QExpr):
    """
        A class for the operator: ** (exponent) for quantum objects.
    """

    def __new__(cls, base, exp):
        base = sympify(base)
        exp = sympify(exp)
        return cls._rules_QPow(base, exp)

    @classmethod
    def _rules_QPow(cls, base, exp):
        from sympy.physics.qadd import QAdd
        from sympy.physics.qmul import QMul
        if not isinstance(base, QExpr):
            if not isinstance(exp, QExpr):
                return Pow(base, exp)
            elif issubclass(exp.acts_like, (Operator, OuterProduct,\
            InnerProduct)):
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = exp.hilbert_space
                ret.acts_like = exp.acts_like
                return ret
        elif not isinstance(exp, QExpr):
            if issubclass(base.acts_like, (Operator, OuterProduct,\
            InnerProduct)):
                if exp == S.Zero:
                    return S.One
                elif exp == S.One:
                    return base
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = base.hilbert_space
                ret.acts_like = base.acts_like
                return ret
        elif issubclass(exp.acts_like, InnerProduct) and issubclass(\
        base.acts_like, (InnerProduct, Operator)):
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = base.hilbert_space
                ret.acts_like = base.acts_like
                return ret
        elif issubclass(base.acts_like, InnerProduct) and issubclass(\
        exp.acts_like, (InnerProduct, Operator)):
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = exp.hilbert_space
                ret.acts_like = exp.acts_like
                return ret

        #make a pretty error message if you have left everything a mess
        if hasattr(exp, 'acts_like'):
            expname = exp.acts_like.__name__
        else:
            expname = exp.__class__.__name__

        if hasattr(base, 'acts_like'):
            basename = base.acts_like.__name__
        else:
            basename = base.__class__.__name__

        raise QuantumError("Can't do (%s)**(%s)" % (basename, expname))

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]

    def _eval_dagger(self):
        return QPow(Dagger(self.base), Dagger(self.exp))

    def _sympystr(self, printer, *args):
        return '(' + sstr(self.base) + ')' '**' + sstr(self.exp)

    def _pretty(self, printer, *args):
        return printer._print(self.args[0], *args)**printer._print(\
        self.args[1], *args)
