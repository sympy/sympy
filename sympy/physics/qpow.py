from sympy import Expr, Pow, S, sympify
from sympy.printing.str import sstr
from sympy.printing.pretty.stringpict import prettyForm

from sympy.physics.qexpr import QuantumError, QExpr


class QPow(QExpr):
    """A quantum version of the Pow class for base**exp operations."""

    def __new__(cls, base, exp):
        base = sympify(base)
        exp = sympify(exp)
        return cls._apply_rules(base, exp)

    @classmethod
    def eval(cls, base, exp):
        return cls._apply_rules(base, exp)

    @classmethod
    def _apply_rules(cls, base, exp):
        """Apply rules to transform the base and exp and validate the expr."""

        from sympy.physics.quantum import InnerProduct, OuterProduct, Operator

        allowed = (Operator, OuterProduct, InnerProduct)
        base_is_qexpr = isinstance(base, QExpr)
        exp_is_qexpr = isinstance(exp, QExpr)

        if not base_is_qexpr:
            if not exp_is_qexpr:
                return Pow(base, exp)
            elif issubclass(exp.acts_like, allowed):
                result = Expr.__new__(cls, base, exp)
                result.hilbert_space = exp.hilbert_space
                result.acts_like = exp.acts_like
                return result
        elif not exp_is_qexpr:
            if issubclass(base.acts_like, allowed):
                if exp == S.Zero:
                    return S.One
                elif exp == S.One:
                    return base
                result = Expr.__new__(cls, base, exp)
                result.hilbert_space = base.hilbert_space
                result.acts_like = base.acts_like
                return result
        elif issubclass(exp.acts_like, InnerProduct) and \
             issubclass(base.acts_like, (InnerProduct, Operator)):
                result = Expr.__new__(cls, base, exp)
                result.hilbert_space = base.hilbert_space
                result.acts_like = base.acts_like
                return result
        elif issubclass(base.acts_like, InnerProduct) and \
             issubclass(exp.acts_like, (InnerProduct, Operator)):
                result = Expr.__new__(cls, base, exp)
                result.hilbert_space = exp.hilbert_space
                result.acts_like = exp.acts_like
                return result

        # Make a pretty error message if the base**exp can't be done.
        if hasattr(exp, 'acts_like'):
            expname = exp.acts_like.__name__
        else:
            expname = exp.__class__.__name__

        if hasattr(base, 'acts_like'):
            basename = base.acts_like.__name__
        else:
            basename = base.__class__.__name__

        raise TypeError("Can't exponentiate: %s**%s" % (basename, expname))

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]

    def _eval_dagger(self):
        from sympy.physics.quantum import Dagger
        return QPow(Dagger(self.base), Dagger(self.exp))

    def _sympystr(self, printer, *args):
        from sympy.physics.qmul import QMul
        from sympy.physics.qadd import QAdd
        base_str = printer._print(self.base, *args)
        exp_str = printer._print(self.exp, *args)
        if isinstance(self.base, (QMul, QAdd)):
            return '(' + base_str + ')' '**' + exp_str
        else:
            return base_str + '**' + exp_str

    def _pretty(self, printer, *args):
        from sympy.physics.qmul import QMul
        from sympy.physics.qadd import QAdd
        base_pform = printer._print(self.base, *args)
        if isinstance(self.base, (QMul, QAdd)):
            base_pform = prettyForm(*base_pform.parens(left='(', right=')'))
        exp_pform = printer._print(self.exp, *args)
        return base_pform**exp_pform
