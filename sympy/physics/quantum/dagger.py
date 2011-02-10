"""Hermitian conjugation."""

from sympy import Expr, sympify, Add, Mul, Matrix, Pow

from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.matrixutils import (
    numpy_ndarray, scipy_sparse_matrix, matrix_dagger
)

__all__ = [
    'Dagger'
]


class Dagger(Expr):
    """General Hermitian conjugate operation.

    For matrices this operation is equivalent to transpose and complex
    conjugate [1].

    Parameters
    ==========
    arg : Expr
        The sympy expression that we want to take the dagger of.

    Examples
    ========

    Daggering various quantum objects:

        >>> from sympy.physics.quantum.dagger import Dagger
        >>> from sympy.physics.quantum.state import Ket, Bra
        >>> from sympy.physics.quantum.operator import Operator
        >>> Dagger(Ket('psi'))
        <psi|
        >>> Dagger(Bra('phi'))
        |phi>
        >>> Dagger(Operator('A'))
        Dagger(A)

    Inner and outer products::

        >>> from sympy.physics.quantum import InnerProduct, OuterProduct
        >>> Dagger(InnerProduct(Bra('a'), Ket('b')))
        <b|a>
        >>> Dagger(OuterProduct(Ket('a'), Bra('b')))
        |b><a|

    Powers, sums and products::

        >>> A = Operator('A')
        >>> B = Operator('B')
        >>> Dagger(A*B)
        Dagger(B)*Dagger(A)
        >>> Dagger(A+B)
        Dagger(A) + Dagger(B)
        >>> Dagger(A**2)
        Dagger(A)**2

    Dagger also seamlessly handles complex numbers and matrices::

        >>> from sympy import Matrix, I
        >>> m = Matrix([[1,I],[2,I]])
        >>> m
        [1, I]
        [2, I]
        >>> Dagger(m)
        [ 1,  2]
        [-I, -I]

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Hermitian_transpose
    """

    def __new__(cls, arg, **old_assumptions):
        # Return the dagger of a sympy Matrix immediately.
        if isinstance(arg, (Matrix, numpy_ndarray, scipy_sparse_matrix)):
            return matrix_dagger(arg)
        arg = sympify(arg)
        r = cls.eval(arg)
        if isinstance(r, Expr):
            return r
        #make unevaluated dagger commutative or non-commutative depending on arg
        if arg.is_commutative:
            obj = Expr.__new__(cls, arg, **{'commutative':True})
        else:
            obj = Expr.__new__(cls, arg, **{'commutative':False})
        if isinstance(obj, QExpr):
            obj.hilbert_space = arg.hilbert_space
        return obj

    @classmethod
    def eval(cls, arg):
        """Evaluates the Dagger instance."""
        from sympy.physics.quantum.operator import Operator
        try:
            d = arg._eval_dagger()
        except (NotImplementedError, AttributeError):
            if isinstance(arg, Expr):
                if isinstance(arg, Operator):
                    # Operator without _eval_dagger
                    return None
                if arg.is_Add:
                    return Add(*[Dagger(i) for i in arg.args])
                if arg.is_Mul:
                    return Mul(*[Dagger(i) for i in reversed(arg.args)])
                if arg.is_Pow:
                    return Pow(Dagger(arg.args[0]),arg.args[1])
                else:
                    if arg.is_Number or arg.is_Function or arg.is_Derivative\
                                     or arg.is_Integer or arg.is_NumberSymbol\
                                     or arg.is_complex or arg.is_integer\
                                     or arg.is_real or arg.is_number:
                        return arg.conjugate()
                    else:
                        return None
            else:
                return None
        else:
            return d

    def _eval_subs(self, old, new):
        r = Dagger(self.args[0].subs(old, new))
        return r

    def _eval_dagger(self):
        return self.args[0]

    def _sympyrepr(self, printer, *args):
        arg0 = printer._print(self.args[0], *args)
        return '%s(%s)' % (self.__class__.__name__, arg0)

    def _sympystr(self, printer, *args):
        arg0 = printer._print(self.args[0], *args)
        return '%s(%s)' % (self.__class__.__name__, arg0)

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = printer._print(self.args[0], *args)
        pform = pform**prettyForm(u'\u2020')
        return pform

    def _latex(self, printer, *args):
        arg = printer._print(self.args[0])
        return '%s^{\\dag}' % arg

