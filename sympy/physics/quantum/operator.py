"""Quantum mechanical operators.

TODO:

* Fix early 0 in apply_operators.
* Debug and test apply_operators.
* Get cse working with classes in this file.
* Doctests and documentation of special methods for InnerProduct, Commutator,
  AntiCommutator, represent, apply_operators.
* Decide how to handle the label of Operators. Currently, if there is only
  1 element in label, we treat it differently.
"""

from sympy import Expr 
from sympy.printing.pretty.stringpict import prettyForm
from sympy.physics.quantum.dagger import Dagger

from sympy.physics.quantum.qexpr import (
    QExpr, dispatch_method
)

__all__ = [
    'Operator',
    'HermitianOperator',
    'UnitaryOperator',
    'OuterProduct'
]

#-----------------------------------------------------------------------------
# Operators and outer products
#-----------------------------------------------------------------------------

class Operator(QExpr):
    """Base class for non-commuting quantum operators.

    An operator maps one ket to another [1]. In quantum mechanics, Hermitian
    operators correspond to observables [2].

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        operator.

    Examples
    ========

    Create an operator and examine its attributes::

        >>> from sympy.physics.quantum import Operator
        >>> from sympy import symbols, I
        >>> A = Operator('A')
        >>> A
        A
        >>> A.hilbert_space
        H
        >>> A.label
        Tuple(A)
        >>> A.is_commutative
        False

    Create another operator and do some arithmetic operations::

        >>> B = Operator('B')
        >>> C = 2*A*A + I*B
        >>> C
        I*B + 2*A**2

    Operators don't commute::

        >>> A.is_commutative
        False
        >>> B.is_commutative
        False
        >>> A*B == B*A
        False

    Polymonials of operators respect the commutation properties::

        >>> e = (A+B)**3
        >>> e.expand()
        A**2*B + B**2*A + A*B**2 + B*A**2 + A**3 + B**3 + A*B*A + B*A*B

    Operator inverses are handle symbolically::

        >>> A.inv()
        1/A
        >>> A*A.inv()
        1

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Operator
    [2] http://en.wikipedia.org/wiki/Observable
    """

    #-------------------------------------------------------------------------
    # Printing
    #-------------------------------------------------------------------------

    _label_separator = ','

    def _print_operator_name(self, printer, *args):
        return printer._print(self.__class__.__name__, *args)

    def _print_operator_name_pretty(self, printer, *args):
        return prettyForm(self.__class__.__name__)

    def _print_contents(self, printer, *args):
        if len(self.label) == 1:
            return self._print_label(printer, *args)
        else:
            return '%s(%s)' % (
                self._print_operator_name(printer, *args),
                self._print_label(printer, *args)
            )

    def _print_contents_pretty(self, printer, *args):
        if len(self.label) == 1:
            return self._print_label_pretty(printer, *args)
        else:
            pform = self._print_operator_name_pretty(printer, *args)
            label_pform = self._print_label_pretty(printer, *args)
            label_pform = prettyForm(
                *label_pform.parens(left='(', right=')')
            )
            pform = prettyForm(*pform.right((label_pform)))
            return pform

    #-------------------------------------------------------------------------
    # _eval_* methods
    #-------------------------------------------------------------------------

    def _eval_commutator(self, other, **options):
        """Evaluate [self, other] if known, return None if not known."""
        return dispatch_method(self, '_eval_commutator', other, **options)

    def _eval_anticommutator(self, other, **options):
        """Evaluate [self, other] if known."""
        return dispatch_method(self, '_eval_anticommutator', other, **options)

    #-------------------------------------------------------------------------
    # Operator application
    #-------------------------------------------------------------------------

    def _apply_operator(self, ket, **options):
        from sympy.physics.quantum.state import KetBase
        if not isinstance(ket, KetBase):
            raise TypeError('KetBase expected, got: %r' % ket)
        return dispatch_method(self, '_apply_operator', ket, **options)

    def matrix_element(self, *args):
        raise NotImplementedError('matrix_elements is not defined')

    #-------------------------------------------------------------------------
    # Printing
    #-------------------------------------------------------------------------

    def inverse(self):
        return self._eval_inverse()

    inv = inverse

    def _eval_inverse(self):
        # TODO: make non-commutative Exprs print powers using A**-1, not 1/A.
        return self**(-1)


class HermitianOperator(Operator):
    """A Hermitian operator that satisfies H == Dagger(H).

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        operator.

    Examples
    ========

    >>> from sympy.physics.quantum import Dagger, HermitianOperator
    >>> H = HermitianOperator('H')
    >>> Dagger(H)
    H
    """

    def _eval_dagger(self):
        return self


class UnitaryOperator(Operator):
    """A unitary operator that satisfies U*Dagger(U) == 1.

    Parameters
    ==========
    label : tuple, sympy.core.containers.Tuple
        The list of numbers or parameters that uniquely specify the
        operator.

    Examples
    ========

    >>> from sympy.physics.quantum import Dagger, UnitaryOperator
    >>> U = UnitaryOperator('U')
    >>> U*Dagger(U)
    1
    """

    def _eval_dagger(self):
        return self._eval_inverse()


class OuterProduct(Operator):
    """An unevaluated outer product between a ket and kra.

    This constructs an outer product between any subclass of KetBase and
    BraBase as |a><b|. An OuterProduct inherits from Operator as they act as
    operators in quantum expressions.  For reference see [1].

    Parameters
    ==========
    ket : KetBase or subclass
        The ket on the left side of the outer product.
    bar : BraBase or subclass
        The bra on the right side of the outer product.

    Examples
    ========

    Create a simple outer product by hand and take its dagger::

        >>> from sympy.physics.quantum import Ket, Bra, OuterProduct, Dagger
        >>> from sympy.physics.quantum import Operator

        >>> k = Ket('k')
        >>> b = Bra('b')
        >>> op = OuterProduct(k, b)
        >>> op
        |k><b|
        >>> op.hilbert_space
        H
        >>> op.ket
        |k>
        >>> op.bra
        <b|
        >>> Dagger(op)
        |b><k|

    In simple products of kets and bras outer products will be automatically
    identified and created::

        >>> k*b
        |k><b|

    But in more complex expressions, outer products are not automatically
    created::

        >>> A = Operator('A')
        >>> A*k*b
        A*|k>*<b|

    A user can force the creation of an outer product in a complex expression
    by using parentheses to group the ket and bra::

        >>> A*(k*b)
        A*|k><b|

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Outer_product
    """

    def __new__(cls, ket, bra, **old_assumptions):
        from sympy.physics.quantum.state import KetBase, BraBase
        if not isinstance(ket, KetBase):
            raise TypeError('KetBase subclass expected, got: %r' % ket)
        if not isinstance(bra, BraBase):
            raise TypeError('BraBase subclass expected, got: %r' % ket)
        if not ket.dual_class == bra.__class__:
            raise TypeError(
                'ket and bra are not dual classes: %r, %r' % \
                (ket.__class__, bra.__class__)
            )
        # TODO: make sure the hilbert spaces of the bra and ket are compatible
        obj = Expr.__new__(cls, *(ket, bra), **{'commutative': False})
        obj.hilbert_space = ket.hilbert_space
        return obj

    @property
    def ket(self):
        """Return the ket on the left side of the outer product."""
        return self.args[0]

    @property
    def bra(self):
        """Return the bra on the right side of the outer product."""
        return self.args[1]

    def _eval_dagger(self):
        return OuterProduct(Dagger(self.bra), Dagger(self.ket))

    def _sympyrepr(self, printer, *args):
        return '%s(%s,%s)' % (self.__class__.__name__, 
            printer._print(self.ket, *args), printer._print(self.bra, *args))

    def _sympystr(self, printer, *args):
        return str(self.ket)+str(self.bra)

    def _pretty(self, printer, *args):
        pform = self.ket._pretty(printer, *args)
        return prettyForm(*pform.right(self.bra._pretty(printer, *args)))

    def _latex(self, printer, *args):
        k = printer._print(self.ket, *args)
        b = printer._print(self.bra, *args)
        return k+b
