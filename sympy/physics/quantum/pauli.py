"""Pauli operators and states"""

from sympy import I, Add, Mul, Pow, Integer, exp, sqrt, conjugate
from sympy.physics.quantum import (Operator, Commutator, AntiCommutator,
                                   Dagger, IdentityOperator, Ket, Bra)
from sympy.physics.quantum import HilbertSpace, ComplexSpace
from sympy.physics.quantum.boson import BosonOp
from sympy.matrices import Matrix
from sympy.functions.special.tensor_functions import KroneckerDelta

class SigmaOpBase(Operator):
    """Pauli sigma x operator"""

    def _eval_commutator_BosonOp(self, other, **hints):
        return Integer(0)


class SigmaX(SigmaOpBase):
    """Pauli sigma x operator"""

    def __new__(cls, *args, **hints):
        return Operator.__new__(cls, *args, **hints)

    def _eval_commutator_SigmaY(self, other, **hints):
        return 2 * I * SigmaZ()

    def _eval_commutator_SigmaZ(self, other, **hints):
        return - 2 * I * SigmaY()

    def _eval_commutator_BosonOp(self, other, **hints):
        return Integer(0)

    def _eval_anticommutator_SigmaY(self, other, **hints):
        return Integer(0)

    def _eval_anticommutator_SigmaZ(self, other, **hints):
        return Integer(0)

    def _eval_adjoint(self):
        return self

    def _print_contents_latex(self, printer, *args):
        return r'{\sigma_x}'

    def _print_contents(self, printer, *args):
        return 'SigmaX()'

    def _eval_power(b, e):
        if e.is_Integer and e.is_positive:
            return SigmaX().__pow__(int(e) % 2)

    def __mul__(self, other):

        if isinstance(other, SigmaX):
            return IdentityOperator(2)

        if isinstance(other, SigmaY):
            return I * SigmaZ()

        if isinstance(other, SigmaZ):
            return - I * SigmaY()

        if isinstance(other, SigmaMinus):
            return (IdentityOperator(2) + SigmaZ())/2

        if isinstance(other, SigmaPlus):
            return (IdentityOperator(2) - SigmaZ())/2

        if isinstance(other, Mul):
            args1 = tuple(arg for arg in other.args if arg.is_commutative)
            args2 = tuple(arg for arg in other.args if not arg.is_commutative)
            x = self
            for y in args2:
                x = x * y
            return Mul(*args1) * x

        return Mul(self, other)

    def _represent_default_basis(self, **options):
        format = options.get('format', 'sympy')
        if format == 'sympy':
            return Matrix([[0, 1], [1, 0]])
        else:
            raise NotImplementedError('Representation in format ' +
                                      format + ' not implemented.')


class SigmaY(SigmaOpBase):
    """Pauli sigma y operator"""

    def __new__(cls, *args, **hints):
        return Operator.__new__(cls, *args)

    def _eval_commutator_SigmaZ(self, other, **hints):
        return 2 * I * SigmaX()

    def _eval_commutator_SigmaX(self, other, **hints):
        return - 2 * I * SigmaZ()

    def _eval_anticommutator_SigmaX(self, other, **hints):
        return Integer(0)

    def _eval_anticommutator_SigmaZ(self, other, **hints):
        return Integer(0)

    def _eval_adjoint(self):
        return self

    def _print_contents_latex(self, printer, *args):
        return r'{\sigma_y}'

    def _print_contents(self, printer, *args):
        return 'SigmaY()'

    def _eval_power(b, e):
        if e.is_Integer and e.is_positive:
            return SigmaY().__pow__(int(e) % 2)

    def __mul__(self, other):

        if isinstance(other, SigmaX):
            return - I * SigmaZ()

        if isinstance(other, SigmaY):
            return IdentityOperator(2)

        if isinstance(other, SigmaZ):
            return I * SigmaX()

        if isinstance(other, SigmaMinus):
            return -I * (IdentityOperator(2) + SigmaZ())/2

        if isinstance(other, SigmaPlus):
            return I * (IdentityOperator(2) - SigmaZ())/2

        if isinstance(other, Mul):
            args1 = tuple(arg for arg in other.args if arg.is_commutative)
            args2 = tuple(arg for arg in other.args if not arg.is_commutative)
            x = self
            for y in args2:
                x = x * y
            return Mul(*args1) * x

        return Mul(self, other)

    def _represent_default_basis(self, **options):
        format = options.get('format', 'sympy')
        if format == 'sympy':
            return Matrix([[0, -I], [I, 0]])
        else:
            raise NotImplementedError('Representation in format ' +
                                      format + ' not implemented.')


class SigmaZ(SigmaOpBase):
    """Pauli sigma z operator"""

    def __new__(cls, *args, **hints):
        return Operator.__new__(cls, *args)

    def _eval_commutator_SigmaX(self, other, **hints):
        return 2 * I * SigmaY()

    def _eval_commutator_SigmaY(self, other, **hints):
        return - 2 * I * SigmaX()

    def _eval_anticommutator_SigmaX(self, other, **hints):
        return Integer(0)

    def _eval_anticommutator_SigmaY(self, other, **hints):
        return Integer(0)

    def _eval_adjoint(self):
        return self

    def _print_contents_latex(self, printer, *args):
        return r'{\sigma_z}'

    def _print_contents(self, printer, *args):
        return 'SigmaZ()'

    def _eval_power(b, e):
        if e.is_Integer and e.is_positive:
            return SigmaZ().__pow__(int(e) % 2)

    def __mul__(self, other):

        if isinstance(other, SigmaX):
            return I * SigmaY()

        if isinstance(other, SigmaY):
            return - I * SigmaX()

        if isinstance(other, SigmaZ):
            return IdentityOperator(2)

        if isinstance(other, SigmaMinus):
            return -(SigmaX() - I * SigmaY())/2  # - SigmaMinus()

        if isinstance(other, SigmaPlus):
            return (SigmaX() + I * SigmaY())/2  # SigmaPlus()

        if isinstance(other, Mul):
            args1 = tuple(arg for arg in other.args if arg.is_commutative)
            args2 = tuple(arg for arg in other.args if not arg.is_commutative)
            x = self
            for y in args2:
                x = x * y
            return Mul(*args1) * x

        return Mul(self, other)

    def _represent_default_basis(self, **options):
        format = options.get('format', 'sympy')
        if format == 'sympy':
            return Matrix([[1, 0], [0, -1]])
        else:
            raise NotImplementedError('Representation in format ' +
                                      format + ' not implemented.')


class SigmaMinus(SigmaOpBase):
    """Pauli sigma minus operator"""

    def __new__(cls, *args, **hints):
        return Operator.__new__(cls, *args)

    def _eval_commutator_SigmaX(self, other, **hints):
        return -SigmaZ()

    def _eval_commutator_SigmaY(self, other, **hints):
        return I * SigmaZ()

    def _eval_commutator_SigmaZ(self, other, **hints):
        return 2 * self

    def _eval_commutator_SigmaMinus(self, other, **hints):
        return SigmaZ()

    def _eval_anticommutator_SigmaZ(self, other, **hints):
        return Integer(0)

    def _eval_anticommutator_SigmaX(self, other, **hints):
        return IdentityOperator(2)

    def _eval_anticommutator_SigmaY(self, other, **hints):
        return - I * IdentityOperator(2)

    def _eval_anticommutator_SigmaPlus(self, other, **hints):
        return IdentityOperator(2)

    def _eval_adjoint(self):
        return SigmaPlus()

    def _eval_power(b, e):
        if e.is_Integer and e.is_positive:
            return Integer(0)

    def __mul__(self, other):
        if isinstance(other, SigmaX):
            return (IdentityOperator(2) - SigmaZ())/2

        if isinstance(other, SigmaY):
            return - I * (IdentityOperator(2) - SigmaZ())/2

        if isinstance(other, SigmaZ):
            return (SigmaX() - I * SigmaY())/2  # SigmaMinus()

        if isinstance(other, SigmaMinus):
            return Integer(0)

        if isinstance(other, SigmaPlus):
            return (IdentityOperator(2) - SigmaZ())/2

        if isinstance(other, Mul):
            args1 = tuple(arg for arg in other.args if arg.is_commutative)
            args2 = tuple(arg for arg in other.args if not arg.is_commutative)
            x = self
            for y in args2:
                x = x * y
            return Mul(*args1) * x

        return Mul(self, other)

    def _print_contents_latex(self, printer, *args):
        return r'{\sigma_-}'

    def _print_contents(self, printer, *args):
        return 'SigmaMinus()'

    def _represent_default_basis(self, **options):
        format = options.get('format', 'sympy')
        if format == 'sympy':
            return Matrix([[0, 0], [1, 0]])
        else:
            raise NotImplementedError('Representation in format ' +
                                      format + ' not implemented.')


class SigmaPlus(SigmaOpBase):
    """Pauli sigma plus operator"""

    def __new__(cls, *args, **hints):
        return Operator.__new__(cls, *args)

    def _eval_commutator_SigmaX(self, other, **hints):
        return SigmaZ()

    def _eval_commutator_SigmaY(self, other, **hints):
        return I * SigmaZ()

    def _eval_commutator_SigmaZ(self, other, **hints):
        return -2 * self

    def _eval_commutator_SigmaMinus(self, other, **hints):
        return SigmaZ()

    def _eval_anticommutator_SigmaZ(self, other, **hints):
        return Integer(0)

    def _eval_anticommutator_SigmaX(self, other, **hints):
        return IdentityOperator(2)

    def _eval_anticommutator_SigmaY(self, other, **hints):
        return I * IdentityOperator(2)

    def _eval_anticommutator_SigmaMinus(self, other, **hints):
        return IdentityOperator(2)

    def _eval_adjoint(self):
        return SigmaMinus()

    def _eval_mul(self, other):
        return self * other

    def _eval_power(b, e):
        if e.is_Integer and e.is_positive:
            return Integer(0)

    def __mul__(self, other):

        if other == IdentityOperator(2):
            return self

        if isinstance(other, SigmaX):
            return (IdentityOperator(2) + SigmaZ())/2

        if isinstance(other, SigmaY):
            return I * (IdentityOperator(2) + SigmaZ())/2

        if isinstance(other, SigmaZ):
            return -(SigmaX() + I * SigmaY())/2  # -SigmaPlus()

        if isinstance(other, SigmaMinus):
            return (IdentityOperator(2) + SigmaZ())/2

        if isinstance(other, SigmaPlus):
            return Integer(0)

        if isinstance(other, Mul):
            args1 = tuple(arg for arg in other.args if arg.is_commutative)
            args2 = tuple(arg for arg in other.args if not arg.is_commutative)
            x = self
            for y in args2:
                x = x * y
            return Mul(*args1) * x

        return Mul(self, other)

    def _print_contents_latex(self, printer, *args):
        return r'{\sigma_+}'

    def _print_contents(self, printer, *args):
        return 'SigmaPlus()'

    def _represent_default_basis(self, **options):
        format = options.get('format', 'sympy')
        if format == 'sympy':
            return Matrix([[0, 1], [0, 0]])
        else:
            raise NotImplementedError('Representation in format ' +
                                      format + ' not implemented.')


class SigmaZKet(Ket):
    """Ket for a two-level system quantum system.

    Parameters
    ==========

    n : Number
        The state number (0 or 1).

    """

    def __new__(cls, n):
        if n not in [0, 1]:
            raise ValueError("n must be 0 or 1")
        return Ket.__new__(cls, n)

    @property
    def n(self):
        return self.label[0]

    @classmethod
    def dual_class(self):
        return SigmaZBra

    @classmethod
    def _eval_hilbert_space(cls, label):
        return ComplexSpace(2)

    def _eval_innerproduct_SigmaZBra(self, bra, **hints):
        return KroneckerDelta(self.n, bra.n)

    def _apply_operator_SigmaZ(self, op, **options):
        if self.n == 0:
            return self
        else:
            return Integer(-1) * self

    def _apply_operator_SigmaX(self, op, **options):
        return SigmaZKet(1) if self.n == 0 else SigmaZKet(0)

    def _apply_operator_SigmaY(self, op, **options):
        return I * SigmaZKet(1) if self.n == 0 else (-I) * SigmaZKet(0)

    def _apply_operator_SigmaMinus(self, op, **options):
        if self.n == 0:
            return SigmaZKet(1)
        else:
            return Integer(0)

    def _apply_operator_SigmaPlus(self, op, **options):
        if self.n == 0:
            return Integer(0)
        else:
            return SigmaZKet(0)

    def _represent_default_basis(self, **options):
        format = options.get('format', 'sympy')
        if format == 'sympy':
            return Matrix([[1], [0]]) if self.n == 0 else Matrix([[0], [1]])
        else:
            raise NotImplementedError('Representation in format ' +
                                      format + ' not implemented.')


class SigmaZBra(Bra):
    """Bra for a two-level quantum system.

    Parameters
    ==========

    n : Number
        The state number (0 or 1).

    """

    def __new__(cls, n):
        if n not in [0, 1]:
            raise ValueError("n must be 0 or 1")
        return Bra.__new__(cls, n)

    @property
    def n(self):
        return self.label[0]

    @classmethod
    def dual_class(self):
        return SigmaZKet
