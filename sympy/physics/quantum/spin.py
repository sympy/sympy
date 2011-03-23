"""Quantum mechanical angular momemtum.

TODO:
* Improve the algorithm for the Wigner D-Function.
* Fix printing of Wigner D-Function.
* Implement inner products using Wigner D-Function.
* Implement rewrite logic for everything.
* Test
"""

from sympy import I, Symbol, S, Integer, Rational, Matrix, sqrt, sympify
from sympy import exp, cos, sin, diff, factorial
from sympy.printing.pretty.stringpict import prettyForm, stringPict
from sympy.matrices.matrices import zeros

from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.operator import (
    HermitianOperator, Operator, UnitaryOperator
)
from sympy.physics.quantum.state import State, Ket, Bra
from sympy.physics.quantum.kronecker import KroneckerDelta
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.hilbert import ComplexSpace


__all__ = [
    'm_values',
    'Jplus',
    'Jminus',
    'Jx',
    'Jy',
    'Jz',
    'J2',
    'JzKet',
    'JzBra',
    'JxKet',
    'JxBra',
    'JyKet',
    'JyBra',
    'Rotation'
]

def m_values(j):
    j = sympify(j)
    size = 2*j + 1
    if not size.is_Integer or not size > 0:
        raise ValueError(
            'Only integer or half-integer values allowed for j, got: : %r' % j
        )
    return size, [j-i for i in range(int(2*j+1))]


#-----------------------------------------------------------------------------
# SpinOperators
#-----------------------------------------------------------------------------


class SpinOpBase(object):
    """Base class for spin operators."""

    @classmethod
    def _eval_hilbert_space(cls, label):
        # We consider all j values so our space is infinite.
        return ComplexSpace(S.Infinity)

    @property
    def name(self):
        return self.args[0]

    def _print_contents(self, printer, *args):
        return '%s%s' % (unicode(self.name), self._coord)

    # def _sympyrepr(self, printer, *args):
    #     return '%s(%s)' % (
    #         self.__class__.__name__, printer._print(self.label,*args)
    #

    def _print_contents_pretty(self, printer, *args):
        a = stringPict(unicode(self.name))
        b = stringPict(self._coord)
        return self._print_subscript_pretty(a, b)

    def _print_contents_latex(self, printer, *args):
        return r'%s_%s' % ((unicode(self.name), self._coord))

    def _represent_base(self, basis, **options):
        j = options.get('j', Rational(1,2))
        size, mvals = m_values(j)
        result = zeros((size, size))
        for p in range(size):
            for q in range(size):
                me = self.matrix_element(j, mvals[p], j, mvals[q])
                result[p, q] = me
        return result


class JplusOp(SpinOpBase, Operator):
    """The J+ operator."""

    _coord = '+'

    def _eval_commutator_JminusOp(self, other):
        return 2*hbar*JzOp(self.name)

    def _apply_operator_JzKet(self, ket, **options):
        j = ket.j
        m = ket.m
        if m.is_Number and j.is_Number:
            if m >= j:
                return S.Zero
        return hbar*sqrt(j*(j+S.One)-m*(m+S.One))*JzKet(j, m+S.One)

    def matrix_element(self, j, m, jp, mp):
        result = hbar*sqrt(j*(j+S.One)-mp*(mp+S.One))
        result *= KroneckerDelta(m, mp+1)
        result *= KroneckerDelta(j, jp)
        return result

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        return self._represent_base(basis, **options)

    def _eval_rewrite_as_xyz(self, *args):
        return JxOp(args[0]) + I*JyOp(args[0])


class JminusOp(SpinOpBase, Operator):
    """The J- operator."""

    _coord = '-'

    def _apply_operator_JzKet(self, ket, **options):
        j = ket.j
        m = ket.m
        if m.is_Number and j.is_Number:
            if m <= -j:
                return S.Zero
        return hbar*sqrt(j*(j+S.One)-m*(m-S.One))*JzKet(j, m-S.One)

    def matrix_element(self, j, m, jp, mp):
        result = hbar*sqrt(j*(j+S.One)-mp*(mp-S.One))
        result *= KroneckerDelta(m, mp-1)
        result *= KroneckerDelta(j, jp)
        return result

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        return self._represent_base(basis, **options)

    def _eval_rewrite_as_xyz(self, *args):
        return JxOp(args[0]) - I*JyOp(args[0])


class JxOp(SpinOpBase, HermitianOperator):
    """The Jx operator."""

    _coord = 'x'

    def _eval_commutator_JyOp(self, other):
        return I*hbar*JzOp(self.name)

    def _eval_commutator_JzOp(self, other):
        return -I*hbar*JyOp(self.name)

    def _apply_operator_JzKet(self, ket, **options):
        jp = JplusOp(self.name)._apply_operator_JzKet(ket, **options)
        jm = JminusOp(self.name)._apply_operator_JzKet(ket, **options)
        return (jp + jm)/Integer(2)

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        jp = JplusOp(self.name)._represent_JzOp(basis, **options)
        jm = JminusOp(self.name)._represent_JzOp(basis, **options)
        return (jp + jm)/Integer(2)

    def _eval_rewrite_as_plusminus(self, *args):
        return (JplusOp(args[0]) + JminusOp(args[0]))/2


class JyOp(SpinOpBase, HermitianOperator):
    """The Jy operator."""

    _coord = 'y'

    def _eval_commutator_JzOp(self, other):
        return I*hbar*JxOp(self.name)

    def _eval_commutator_JxOp(self, other):
        return -I*hbar*J2Op(self.name)

    def _apply_operator_JzKet(self, ket, **options):
        jp = JplusOp(self.name)._apply_operator_JzKet(ket, **options)
        jm = JminusOp(self.name)._apply_operator_JzKet(ket, **options)
        return (jp - jm)/(Integer(2)*I)

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        jp = JplusOp(self.name)._represent_JzOp(basis, **options)
        jm = JminusOp(self.name)._represent_JzOp(basis, **options)
        return (jp - jm)/(Integer(2)*I)

    def _eval_rewrite_as_plusminus(self, *args):
        return (JplusOp(args[0]) - JminusOp(args[0]))/(2*I)


class JzOp(SpinOpBase, HermitianOperator):
    """The Jz operator."""

    _coord = 'z'

    def _eval_commutator_JxOp(self, other):
        return I*hbar*JyOp(self.name)

    def _eval_commutator_JyOp(self, other):
        return -I*hbar*JxOp(self.name)

    def _eval_commutator_JplusOp(self, other):
        return hbar*JplusOp(self.name)

    def _eval_commutator_JminusOp(self, other):
        return -hbar*JminusOp(self.name)

    def _apply_operator_JzKet(self, ket, **options):
        return (hbar*ket.m)*ket

    def matrix_element(self, j, m, jp, mp):
        result = hbar*mp
        result *= KroneckerDelta(m, mp)
        result *= KroneckerDelta(j, jp)
        return result

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        return self._represent_base(basis, **options)


class J2Op(SpinOpBase, HermitianOperator):
    """The J^2 operator."""

    _coord = '2'

    def _eval_commutator_JxOp(self, other):
        return S.Zero

    def _eval_commutator_JyOp(self, other):
        return S.Zero

    def _eval_commutator_JzOp(self, other):
        return S.Zero

    def _eval_commutator_JplusOp(self, other):
        return S.Zero

    def _eval_commutator_JminusOp(self, other):
        return S.Zero

    def _apply_operator_JzKet(self, ket, **options):
        j = ket.j
        return hbar**2**j*(j+1)*ket

    def matrix_element(self, j, m, jp, mp):
        result = (hbar**2)*j*(j+1)
        result *= KroneckerDelta(m, mp)
        result *= KroneckerDelta(j, jp)
        return result

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        return self._represent_base(basis, **options)

    def _pretty(self, printer, *args):
        a = stringPict('J')
        b = stringPict('2')
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.above(top))

    def _latex(self, printer, *args):
        return r'%s^2' % str(self.name)

    def _eval_rewrite_as_xyz(self, *args):
        return JxOp(args[0])**2 + JyOp(args[0])**2 + JzOp(args[0])**2

    def _eval_rewrite_as_plusminus(self, *args):
        a = args[0]
        return JzOp(a)**2 +\
            Rational(1,2)*(JplusOp(a)*JminusOp(a) + JminusOp(a)*JplusOp(a))


class Rotation(UnitaryOperator):
    """Wigner D operator in terms of Euler angles."""

    @classmethod
    def _eval_args(cls, args):
        args = QExpr._eval_args(args)
        if len(args) != 3:
            raise ValueError('3 Euler angles required, got: %r' % args)
        return args

    @classmethod
    def _eval_hilbert_space(cls, label):
        # We consider all j values so our space is infinite.
        return ComplexSpace(S.Infinity)

    @property
    def alpha(self):
        return self.label[0]

    @property
    def beta(self):
        return self.label[1]

    @property
    def gamma(self):
        return self.label[2]

    def _print_operator_name(self, printer, *args):
        return printer._print('R', *args)

    def _print_operator_name_pretty(self, printer, *args):
        return prettyForm(u"\u211B" + u" ")

    def _eval_inverse(self):
        return Rotation((-self.gamma, -self.beta, -self.alpha))

    @classmethod
    def D(cls, j, m, mp, alpha, beta, gamma):
        """Wigner-D function."""
        result = exp(-I*m*alpha)*exp(-I*mp*gamma)
        result *= cls.d(j, m, mp, beta)
        return result

    @classmethod
    def d(cls, j, m, mp, beta):
        """Wigner's lowercase d function."""
        # TODO: This does not do a good job of simplifying the trig functions.
        # The Jacobi Polynomial expansion is probably best as it is a simple
        # sum. But, this version does give correct answers and uses
        # Eq. 7 in Section 4.3.2 of Varshalovich.
        cosbeta = Symbol('cosbeta')
        x = 1-cosbeta
        y = 1+cosbeta
        jmmp = j-mp
        jpm = j+m
        jpmp = j+mp
        jmm = j-m
        mpmm2 = (mp-m)/2
        mpmp2 = (m+mp)/2
        result = (-1)**jmmp
        result *= 2**(-j)
        result *= sqrt(factorial(jpm)/(factorial(jmm)*factorial(jpmp)*factorial(jmmp)))
        result *= (x**mpmm2)*(y**(-mpmp2))
        result *= diff((x**jmmp)*(y**jpmp), cosbeta, jmm)
        if (2*j).is_Integer and not j.is_Integer:
            result = result.subs(1+cosbeta, 2*cos(beta/2)**2)
            result = result.subs(1-cosbeta, 2*sin(beta/2)**2)
        else:
            result = result.subs(cosbeta, cos(beta))
        return result

    def matrix_element(self, j, m, jp, mp):
        result = self.__class__.D(
            jp, m, mp, self.alpha, self.beta, self.gamma
        )
        result *= KroneckerDelta(j,jp)
        return result

    def _represent_base(self, basis, **options):
        j = sympify(options.get('j', Rational(1,2)))
        size, mvals = m_values(j)
        result = zeros((size, size))
        for p in range(size):
            for q in range(size):
                me = self.matrix_element(j, mvals[p], j, mvals[q])
                result[p, q] = me
        return result

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        return self._represent_base(basis, **options)


Jx = JxOp('J')
Jy = JyOp('J')
Jz = JzOp('J')
J2 = J2Op('J')
Jplus = JplusOp('J')
Jminus = JminusOp('J')


#-----------------------------------------------------------------------------
# Spin States
#-----------------------------------------------------------------------------


class SpinState(State):
    """Base class for angular momentum states."""

    _label_separator = ','

    @property
    def j(self):
        return self.label[0]

    @property
    def m(self):
        return self.label[1]

    @classmethod
    def _eval_hilbert_space(cls, label):
        return ComplexSpace(2*label[0]+1)


class JzKet(SpinState, Ket):
    """Eigenket of Jz."""

    def _eval_innerproduct_JzBra(self, bra, **hints):
        d1 = KroneckerDelta(self.j, bra.j)
        d2 = KroneckerDelta(self.m, bra.m)
        return d1*d2

    @property
    def dual_class(self):
        return JzBra

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        if self.j == Rational(1,2):
            if self.m == Rational(1,2):
                return Matrix([1,0])
            elif self.m == -Rational(1,2):
                return Matrix([0,1])


class JzBra(SpinState, Bra):
    """Eigenbra of Jz."""

    @property
    def dual_class(self):
        return JzKet


class JxKet(SpinState, Ket):
    """Eigenket of Jx."""

    @property
    def dual_class(self):
        return JxBra

    def _eval_innerproduct_JxBra(self, bra, **hints):
        d1 = KroneckerDelta(self.j, bra.j)
        d2 = KroneckerDelta(self.m, bra.m)
        return d1*d2

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        if self.j == Rational(1,2):
            if self.m == Rational(1,2):
                return Matrix([1,1])/sqrt(2)
            elif self.m == -Rational(1,2):
                return Matrix([1,-1])/sqrt(2)


class JxBra(SpinState, Bra):
    """Eigenbra of Jx."""

    @property
    def dual_class(self):
        return JxKet


class JyKet(SpinState, Ket):
    """Eigenket of Jy."""

    @property
    def dual_class(self):
        return JyBra

    def _eval_innerproduct_JyBar(self, bra, **hints):
        d1 = KroneckerDelta(self.s, bra.s)
        d2 = KroneckerDelta(self.ms, bra.ms)
        return d1*d2

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JzOp(self, basis, **options):
        if self.j == Rational(1,2):
            if self.m == Rational(1,2):
                return Matrix([1,I])/sqrt(2)
            elif self.m == -Rational(1,2):
                return Matrix([1,-I])/sqrt(2)


class JyBra(SpinState, Bra):
    """Eigenbra of Jy."""

    @property
    def dual_class(self):
        return JyKet
