"""Quantum mechanical angular momemtum."""

from sympy import (
        AtomicExpr, binomial, cos, diff, exp, factorial, I, Integer, Matrix, N, pi,
        Rational, S, sin, sqrt, Symbol, sympify, simplify
)
from sympy.matrices.matrices import zeros
from sympy.printing.pretty.stringpict import prettyForm, stringPict

from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.operator import (
    HermitianOperator, Operator, UnitaryOperator
)
from sympy.physics.quantum.state import Bra, Ket, State
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
    'Rotation',
    'WignerD'
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
        return hbar**2*j*(j+1)*ket

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
    """Wigner D operator in terms of Euler angles.

    Defines the rotation operator in terms of the Euler angles defined by
    the z-y-z convention for a passive transformation. That is the coordinate
    axes are rotated first about the z-axis, giving the new x'-y'-z' axes. Then
    this new coordinate system is rotated about the new y'-axis, giving new
    x''-y''-z'' axes. Then this new coordinate system is rotated about the
    z''-axis. Conventions follow those laid out in [1].

    See the Wigner D-function, Rotation.D, and the Wigner small-d matrix for
    the evaluation of the rotation operator on spin states.

    Parameters
    ==========

    alpha : Number, Symbol
        First Euler Angle
    beta : Number, Symbol
        Second Euler angle
    gamma : Number, Symbol
        Third Euler angle

    Examples
    ========

    A simple example rotation operator:

        >>> from sympy import pi
        >>> from sympy.physics.quantum.spin import Rotation
        >>> Rotation(pi, 0, pi/2)
        'R'(pi,0,pi/2)

    With symbolic Euler angles and calculating the inverse rotation operator:

        >>> from sympy import symbols
        >>> a, b, c = symbols('a b c')
        >>> Rotation(a, b, c)
        'R'(a,b,c)
        >>> Rotation(a, b, c).inverse()
        'R'(-c,-b,-a)


    References
    ==========

    [1] Varshalovich, D A, Quantum Theory of Angular Momentum. 1988.
    """

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
        return Rotation(-self.gamma, -self.beta, -self.alpha)

    @classmethod
    def D(cls, j, m, mp, alpha, beta, gamma):
        """Wigner D-function.

        Returns an instance of the WignerD class. See the corresponding
        docstring for more information on the Wigner-D matrix.

        Parameters
        ===========

        j : Number
            Total angular momentum
        m : Number
            Eigenvalue of angular momentum along axis after rotation
        mp : Number
            Eigenvalue of angular momentum along rotated axis
        alpha : Number, Symbol
            First Euler angle of rotation
        beta : Number, Symbol
            Second Euler angle of rotation
        gamma : Number, Symbol
            Third Euler angle of rotation

        Examples
        ========

        Return the Wigner-D matrix element for a defined rotation, both
        numerical and symbolic:

            >>> from sympy.physics.quantum.spin import Rotation
            >>> from sympy import pi, symbols
            >>> alpha, beta, gamma = symbols('alpha beta gamma')
            >>> Rotation.D(1, 1, 0,pi, pi/2,-pi)
            WignerD(1, 1, 0, pi, pi/2, -pi)

        """
        return WignerD(j,m,mp,alpha,beta,gamma)

    @classmethod
    def d(cls, j, m, mp, beta):
        """Wigner small-d function.

        Roturns an instance of cthe WignerD class with the alpha and gamma
        angles given as 0. See the corresponding docstring for more
        information on the Wigner small-d matrix.

        Parameters
        ===========

        j : Number
            Total angular momentum
        m : Number
            Eigenvalue of angular momentum along axis after rotation
        mp : Number
            Eigenvalue of angular momentum along rotated axis
        beta : Number, Symbol
            Second Euler angle of rotation

        Examples
        ========

        Return the Wigner-D matrix element for a defined rotation, both
        numerical and symbolic:

            >>> from sympy.physics.quantum.spin import Rotation
            >>> from sympy import pi, symbols
            >>> beta = symbols('beta')
            >>> Rotation.d(1, 1, 0, pi/2)
            WignerD(1, 1, 0, 0, pi/2, 0)

        """
        return WignerD(j,m,mp,0,beta,0)

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


class WignerD(AtomicExpr):
    """Wigner-D function

    The Wigner D-function gives the matrix elements of the rotation
    operator in the jm-representation. For the Euler angles alpha, beta,
    gamma, the the D-function is defined such that:
    <j,m| R(alpha,beta,gamma) |j',m'> = delta_jj' * D(j, m, m', alpha, beta, gamma)
    Where the rotation operator is as defined by the Rotation class.

    The Wigner D-function defined in this way gives:
    D(j, m, m', alpha, beta, gamma) = exp(-i*m*alpha) * d(j, m, m', beta) * exp(-i*m'*gamma)
    Where d is the Wigner small-d function, which is given by Rotation.d.

    The Wigner small-d function gives the component of the Wigner
    D-function that is determined by the second Euler angle. That is the
    Wigner D-function is:
    D(j, m, m', alpha, beta, gamma) = exp(-i*m*alpha) * d(j, m, m', beta) * exp(-i*m'*gamma)
    Where d is the small-d function. The Wigner D-function is given by
    Rotation.D.

    Note that to evaluate the D-function, the j, m and mp parameters must
    be integer or half integer numbers.

    Parameters
    ==========

    j : Number
        Total angular momentum
    m : Number
        Eigenvalue of angular momentum along axis after rotation
    mp : Number
        Eigenvalue of angular momentum along rotated axis
    alpha : Number, Symbol
        First Euler angle of rotation
    beta : Number, Symbol
        Second Euler angle of rotation
    gamma : Number, Symbol
        Third Euler angle of rotation

    Examples
    ========

    Evaluate the Wigner-D matrix elements of a simple rotation:

        >>> from sympy.physics.quantum.spin import Rotation
        >>> from sympy import pi
        >>> rot = Rotation.D(1, 1, 0, pi, pi/2, 0)
        >>> rot
        WignerD(1, 1, 0, pi, pi/2, 0)
        >>> rot.doit()
        2**(1/2)/2

    Evaluate the Wigner-d matrix elements of a simple rotation

        >>> rot = Rotation.d(1, 1, 0, pi/2)
        >>> rot
        WignerD(1, 1, 0, 0, pi/2, 0)
        >>> rot.doit()
        -2**(1/2)/2

    References
    ==========

    [1] Varshalovich, D A, Quantum Theory of Angular Momentum. 1988.
    """

    def __new__(cls, *args, **hints):
        if not len(args) == 6:
            raise ValueError('6 parameters expected, got %s' % args)
        evaluate = hints.get('evaluate', False)
        if evaluate:
            return AtomicExpr.__new__(cls, *args)._eval_wignerd()
        return AtomicExpr.__new__(cls, *args, **{'evaluate': False})

    @property
    def j(self):
        return self.args[0]

    @property
    def m(self):
        return self.args[1]

    @property
    def mp(self):
        return self.args[2]

    @property
    def alpha(self):
        return self.args[3]

    @property
    def beta(self):
        return self.args[4]

    @property
    def gamma(self):
        return self.args[5]

    def _latex(self, printer, *args):
        if self.alpha == 0 and self.gamma == 0:
            return r'd^{%s}_{%s,%s}\left(%s\right)' % \
                ( printer._print(self.j), printer._print(self.m), printer._print(self.mp),
                printer._print(self.beta) )
        return r'D^{%s}_{%s,%s}\left(%s,%s,%s\right)' % \
            ( printer._print(self.j), printer._print(self.m), printer._print(self.mp),
            printer._print(self.alpha), printer._print(self.beta), printer._print(self.gamma) )

    def _pretty(self, printer, *args):
        top = printer._print(self.j)

        bot = printer._print(self.m)
        bot = prettyForm(*bot.right(','))
        bot = prettyForm(*bot.right(printer._print(self.mp)))

        pad = max(top.width(), bot.width())
        top = prettyForm(*top.left(' '))
        bot = prettyForm(*bot.left(' '))
        if pad > top.width():
            top = prettyForm(*top.right(' ' * (pad-top.width())))
        if pad > bot.width():
            bot = prettyForm(*bot.right(' ' * (pad-bot.width())))

        if self.alpha == 0 and self.gamma == 0:
            args = printer._print(self.beta)

            s = stringPict('d' + ' '*pad)
        else:
            args = printer._print(self.alpha)
            args = prettyForm(*args.right(','))
            args = prettyForm(*args.right(printer._print(self.beta)))
            args = prettyForm(*args.right(','))
            args = prettyForm(*args.right(printer._print(self.gamma)))

            s = stringPict('D' + ' '*pad)

        args = prettyForm(*args.parens())
        s = prettyForm(*s.above(top))
        s = prettyForm(*s.below(bot))
        s = prettyForm(*s.right(args))
        return s

    def doit(self, **hints):
        hints['evaluate'] = True
        return WignerD(*self.args, **hints)

    def _eval_wignerd(self):
        j = sympify(self.j)
        m = sympify(self.m)
        mp = sympify(self.mp)
        alpha = sympify(self.alpha)
        beta = sympify(self.beta)
        gamma = sympify(self.gamma)
        if not j.is_number:
            raise ValueError("j parameter must be numerical to evaluate, got %s", j)
        r = 0
        if beta == pi/2:
            # Varshalovich Equation (5), Section 4.16, page 113, setting
            # alpha=gamma=0.
            for k in range(2*j+1):
                if k > j+mp or k > j-m or k < mp-m:
                    continue
                r += (-S(1))**k * binomial(j+mp, k) * binomial(j-mp, k+m-mp)
            r *= (-S(1))**(m-mp) / 2**j * sqrt(factorial(j+m) * \
                    factorial(j-m) / (factorial(j+mp) * factorial(j-mp)))
        else:
            # Varshalovich Equation(5), Section 4.7.2, page 87, where we set
            # beta1=beta2=pi/2, and we get alpha=gamma=pi/2 and beta=phi+pi,
            # then we use the Eq. (1), Section 4.4. page 79, to simplify:
            # d(j, m, mp, beta+pi) = (-1)**(j-mp) * d(j, m, -mp, beta)
            # This happens to be almost the same as in Eq.(10), Section 4.16,
            # except that we need to substitute -mp for mp.
            size, mvals = m_values(j)
            for mpp in mvals:
                r += Rotation.d(j, m, mpp, pi/2).doit() * (cos(-mpp*beta)+I*sin(-mpp*beta)) * \
                    Rotation.d(j, mpp, -mp, pi/2).doit()
            # Empirical normalization factor so results match Varshalovich
            # Tables 4.3-4.12
            # Note that this exact normalization does not follow from the
            # above equations
            r = r * I**(2*j-m-mp) * (-1)**(2*m)
            # Finally, simplify the whole expression
            r = simplify(r)
        r *= exp(-I*m*alpha)*exp(-I*mp*gamma)
        return r


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

    def _represent_base(self, **options):
        j = self.j
        alpha = options.get('alpha', 0)
        beta = options.get('beta', 0)
        gamma = options.get('gamma', 0)
        size, mvals = m_values(j)
        result = zeros((size,1))
        for p in range(size):
            result[p,0] = Rotation.D(self.j, mvals[p], self.m, alpha, beta, gamma).doit()
        return result

    def _eval_rewrite_as_Jx(self, *args):
        if isinstance(self,Bra):
            return self._rewrite_basis(Jx, JxBra, *args)
        return self._rewrite_basis(Jx, JxKet, *args)

    def _eval_rewrite_as_Jy(self, *args):
        if isinstance(self,Bra):
            return self._rewrite_basis(Jy, JyBra, *args)
        return self._rewrite_basis(Jy, JyKet, *args)

    def _eval_rewrite_as_Jz(self, *args):
        return self._rewrite_basis(Jz, JzKet, *args)

    def _rewrite_basis(self, basis, evect_cls, *args):
        j = self.j
        size, mvals = m_values(j)
        result = 0
        vect = represent(self, basis=basis)
        for p in range(size):
            me = vect[p]
            me *= evect_cls.__new__(evect_cls, j, mvals[p])
            result += me
        return result

    def _eval_innerproduct_JxBra(self, bra, **hints):
        result = KroneckerDelta(self.j, bra.j)
        if not isinstance(bra, self.__class__):
            result *= self._represent_JxOp(None)[bra.j-bra.m]
        else:
            result *= KroneckerDelta(self.m, bra.m)
        return result

    def _eval_innerproduct_JyBra(self, bra, **hints):
        result = KroneckerDelta(self.j, bra.j)
        if not isinstance(bra, self.__class__):
            result *= self._represent_JyOp(None)[bra.j-bra.m]
        else:
            result *= KroneckerDelta(self.m, bra.m)
        return result

    def _eval_innerproduct_JzBra(self, bra, **hints):
        result = KroneckerDelta(self.j, bra.j)
        if not issubclass(bra.dual_class, self.__class__):
            result *= self._represent_JzOp(None)[bra.j-bra.m]
        else:
            result *= KroneckerDelta(self.m, bra.m)
        return result


class JxKet(SpinState, Ket):
    """Eigenket of Jx.

    See JzKet for the usage of spin eigenstates.
    """

    @property
    def dual_class(self):
        return JxBra

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JxOp(self, basis, **options):
        return self._represent_base()

    def _represent_JyOp(self, basis, **options):
        return self._represent_base(alpha=3*pi/2)

    def _represent_JzOp(self, basis, **options):
        return self._represent_base(beta=pi/2)


class JxBra(SpinState, Bra):
    """Eigenbra of Jx.

    See JzKet for the usage of spin eigenstates.
    """

    @property
    def dual_class(self):
        return JxKet


class JyKet(SpinState, Ket):
    """Eigenket of Jy.

    See JzKet for the usage of spin eigenstates.
    """

    @property
    def dual_class(self):
        return JyBra

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JxOp(self, basis, **options):
        return self._represent_base(gamma=pi/2)

    def _represent_JyOp(self, basis, **options):
        return self._represent_base()

    def _represent_JzOp(self, basis, **options):
        return self._represent_base(alpha=3*pi/2,beta=-pi/2,gamma=pi/2)


class JyBra(SpinState, Bra):
    """Eigenbra of Jy.

    See JzKet for the usage of spin eigenstates.
    """

    @property
    def dual_class(self):
        return JyKet


class JzKet(SpinState, Ket):
    """Eigenket of Jz.

    Spin state which is an eigenstate of the Jz operator

    Parameters
    ==========

    j : Number, Symbol
        Total spin angular momentum
    m : Number, Symbol
        Eigenvalue of the Jz spin operator

    Examples
    ========

    Defining simple spin states, both numerical and symbolic:

        >>> from sympy.physics.quantum.spin import JzKet
        >>> from sympy import symbols
        >>> JzKet(1, 0)
        |1,0>
        >>> j, m = symbols('j m')
        >>> JzKet(j, m)
        |j,m>

    Rewriting the JzKet in terms of eigenkets of the Jx operator:
    Note: that the resulting eigenstates are JxKet's

        >>> JzKet(1,1).rewrite("Jx")
        |1,-1>/2 - 2**(1/2)*|1,0>/2 + |1,1>/2

    Get the vector representation of a state in terms of the basis elements
    of the Jx operator:

        >>> from sympy.physics.quantum.represent import represent
        >>> from sympy.physics.quantum.spin import Jx
        >>> represent(JzKet(1,-1), basis=Jx)
        [       1/2]
        [2**(1/2)/2]
        [       1/2]

    Apply innerproducts between states:

        >>> from sympy.physics.quantum.innerproduct import InnerProduct
        >>> from sympy.physics.quantum.spin import JxBra
        >>> i = InnerProduct(JxBra(1,1), JzKet(1,1))
        >>> i
        <1,1|1,1>
        >>> i.doit()
        1/2

    """

    @property
    def dual_class(self):
        return JzBra

    def _represent_default_basis(self, **options):
        return self._represent_JzOp(None, **options)

    def _represent_JxOp(self, basis, **options):
        return self._represent_base(beta=3*pi/2)

    def _represent_JyOp(self, basis, **options):
        return self._represent_base(alpha=3*pi/2,beta=pi/2,gamma=pi/2)

    def _represent_JzOp(self, basis, **options):
        return self._represent_base()


class JzBra(SpinState, Bra):
    """Eigenbra of Jz.

    See the JzKet for the usage of spin eigenstates.
    """

    @property
    def dual_class(self):
        return JzKet
