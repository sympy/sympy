from sympy.concrete.summations import Sum
from sympy.core import oo
from sympy.core.mul import Mul
from sympy.core.numbers import I, Rational
from sympy.core.singleton import S
from sympy.simplify import simplify
from sympy.core.symbol import symbols, Symbol
from sympy.core.sympify import sympify
from sympy.integrals import Integral
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.special.delta_functions import DiracDelta
from sympy.functions.special.tensor_functions import KroneckerDelta


from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.cartesian import X, XKet, XBra, PxOp
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.gate import H, XGate, IdentityGate
from sympy.physics.quantum.operator import Operator, IdentityOperator, OuterProduct
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.spin import Jx, Jy, Jz, Jplus, Jminus, J2, JzKet
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.state import Ket, Bra
from sympy.physics.quantum.density import Density
from sympy.physics.quantum.qubit import Qubit, QubitBra
from sympy.physics.quantum.boson import BosonOp, BosonFockKet, BosonFockBra
from sympy.physics.quantum.fermion import FermionOp, FermionFockKet, FermionFockBra
from sympy.testing.pytest import warns_deprecated_sympy


# -----------------------------------------------------------------------------
# Variables and utilities
# -----------------------------------------------------------------------------


j, jp, m, mp = symbols("j j' m m'")
p, q, r, st = symbols("p q r s", integer=True, nonnegative=True)
x, y, z = symbols("x y z", real=True)
omega = Symbol("omega", real=True)
alpha = Symbol("alpha", complex=True)

z = JzKet(1, 0)
po = JzKet(1, 1)
mo = JzKet(1, -1)

a = BosonOp("a")
f = FermionOp("f")

A = Operator("A")
B = Operator("B")

class SimpleBra(Bra):
    pass


class SimpleKet(Ket):
    def _eval_innerproduct_SimpleBra(self, bra, **hints):
        return alpha * exp(I * omega)


TP = TensorProduct

A = Operator("A")


# -----------------------------------------------------------------------------
# Main qapply Unary Handler Tests
# -----------------------------------------------------------------------------


def test_unary_dagger():
    """Test Dagger expression handling"""
    assert qapply(Dagger(a * BosonFockKet(p)), dagger=True) == sqrt(p) * BosonFockBra(
        p - 1
    )


def test_unary_density():
    """Test Density matrix handling"""
    d = Density([a * BosonFockKet(p), 0.5], [a * BosonFockKet(q), 0.5])
    assert qapply(d) == Density(
        [sqrt(p) * BosonFockKet(p - 1), 0.5], [sqrt(q) * BosonFockKet(q - 1), 0.5]
    )
    d = Density([Jz * mo, 0.5], [Jz * po, 0.5])
    assert qapply(d) == Density([-hbar * mo, 0.5], [hbar * po, 0.5])


def test_unary_tensorproduct():
    """Test TensorProduct expression handling"""
    assert qapply(TP(a * BosonFockKet(p), a * BosonFockKet(q))) == TP(
        sqrt(p) * BosonFockKet(p - 1), sqrt(q) * BosonFockKet(q - 1)
    )


def test_unary_add():
    """Test distribution over addition"""
    assert qapply(a * BosonFockKet(p) + Dagger(a) * BosonFockKet(q)) == sqrt(
        p
    ) * BosonFockKet(p - 1) + sqrt(q + 1) * BosonFockKet(q + 1)


def test_unary_pow():
    """Test power expression handling"""
    assert qapply((SimpleBra() * SimpleKet()) ** 2) == alpha**2 * exp(2 * I * omega)
    assert qapply(sqrt(SimpleBra() * SimpleKet())) == sqrt(alpha * exp(I * omega))
    assert qapply((po.dual * Jz * po) ** 2) == hbar**2
    assert qapply((BosonFockBra(1) * a * BosonFockKet(1)) ** 2) == S.Zero


def test_unary_sum():
    """Test Sum expression handling"""
    assert qapply(Sum(a * BosonFockKet(p), (p, 0, 2))) == Sum(
        sqrt(p) * BosonFockKet(p - 1), (p, 0, 2)
    )
    assert qapply(Sum(a * BosonFockKet(p), (p, 0, 2)), sum_doit=True) == sqrt(
        1
    ) * BosonFockKet(0) + sqrt(2) * BosonFockKet(1)
    assert qapply(Jz * po + Jz * mo) == hbar * po - hbar * mo


def test_unary_integral():
    """Test Integral expression handling"""
    assert qapply(Integral(XBra(x) * XKet(y), (x, -oo, oo))).doit() == S.One
    assert qapply(Integral(XBra(x) * XKet(y), (x, -oo, oo))) == Integral(
        DiracDelta(y-x), (x, -oo, oo)
    )


def test_unary_abs():
    """Test absolute value expressions"""
    assert qapply(Abs(alpha)) == Abs(alpha)


def test_unary_number():
    """Test numeric type conversion"""
    assert qapply(0) == S.Zero
    assert qapply(1.0) == sympify(1.0)


# -----------------------------------------------------------------------------
# qapply_Mul: SlidingTransform unary handler tests
# -----------------------------------------------------------------------------


def test_mul_unary_outerproduct():
    """Test OuterProduct decomposition"""
    assert qapply(
        OuterProduct(BosonFockKet(p), BosonFockBra(q)) * a * BosonFockKet(q + 1)
    ) == sqrt(q + 1) * BosonFockKet(p)


def test_mul_unary_commutator():
    """Test commutator evaluation in products"""
    assert qapply(Commutator(Jx, Jy) * Jz * po) == I * hbar**3 * po
    assert qapply(Commutator(J2, Jz) * Jz * po) == 0
    assert qapply(Commutator(a, Dagger(a)) * BosonFockKet(1)) == BosonFockKet(1)
    assert qapply(Commutator(X, PxOp()) * XKet(x)) == I * hbar * XKet(x)


def test_mul_unary_anticommutator():
    """Test anticommutator evaluation in products"""
    assert qapply(AntiCommutator(f, Dagger(f)) * FermionFockKet(0)) == FermionFockKet(0)
    assert qapply(AntiCommutator(Jx, Jy) * Jz * po).expand() == hbar**3*I*mo

def test_mul_unary_tensorproduct():
    """Test tensor product in products."""
    assert qapply(A*TensorProduct(Jz*po, Jz*mo)*B) == -hbar**2*A*TensorProduct(po, mo)*B


def test_mul_unary_density():
    """Test density matrix in products."""
    d = Density([a * BosonFockKet(p), 0.5], [a * BosonFockKet(q), 0.5])
    assert qapply(A*d*B) == A*Density(
        [sqrt(p) * BosonFockKet(p - 1), 0.5], [sqrt(q) * BosonFockKet(q - 1), 0.5]
    )*B
    d = Density([Jz * mo, 0.5], [Jz * po, 0.5])
    assert qapply(A*d*B) == A*Density([-hbar * mo, 0.5], [hbar * po, 0.5])*B


def test_mul_unary_pow():
    """Test power expansion in multiplications"""
    assert qapply(a**2 * BosonFockKet(p)) == sqrt(p) * sqrt(p - 1) * BosonFockKet(p - 2)
    assert qapply(BosonFockBra(q) * a**2 * BosonFockKet(p)) == KroneckerDelta(
        q, p - 2
    ) * sqrt(p) * sqrt(p - 1)
    # Here we are testing the path for a Mul that has a Pow with a non-integer
    # power
    assert qapply(
        alpha
        * sqrt(
            BosonFockBra(p)
            * a
            * BosonFockKet(p + 1)
            * BosonFockBra(p + 1)
            * Dagger(a)
            * BosonFockKet(p)
        )
    ) == alpha * sqrt(p + 1)


# -----------------------------------------------------------------------------
# qapply_Mul: SlidingTransform binary handler tests
# -----------------------------------------------------------------------------

def test_mul_binary_op_ket():
    """Test operator application to kets"""
    assert qapply(a * BosonFockKet(q)) == sqrt(q) * BosonFockKet(q - 1)
    assert qapply(A * BosonFockKet(q)) == A * BosonFockKet(q)
    assert qapply(Jplus * Jplus * mo) == 2 * hbar**2 * po


def test_mul_binary_op_tensorproduct():
    pass


def test_mul_binary_op_wavefunction():
    pass


def test_mul_binary_bra_op():
    """Test bra-operator interactions"""
    assert qapply(BosonFockBra(q) * a, dagger=True) == sqrt(q + 1) * BosonFockBra(q + 1)


def test_mul_binary_sum_sum():
    pass


def test_mul_binary_integral_integral():
    """Test integral-expression interactions"""
    e = qapply(XBra(x) * Integral(X * XKet(y), (y, -oo, oo)))
    assert e == Integral(y * DiracDelta(y - x), (y, -oo, oo))
    assert e.doit() == x

    e = qapply(Integral(XBra(y) * X, (y, -oo, oo)) * XKet(x))
    assert e == Integral(x * DiracDelta(x - y), (y, -oo, oo))
    assert e.doit() == x

    e = qapply(Integral(XBra(x), (x, -oo, oo)) * Integral(XKet(y), (y, -oo, oo)))
    assert e == Integral(DiracDelta(y - x), (x, -oo, oo), (y, -oo, oo))


def test_mul_binary_add_add():
    pass


def test_mul_binary_sum_add():
    pass


def test_mul_binary_add_sum():
    pass


def test_mul_binary_integral_add():
    pass


def test_mul_binary_add_integral():
    pass


def test_mul_binary_sum_integral():
    pass


def test_mul_binary_integral_sum():
    pass


def test_mul_binary_add_expr():
    pass


def test_mul_binary_expr_add():
    pass


def test_mul_binary_sum_expr():
    pass


def test_mul_binary_expr_sum():
    pass


def test_mul_binary_integral_expr():
    pass


def test_mul_binary_expr_integral():
    pass


def test_mul_binary_add():
    """Test distribution over addition"""
    assert (
        qapply((Jplus + Jminus) * z / sqrt(2)).expand() == (hbar * (po + mo)).expand()
    )
    assert qapply(Jz * (po + mo)).expand() == (hbar * (po - mo)).expand()
    assert qapply((a + Dagger(a)) * BosonFockKet(1)).expand() == (BosonFockKet(0) + sqrt(2) * BosonFockKet(2)).expand()


def test_mul_binary_sum():
    """Test sum-expression interactions"""
    e = qapply(Jz * Sum(JzKet(j, m), (m, -1, 1)))
    assert e == Sum(hbar * m * JzKet(j, m), (m, -1, 1))
    assert e.doit() == hbar * (JzKet(j, 1) - JzKet(j, -1))
    assert qapply(
        Sum(BosonFockBra(p), (p, 0, 3)) * Sum(BosonFockKet(q), (q, 0, 3))
    ) == Sum(KroneckerDelta(p, q), (p, 0, 3), (q, 0, 3))


# -----------------------------------------------------------------------------
# Integration Tests
# -----------------------------------------------------------------------------


def test_basic():
    """Test basic operator applications"""
    assert qapply(Jz * po) == hbar * po
    assert qapply(Jx * z).expand() == (sqrt(2) * hbar * (po + mo) / sympify(2)).expand()
    assert qapply(Jminus * Jminus * po) == 2 * hbar**2 * mo
    assert qapply(Jplus**2 * mo) == 2 * hbar**2 * po
    assert qapply(Jplus**2 * Jminus**2 * po) == 4 * hbar**4 * po


def test_extra():
    """Test complex expressions with non-quantum factors"""
    extra = z.dual * A * z
    assert qapply(Jz * po * extra) == hbar * po * extra
    assert (
        qapply(Jx * z * extra).expand()
        == ((hbar * po / sqrt(2) + hbar * mo / sqrt(2)) * extra).expand()
    )
    assert (
        qapply((Jplus + Jminus) * z / sqrt(2) * extra).expand()
        == (hbar * (po + mo) * extra).expand()
    )
    assert (
        qapply(Jz * (po + mo) * extra).expand() == (hbar * (po - mo) * extra).expand()
    )
    assert (
        qapply(Jz * po * extra + Jz * mo * extra).expand()
        == (hbar * (po - mo) * extra).expand()
    )
    assert qapply(Jminus * Jminus * po * extra) == 2 * hbar**2 * mo * extra
    assert qapply(Jplus**2 * mo * extra) == 2 * hbar**2 * po * extra
    assert qapply(Jplus**2 * Jminus**2 * po * extra) == 4 * hbar**4 * po * extra


def test_innerproduct():
    """Test inner product handling"""
    assert qapply(po.dual * Jz * po, ip_doit=False) == hbar * (po.dual * po)
    assert qapply(po.dual * Jz * po) == hbar


def test_outerproduct():
    """Test outer product handling"""
    e = Jz * (mo * po.dual) * Jz * po
    assert qapply(e) == -(hbar**2) * mo
    assert qapply(e, ip_doit=False) == -(hbar**2) * (po.dual * po) * mo
    assert qapply(e).doit() == -(hbar**2) * mo


def test_tensorproduct():
    """Test tensor product applications"""
    a = BosonOp("a")
    b = BosonOp("b")
    ket1 = TensorProduct(BosonFockKet(1), BosonFockKet(2))
    ket2 = TensorProduct(BosonFockKet(0), BosonFockKet(0))
    ket3 = TensorProduct(BosonFockKet(0), BosonFockKet(2))
    bra1 = TensorProduct(BosonFockBra(0), BosonFockBra(0))
    bra2 = TensorProduct(BosonFockBra(1), BosonFockBra(2))
    assert qapply(TensorProduct(a, b**2) * ket1) == sqrt(2) * ket2
    assert qapply(TensorProduct(a, Dagger(b) * b) * ket1) == 2 * ket3
    assert qapply(bra1 * TensorProduct(a, b * b), dagger=True) == sqrt(2) * bra2
    assert qapply(bra2 * ket1).doit() == S.One
    assert qapply(TensorProduct(a, b * b) * ket1) == sqrt(2) * ket2
    assert qapply(Dagger(TensorProduct(a, b * b) * ket1), dagger=True) == sqrt(
        2
    ) * Dagger(ket2)


def test_dagger():
    """Test dagger option handling"""
    lhs = Dagger(Qubit(0)) * Dagger(H(0))
    rhs = Dagger(Qubit(1)) / sqrt(2) + Dagger(Qubit(0)) / sqrt(2)
    assert qapply(lhs, dagger=True) == rhs


# -----------------------------------------------------------------------------
# Historical Issues
# -----------------------------------------------------------------------------


def test_issue_6073():
    """Test handling of non-commutative symbols in states"""
    x, y = symbols("x y", commutative=False)
    A = Ket(x, y)
    B = Operator("B")
    assert qapply(A) == A
    assert qapply(A.dual * B) == A.dual * B


def test_issue_3044():
    """Test tensor product of spin operators"""
    expr1 = TensorProduct(
        Jz * JzKet(S(2), S.NegativeOne) / sqrt(2), Jz * JzKet(S.Half, S.Half)
    )
    result = Mul(S.NegativeOne, Rational(1, 4), 2**S.Half, hbar**2)
    result *= TensorProduct(JzKet(2, -1), JzKet(S.Half, S.Half))
    assert qapply(expr1) == result


def test_issue_24158_ket_times_op():
    """Test that qapply doesn't incorrectly evaluate ket*op as op*ket"""
    P = BosonFockKet(0) * BosonOp("a")  # undefined term
    # Does lhs._apply_operator_BosonOp(rhs) still evaluate ket*op as op*ket?
    assert qapply(P) == P  # qapply(P) -> BosonOp("a")*BosonFockKet(0) = 0 before fix
    P = Qubit(1) * XGate(0)  # undefined term
    # Does rhs._apply_operator_Qubit(lhs) still evaluate ket*op as op*ket?
    assert qapply(P) == P  # qapply(P) -> Qubit(0) before fix
    P1 = Mul(
        QubitBra(0), Mul(QubitBra(0), Qubit(0)), XGate(0)
    )  # legal expr <0| * (<1|*|1>) * X
    assert qapply(P1) == QubitBra(0) * XGate(0)  # qapply(P1) -> 0 before fix
    P1 = qapply(
        P1, dagger=True
    )  # unsatisfactorily -> <0|*X(0), expect <1| since dagger=True
    assert qapply(P1, dagger=True) == QubitBra(
        1
    )  # qapply(P1, dagger=True) -> 0 before fix
    P2 = QubitBra(0) * (QubitBra(0) * Qubit(0)) * XGate(0)  # 'forgot' to set brackets
    P2 = qapply(
        P2, dagger=True
    )  # unsatisfactorily -> <0|*X(0), expect <1| since dagger=True
    assert P2 == QubitBra(1)  # qapply(P1) -> 0 before fix
    # Pull Request 24237: IdentityOperator from the right without dagger=True option
    with warns_deprecated_sympy():
        assert qapply(QubitBra(1) * IdentityOperator()) == QubitBra(1)
        assert qapply(IdentityGate(0) * (Qubit(0) + Qubit(1))) == Qubit(0) + Qubit(1)
