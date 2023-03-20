from sympy.core.mul import Mul
from sympy.core.numbers import (I, E, Integer, Rational)
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.core import Pow
from sympy.core.expr import unchanged
from sympy.functions.elementary.trigonometric import sin, cos, pi
from sympy.matrices import ImmutableMatrix

from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.gate import H, XGate, IdentityGate, TGate, UGate
from sympy.physics.quantum.operator import (Operator, IdentityOperator,
                                            UnitaryOperator, OuterProduct)
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.spin import Jx, Jy, Jz, Jplus, Jminus, J2, JzKet
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.state import Ket, Bra, OrthogonalKet
from sympy.physics.quantum.density import Density
from sympy.physics.quantum.qubit import Qubit, QubitBra
from sympy.physics.quantum.boson import BosonOp, BosonFockKet, BosonFockBra


j, jp, m, mp = symbols("j j' m m'")

z = JzKet(1, 0)
po = JzKet(1, 1)
mo = JzKet(1, -1)

A = Operator('A')


class Foo(Operator):
    def _apply_operator_JzKet(self, ket, **options):
        return ket


def test_basic():
    assert qapply(Jz*po) == hbar*po
    assert qapply(Jx*z) == hbar*po/sqrt(2) + hbar*mo/sqrt(2)
    assert qapply((Jplus + Jminus)*z/sqrt(2)) == hbar*po + hbar*mo
    assert qapply(Jz*(po + mo)) == hbar*po - hbar*mo
    assert qapply(Jz*po + Jz*mo) == hbar*po - hbar*mo
    assert qapply(Jminus*Jminus*po) == 2*hbar**2*mo
    assert qapply(Jplus**2*mo) == 2*hbar**2*po
    assert qapply(Jplus**2*Jminus**2*po) == 4*hbar**4*po


def test_extra():
    extra = z.dual*A*z
    assert qapply(Jz*po*extra, op_join=False) == hbar*po*extra # op_join
    assert qapply(Jx*z*extra, op_join=False) == \
        hbar*po/sqrt(2)*extra + hbar*mo/sqrt(2)*extra #2 op_join
    assert qapply(
        (Jplus + Jminus)*z/sqrt(2)*extra, op_join=False) == hbar*po*extra + hbar*mo*extra # op_join
    assert qapply(Jz*(po + mo)*extra, op_join=False) == hbar*po*extra - hbar*mo*extra # op_join
    assert qapply(Jz*po*extra + Jz*mo*extra, op_join=False) == hbar*po*extra - hbar*mo*extra # op_join
    assert qapply(Jminus*Jminus*po*extra, op_join=False) == 2*hbar**2*mo*extra # op_join
    assert qapply(Jplus**2*mo*extra, op_join=False) == 2*hbar**2*po*extra # op_join
    assert qapply(Jplus**2*Jminus**2*po*extra, op_join=False) == 4*hbar**4*po*extra # op_join


def test_innerproduct():
    assert qapply(po.dual*Jz*po, ip_doit=False) == hbar*(po.dual*po)
    assert qapply(po.dual*Jz*po) == hbar


def test_zero():
    assert qapply(0) == 0
    assert qapply(Integer(0)) == 0


def test_commutator():
    assert qapply(Commutator(Jx, Jy)*Jz*po) == I*hbar**3*po
    assert qapply(Commutator(J2, Jz)*Jz*po) == 0
    assert qapply(Commutator(Jz, Foo('F'))*po) == 0
    assert qapply(Commutator(Foo('F'), Jz)*po) == 0


def test_anticommutator():
    assert qapply(AntiCommutator(Jz, Foo('F'))*po) == 2*hbar*po
    assert qapply(AntiCommutator(Foo('F'), Jz)*po) == 2*hbar*po


def test_outerproduct():
    e = Jz*(mo*po.dual)*Jz*po
    assert qapply(e) == -hbar**2*mo
    assert qapply(e, ip_doit=False) == -hbar**2*(po.dual*po)*mo
    assert qapply(e).doit() == -hbar**2*mo


def test_tensorproduct():
    a = BosonOp("a")
    b = BosonOp("b")
    ket1 = TensorProduct(BosonFockKet(1), BosonFockKet(2))
    ket2 = TensorProduct(BosonFockKet(0), BosonFockKet(0))
    ket3 = TensorProduct(BosonFockKet(0), BosonFockKet(2))
    bra1 = TensorProduct(BosonFockBra(0), BosonFockBra(0))
    bra2 = TensorProduct(BosonFockBra(1), BosonFockBra(2))
    assert qapply(TensorProduct(a, b ** 2) * ket1) == sqrt(2) * ket2
    assert qapply(TensorProduct(a, Dagger(b) * b) * ket1) == 2 * ket3
    assert qapply(bra1 * TensorProduct(a, b * b),
                  dagger=True) == sqrt(2) * bra2
    assert qapply(bra2 * ket1).doit() == TensorProduct(1, 1)
    assert qapply(TensorProduct(a, b * b) * ket1) == sqrt(2) * ket2
    assert qapply(Dagger(TensorProduct(a, b * b) * ket1),
                  dagger=True) == sqrt(2) * Dagger(ket2)


def test_dagger():
    lhs = Dagger(Qubit(0))*Dagger(H(0))
    rhs = Dagger(Qubit(1))/sqrt(2) + Dagger(Qubit(0))/sqrt(2)
    assert qapply(lhs, dagger=True) == rhs


def test_issue_6073():
    x, y = symbols('x y', commutative=False)
    A = Ket(x, y)
    B = Operator('B')
    assert qapply(A) == A
    assert qapply(A.dual*B) == A.dual*B


def test_density():
    d = Density([Jz*mo, 0.5], [Jz*po, 0.5])
    assert qapply(d) == Density([-hbar*mo, 0.5], [hbar*po, 0.5])


def test_issue3044():
    expr1 = TensorProduct(Jz*JzKet(S(2),S.NegativeOne)/sqrt(2), \
                                                    Jz*JzKet(S.Half,S.Half))
    result = Mul(S.NegativeOne, Rational(1, 4), 2**S.Half, hbar**2)
    result *= TensorProduct(JzKet(2,-1), JzKet(S.Half,S.Half))
    assert qapply(expr1) == result


# Issue 24158: Tests whether qapply incorrectly evaluates some ket*op as op*ket
def test_issue24158_ket_times_op():
    P = BosonFockKet(0) * BosonOp("a") # undefined term
    # Does lhs._apply_operator_BosonOp(rhs) still evaluate ket*op as op*ket?
    # qapply(P) -> BosonOp("a")*BosonFockKet(0) = 0 before fix
    assert qapply(P) == P
    P = Qubit(1) * XGate(0) # undefined term
    # Does rhs._apply_operator_Qubit(lhs) still evaluate ket*op as op*ket?
    assert qapply(P) == P   # qapply(P) -> Qubit(0) before fix
    # check legal expr <0| * (<1|*|1>) * X
    P1 = Mul(QubitBra(0), Mul(QubitBra(0), Qubit(0)), XGate(0))
    assert qapply(P1) == QubitBra(0) * XGate(0)     # qapply(P1) -> 0 before fix
    P1 = qapply(P1, dagger = True)  # classic: <0|*X(0), -> <1| w/ dagger=True
    assert qapply(P1, dagger = True) == QubitBra(1) # 0 before fix
    P2 = QubitBra(0) * QubitBra(0) * Qubit(0) * XGate(0) # 'forgot' to set brackets
    P2 = qapply(P2, dagger = True) # classic: <0|*X(0), -> <1| w/ dagger=True
    assert qapply(P2, dagger = True) == QubitBra(1) # 0 before fix
    # Pull Request 24237: IdentityOperator from the right and no dagger=True
    assert qapply(QubitBra(1)*IdentityOperator()) == QubitBra(1)
    assert qapply(IdentityGate(0)*(Qubit(0) + Qubit(1))) == Qubit(0) + Qubit(1)

#==============================================================================
#= Additional test cases for re-write
#==============================================================================

def test_suite2022_10():
    qb0 = QubitBra(0) # <0|
    qb1 = QubitBra(1) # <1|
    q0  = Qubit(0)    # |0>
    q1  = Qubit(1)    # |1>
    dyad11 = q1*qb1   # |1><1|
    dyad01 = q0*qb1   # |0><1|
    dyad10 = q1*qb0   # |1><0|

    assert qapply(qb1*q0) == 0 # honor ip_doit=True and evaluate InnerProduct
    assert qapply(dyad11**2 * dyad01**2) == 0  # evaluate dyads in powers
    assert qapply(qb1 * qb1 * q1 * q1 * q1) == q1 # 'missing brackets', right
    assert qapply(qb1 * qb1 * qb1 * q1 * q1) == qb1 # 'missing brackets' , left
    # deal with missing brackets that Mul simplifies to Pow objects
    assert qapply(qb1**2 * q1**3) == q1
    assert qapply(qb1**3 * q1**2) == qb1
    # legal expr <0| * (<1|*|1>) * X
    P1 = Mul(QubitBra(0), Mul(QubitBra(0), Qubit(0)), XGate(0))
    assert qapply(P1, dagger = True) == qb1  # classic: was 0, then <0|*X(0).


    assert qapply(TensorProduct(dyad11,dyad11) * \
                TensorProduct(dyad11,dyad01)) == TensorProduct(dyad11, 0)
    assert qapply(TensorProduct(dyad11, Mul(q1, qb1)) * \
                TensorProduct(dyad10, dyad10)) == TensorProduct(dyad10,dyad10)
    U = UnitaryOperator('U')
    assert qapply(TensorProduct(dyad11*dyad11, U)) == TensorProduct(dyad11, U)
    assert qapply(TensorProduct((dyad11+dyad01)**2, U), tensorproduct=False) \
                                            == TensorProduct(dyad11+dyad01, U)

    N = symbols("N", integer=True)
    IN = IdentityOperator(N) #N-dim Identity
    T1 = TensorProduct(dyad11,IN)
    T2 = TensorProduct(dyad11, U)
    P = T1 * T2
    assert qapply(P) == TensorProduct(dyad11, U)

    P = qb1 * qb1*q1 * XGate(0) # <0|
    assert qapply(P) == qb1*XGate(0)
    assert qapply(P, dagger = True) == qb0

    P = qb1 * qb1*q1 * q1 # 1
    assert qapply(P) == 1
    assert qapply(P, dagger = True) == 1

    P = dyad11 * dyad01 # 0
    assert qapply(P, ip_doit = False) == InnerProduct(qb1, q0)*dyad11
    assert qapply(P) == 0
    assert qapply(P, dagger = True) == 0

    P = dyad01 * dyad01 * dyad01 # 0
    assert qapply(P, ip_doit = False) == InnerProduct(qb1, q0)**2 * dyad01
    assert qapply(P) == 0
    assert qapply(P, dagger = True) == 0

    P = dyad11 * dyad11 * dyad11 # |1><1|
    assert qapply(P, ip_doit = False) == InnerProduct(qb1, q1)**2 * dyad11
    assert qapply(P) == dyad11
    assert qapply(P, dagger = True) == dyad11

    dyadphipsi = Ket('phi') * Bra('psi')
    P = qapply(dyadphipsi * dyadphipsi * dyadphipsi)
    assert qapply(P, ip_doit = False) == \
                        InnerProduct(Bra('psi'), Ket('phi'))**2 * dyadphipsi
    assert qapply(P) == InnerProduct(Bra('psi'), Ket('phi'))**2 * dyadphipsi
    assert qapply(P, dagger = True) == \
                        InnerProduct(Bra('psi'), Ket('phi'))**2 * dyadphipsi

    P = (TGate(0) ** 8) * Qubit(1) # Expect |1>
    assert qapply(P) == Qubit(1)   # classic qapply: -i*i*|1>
    P = (TGate(0) ** 8) * Qubit(0) # Expect |0>
    assert qapply(P) == Qubit(0)
    P = TensorProduct(XGate(0), XGate(0).inv()) * \
                                    TensorProduct(XGate(0).inv(), XGate(0))
    assert qapply(P) == TensorProduct(1, 1) # 1 x 1  classic: (XxX)**2

    assert qapply(Operator('A')*Ket(1)*Bra(1)) == \
            Operator('A') * OuterProduct(Ket(1), Bra(1)) # classic: Ket*Bra


def test_U2():
    # An extended example for qapply
    # Expression taken from Braun et al., Error Resilient Quantum Amplitude
    # Extimation, 2022, arXiv:2204.01337, 2022
    #
    # The following qubit registers are used:
    #   c, d : 1-qubit each, c is the master control register, d is its copy
    #   x, y : n qubit each, operator U operates on them
    # Their order is: c,x,d,y
    #
    # We want to compute the product c_Xd c_Ux d_Uy c_Xd
    # where d_Uy means U on register y controlled by d=1.

    U = UnitaryOperator('U')
    N = symbols("N", integer=True)

    I1 = IdentityOperator(1) #1-dim Identity
    IN = IdentityOperator(N) #N-dim Identity
    Iall = TensorProduct(I1, IN, I1, IN)
    #Dyads |1><1|, |0><0|
    dyad11 = Qubit(1) * QubitBra(1) # OuterProduct(|1>,<1|)
    dyad00 = Qubit(0) * QubitBra(0) # OuterProduct(|0>,<0|)
    #XGate build from Dyads
    XG = Qubit(1) * QubitBra(0) + Qubit(0) * QubitBra(1)

    c_Xd = TensorProduct(dyad11, IN, XG - I1, IN) + Iall
    # print("c_Xd:", c_Xd)
    c_Ux = TensorProduct(dyad11, U - IN, I1, IN) + Iall
    # print("c_Ux:", c_Ux)
    d_Uy = TensorProduct(I1, IN, dyad11, U - IN) + Iall
    # print("d_Uy:", d_Uy)

    # Compute the product c_Xd c_Ux d_Uy c_Xd
    P = c_Xd * c_Ux * d_Uy * c_Xd
    R = qapply(P, dagger=True) # option tensorproduct=True is on by default
    # Expected: IxIxIxI - IxIx|1><1|xI + IxIx|1><1|xU - |1><1|xIxIxI +
    # |1><1|xIx|1><1|xI - |1><1|xIx|1><1|xU + |1><1|xUx|0><0|xU + |1><1|xUx|1><1|xI
    assert R == \
             TensorProduct(I1,     IN,     I1, IN) \
           - TensorProduct(I1,     IN, dyad11, IN) \
           + TensorProduct(I1,     IN, dyad11,  U) \
           - TensorProduct(dyad11, IN,     I1, IN) \
           + TensorProduct(dyad11, IN, dyad11, IN) \
           - TensorProduct(dyad11, IN, dyad11,  U) \
           + TensorProduct(dyad11,  U, dyad00,  U) \
           + TensorProduct(dyad11,  U, dyad11, IN)


def test_doc_example(): # the example from qapply doc string
    b = Bra('b')
    k = Ket('k')
    A = k * b
    assert qapply(A * b.dual / (b * b.dual)) == k
    assert qapply(k.dual * A / (k.dual * k)) == b
    assert qapply(QubitBra(1) * XGate(0), dagger = True) == QubitBra(0)
    assert qapply(A * k * b * A, op_join = False) == \
                                        InnerProduct(b, k)**2 * k * b
    assert qapply(A * k * b * A) == InnerProduct(b, k)**2 * A
    n = symbols("n", integer=True, nonnegative=True)
    assert qapply(A ** (n + 2), power_exp=False) == \
                                        InnerProduct(b, k) ** (n + 1) * A


def test_issue19540():  # Issue #19540 https://github.com/sympy/sympy/issues/19540
    k1 = OrthogonalKet(1)
    k2 = OrthogonalKet(2)
    lh = TensorProduct(k1, k2)           # |1>|2>
    rh = TensorProduct(k1 + k2, k1 - k2) # |1>|1> - |1>|2> + |2>|1> - |2>|2>

    # InnerProduct(Dagger(lh), rh) throws a TypeError (as described in
    # issue #19540)
    # Workaround: Just use Dagger(lh)*rh. It's a completely legal product
    # and qapply(Dagger(lh)*rh) computes it. What is missing is that a
    # TensorProduct of all scalars is reduced to a scalar (could that be added
    # to TensorProduct.doit()?)
    assert qapply(Dagger(lh)*rh) == -TensorProduct(1,1)


def test_suite2022_11():

    # Handle negative powers sensibly, if the operator knows its inverse better
    # than just op.inv()==Pow(op,-1)
    assert qapply(TGate(0)**(8) * Qubit(1))  == Qubit(1)  # classic: -I*I*|1>
    # classic: recursion exponent to -infinity:
    assert qapply(TGate(0)**(-8) * Qubit(1)) == TGate(0)**(-8) * Qubit(1)
    assert qapply(TGate(0)**(-1) * Qubit(1)) == TGate(0)**(-1) * Qubit(1)
    assert qapply(XGate(0)**(-1) * Qubit(1)) == Qubit(0)

    In = IdentityOperator()
    P = TensorProduct(XGate(0), TGate(0)) * TensorProduct(XGate(0), In)
    U = Operator('U')
    n = symbols("n")
    # PR #24178:
    assert qapply(P) == TensorProduct(1, TGate(0)) # instead of I x T
    assert TensorProduct(U * U ** -1, pow(U, n) * U ** -n) == TensorProduct(1, 1)

    n = 7
    U = UGate(0, ImmutableMatrix([[cos(2 * pi / n), -sin(2 * pi / n)], \
                                  [sin(2 * pi / n), cos(2 * pi / n)]]))
    assert qapply(U * Qubit(0)) == \
                        cos(2 * pi / n)*Qubit(0) + sin(2 * pi / n)*Qubit(1)
    # qapply relies on Pow(U,2) to compute U*U!
    assert qapply(U * IdentityOperator(2) * U) == Pow(U, 2)
    assert qapply(U * IdentityOperator(2) * U ** -3) == Pow(U, -2) # addition of exps
    assert qapply((U ** 2 * Qubit(0))) == - sin(2 * pi / n)**2 * \
                        Qubit(0) + cos(2 * pi / n)**2 * Qubit(0) + \
                               2 * sin(2 * pi / n) * cos(2 * pi / n) * Qubit(1)
    m = 2 # simplification with trigsimp doesn't work for m>=3
    P2b = qapply(U ** m * Qubit(0), mul=False) # non-expanded handled by trigsimp()
    Ptsb = P2b.trigsimp() # won't-> cos(m*2*pi/n)*|0>+sin(m*2*pi/n)*|1> for m >=3
    assert Ptsb == (cos(m * 2 * pi / n) * Qubit(0) + \
                    sin(m * 2 * pi / n) * Qubit(1)).trigsimp()

    n = symbols("n", Integer = True)
    U = UGate(0, ImmutableMatrix([[cos(2 * pi / n), -sin(2 * pi / n)], \
                                  [sin( 2 * pi / n), cos(2 * pi / n)]]))
    assert qapply(U * Qubit(0)) == \
                    cos(2 * pi / n) * Qubit(0) + sin(2 * pi / n) * Qubit(1)
    # qapply relies on Pow(U,2) to compute U*U!
    assert qapply(Mul(U, IdentityOperator(2), U)) == Pow(U, 2)
    # addition of numeric exponents and of symbolic exponents
    assert qapply(Mul(U, IdentityOperator(2), U ** -3)) == Pow(U, -2)
    assert qapply(Mul(U ** n, IdentityOperator(2), U ** (-n + 2))) == Pow(U, 2)
    assert qapply(Mul(U ** n, IdentityOperator(2), U ** (-n + 2), Qubit(0)), \
        mul=False) == sin(2 * pi / n)*(-sin(2 * pi / n) * Qubit(0) + \
                                        cos(2 * pi / n) * Qubit(1)) + \
                      cos(2 * pi / n)*( cos(2 * pi / n) * Qubit(0) + \
                                        sin(2 * pi / n) * Qubit(1))
    assert qapply(Mul(U ** n, IdentityOperator(2), U ** (-n + 2), Qubit(0))) == \
                                  -sin(2 * pi / n)**2 * Qubit(0) + \
                2 * sin(2 * pi / n) * cos(2 * pi / n) * Qubit(1) + \
                                   cos(2 * pi / n)**2 * Qubit(0)

    m = 1 # trigsimp wont -> cos(m*2*pi/n)*|0>+sin(m*2*pi/n)*|1> for m == 2
    P2b = qapply(U ** m * Qubit(0), mul=False)
    Ptsb = P2b.trigsimp()
    assert Ptsb == (cos(m * 2 * pi / n) * Qubit(0) + \
                    sin(m * 2 * pi / n) * Qubit(1)).trigsimp()

    U = UGate(0, ImmutableMatrix([[E ** (2 * pi * I / n),  0], \
                                  [0, E ** -(2 * pi * I / n)]]))
    m = 5 # qapply relies on Pow, pow etc. to compute powers!
    assert qapply(U ** m) == U ** m
    assert qapply(U ** m * Qubit(0)) == E ** (m * 2 * I * pi / n) * Qubit(0)

    # simple addition of exponents
    assert qapply(TGate(0) ** (m - n) * IdentityOperator(2) * TGate(0) ** n) == \
        TGate(0) ** m # qapply relies on Pow for powers!
    assert qapply(TGate(0) ** (m - n) * IdentityOperator(2) * \
        TGate(0) ** n * Qubit(1)) == E ** -((8-m) * 2 * I * pi / 8) * Qubit(1)

    # Dagger of unitary matrices U is U**-1, test qapply behaviour
    U12 = UGate(0, ImmutableMatrix([[2, 0], [0, 2]]))
    assert qapply(QubitBra(1)*U12) == QubitBra(1)*U12 # no method for Bra*UGate
    # (Dagger(UGate)==UGate**-1) which doesn't know about |1>
    assert qapply(Dagger(U12)*Qubit(1), dagger=True) == U12**(-1) * Qubit(1)
    assert Dagger(Dagger(U12)) == U12
    # with dagger=True it boils down to Dagger(U12*Qubit(1))
    assert qapply(QubitBra(1)*Dagger(U12), dagger=True) == 2 * QubitBra(1)
    assert qapply(Dagger(Dagger(U12))*Dagger(QubitBra(1))) == 2 * Qubit(1)


def test_powers2022_11():
    from sympy import exp, I, pi
    from sympy.matrices import ImmutableMatrix, eye

    # The quantum package assumes that the user adds or modifies the methods
    # _apply_operator_XXX / _apply_from_the_right_to_XXX of operators if not in
    # source code then at runtime. In a single user environment the user
    # themselves may handle concurrency issues that arise from modifying the
    # global operator class at runtime. A safer approach is to clone the operator
    # class and locally modify the clone. As SymPy objects are supposed to be
    # non-mutable, class inheritance does suffice for this.
    from sympy.physics.quantum.gate import UGate as UGate_Org # global UGate
    class UGate(UGate_Org):     # UGate is now a local clone of the global UGate
        pass

    # Use an extended definition of UGate*UGate with diverse result types in order
    # to build test cases with commutative factors etc.
    def UGate_apply_operator_UGate(u1, u2, **options):
        if u1.args[0] != u2.args[0]: # This case is not implemented
            raise NotImplementedError("UGate*UGate on different target qubits")
        u1id = (u1.args[1] == u1.args[1][0,0] * eye(u1.args[1].shape[0]))
        u2id = (u2.args[1] == u2.args[1][0,0] * eye(u2.args[1].shape[0]))
        if u1id and u2id:
            return u1.args[1][0,0] * u2.args[1][0,0] # scalar result
        elif u1id:
            return u1.args[1][0,0] * u2 # scalar * UGate
        elif u2id:
            return u2.args[1][0,0] * u1 # scalar * UGate
        else:
            res = u1.args[1] * u2.args[1] # assuming res again immutableMatrix
            if res == res[0,0] * eye(res.shape[0]):
                return res[0,0]  # scalar result
            else:
                return UGate(u1.args[0], res)

    # Set the method on the local clone UGate of the global class UGate_Org
    prev_op = getattr(UGate, "_apply_operator_UGate", None) # cannot fail
    setattr(UGate, "_apply_operator_UGate", UGate_apply_operator_UGate) # no fail

    n = 7
    U71   = UGate(0, ImmutableMatrix([[1*exp(2*pi*I/n),  0],
                                      [0, 1*exp(-2*pi*I/n)]]))
    U71p2 = UGate(0, ImmutableMatrix([[1*exp(4*pi*I/n),  0],
                                      [0, 1*exp(-4*pi*I/n)]])) # == U71 ** 2
    U71p3 = UGate(0, ImmutableMatrix([[1*exp(6*pi*I/n),  0],
                                      [0, 1*exp(-6*pi*I/n)]])) # == U71 ** 3
    U71p5 = UGate(0, ImmutableMatrix([[1*exp(-4*pi*I/n), 0],   # 2023-01-31
                                      [0, 1*exp(4*pi*I/n) ]])) # == U71 ** 5
    U71p6 = UGate(0, ImmutableMatrix([[1*exp(-2*pi*I/n), 0],
                                      [0, 1*exp(2*pi*I/n) ]])) # == U71 ** 6
    U72   = UGate(0, ImmutableMatrix([[2*exp(2*pi*I/n), 0 ],
                                      [0, 2*exp(-2*pi*I/n)]]))
    U72p2 = UGate(0, ImmutableMatrix([[4*exp(4*pi*I/n), 0 ],
                                      [0, 4*exp(-4*pi*I/n)]])) # == U72 ** 2
    #U73 = UGate(0, ImmutableMatrix([[3*exp(2*pi*I/n), 0], [0, 3*exp(-2*pi*I/n)]]))
    #U11 = UGate(0, ImmutableMatrix([[1, 0], [0, 1]]))
    U12 = UGate(0, ImmutableMatrix([[2, 0], [0, 2]]))
    U13 = UGate(0, ImmutableMatrix([[3, 0], [0, 3]]))

    n = symbols("n", commutative=True)
    m = symbols("m", integer=True, nonnegative=True)
    U1n = UGate(0, ImmutableMatrix([[n, 0], [0, n]]))

    assert qapply((m * XGate(0)*U1n*XGate(0))**7) == \
            m**7 * n**6 * XGate(0) * U1n * XGate(0) # case 5a: with expi uneven
    assert qapply((m * XGate(0)*U1n*XGate(0))**8) == \
            m**8 * n**8 # case 5a: with expi even
    assert qapply((n * U72) ** Rational(3,2)) == \
            n * sqrt(n * U72) * U72 # case 7: expi < 2
    # case 6: base atomic but not "idempotent"
    assert qapply((n * U72) ** Rational(5,2)) == \
            n**2 * sqrt(n * U72) * U72p2
    assert qapply((n * U12 * U71) ** Rational(5,2)) == \
            n**2 * 2**Rational(5,2) * sqrt(n*U71) * U71p2 # Case 9, 6
    # case 5a. base atomic + "involutoric", expi uneven, fraction
    assert qapply((n * U12) ** Rational(7,2)) == n**3 * 2**2 * sqrt(n * U12) * U12
    # case 5a. base atomic + "involutoric", expi even, fraction
    assert qapply((n * U12) ** Rational(5,2)) == n**2 * 2**2 * sqrt(n * U12)
    # case 5a: expi uneven, no fraction
    assert qapply((n * U12) ** 3) == n**3 * 2**2 * U12
    # 5a. case base atomic + involutoric, expi even, no fraction
    assert qapply((n * U12) ** 4) == n**4 * 2**4

    # case 5b: base atomic + "idempotent"
    k1 = Qubit(1); b1 = QubitBra(1)
    assert qapply((n * U12 * OuterProduct(k1, b1)) ** Rational(7/2)) == \
        n ** 3 * 2 ** Rational(7,2) * InnerProduct(b1, k1).doit() ** 2 * \
        sqrt(n * OuterProduct(k1, b1)) * OuterProduct(k1, b1)
    assert qapply((n * U12 * OuterProduct(k1, b1)) ** Rational(7/2), \
        ip_doit = False) == \
                n ** 3 * 2 ** Rational(7,2) * InnerProduct(b1, k1) ** 2 * \
                sqrt(n * OuterProduct(Qubit(1), QubitBra(1))) * \
                OuterProduct(Qubit(1), QubitBra(1))
    assert qapply((n * U12 * OuterProduct(Qubit(1), QubitBra(1))) ** 10) == \
             n**10 * 2**10 * OuterProduct(Qubit(1), QubitBra(1))
    assert qapply((m * U1n * OuterProduct(Qubit(1), QubitBra(1))) ** 10) == \
             m**10 * n**10 * OuterProduct(Qubit(1), QubitBra(1))
    assert qapply((m * U1n * OuterProduct(Qubit(1), QubitBra(1))) ** (m+2), \
        power_exp=False) == m**(m+2) * n**(m+2) * OuterProduct(k1, b1)
    assert qapply((m * U1n * OuterProduct(k1, b1)) ** (m+2), \
        ip_doit = False, power_exp=False) == \
        m**(m+2) * n**(m+2) * InnerProduct(b1, k1)**(m+1) * OuterProduct(k1, b1)
    assert qapply((n * OuterProduct(Ket(1), Bra(1))) ** (m+2), \
        ip_doit = False, power_exp=False) == n**(m+2) * \
        InnerProduct(Bra(1), Ket(1))**(m+1) * OuterProduct(Ket(1), Bra(1))
    assert qapply((n * OuterProduct(Ket(1), Bra(1))) ** (m+Rational(5,2)), \
        ip_doit = False, power_exp=False) == \
        n**(m+2) * InnerProduct(Bra(1), Ket(1))**(m+1) * \
        sqrt(n * OuterProduct(Ket(1), Bra(1))) * OuterProduct(Ket(1), Bra(1))

    # Cases 7, 9, 6:
    assert qapply((m*U12*U71) ** Rational(3,2) * (n*U13*U71) ** Rational(11,2)) \
        ==  n**5 * m**Rational(3,2) * 2**Rational(3,2) * 3**Rational(11,2) * \
            U71**Rational(3,2) * sqrt(n*U71) * U71p5 # 2023-01-31 7.Case expi < 2
    assert qapply((n*U13*U71) ** Rational(11,2) * (m*U12*U71) ** Rational(3,2)) \
        ==  n**5 * m**Rational(3,2) * 2**Rational(3,2) * 3**Rational(11,2) * \
            sqrt(n*U71) * sqrt(U71) * U71p6 # 2023-01-31 7.Case expi  < 2
    assert qapply((n * U13*U71) ** 5 * (m * U12*U71) ** 2) == \
            n**5 * m**2 * 2 ** 2 * 3 ** 5 # Case 9, 6

    # case 3: component base with interaction in base*base
    #print(qapply((U71*XGate(0)*U71*U13)))
    assert qapply((n*U71*XGate(0)*U71*U13)**3) == \
        n**3 * 3**3 * U71 * XGate(0) * U71p2 * XGate(0) * U71p2 * XGate(0) * U71
    # case 1: component base with no interaction in base*base
    assert qapply((n * U71*XGate(0))**3) == n**3 * (U71*XGate(0)) ** 3
    assert qapply((n * U71*XGate(0))**3 * Qubit(0)) == n**3 * exp(-2*I*pi/7) * k1
    # case 2: component base, with interaction in base*base, commutative
    assert qapply((n * XGate(0)*U71*U13*XGate(0))**3) == \
                                    n**3 * 3**3 * XGate(0) * U71p3 * XGate(0)
    assert qapply((m * XGate(0)*U71*U13*XGate(0))**Rational(7,2)) == \
        m**Rational(7,2) * 3**Rational(7,2) * \
        (XGate(0) * U71 * XGate(0))**Rational(1,2) * XGate(0) * U71p3 * XGate(0)

    # no need to restore the previous operator _apply_operator_UGate on the
    # local clone UGate. Done anyway as code demo.
    if prev_op: # restore previous operator
        setattr(UGate, "_apply_operator_UGate", prev_op)
    else:
        delattr(UGate, "_apply_operator_UGate")


def test_expansion():
    # check expansion behaviour on options mul and power_base,
    # power_exp on commutative factors. qapply() applies hint
    # functions for options  mul, power_base and power_exp of expand()

    A1 = Operator("A1"); A2 = Operator("A2"); A3 = Operator("A3")

    p = (2*(A1 + A2) + 2*(A1 - A2))*A3  # distribute factor A3 over sums
    assert qapply(p, mul=False) == 2*(A1*A3 - A2*A3) + 2*(A1*A3 + A2*A3)
    # distribute commutative factors over sums, too  (default)
    assert qapply(p, mul=True) == 4*A1*A3

    u, v = symbols("u v", commutative=True)
    w, z = symbols("w z", integer=True, positive=True)

    p = (u*(A1 + A2) + u*(A1 - A2))*A3
    assert qapply(p, mul=False) == u*(A1*A3 - A2*A3) + u*(A1*A3 + A2*A3)
    assert qapply(p, mul=True) == 2*u*A1*A3 # mul=True is default

    p = (u*(A1 + A2) + v*(A1 - A2))*A3
    assert qapply(p, mul=False) == u*(A1*A3 + A2*A3) + v*(A1*A3 - A2*A3)
    assert qapply(p, mul=True) == u*A1*A3 + u*A2*A3 + v*A1*A3 - v*A2*A3

    p = (u*(A1 + A2) + u*(A1 - A2))**2 * A3
    r = u*(u*(A1**2*A3 - A2*A1*A3) + u*(A1**2*A3 + A2*A1*A3) + \
        Mul(-1, (u*(A1*A2*A3 - A2**2*A3) + u*(A1*A2*A3 + A2**2*A3)), evaluate=False)) +\
        u*(u*(A1**2*A3 - A2*A1*A3) + u*(A1**2*A3 + A2*A1*A3) +  \
        u*(A1*A2*A3 - A2**2*A3) + u*(A1*A2*A3 + A2**2*A3) )
    assert qapply(p, mul=False) == r
    assert qapply(p) == 4 * u**2 * A1**2 * A3 # mul=True is default

    # Check  effects of expand options power_base and power_exp on
    # communicative factors. Note that sqrt(u*(u+1) + w*(w+1)) is unmodified.
    p = u*(w + sqrt(u*(u + 1) + w*(w + 1)))**(w + z)
    # commutative factors left untouched
    assert qapply(p, mul=False, power_base=False, power_exp=False) == p
    # mul/add expansion has no effect on p
    assert qapply(p, mul=True, power_base=False, power_exp=False) == p
    # but power_exp has
    assert qapply(p, mul=False, power_base=False, power_exp=True) == \
        u*(w + sqrt(u*(u + 1) + w*(w + 1)))**w*(w + sqrt(u*(u + 1) + w*(w + 1)))**z
    assert qapply(p, mul=True, power_base=False, power_exp=True) == \
        u*(w + sqrt(u*(u + 1) + w*(w + 1)))**w*(w + sqrt(u*(u + 1) + w*(w + 1)))**z

    p = (u*(w + sqrt(u*(u + 1) + w*(w + 1))))**(w + z)
    # commutative factors left untouched
    assert qapply(p, mul=False, power_base=False, power_exp=False) == p
    # power_base has effect
    assert qapply(p, mul=False, power_base=True, power_exp=False) == \
        u**(w + z) * (w + sqrt(u*(u + 1) + w*(w + 1)))**(w + z)
    # As exception mul has effect in commutative in basis IF Pow is the only factor
    assert qapply(p, mul=True, power_base=False, power_exp=False) == \
        (u*w + u*sqrt(u*(u + 1) + w*(w + 1)))**(w + z)
    assert qapply(p, mul=False, power_base=True, power_exp=True) == \
        u**w * u**z * (w + sqrt(u*(u + 1) + w*(w + 1)))**w * \
        (w + sqrt(u*(u + 1) + w*(w + 1)))**z # power_xxx has effect


def test_powers2023_01():

    A, B = symbols('A B',commutative=False)
    c, d = symbols("c d", commutative=True)
    f, g = symbols("f g", commutative=True, nonnegative=True)
    o = symbols("o", commutative=True, integer=True, positive=True)
    r = symbols("r", real=True)
    In = IdentityOperator()

    # Results verified using Pow() resp. expand().
    # Also shows that qapply essentially rebuilds expand() without deep=True
    p = Pow(c*d*f*In, c)  # all commutative, f positive, unknown exponent
    # IdentityOperator I remains because c, d are scalars. f is pulled out.
    assert qapply(p) == f**c*(c*d*In)**c
    p = Pow(c*d*A*In, -o)  # negative integer exponent. Fully expands.
    assert qapply(p) == c**(-o)*d**(-o)*A**(-o)
    p = Pow(c*d*A*In, Rational(1/2)) # fractional exponent. No extraction.
    assert qapply(p) == sqrt(c*d*A)
    p = Pow(c*d*f*A*In, Rational(1/2))
    # fractional exponent, f non-negative. expands() same
    assert qapply(p) == sqrt(f)*sqrt(c*d*A)
    p = Pow(c*d*f*A*In, c)  # arbitrary exponent. expand() same
    assert qapply(p) == f**c*(c*d*A)**c

    # numeric positive integer exponents; Pow() extracts same
    p = Pow(c*d*f*A*In, 4)*Pow(d*g*A*In, 3)
    assert qapply(p) == c**4*d**7*f**4*g**3*A**7
    # symbolic positive integer exponents; expand() extracts same
    p = Pow(c*d*f*A*In, o+4)*Pow(d*g*In*A, o)
    assert qapply(p) == c**(o + 4)*d**o*d**(o + 4)* \
                                                f**(o + 4)*g**o*A**(2*o + 4)
    # symbolic positive and negative integer exponents; expand extracts same
    p = Pow(c*d*f*A*In, o+4)*Pow(d*g*In*A, o-3)
    assert qapply(p) == c**(o + 4)*d**(o - 3)*d**(o + 4)* \
                                          f**(o + 4)*g**(o - 3)*A**(2*o + 1)
    # symbolic real exponents non-integer; expand extracts same
    p = Pow(d*f*A*In, o+r)*Pow(d*g*In*A, o-r)
    assert qapply(p) == d**(2*o) * f**(o + r) * g**(o - r) * A**(2*o) #2023-01-31
    # symbolic real exponents; expand extracts same
    p = Pow(d*f*A*In, r)*Pow(d*g*A*In, -r)
    assert qapply(p) == f**r*g**(-r)
    # numeric rational exp; expand extracts same
    p = Pow(c*A*In, 1+Rational(1,2))*Pow(d*A, Rational(1,2))
    assert qapply(p) == c*sqrt(c*A)*A*sqrt(d*A)

    # positive integer exponents, nested pos integer exp; expand extracts same
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*B*In, 3), 7)
    assert qapply(p) == c**(7*o + 21)*d**21*f**(7*o + 21)*g**21* \
                                                        (A**(o + 3)*B**3)**7
    # positive integer exponents, nested pos integer exp
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*A, 3), 7)
    assert qapply(p) == c**(7*o + 21)*d**21*f**(7*o + 21)*g**21*A**(7*o + 42)
    # pos int exp, unknown int exp, nested pos int exp: expand extracts same
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*A, 3-o), o+2)
    assert qapply(p) == c**(o**2 + 5*o + 6)*d**(-o**2 + o + 6)* \
                        f**(o**2 + 5*o + 6)*g**(-o**2 + o + 6)*A**(6*o + 12)
    assert qapply(p, nested_exp=False) == c**((o + 2)*(o + 3))* \
      d**((3 - o)*(o + 2))*f**((o + 2)*(o + 3))*g**((3 - o)*(o + 2))*A**(6*o + 12)

    # pos int exp, real exp, nested pos int exp: expand extracts same
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*A, r), o+2)
    assert qapply(p) == c**(o**2 + 5*o + 6)*f**(o**2 + 5*o + 6)* \
                        g**(o*r + 2*r)*(A**(o + 3)*(d*A)**r)**(o + 2)
    assert qapply(p, nested_exp=False) == c**((o + 2)*(o + 3))* \
        f**((o + 2)*(o + 3))*g**(r*(o + 2))*(A**(o + 3)*(d*A)**r)**(o + 2)

    # pos int exp, unknown int exp, nested real exp: expand extracts same
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*A, 3-o), r)
    assert qapply(p) == f**(o*r + 3*r)*(c**(o + 3)*d**(3 - o)*g**(3 - o)*A**6)**r
    assert qapply(p, nested_exp=False) == f**(r*(o + 3))* \
                                 (c**(o + 3)*d**(3 - o)*g**(3 - o)*A**6)**r

    # pos exp, unknown int exp, nested unknown exp: expand extracts same
    p = Pow(Pow(c*f*A*In, f)*Pow(d*g*A, 3-o), c)
    assert qapply(p) == f**(c*f)*(d**(3 - o)*g**(3 - o)*(c*A)**f*A**(3 - o))**c
    # unknown exp, unknown integer exp, nested unknown exp: no extraction
    p = Pow(Pow(c*f*A*In, d)*Pow(d*g*A, 3-o), c)
    assert qapply(p) == (d**(3 - o)*f**d*g**(3 - o)*(c*A)**d*A**(3 - o))**c
    # qapply won't expand the exponent, expand does
    p = Pow(c*f*A*In, f*(f+c)*(d-g))
    assert qapply(p) == f**(f*(c + f)*(d - g))*(c*A)**(f*(c + f)*(d - g))

    from sympy import evaluate, expand
    with evaluate(False): # creates objects in non-canonical form
        p0 = Pow(Pow(Pow(Pow(2*A, c), d), c), d)
        p1 = c**A * c**B
        p2 = c**A * c**d # c**A * c**d unsorted
        p3 = Pow(p2, d*(A+B))
        p4 = Pow(c, A*(A*c + d)*B)

    # prove coincidence with expand()
    assert qapply(p0) == expand(p0)
    assert qapply(p1) == expand(p1) # c**A * c**B, not a**(A+B)
    assert qapply(p2) == expand(p2.doit())
    assert qapply(p3) == p3.doit()  # exponent d*(A+B) un-expanded
    assert qapply(p3, apply_exp=True) == expand(p3.doit())
    assert qapply(p4, apply_exp=True) == c**(c*A**2*B + d*A*B)

def test_density2023_01():
    # class Density only capable of (unitary) operator from the left
    d1 = Density((Qubit(0), 0.3), (Qubit(1), 0.7))
    assert qapply(XGate(0)*d1) == Density((Qubit(1), 0.3), (Qubit(0), 0.7))

def test_powers2023_02_11():
    # Use a diversified definition of UGate*UGate in order
    # to build test cases with commutative factors etc.
    from sympy.physics.quantum.gate import UGate as UGate_Org # global UGate
    class UGate(UGate_Org):     # UGate is now a local clone of the global UGate
        def _apply_operator_UGate(u1,u2, **options):
            from sympy.matrices import eye
            if u1.args[0] != u2.args[0]: # This case is not implemented
                raise NotImplementedError("UGate*UGate on different target qubits")
            u1id = (u1.args[1] == u1.args[1][0,0] * eye(u1.args[1].shape[0]))
            u2id = (u2.args[1] == u2.args[1][0,0] * eye(u2.args[1].shape[0]))
            if u1id and u2id:
                return u1.args[1][0,0] * u2.args[1][0,0] # scalar result
            elif u1id:
                return u1.args[1][0,0] * u2 # scalar * UGate
            elif u2id:
                return u2.args[1][0,0] * u1 # scalar * UGate
            else:
                res = u1.args[1] * u2.args[1] # assuming res again immutableMatrix
                if res == res[0,0] * eye(res.shape[0]):
                    return res[0,0]  # scalar result
                else:
                    return UGate(u1.args[0], res)

    c, d = symbols("c d", commutative=True)
    # expansion behaviour of x**(a+b) changed with 1.12, PR #23962
    f, g = symbols("f g", commutative=True, positive=True)
    o, u = symbols("o u", commutative=True, integer=True, positive=True)
    r = symbols("r", real=True)

    X = XGate(0)
    Ud = UGate(0, ImmutableMatrix([[1, 2], [0, 2]]))
    Ud4 = qapply(Ud ** 4)
    # Expansion of exp ((c+d)*(f+g))! when apply_exp=True
    p = Pow(Ud, ((c+d)*(f+g)), evaluate=False) * Ud**o * Ud**3
    q = qapply(p, apply_exp=True)
    assert q == Ud**(c*f + c*g + d*f + d*g + o - 1) * Ud4

    # test nested powers
    p = Pow(Ud, o+1) * Pow(Pow(Ud, o+1), o)
    q = qapply(p, nexted_exp=True, apply_exp=True)
    assert q == Ud**(o**2 + 2*o - 3) * Ud4
    p = Pow(Ud, r) * Pow(Pow(Ud, r), o)
    q = qapply(p, nexted_exp=True)
    assert q == Ud**(o*r +r)

    # test involutoric operator
    p = Pow(X, o+10) * Pow(X, o+8)
    assert unchanged(Pow, X, o+10) and unchanged(Pow, X, o+8)
    assert unchanged(Mul, Pow(X, o+10), Pow(X, o+8))
    q = qapply(p, apply_exp=True)
    assert q == 1 # involutoric operator simplified

    p = Pow(Ud, o+17)
    q = qapply(p, apply_exp=True)
    assert q == Ud**(o-1) * \
        UGate(0, ImmutableMatrix([[1, 524286],[0, 2**18]]))

    p = Pow(Ud, 10*o+8) # sympy knows that this is >= 18!
    q = qapply(p, apply_exp=True)
    assert q == Ud**(10*o-10) * \
        UGate(0, ImmutableMatrix([[1, 524286],[0, 2**18]]))

    p = Pow(Ud, 5*(o+1)*(u+1))
    q = qapply(p, apply_exp=True) # exp=5*o*u+5*o+5*u+5, sympy verifies >=5 only
    assert q == Ud**(5*o*u+5*o+5*u) * \
        UGate(0, ImmutableMatrix([[1, 62],[0, 2**5]]))

    p = Pow(Ud, 5*(o+1)*(u+1))
    q = qapply(p) # exp = (5*o+5)*(u+1), sympy verifies >= 4 only
    assert q == Ud**((5*o+5)*(u+1)-4) * \
        UGate(0, ImmutableMatrix([[1, 30],[0, 2**4]]))

    p = Pow(Ud, 10*o-3*u) # this is not verifiable >= 0
    q = qapply(p, apply_exp=True)
    assert q == Ud**(10*o-3*u)

    e = Pow(Pow(o*f*g, o+3)*Pow(f*g*o, 3-o)+17, o+2)
    ep = qapply(e, power_exp=True) # power_exp simplifies
    p = Pow(Ud, e+18)
    q = qapply(p, apply_exp=True, power_exp=True)
    # expansion behaviour of x**(a+b) changed with 1.12, PR #23962
    assert q == Ud**ep * \
        UGate(0, ImmutableMatrix([[1, 524286], [0, 262144]]))
