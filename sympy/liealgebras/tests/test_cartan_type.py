from sympy.core.numbers import Rational
from sympy.liealgebras.cartan_type import CartanType
from sympy import eye, Matrix, S

def test_Standard_Cartan():
    c = CartanType("A4")
    assert c.rank == 4
    assert c.series == "A"

    b = CartanType("B12")
    assert b.rank == 12
    assert b.series == "B"

def test_CartanType():
    c = CartanType("A4")

    omega = c.omega_matrix()
    cocar = c.cocartan_matrix()

    assert cocar * omega.T == eye(4)

    c2 = CartanType("A2")
    assert c2.omega_matrix() == Matrix([
                        [Rational(2,3), Rational(-1,3), Rational(-1,3)],
                        [Rational(1,3),  Rational(1,3), Rational(-2,3)]])

    assert c2.fundamental_weight(1) == Matrix([[Rational(2,3), Rational(-1,3), Rational(-1,3)]])

def test_CartanTypeBasisOverride():
    c = CartanType("F4")
    omega_orig = c.omega_matrix()
    cartan_orig = c.cartan_matrix()

    my_root_basis = [
        Matrix([[1,-1,0,0]]),
        Matrix([[0,1,-1,0]]),
        Matrix([[0,0,1,0]]),
        Matrix([[-S.Half, -S.Half, -S.Half, -S.Half]]),
    ]

    c.simple_roots = my_root_basis
    omega_should_changed = c.omega_matrix()
    new_cartan = c.cartan_matrix()
    # cartan is basis independent while omega
    # is basis dependent
    assert new_cartan == cartan_orig
    assert omega_orig != omega_should_changed

    prnew = c.positive_roots()

    assert prnew == [
        Matrix([[1, 0, 0, -1]]),
        Matrix([[0, 1, 0, -1]]),
        Matrix([[0, 0, 1, -1]]),
        Matrix([[0, 0, 0, -1]]),
        Matrix([[0, 0, -1, -1]]),
        Matrix([[1/2, 1/2, 1/2, -1/2]]),
        Matrix([[0, -1, 0, -1]]),
        Matrix([[1/2, 1/2, -1/2, -1/2]]),
        Matrix([[-1, 0, 0, -1]]),
        Matrix([[1/2, -1/2, 1/2, -1/2]]),
        Matrix([[1, 1, 0, 0]]),
        Matrix([[-1/2, 1/2, 1/2, -1/2]]),
        Matrix([[1/2, -1/2, -1/2, -1/2]]),
        Matrix([[1, 0, 1, 0]]),
        Matrix([[-1/2, 1/2, -1/2, -1/2]]),
        Matrix([[0, 1, 1, 0]]),
        Matrix([[1, 0, 0, 0]]),
        Matrix([[-1/2, -1/2, 1/2, -1/2]]),
        Matrix([[0, 1, 0, 0]]),
        Matrix([[1, 0, -1, 0]]),
        Matrix([[-1/2, -1/2, -1/2, -1/2]]),
        Matrix([[0, 0, 1, 0]]),
        Matrix([[0, 1, -1, 0]]),
        Matrix([[1, -1, 0, 0]])]
