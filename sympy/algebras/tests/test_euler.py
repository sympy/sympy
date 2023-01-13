from sympy.core.symbol import Symbol
from sympy.core.numbers import pi
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.algebras.quaternion import Quaternion
from sympy.algebras.euler import Euler
from sympy.testing.pytest import raises

from sympy.abc import x, y, z, alpha, beta, gamma


def test_euler_construction():
    assert Euler(alpha, beta, gamma).alpha == alpha

    raises(ValueError, lambda: Euler(alpha, beta, gamma, seq='ZYZY'))
    raises(ValueError, lambda: Euler(alpha, beta, gamma, seq='ZyZ'))
    raises(ValueError, lambda: Euler(alpha, beta, gamma, seq='ZAZ'))
    raises(ValueError, lambda: Euler(alpha, beta, gamma, seq='ZA'))
    raises(ValueError, lambda: Euler(alpha, beta, gamma, seq='zzx'))
    raises(ValueError, lambda: Euler(alpha, beta, gamma, seq='zxx'))

    nc = Symbol('nc', commutative=False)
    raises(ValueError, lambda: Euler(alpha, nc, gamma))

    angles = (3*pi, 4*pi, 5*pi)
    assert Euler(*angles).args == (-pi, 0, -pi)

    euler = Euler(x, y, z)
    euler = euler.subs([(i, angle) for i, angle in zip(euler.args, angles)])
    assert euler.args == (-pi, 0, -pi)


def test_euler_to_quaternion():
    # since Euler is just a convenience class, tests on Quaternion are already
    # enough
    assert Euler(0, 0, 0).to_quaternion() == Quaternion(1, 0, 0, 0)
    roll = Euler(pi/2, 0, 0, seq='ZYX').to_quaternion()
    pitch = Euler(0, pi/2, 0, seq='ZYX').to_quaternion()
    yaw = Euler(0, 0, pi/2, seq='ZYX').to_quaternion()
    assert roll == Quaternion(sqrt(2)/2, 0, 0, sqrt(2)/2)
    assert pitch == Quaternion(sqrt(2)/2, 0, sqrt(2)/2, 0)
    assert yaw == Quaternion(sqrt(2)/2, sqrt(2)/2, 0, 0)


def test_euler_inverse():
    angles = [
            Euler(1, 2, 3, 'ZYX'),
            Euler(1, 2, 3, 'ZXZ'),
            Euler(1, 2, 3, 'zyx'),
            Euler(1, 2, 3, 'xyx')]
    for angle in angles:
        test1 = angle.to_quaternion().inverse()
        test2 = angle.inverse().to_quaternion()
        assert test1.evalf() == test2.evalf()
