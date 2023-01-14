from sympy.core.symbol import Symbol
from sympy.core.numbers import pi
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.algebras.quaternion import Quaternion
from sympy.algebras.euler import Euler
from sympy.testing.pytest import raises, warns
from sympy.core.sympify import sympify

from sympy.abc import x, y, z, alpha, beta, gamma
from sympy.abc import a, b, c, d


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


def test_euler_is_zero_rotation():
    assert Euler(0, 0, 0).is_zero_rotation
    assert not Euler(1, 0, 0).is_zero_rotation


def test_euler_basic_operations():
    with warns(UserWarning, match='Adding', test_stacklevel=False):
        euler1 = Euler(0.1, 0.2, 0.3)
        euler2 = Euler(0.2, 0.4, 0.6)
        assert -euler1 == Euler(-0.1, -0.2, -0.3)
        assert euler1 + euler1 == euler2
        assert 2 * euler1 == euler2
        assert euler2 - euler1 == euler1
        assert euler2 / 2 == euler1
        assert euler2 / 2 == euler1
        assert 3 * Euler(pi/2, 0, pi / 4) == Euler(- pi/2, 0, 3 * pi / 4)
        raises(ValueError, lambda: 2 / euler1)


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


def test_euler_from_quaternion():
    # since Euler is just a convenience class, tests on Quaternion are already
    # enough
    q = Quaternion(a, b, c, d)
    seq = 'zyz'
    euler1 = q.to_euler(seq)
    euler2 = Euler.from_quaternion(q, seq)
    euler3 = Euler.from_quaternion(q.args, seq)
    assert euler1 == euler2.args
    assert euler2 == euler3


def test_euler_mul_with_euler():
    euler1 = Euler(pi/2, 0, 0, 'ZYX')
    euler2 = Euler(pi/2, 0, 0, 'ZXY')
    with warns(UserWarning, match='Euler', test_stacklevel=False):
        assert euler1 * euler1 == 2 * euler1  # simple case, only one angle

    with warns(UserWarning, match='Different', test_stacklevel=False):
        assert euler1 * euler2 == 2 * euler1
        assert euler2 * euler1 == 2 * euler2


def test_euler_evalf():
    euler1 = Euler(1 / 2, 0, 0)
    euler2 = Euler(1 / (sympify(2)), 0, 0)
    assert euler1 != euler2
    assert euler1 == euler2.evalf()


def test_euler_set_sequence():
    q = Quaternion(1, 0, 0, 1)
    euler1 = Euler.from_quaternion(q, 'xyz')
    euler2 = Euler.from_quaternion(q, 'zyx')
    assert euler1 != euler2
    assert euler1 == euler2.set_sequence(euler1.seq)
