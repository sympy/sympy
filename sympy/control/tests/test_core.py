"""
sympy.control.core unit tests
"""

from sympy.control.core import tf, ss, tf2ss, ss2tf
from sympy.matrices import Matrix


def test_tf():
    """
    tranfer function creation
    """
    tf1 = tf([1, 2], [1, 2, 3])
    assert tf1 is not None


def test_ss():
    """
    state space creation
    """
    ss1 = ss([1], [1], [1], [1], 0)
    assert ss1 is not None


def test_tf_ss_conv():
    """
    transfer function to space space conversion
    """
    tf1 = Matrix([tf([1, 2], [1, 2, 3])])
    ss1 = tf2ss(tf1)
    tf1_new = ss2tf(ss1)
    assert tf1 == tf1_new
