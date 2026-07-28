from __future__ import annotations
from sympy.physics.mechanics.method import MethodBase
from sympy.testing.pytest import raises

def test_method():
    raises(TypeError, lambda: MethodBase())
