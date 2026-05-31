from __future__ import annotations
from sympy.physics.mechanics.method import Method
from sympy.testing.pytest import raises

def test_method():
    raises(TypeError, lambda: Method())
