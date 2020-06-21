"""Tests for Finite Fields."""

from sympy.polys.domains.finitefield import FiniteField
from sympy.testing.pytest import raises

def test_FiniteField():
    raises(ValueError, lambda: FiniteField(42))
