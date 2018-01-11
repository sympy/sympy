"""Tests for Finite Fields."""

from sympy.polys.domains.finitefield import FiniteField

from sympy.utilities.pytest import raises

def test_finitefield():
    raises(ValueError, lambda: FiniteField(42))
