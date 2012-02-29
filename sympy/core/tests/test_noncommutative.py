"""Tests for noncommutative symbols and expressions."""

from sympy import (
    conjugate,
    expand,
    factor,
    radsimp,
    simplify,
    symbols,
    I,
)
from sympy.abc import x, y, z

A, B, C = symbols("A B C", commutative=False)
X, Y = symbols("X Y", commutative=False, real=True)
Z = X + I*Y

def test_complex():
    assert Z.conjugate() == X - I*Y
    assert (Z*Z.conjugate()).expand() == X**2 + Y**2 + I*Y*X - I*X*Y

def test_conjugate():
    assert conjugate(A).is_commutative == False
    assert (A*A).conjugate() == conjugate(A)**2
    assert (A*B).conjugate() == conjugate(A)*conjugate(B)
    assert (A*B**2).conjugate() == conjugate(A)*conjugate(B)**2
    assert (A*B - B*A).conjugate() == conjugate(A)*conjugate(B) - conjugate(B)*conjugate(A)
    assert (A*B).conjugate() - (B*A).conjugate() == conjugate(A)*conjugate(B) - conjugate(B)*conjugate(A)

def test_expand():
    assert expand((A*B)**2) == A*B*A*B
    assert expand(A*B - B*A) == A*B - B*A
    assert expand((A*B/A)**2) == A*B*B/A
    assert expand(B*A*(A + B)*B) == B*A**2*B + B*A*B**2
    assert expand(B*A*(A + C)*B) == B*A**2*B + B*A*C*B

def test_factor():
    assert factor(A*B - B*A) == A*B - B*A

def test_radsimp():
    assert radsimp(A*B - B*A) == A*B - B*A

def test_simplify():
    assert simplify(A*B - B*A) == A*B - B*A

def test_subs():
    assert (x*y*A).subs(x*y, z) == A*z
    assert (x*A*B).subs(x*A, C) == C*B
    assert (x*A*x*x).subs(x**2*A, C) == x*C
    assert (x*A*x*B).subs(x**2*A, C) == C*B
    assert (A**2*B**2).subs(A*B**2, C) == A*C
    assert (A*A*A + A*B*A).subs(A*A*A, C) == C + A*B*A
