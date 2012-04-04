"""Tests for noncommutative symbols and expressions."""

from sympy import (
    cancel,
    collect,
    combsimp,
    conjugate,
    expand,
    factor,
    posify,
    radsimp,
    ratsimp,
    rcollect,
    simplify,
    symbols,
    trigsimp,
    I,
)
from sympy.abc import x, y, z
from sympy.utilities.pytest import XFAIL

A, B, C = symbols("A B C", commutative=False)

@XFAIL
def test_cancel():
    assert cancel(A*B - B*A) == A*B - B*A

@XFAIL
def test_collect():
    assert collect(A*B - B*A, A) == A*B - B*A
    assert collect(A*B - B*A, B) == A*B - B*A
    assert collect(A*B - B*A, x) == A*B - B*A

def test_combsimp():
    assert combsimp(A*B - B*A) == A*B - B*A

def test_conjugate():
    assert conjugate(A).is_commutative == False
    assert (A*A).conjugate() == conjugate(A)**2
    assert (A*B).conjugate() == conjugate(A)*conjugate(B)
    assert (A*B**2).conjugate() == conjugate(A)*conjugate(B)**2
    assert (A*B - B*A).conjugate() == conjugate(A)*conjugate(B) - conjugate(B)*conjugate(A)
    assert (A*B).conjugate() - (B*A).conjugate() == conjugate(A)*conjugate(B) - conjugate(B)*conjugate(A)
    assert (A + I*B).conjugate() == conjugate(A) - I*conjugate(B)

def test_expand():
    assert expand((A*B)**2) == A*B*A*B
    assert expand(A*B - B*A) == A*B - B*A
    assert expand((A*B/A)**2) == A*B*B/A
    assert expand(B*A*(A + B)*B) == B*A**2*B + B*A*B**2
    assert expand(B*A*(A + C)*B) == B*A**2*B + B*A*C*B

def test_factor():
    assert factor(A*B - B*A) == A*B - B*A

def test_posify():
    assert posify(A)[0].is_commutative == False
    for q in (A*B/A, (A*B/A)**2, (A*B)**2, A*B - B*A):
        p = posify(q)
        assert p[0].subs(p[1]) == q

def test_radsimp():
    assert radsimp(A*B - B*A) == A*B - B*A

@XFAIL
def test_ratsimp():
    assert ratsimp(A*B - B*A) == A*B - B*A

@XFAIL
def test_rcollect():
    assert rcollect(A*B - B*A, A) == A*B - B*A
    assert rcollect(A*B - B*A, B) == A*B - B*A
    assert rcollect(A*B - B*A, x) == A*B - B*A

def test_simplify():
    assert simplify(A*B - B*A) == A*B - B*A

def test_subs():
    assert (x*y*A).subs(x*y, z) == A*z
    assert (x*A*B).subs(x*A, C) == C*B
    assert (x*A*x*x).subs(x**2*A, C) == x*C
    assert (x*A*x*B).subs(x**2*A, C) == C*B
    assert (A**2*B**2).subs(A*B**2, C) == A*C
    assert (A*A*A + A*B*A).subs(A*A*A, C) == C + A*B*A

@XFAIL
def test_trigsimp():
    assert trigsimp(A*sin(x)**2 + A*cos(x)**2) == A
