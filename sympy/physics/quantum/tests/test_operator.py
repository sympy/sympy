# -*- encoding: utf-8 -*-

from sympy import Symbol, Integer, Mul, Pow, pretty, latex, srepr

from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.hilbert import HilbertSpace
from sympy.physics.quantum.operator import (
    Operator, UnitaryOperator, HermitianOperator, OuterProduct
)
from sympy.physics.quantum.state import Ket, Bra


def test_operator():
    A = Operator('A')
    B = Operator('B')
    C = Operator('C')

    assert isinstance(A, Operator)
    assert isinstance(A, QExpr)

    assert A.label == (Symbol('A'),)
    assert A.is_commutative == False
    assert A.hilbert_space == HilbertSpace()

    assert A*B != B*A

    assert (A*(B+C)).expand() == A*B + A*C
    assert ((A+B)**2).expand() == A**2 + A*B + B*A + B**2

    assert pretty(A) == u'A'
    assert latex(A) == 'A'
    assert eval(srepr(A)) == A

    D = Operator('D', Symbol('t'), '1/2')
    assert pretty(D) == u'Operator(D,t,1/2)'
    assert latex(D) == 'Operator(D,t,1/2)'
    assert eval(srepr(D)) == D

def test_operator_inv():
    A = Operator('A')
    assert A*A.inv() == 1
    assert A.inv()*A == 1
    assert pretty(A.inv()) == u'1\n─\nA'
    assert latex(A.inv()) == r'\frac{1}{A}'
    assert eval(srepr(A.inv())) == A.inv()

def test_hermitian():
    H = HermitianOperator('H')

    assert isinstance(H, HermitianOperator)
    assert isinstance(H, Operator)

    assert Dagger(H) == H
    assert H.inv() != H
    assert H.is_commutative == False
    assert Dagger(H).is_commutative == False

def test_unitary():
    U = UnitaryOperator('U')

    assert isinstance(U, UnitaryOperator)
    assert isinstance(U, Operator)

    assert U.inv() == Dagger(U)
    assert U*Dagger(U) == 1
    assert Dagger(U)*U == 1
    assert U.is_commutative == False
    assert Dagger(U).is_commutative == False


def test_outer_product():
    k = Ket('k')
    b = Bra('b')
    op = OuterProduct(k, b)

    assert isinstance(op, OuterProduct)
    assert isinstance(op, Operator)

    assert op.ket == k
    assert op.bra == b
    assert op.label == (k, b)
    assert op.is_commutative == False

    op = k*b

    assert isinstance(op, OuterProduct)
    assert isinstance(op, Operator)

    assert op.ket == k
    assert op.bra == b
    assert op.label == (k, b)
    assert op.is_commutative == False

    op = 2*k*b

    assert op == Mul(Integer(2), k, b)

    op = 2*(k*b)

    assert op == Mul(Integer(2), OuterProduct(k, b))

    assert Dagger(k*b) == OuterProduct(Dagger(b),Dagger(k))
    assert Dagger(k*b).is_commutative == False

    op = k*b
    assert pretty(op) == u'\u2758k\u27e9\u27e8b\u2758'
    assert latex(op) == r'{\left| k \right\rangle }{\left\langle b \right| }'
    assert eval(srepr(op)) == op

def test_operator_dagger():
    A = Operator('A')
    B = Operator('B')
    assert Dagger(A*B) == Dagger(B)*Dagger(A)
    assert Dagger(A+B) == Dagger(A) + Dagger(B)
    assert Dagger(A**2) == Dagger(A)**2

    assert pretty(Dagger(A)) == u" †\nA "
    assert latex(Dagger(A)) == r'A^{\dag}'
    assert eval(srepr(Dagger(A))) == Dagger(A)
