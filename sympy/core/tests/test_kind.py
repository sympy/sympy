from sympy.core.add import Add
from sympy.core.kind import NumberKind, UndefinedKind
from sympy.core.numbers import pi, zoo, I, AlgebraicNumber
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.integrals.integrals import Integral
from sympy.matrices import (Matrix, SparseMatrix, ImmutableMatrix,
    ImmutableSparseMatrix, MatrixSymbol, MatrixKind)

comm_x = Symbol('x')
noncomm_x = Symbol('x', commutative=False)

def test_NumberKind():
    assert S.One.kind is NumberKind
    assert pi.kind is NumberKind
    assert S.NaN.kind is NumberKind
    assert zoo.kind is NumberKind
    assert I.kind is NumberKind
    assert AlgebraicNumber(1).kind is NumberKind

def test_Add_kind():
    assert Add(2, 3, evaluate=False).kind is NumberKind
    assert Add(2,comm_x).kind is NumberKind
    assert Add(2,noncomm_x).kind is UndefinedKind

def test_Symbol_kind():
    assert comm_x.kind is NumberKind
    assert noncomm_x.kind is UndefinedKind

def test_Integral_kind():
    A = MatrixSymbol('A', 2,2)
    assert Integral(comm_x, comm_x).kind is NumberKind
    assert Integral(A, comm_x).kind is MatrixKind(NumberKind)

def test_Matrix_kind():
    classes = (Matrix, SparseMatrix, ImmutableMatrix, ImmutableSparseMatrix)
    for cls in classes:
        m = cls.zeros(3, 2)
        assert m.kind is MatrixKind(NumberKind)
