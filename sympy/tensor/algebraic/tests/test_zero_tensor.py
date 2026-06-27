from __future__ import annotations

from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol

from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor, algebraic_zero_tensor

# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
I3 = MatrixSymbol("I3", 3, 3)


def test_zero_tensor_construction():
    z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert z.shape == ((3, 4), (4, 5))
    assert z.tensor_shape == ((3, 4), (4, 5))
    assert not z
    assert z == AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert z != AlgebraicZeroTensor(((3, 5),))
    assert -z is z


def test_zero_tensor_bare_pair_wrapping():
    z = AlgebraicZeroTensor((3, 4))
    assert z.shape == ((3, 4),)


def test_zero_tensor_hash():
    z1 = AlgebraicZeroTensor(((3, 4),))
    z2 = AlgebraicZeroTensor(((3, 4),))
    assert hash(z1) == hash(z2)
    assert z1 == z2
    assert {z1} == {z2}


def test_zero_tensor_repr_str():
    z = AlgebraicZeroTensor(((3, 4),))
    assert repr(z) == "AlgebraicZeroTensor((3, 4),)"
    assert str(z) == "0_((3, 4),)"


def test_zero_tensor_add():
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    z = AlgebraicZeroTensor(((3, 4), (3, 3)))
    pt = AlgebraicPureTensor(A, I3)
    result = z + pt
    assert isinstance(result, (AlgebraicPureTensor, AlgebraicTensor))


def test_zero_tensor_radd():
    z = AlgebraicZeroTensor(((3, 4), (3, 3)))
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    pt = AlgebraicPureTensor(A, I3)
    result = pt + z
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_sub():
    z = AlgebraicZeroTensor(((3, 4), (3, 3)))
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    pt = AlgebraicPureTensor(A, I3)
    result = z - pt
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_rsub():
    z = AlgebraicZeroTensor(((3, 4), (3, 3)))
    from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    pt = AlgebraicPureTensor(A, I3)
    result = pt - z
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_convenience():
    z = algebraic_zero_tensor(((3, 4), (4, 5)))
    assert isinstance(z, AlgebraicZeroTensor)
    assert z.shape == ((3, 4), (4, 5))
