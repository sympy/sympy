from __future__ import annotations

from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol

from sympy.tensor.algebraic.zero_tensor import ZeroTensor, zero_tensor

# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
I3 = MatrixSymbol("I3", 3, 3)


def test_zero_tensor_construction():
    z = ZeroTensor(((3, 4), (4, 5)))
    assert z.shape == ((3, 4), (4, 5))
    assert z.tensor_shape == ((3, 4), (4, 5))
    assert not z
    assert z == ZeroTensor(((3, 4), (4, 5)))
    assert z != ZeroTensor(((3, 5),))
    assert -z is z


def test_zero_tensor_bare_pair_wrapping():
    z = ZeroTensor((3, 4))
    assert z.shape == ((3, 4),)


def test_zero_tensor_hash():
    z1 = ZeroTensor(((3, 4),))
    z2 = ZeroTensor(((3, 4),))
    assert hash(z1) == hash(z2)
    assert z1 == z2
    assert {z1} == {z2}


def test_zero_tensor_repr_str():
    z = ZeroTensor(((3, 4),))
    assert repr(z) == "ZeroTensor((3, 4),)"
    assert str(z) == "0_((3, 4),)"


def test_zero_tensor_add():
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = z + pt
    assert isinstance(result, (PureTensor, AlgebraicTensor))


def test_zero_tensor_radd():
    z = ZeroTensor(((3, 4), (3, 3)))
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    pt = PureTensor(A, I3)
    result = pt + z
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_sub():
    z = ZeroTensor(((3, 4), (3, 3)))
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    pt = PureTensor(A, I3)
    result = z - pt
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_rsub():
    z = ZeroTensor(((3, 4), (3, 3)))
    from sympy.tensor.algebraic.pure_tensor import PureTensor
    from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
    pt = PureTensor(A, I3)
    result = pt - z
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_convenience():
    z = zero_tensor(((3, 4), (4, 5)))
    assert isinstance(z, ZeroTensor)
    assert z.shape == ((3, 4), (4, 5))
