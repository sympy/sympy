from __future__ import annotations

from sympy import Array, MatrixSymbol, symbols
from sympy.matrices.expressions.special import Identity
from sympy.tensor.array.expressions import ArrayAdd, ArrayContraction, ArrayTensorProduct, \
    ArraySymbol, PermuteDims, collect_tensor_products
from sympy.tensor.array.expressions.array_expressions import ZeroArray

x, y, k = symbols("x y k")

A = MatrixSymbol("A", 2, 2)
B = MatrixSymbol("B", 2, 2)
C = MatrixSymbol("C", 2, 2)
D = MatrixSymbol("D", 2, 2)

Ra = ArraySymbol("Ra", (2, 3))
Rb = ArraySymbol("Rb", (3, 4))
Rc = ArraySymbol("Rc", (2, 3))


def test_array_add_collect_coefficients():
    # Automatic collection of proportional terms in ArrayAdd.doit():
    T = ArrayTensorProduct(Ra, Rb)
    assert ArrayAdd(T, T).doit() == ArrayTensorProduct(2, Ra, Rb)
    assert ArrayAdd(ArrayTensorProduct(2, Ra, Rb),
                    ArrayTensorProduct(3, Ra, Rb)).doit() == ArrayTensorProduct(5, Ra, Rb)
    assert ArrayAdd(T, ArrayTensorProduct(-1, Ra, Rb)).doit() == ZeroArray(2, 3, 3, 4)
    assert ArrayAdd(ArrayTensorProduct(x, Ra, Rb),
                    ArrayTensorProduct(y, Ra, Rb)).doit() == ArrayTensorProduct(x + y, Ra, Rb)

    # Scalar coefficients inside MatrixExpr factors are also recognized:
    assert ArrayAdd(ArrayTensorProduct(2*A, B),
                    ArrayTensorProduct(3, A, B)).doit() == ArrayTensorProduct(5, A, B)

    # Single matrix leaves merge to a MatMul:
    assert ArrayAdd(A, A).doit() == 2*A

    # Terms that are not proportional are not merged:
    e = ArrayAdd(ArrayTensorProduct(Ra, Rb), ArrayTensorProduct(Rc, Rb))
    assert e.doit() == e

    # A Mul argument containing an array factor is normalized:
    assert ArrayAdd(2*Ra, Rc).doit() == ArrayAdd(ArrayTensorProduct(2, Ra), Rc)

    # Sums of scalars remain ArrayAdd objects:
    M3 = MatrixSymbol("M3", 3, 3)
    e = ArrayAdd(M3[0, 0], M3[1, 1])
    assert e.doit() == e


def test_collect_tensor_products():
    # Single differing slot merges by linearity:
    assert collect_tensor_products(
        ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(A, C))) == \
        ArrayTensorProduct(A, B + C)
    assert collect_tensor_products(
        ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(C, B))) == \
        ArrayTensorProduct(A + C, B)

    # Scalar coefficients are absorbed into the merged slot:
    assert collect_tensor_products(
        ArrayAdd(ArrayTensorProduct(2, A, B), ArrayTensorProduct(3, C, B))) == \
        ArrayTensorProduct(2*A + 3*C, B)

    # More than one differing slot: no merge
    e = ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(C, D))
    assert collect_tensor_products(e) == e

    # A merge can enable further merges:
    e = ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(A, C),
                 ArrayTensorProduct(C, B + C))
    assert collect_tensor_products(e) == ArrayTensorProduct(A + C, B + C)

    # Three factors, middle slot:
    e = ArrayAdd(ArrayTensorProduct(A, B, C), ArrayTensorProduct(A, D, C))
    assert collect_tensor_products(e) == ArrayTensorProduct(A, B + D, C)

    # Applies recursively inside other array operators:
    e = ArrayContraction(ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(A, C)), (1, 2))
    assert collect_tensor_products(e) == ArrayContraction(ArrayTensorProduct(A, B + C), (1, 2))

    # Non-matrix array factors:
    e = ArrayAdd(ArrayTensorProduct(Ra, Rb), ArrayTensorProduct(Rc, Rb))
    assert collect_tensor_products(e) == ArrayTensorProduct(ArrayAdd(Ra, Rc), Rb)

    # Non-ArrayAdd expressions pass through unchanged:
    assert collect_tensor_products(ArrayTensorProduct(A, B)) == ArrayTensorProduct(A, B)
    assert collect_tensor_products(x + y) == x + y


def test_collect_tensor_products_numeric():
    a = Array([[1, 2], [3, 4]])
    b = Array([[0, 1], [1, 0]])
    c = Array([[2, 0], [0, 2]])

    for e in [
        ArrayAdd(ArrayTensorProduct(a, b), ArrayTensorProduct(a, c)),
        ArrayAdd(ArrayTensorProduct(2, a, b), ArrayTensorProduct(3, c, b)),
        ArrayAdd(ArrayTensorProduct(a, b, c), ArrayTensorProduct(a, c, c)),
    ]:
        collected = collect_tensor_products(e)
        assert collected.as_explicit() == e.as_explicit()
        # something was actually collected:
        assert not isinstance(collected, ArrayAdd)


def test_array_add_collect_derivative_regression():
    # Collection makes derivatives of sums with repeated addends compact:
    from sympy.combinatorics import Permutation
    from sympy.matrices.expressions.matadd import MatAdd
    Ak = MatrixSymbol("Ak", k, k)
    I = Identity(k)
    assert MatAdd(Ak, Ak).diff(Ak) == \
        PermuteDims(ArrayTensorProduct(2, I, I), Permutation(3)(1, 2))
