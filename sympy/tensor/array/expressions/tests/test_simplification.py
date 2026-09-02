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
        PermuteDims(ArrayTensorProduct(2*I, I), Permutation(3)(1, 2))


def test_commutativity_simplify():
    from sympy import Symbol
    from sympy.tensor.array.expressions import commutativity_simplify

    z, w = symbols("z w")
    u = Symbol("u", commutative=False)
    v = Symbol("v", commutative=False)

    E00 = Array([[1, 0], [0, 0]])
    E11 = Array([[0, 0], [0, 1]])
    X = Array([[u, 0], [0, v]])
    zE00 = Array([[z, 0], [0, 0]])

    # The scalar z can be moved between commutative factors; the two
    # addends differ in two slots, so collect_tensor_products alone
    # cannot merge them, but the basis decomposition can:
    expr = ArrayAdd(ArrayTensorProduct(zE00, E00, X),
                    ArrayTensorProduct(E00, zE00, X))
    collected = collect_tensor_products(expr)
    assert isinstance(collected, ArrayAdd)
    assert commutativity_simplify(expr) == ArrayTensorProduct(2*z, E00, E00, X)

    # Cancellation of terms with redistributed entries:
    expr = ArrayAdd(ArrayTensorProduct(zE00, E00, X),
                    ArrayTensorProduct(-1, E00, zE00, X))
    assert commutativity_simplify(expr) == ZeroArray(2, 2, 2, 2, 2, 2)

    # Decomposing a full matrix regroups its entries; the result is
    # value-equivalent to the input:
    M = Array([[z, w], [w, z]])
    expr = ArrayAdd(ArrayTensorProduct(M, X), ArrayTensorProduct(zE00, X))
    simplified = commutativity_simplify(expr)
    assert simplified.as_explicit() == expr.as_explicit()

    # Non-decomposable sums fall back to plain collection:
    expr = ArrayAdd(ArrayTensorProduct(Ra, Rb), ArrayTensorProduct(Rc, Rb))
    assert commutativity_simplify(expr) == ArrayTensorProduct(ArrayAdd(Ra, Rc), Rb)

    # Terms whose noncommutative factors coincide merge even when every
    # commutative factor differs, when the basis components match:
    Y = Array([[u, v], [v, u]])
    expr = ArrayAdd(ArrayTensorProduct(2*E00, Y), ArrayTensorProduct(E00, Y),
                    ArrayTensorProduct(E11, Y))
    simplified = commutativity_simplify(expr)
    assert simplified.as_explicit() == expr.as_explicit()


def test_commutativity_simplify_commutator_example():
    # A miniature version of the Dirac-commutator computation of
    # PR #30131: expand [D, a] = D.a - a.D as a sum of tensor products
    # of explicit matrices with noncommutative entries, then verify that
    # commutativity_simplify preserves the value while reducing the
    # number of addends:
    from sympy import Symbol
    from sympy.tensor.array.expressions import commutativity_simplify

    z, w = symbols("z w")
    ups = Symbol("Upsilon", commutative=False)
    upsc = Symbol("Upsilon_c", commutative=False)

    pL = Array([[1, 0], [0, 0]])
    pR = Array([[0, 0], [0, 1]])
    m = Array([[0, ups], [upsc, 0]])
    a1 = Array([[z, 0], [0, w]])
    a2 = Array([[w, 0], [0, z]])
    id2 = Array([[1, 0], [0, 1]])

    # Terms of the expanded commutator (schematically).  With this
    # symmetric choice of algebra elements the commutator vanishes
    # identically, which the basis decomposition detects:
    terms = [
        ArrayTensorProduct(pL, m, a1),
        ArrayTensorProduct(pR, m, a2),
        ArrayTensorProduct(-1, a1, m, pL),
        ArrayTensorProduct(-1, a2, m, pR),
        ArrayTensorProduct(z, pL, m, id2),
        ArrayTensorProduct(-z, pL, m, id2),
    ]
    expr = ArrayAdd(*terms)
    assert expr.as_explicit() == ZeroArray(2, 2, 2, 2, 2, 2).as_explicit()
    assert commutativity_simplify(expr) == ZeroArray(2, 2, 2, 2, 2, 2)

    # An asymmetric variant does not vanish; the six addends collapse to
    # two basis terms:
    a1b = Array([[z, 0], [0, 3*w]])
    terms = [
        ArrayTensorProduct(pL, m, a1b),
        ArrayTensorProduct(pR, m, a2),
        ArrayTensorProduct(-1, a1b, m, pL),
        ArrayTensorProduct(-1, a2, m, pR),
        ArrayTensorProduct(z, pL, m, id2),
        ArrayTensorProduct(-z, pL, m, id2),
    ]
    expr = ArrayAdd(*terms)
    simplified = commutativity_simplify(expr)
    assert simplified.as_explicit() == expr.as_explicit()
    assert isinstance(simplified, ArrayAdd)
    assert len(simplified.args) == 2
