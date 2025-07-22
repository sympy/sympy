import random
import string
from collections import Counter

from sympy import Array
from sympy.external import import_module
from sympy.tensor.array.expressions import ArraySymbol, ArrayAdd, ArrayContraction, ArrayTensorProduct, PermuteDims, \
    ArrayDiagonal
from sympy.tensor.array.expressions.einsum_sympy import _convert_einsum_to_sympy_array, Einsum, convert_array_to_einsum
from sympy.testing.pytest import raises

numpy = import_module("numpy")


def _tolist(expr):
    if hasattr(expr, "tolist"):
        return expr.tolist()
    return expr


def _compare_einsum_numpy(path, *a):
    sa = [Array(i) for i in a]
    result_numpy = numpy.einsum(path, *a)
    expr_array = _convert_einsum_to_sympy_array(path, *sa)
    assert not isinstance(expr_array, Einsum)
    expr_einsum = Einsum(path, *sa)
    assert (result_numpy == numpy.array(_tolist(expr_array.as_explicit()))).all()
    assert (result_numpy == numpy.array(_tolist(expr_einsum.as_explicit()))).all()


def _get_random_path(*a):
    dim = sum([len(i.shape) for i in a])
    letters = string.ascii_letters[:dim]
    indices_src = numpy.random.choice(list(letters), size=dim, replace=True)
    indices_src_by_arg = []
    counter = 0
    for i in a:
        indices_src_by_arg.append(list(indices_src[counter:(counter+len(i.shape))]))
        counter += len(i.shape)
    indices_counted = Counter(indices_src)
    indices_sing = [k for k, v in indices_counted.items() if v == 1]
    indices_mult = [k for k, v in indices_counted.items() if v > 1]
    min0 = 0 if indices_sing else 1
    if random.random() > 0.5:
        indices_sing = [i for i in indices_sing if random.random() > 0.5]
    indices_dst = indices_sing + list(numpy.random.choice(indices_mult, size=numpy.random.randint(min0, len(indices_mult)), replace=False))
    random.shuffle(indices_dst)
    path = f"{','.join([''.join(i) for i in indices_src_by_arg])}->{''.join(indices_dst)}"
    return path


def test_einsum_in_sympy():
    Ra = ArraySymbol("Ra", (1, 2))
    Rb = ArraySymbol("Rb", (2, 3, 4))
    Rc = ArraySymbol("Rc", (2, 3, 5, 6, 7))
    Rd = ArraySymbol("Rd", (8, 9, 2, 4, 6, 10))

    raises(ValueError, lambda: Einsum("a->a", Ra))
    raises(ValueError, lambda: Einsum("abc->cab", Ra))

    ei = Einsum("ab,cde,fghij,klmnop->abcdefghijklmnop", Ra, Rb, Rc, Rd)
    assert ei.shape == (1, 2, 2, 3, 4, 2, 3, 5, 6, 7, 8, 9, 2, 4, 6, 10)
    ei = Einsum("ab,cde,fghij,klmnop->ponmlkjihgfedcba", Ra, Rb, Rc, Rd)
    assert ei.shape == (10, 6, 4, 2, 9, 8, 7, 6, 5, 3, 2, 4, 3, 2, 2, 1)
    ei = Einsum("ab,cde,fghij,klmnop->jlbophmicgdfekan", Ra, Rb, Rc, Rd)
    assert ei.shape == (7, 9, 2, 6, 10, 5, 2, 6, 2, 3, 3, 2, 4, 8, 1, 4)
    ei = Einsum("ab,cde,fghij,klmnop->bfgmkijlaoecdhpn", Ra, Rb, Rc, Rd)
    assert ei.shape == (2, 2, 3, 2, 8, 6, 7, 9, 1, 6, 4, 2, 3, 5, 10, 4)

    ei = Einsum("ab,cde,fghij,klmnop->cimdb", Ra, Rb, Rc, Rd)
    assert ei.shape == (2, 6, 2, 3, 2)

    A2 = ArraySymbol("A2", (2, 2))

    ei = Einsum("ij->ij", A2)
    assert ei.shape == (2, 2)
    assert ei.as_explicit() == A2.as_explicit()

    ei = Einsum("ij->", A2)
    assert ei.shape == ()
    assert ei.as_explicit() == A2[0, 0] + A2[0, 1] + A2[1, 0] + A2[1, 1]

    ei = Einsum("ij->j", A2)
    assert ei.shape == (2,)
    assert ei.as_explicit() == Array([A2[0, 0] + A2[1, 0], A2[0, 1] + A2[1, 1]])

    ei = Einsum("ii->", A2)
    assert ei.shape == ()
    assert ei.as_explicit() == A2[0, 0] + A2[1, 1]

    ei = Einsum("ii->i", A2)
    assert ei.shape == (2,)
    assert ei.as_explicit() == Array([A2[0, 0], A2[1, 1]])

    A = ArraySymbol("A", (3, 3))
    B = ArraySymbol("B", (3, 3))
    C = ArraySymbol("C", (3, 3))
    D = ArraySymbol("D", (3, 3))

    assert convert_array_to_einsum(_convert_einsum_to_sympy_array("ij,jk->ik", A, B)) == Einsum("ab,bc->ac", A, B)
    assert convert_array_to_einsum(_convert_einsum_to_sympy_array("ij,jk->ijk", A, B)) == Einsum("ab,bc->abc", A, B)
    assert convert_array_to_einsum(_convert_einsum_to_sympy_array("ij->ji", A)) == Einsum("ab->ba", A)
    assert convert_array_to_einsum(ArrayAdd(A, B)) == ArrayAdd(A, B)

    expr_array = ArrayContraction(ArrayTensorProduct(A, ArrayAdd(A, B)), (1, 2))
    expr_einsum = Einsum("ab,bc->ac", A, ArrayAdd(A, B))
    assert convert_array_to_einsum(expr_array) == expr_einsum
    assert expr_einsum.as_array_expression() == expr_array

    expr_array = ArrayContraction(ArrayTensorProduct(A, ArrayAdd(B, ArrayContraction(ArrayTensorProduct(C, D), (0, 2)))), (1, 2))
    expr_einsum = Einsum("ab,bc->ac", A, ArrayAdd(B, Einsum("ab,ac->bc", C, D)))
    assert convert_array_to_einsum(expr_array) == expr_einsum
    # .as_array_expression() does not apply recursively:
    assert expr_einsum.as_array_expression() == ArrayContraction(ArrayTensorProduct(A, ArrayAdd(B, Einsum("ab,ac->bc", C, D))), (1, 2))

    tp = ArrayTensorProduct(A, B, C, D)

    expr_array = PermuteDims(
        ArrayContraction(tp, (2, 5, 7)),
        index_order_old="abdef",
        index_order_new="fdeab")
    expr_einsum = Einsum('ab,cd,ec,fc->fdeab', A, B, C, D)
    assert convert_array_to_einsum(expr_array) == expr_einsum
    assert expr_einsum.as_array_expression() == expr_array

    expr_array = ArrayContraction(tp, (0,), (1, 4), (2, 5), (3,), (6,), (7,))
    expr_einsum = Einsum('ab,cd,bc,ef->', A, B, C, D)
    assert convert_array_to_einsum(expr_array) == expr_einsum
    assert expr_einsum.as_array_expression() == expr_array

    expr_array = PermuteDims(ArrayDiagonal(
            ArrayContraction(tp, (3, 6, 7), (5,)),
            (0, 1)
        ),
        index_order_old="bda",
        index_order_new="dab")
    expr_einsum = Einsum('aa,bc,de,cc->dab', A, B, C, D)
    assert convert_array_to_einsum(expr_array) == expr_einsum
    assert expr_einsum.as_array_expression() == expr_array

    expr_array = PermuteDims(
        ArrayContraction(tp, (3,), (4, 6, 7)),
        index_order_old="abcf",
        index_order_new="bcfa")
    expr_einsum = Einsum('ab,cd,ef,ee->bcfa', A, B, C, D)
    assert convert_array_to_einsum(expr_array) == expr_einsum
    assert expr_einsum.as_array_expression() == expr_array


def test_from_to_einsum_numpy():
    if numpy is None:
        return

    base_dim = 2

    a = numpy.random.randint(0, 100, base_dim**4).reshape(base_dim, base_dim, base_dim, base_dim)
    b = numpy.random.randint(0, 100, base_dim**4).reshape(base_dim, base_dim, base_dim, base_dim)
    c = numpy.random.randint(0, 100, base_dim**2).reshape(base_dim, base_dim)

    _compare_einsum_numpy("ijkl->jkil", a)
    _compare_einsum_numpy("iijk->kj", a)
    _compare_einsum_numpy("iijk->kji", a)

    for i in range(5):
        path = _get_random_path(a, b)
        _compare_einsum_numpy(path, a, b)

    for i in range(2):
        path = _get_random_path(c, c, c)
        _compare_einsum_numpy(path, c % 7, c % 11, c % 13)
