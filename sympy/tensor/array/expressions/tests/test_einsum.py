from __future__ import annotations

import random
import string
from collections import Counter

from sympy import Array, Integer, Symbol
from sympy.core.symbol import Str
from sympy.external import import_module
from sympy.printing.str import sstr
from sympy.tensor.array.expressions import ArraySymbol, ArrayAdd, ArrayContraction, ArrayTensorProduct, \
    PermuteDims, ArrayDiagonal, Reshape, Einsum, convert_array_to_einsum, convert_einsum_to_array
from sympy.testing.pytest import raises

numpy = import_module("numpy")


def _tolist(expr):
    if hasattr(expr, "tolist"):
        return expr.tolist()
    return expr


def _compare_einsum_numpy(path, *a):
    sa = [Array(i) for i in a]
    result_numpy = numpy.einsum(path, *a)
    expr_array = convert_einsum_to_array(path, *sa)
    assert not isinstance(expr_array, Einsum)
    expr_einsum = Einsum(path, *sa)
    assert tuple(result_numpy.shape) == tuple(int(i) for i in expr_einsum.shape)
    assert (result_numpy == numpy.array(_tolist(expr_array.as_explicit()))).all()
    assert (result_numpy == numpy.array(_tolist(expr_einsum.as_explicit()))).all()


def _get_random_path(*a):
    # Uses the ``random`` module only, which the test runner seeds, so that
    # a failure can be reproduced from the reported seed:
    dim = sum([len(i.shape) for i in a])
    letters = string.ascii_letters[:dim]
    indices_src = [random.choice(letters) for _ in range(dim)]
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
    indices_dst = indices_sing + random.sample(indices_mult, random.randint(min0, len(indices_mult)))
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

    assert convert_array_to_einsum(convert_einsum_to_array("ij,jk->ik", A, B)) == Einsum("ab,bc->ac", A, B)
    assert convert_array_to_einsum(convert_einsum_to_array("ij,jk->ijk", A, B)) == Einsum("ab,bc->abc", A, B)
    assert convert_array_to_einsum(convert_einsum_to_array("ij->ji", A)) == Einsum("ab->ba", A)
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


def test_einsum_basic_properties():
    A = ArraySymbol("A", (2, 3))
    B = ArraySymbol("B", (3, 4))

    ei = Einsum("ij,jk->ik", A, B)
    assert ei.path_string == "ij,jk->ik"
    assert ei.path_src == (("i", "j"), ("j", "k"))
    assert ei.path_dst == ("i", "k")
    assert ei.array_args == (A, B)
    assert ei.shape == (2, 4)

    # Rebuildable from args:
    assert ei.func(*ei.args) == ei
    # The subscript string may also be given as a Str object:
    assert Einsum(Str("ij,jk->ik"), A, B) == ei

    # Printing quotes the subscript string:
    assert sstr(ei) == "Einsum('ij,jk->ik', A, B)"

    # Symbolic shapes are supported:
    k, l, m = Symbol("k"), Symbol("l"), Symbol("m")
    Sa = ArraySymbol("Sa", (k, l))
    Sb = ArraySymbol("Sb", (l, m))
    assert Einsum("ij,jk->ik", Sa, Sb).shape == (k, m)
    # Mismatched symbolic dimensions are rejected:
    Sc = ArraySymbol("Sc", (m, m))
    raises(ValueError, lambda: Einsum("ij,jk->ik", Sa, Sc))


def test_einsum_argument_coercion():
    # Nested lists are converted to arrays:
    ei = Einsum("ij->ji", [[1, 2], [3, 4]])
    assert ei.array_args[0] == Array([[1, 2], [3, 4]])
    assert ei.as_explicit() == Array([[1, 3], [2, 4]])

    # Matrices are accepted as operands:
    from sympy import Matrix
    ei = Einsum("ij,jk->ik", Matrix([[1, 2], [3, 4]]), Matrix([[5, 6], [7, 8]]))
    assert ei.as_explicit() == Array(Matrix([[1, 2], [3, 4]]) * Matrix([[5, 6], [7, 8]]))

    if numpy:
        # NumPy arrays are converted to arrays:
        ei = Einsum("ij->ji", numpy.array([[1, 2], [3, 4]]))
        assert ei.array_args[0] == Array([[1, 2], [3, 4]])


def test_einsum_implicit_mode():
    A = ArraySymbol("A", (2, 3))
    B = ArraySymbol("B", (3, 4))
    V = ArraySymbol("V", (2, 3, 4))

    # Implicit mode: output indices are the ones appearing exactly once,
    # sorted alphabetically:
    assert Einsum("ij,jk", A, B) == Einsum("ij,jk->ik", A, B)
    assert Einsum("ji", A) == Einsum("ji->ij", A)
    assert Einsum("ji", A).shape == (3, 2)
    assert Einsum("kji", V) == Einsum("kji->ijk", V)
    A2 = ArraySymbol("A2", (2, 2))
    assert Einsum("ii", A2) == Einsum("ii->", A2)
    # Uppercase letters sort before lowercase ones (ASCII order, as in NumPy):
    assert Einsum("bA", A).shape == (3, 2)


def test_einsum_ellipsis():
    D = ArraySymbol("D", (5, 2, 3))
    E = ArraySymbol("E", (5, 3, 4))
    V = ArraySymbol("V", (2, 3, 4))

    # Batched matrix multiplication:
    ei = Einsum("...ij,...jk->...ik", D, E)
    assert ei.shape == (5, 2, 4)
    assert ei == Einsum("aij,ajk->aik", D, E)

    # Ellipsis dimensions come first in the implicit-mode output:
    assert Einsum("k...i", V).shape == (3, 4, 2)
    assert Einsum("k...i", V) == Einsum("kai->aik", V)

    # Sum over a named index, keeping the ellipsis dimensions:
    assert Einsum("...i->...", V).shape == (2, 3)
    assert Einsum("i...->...", V).shape == (3, 4)
    assert Einsum("...i->...i", V) == Einsum("abi->abi", V)

    # Diagonal across the ellipsis boundary:
    W = ArraySymbol("W", (2, 3, 2))
    assert Einsum("a...a->...a", W).shape == (3, 2)

    # An ellipsis matching zero dimensions is allowed:
    A = ArraySymbol("A", (2, 3))
    assert Einsum("...ij->...ji", A) == Einsum("ij->ji", A)

    # Broadcasting of ellipsis dimensions is not supported:
    u = ArraySymbol("u", (4,))
    raises(NotImplementedError, lambda: Einsum("...i,...i->...", V, u))


def test_einsum_scalar_operands():
    x = Symbol("x")
    A = ArraySymbol("A", (2, 3))

    ei = Einsum(",ab->ab", x, A)
    assert ei.shape == (2, 3)
    ei = Einsum("ab,->ab", A, x)
    assert ei.shape == (2, 3)
    ei = Einsum(",,->", x, Integer(3), Integer(5))
    assert ei.shape == ()
    assert Einsum(",ab->ab", 2, Array([[1, 2], [3, 4]])).as_explicit() == Array([[2, 4], [6, 8]])


def test_einsum_validation():
    A = ArraySymbol("A", (2, 3))
    B = ArraySymbol("B", (5, 4))
    V = ArraySymbol("V", (2, 3, 4))

    raises(ValueError, lambda: Einsum())
    raises(ValueError, lambda: Einsum(A))
    # Wrong number of operands:
    raises(ValueError, lambda: Einsum("ij,jk->ik", A))
    # Wrong number of indices:
    raises(ValueError, lambda: Einsum("ij", V))
    raises(ValueError, lambda: Einsum("ijk", A))
    # Invalid characters and malformed ellipses:
    raises(ValueError, lambda: Einsum("i j+k", A))
    raises(ValueError, lambda: Einsum("i..", A))
    raises(ValueError, lambda: Einsum(".i...", V))
    raises(ValueError, lambda: Einsum("ij->j.", A))
    # More than one '->':
    raises(ValueError, lambda: Einsum("ij->j->i", A))
    # Output indices must appear in the input:
    raises(ValueError, lambda: Einsum("ij->ik", A))
    # Output indices must not repeat:
    raises(ValueError, lambda: Einsum("ij->iij", A))
    # Missing output ellipsis:
    raises(ValueError, lambda: Einsum("...i->i", V))
    # Inconsistent dimensions for the same index:
    raises(ValueError, lambda: Einsum("ij,jk->ik", A, B))
    raises(ValueError, lambda: Einsum("ii->i", A))


def test_convert_array_to_einsum_positions():
    # Contraction and diagonal indices are axis positions of the nested
    # expression: converting them to einsum index letters has to take the
    # nested expression's axis order into account (this used to be wrong
    # when the nested expression was itself a PermuteDims or ArrayDiagonal):
    A = ArraySymbol("A", (2, 2))
    B = ArraySymbol("B", (2, 2))
    tp = ArrayTensorProduct(A, B)

    expr = ArrayContraction(PermuteDims(tp, [1, 2, 3, 0]), (0, 1))
    ei = convert_array_to_einsum(expr)
    assert ei == Einsum("ab,bc->ca", A, B)
    assert ei.as_explicit() == expr.as_explicit()

    expr = ArrayDiagonal(PermuteDims(tp, [1, 2, 3, 0]), (0, 3))
    ei = convert_array_to_einsum(expr)
    assert ei.as_explicit() == expr.as_explicit()

    expr = ArrayContraction(ArrayDiagonal(tp, (0, 3)), (0, 1))
    ei = convert_array_to_einsum(expr)
    assert ei.as_explicit() == expr.as_explicit()

    expr = ArrayDiagonal(ArrayDiagonal(tp, (0, 3)), (0, 1))
    ei = convert_array_to_einsum(expr)
    assert ei.as_explicit() == expr.as_explicit()


def test_convert_array_to_einsum_leaves():
    A = ArraySymbol("A", (2, 2))

    # Explicit arrays are valid leaves:
    a = Array([[1, 2], [3, 4]])
    ei = convert_array_to_einsum(ArrayContraction(a, (0, 1)))
    assert ei == Einsum("aa->", a)
    assert ei.as_explicit() == 5

    # Matrix expressions are valid leaves:
    from sympy import MatrixSymbol
    M = MatrixSymbol("M", 2, 2)
    assert convert_array_to_einsum(ArrayContraction(M, (0, 1))) == Einsum("aa->", M)

    # Objects that cannot be expressed in the einsum syntax are kept as
    # opaque operands:
    r = Reshape(A, (4,))
    assert convert_array_to_einsum(r) == r
    assert convert_array_to_einsum(ArrayTensorProduct(r, r)) == Einsum("a,b->ab", r, r)

    # An Einsum operand is kept as-is:
    ei = Einsum("ab->ba", A)
    assert convert_array_to_einsum(ei) == ei


def test_einsum_doit():
    A = ArraySymbol("A", (3, 3))
    B = ArraySymbol("B", (3, 3))

    expr = Einsum("ij,jk->ik", A, Einsum("ab->ba", B))
    assert expr.doit() == ArrayContraction(ArrayTensorProduct(A, B), (1, 3))
    # doit(deep=False) does not convert the nested Einsum:
    assert expr.doit(deep=False) == ArrayContraction(
        ArrayTensorProduct(A, Einsum("ab->ba", B)), (1, 2))


# The following subscript strings are taken from the test suite of the
# opt_einsum package (opt_einsum/tests/test_contract.py), BSD-3-Clause,
# https://github.com/dgasmith/opt_einsum
_opt_einsum_tests = [
    # Test scalar-like operations
    "a,->a",
    "ab,->ab",
    ",ab,->ab",
    ",,->",
    # Test hadamard-like products
    "a,ab,abc->abc",
    "a,b,ab->ab",
    # Test index-transformations
    "ea,fb,gc,hd,abcd->efgh",
    "ea,fb,abcd,gc,hd->efgh",
    "abcd,ea,fb,gc,hd->efgh",
    # Test complex contractions
    "acdf,jbje,gihb,hfac,gfac,gifabc,hfac",
    "cd,bdhe,aidb,hgca,gc,hgibcd,hgac",
    "abhe,hidj,jgba,hiab,gab",
    "bde,cdh,agdb,hica,ibd,hgicd,hiac",
    "chd,bde,agbc,hiad,hgc,hgi,hiad",
    "chd,bde,agbc,hiad,bdi,cgh,agdb",
    "bdhe,acad,hiab,agac,hibd",
    # Test collapse
    "ab,ab,c->",
    "ab,ab,c->c",
    "ab,ab,cd,cd->",
    "ab,ab,cd,cd->ac",
    "ab,ab,cd,cd->cd",
    "ab,ab,cd,cd,ef,ef->",
    # Test outer products
    "ab,cd,ef->abcdef",
    "ab,cd,ef->acdf",
    "ab,cd,de->abcde",
    "ab,cd,de->be",
    "ab,bcd,cd->abcd",
    "ab,bcd,cd->abd",
    # Random test cases that have previously failed
    "eb,cb,fb->cef",
    "dd,fb,be,cdb->cef",
    "bca,cdb,dbf,afc->",
    "dcc,fce,ea,dbf->ab",
    "fdf,cdd,ccd,afe->ae",
    "abcd,ad",
    "ed,fcd,ff,bcf->be",
    "baa,dcf,af,cde->be",
    "bd,db,eac->ace",
    "fff,fae,bef,def->abd",
    "efc,dbc,acf,fd->abe",
    # Inner products
    "ab,ab",
    "ab,ba",
    "abc,abc",
    "abc,bac",
    "abc,cba",
    # GEMM test cases
    "ab,bc",
    "ab,cb",
    "ba,bc",
    "ba,cb",
    "abcd,cd",
    "abcd,ab",
    "abcd,cdef",
    "abcd,cdef->feba",
    "abcd,efdc",
    # Inner than dot
    "aab,bc->ac",
    "ab,bcc->ac",
    "aab,bcc->ac",
    "baa,bcc->ac",
    "aab,ccb->ac",
    # Randomly built test cases
    "aab,fa,df,ecc->bde",
    "ecb,fef,bad,ed->ac",
    "bcf,bbb,fbf,fc->",
    "bb,ff,be->e",
    "bcb,bb,fc,fff->",
    "fbb,dfd,fc,fc->",
    "afd,ba,cc,dc->bf",
    "adb,bc,fa,cfc->d",
    "bbd,bda,fc,db->acf",
    "dba,ead,cad->bce",
    "aef,fbc,dca->bde",
]

# Ellipsis test cases from opt_einsum/tests/test_input.py:
_opt_einsum_ellipsis_tests = [
    "...a",
    "a...",
    "...a->...",
    "a...->...",
    "...a->...a",
    "a...a->...a",
    "...,...",
    "...a,...b",
    "...b,...a->...",
    "...ij,...jk->...ik",
]


def _build_numpy_views(path, dim=2):
    """Build small random integer numpy arrays matching *path* (ellipsis
    dimensions get an extra dimension of size *dim*)."""
    lhs = path.split("->")[0]
    views = []
    rng = numpy.random.default_rng(12345 + len(path))
    for spec in lhs.split(","):
        ndim = len(spec.replace("...", "", 1)) + (2 if "..." in spec else 0)
        if ndim == 0:
            views.append(rng.integers(0, 10))
        else:
            views.append(rng.integers(0, 10, dim**ndim).reshape((dim,)*ndim))
    return views


def test_opt_einsum_test_cases():
    if not numpy:
        return

    for path in _opt_einsum_tests:
        views = _build_numpy_views(path)
        total_ndim = sum(v.ndim if hasattr(v, "ndim") else 0 for v in views)
        expected = numpy.einsum(path, *views)
        sa = [Array(v.tolist()) if hasattr(v, "ndim") and v.ndim else Integer(int(v)) for v in views]
        ei = Einsum(path, *sa)
        assert tuple(int(i) for i in ei.shape) == expected.shape, path
        if total_ndim <= 12:
            result = numpy.array(_tolist(ei.as_explicit()))
            assert (result == expected).all(), path


def test_opt_einsum_ellipsis_test_cases():
    if not numpy:
        return

    for path in _opt_einsum_ellipsis_tests:
        views = _build_numpy_views(path)
        expected = numpy.einsum(path, *views)
        sa = [Array(v.tolist()) for v in views]
        ei = Einsum(path, *sa)
        assert tuple(int(i) for i in ei.shape) == expected.shape, path
        result = numpy.array(_tolist(ei.as_explicit()))
        assert (result == expected).all(), path


def test_from_to_einsum_numpy():
    if numpy is None:
        return

    base_dim = 2

    # Seed NumPy's generator from the (seeded) ``random`` module, so that
    # the arrays depend on the test seed reported by the test runner:
    rng = numpy.random.default_rng(random.getrandbits(32))
    a = rng.integers(0, 100, base_dim**4).reshape(base_dim, base_dim, base_dim, base_dim)
    b = rng.integers(0, 100, base_dim**4).reshape(base_dim, base_dim, base_dim, base_dim)
    c = rng.integers(0, 100, base_dim**2).reshape(base_dim, base_dim)

    _compare_einsum_numpy("ijkl->jkil", a)
    _compare_einsum_numpy("iijk->kj", a)
    _compare_einsum_numpy("iijk->kji", a)

    for i in range(5):
        path = _get_random_path(a, b)
        _compare_einsum_numpy(path, a, b)

    for i in range(2):
        path = _get_random_path(c, c, c)
        _compare_einsum_numpy(path, c % 7, c % 11, c % 13)


def test_convert_array_to_einsum_random_numpy():
    if numpy is None:
        return

    # Round-trip random array expressions through Einsum, checking the
    # explicit values against the original expression:
    rng = numpy.random.default_rng(random.getrandbits(32))
    a = rng.integers(0, 10, 16).reshape(2, 2, 2, 2)
    b = rng.integers(0, 10, 4).reshape(2, 2)
    sa = Array(a.tolist())
    sb = Array(b.tolist())

    exprs = [
        ArrayContraction(ArrayTensorProduct(sa, sb), (1, 4), (2, 5)),
        ArrayDiagonal(ArrayTensorProduct(sa, sb), (0, 4), (1, 5)),
        PermuteDims(ArrayContraction(ArrayTensorProduct(sa, sb), (0, 4)), [2, 0, 1, 3]),
        ArrayContraction(PermuteDims(ArrayTensorProduct(sa, sb), [5, 0, 1, 2, 3, 4]), (0, 2)),
        ArrayContraction(ArrayDiagonal(ArrayTensorProduct(sa, sb), (0, 4)), (0, 3)),
        ArrayDiagonal(PermuteDims(ArrayTensorProduct(sa, sb), [5, 0, 1, 2, 3, 4]), (0, 2), (1, 3)),
    ]
    for expr in exprs:
        ei = convert_array_to_einsum(expr)
        assert ei.as_explicit() == expr.as_explicit(), expr
