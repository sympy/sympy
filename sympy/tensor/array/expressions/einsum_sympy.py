from collections import Counter, defaultdict
from functools import singledispatch
from string import ascii_letters
from typing import List, Tuple

from sympy import Expr
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from sympy.tensor.array.expressions import ArrayTensorProduct, ArrayContraction, ArrayDiagonal, PermuteDims
from sympy.tensor.array.expressions.array_expressions import _ArrayExpr, get_rank, ArrayAdd, get_shape


def _convert_einsum_to_sympy_array(path, *args):
    path_src, path_dst = _parse_path_from_string(path)
    _check_path_compatibility_with_array_args(path_src, path_dst, args)

    counted_src = Counter([j for i in path_src for j in i])
    counted_dst = Counter([i for i in path_dst])

    indices_contraction = {i for i, count in counted_src.items() if (count > 1 or i not in path_dst) and i not in counted_dst}
    indices_diagonalization = [i for i, count in counted_src.items() if count > 1 and i in counted_dst]

    path_src_flat = [j for i in path_src for j in i]

    order_after_contraction = [
        i for i in path_src_flat if i not in indices_contraction]
    order_after_diagonalization = [
        i for i in order_after_contraction if i not in indices_diagonalization]
    order_after_diagonalization_part2 = [
        i for i in order_after_contraction if i in indices_diagonalization]
    [
        order_after_diagonalization.append(i)
        for i in order_after_diagonalization_part2
        if i not in order_after_diagonalization
    ]

    src_indices_pos = defaultdict(list)
    for i, ind in enumerate(path_src_flat):
        src_indices_pos[ind].append(i)

    contraction_indices = []
    for i in indices_contraction:
        contraction_indices.append(src_indices_pos[i])

    diagonalization_indices = []
    for i in indices_diagonalization:
        diagonalization_indices.append(src_indices_pos[i])

    # Wrap in an ArrayTensorProduct object:
    if len(args) > 1:
        base = ArrayTensorProduct(*args)
    else:
        base = args[0]

    # Wrap in an ArrayContraction object:
    if len(contraction_indices) > 0:
        array_contraction = ArrayContraction(base, *contraction_indices)
        diagonalization_indices = ArrayContraction._push_indices_up(
            array_contraction.contraction_indices, diagonalization_indices)
    else:
        array_contraction = base

    # Wrap in an ArrayDiagonal object:
    array_diagonal = (
        ArrayDiagonal(array_contraction, *diagonalization_indices)
        if len(diagonalization_indices) > 0 else
        array_contraction
    )

    # Wrap in a PermuteDims object:
    permu = (
        PermuteDims(array_diagonal, index_order_old=order_after_diagonalization, index_order_new=path_dst)
        if order_after_diagonalization != list(path_dst) else
        array_diagonal
    )
    return permu


def _parse_path_from_string(path):
    path = str(path)
    path = path.replace(" ", "").replace("\t", "").strip()
    indices1, indices2 = path.split("->")
    indices1_by_arg = indices1.split(",")
    return indices1_by_arg, indices2


def _check_path_compatibility_with_array_args(path_src, path_dst, array_args):
    if len(path_src) != len(array_args):
        raise ValueError("mismatch between indices and number of arrays")
    for i, arg in zip(path_src, array_args):
        if len(i) != get_rank(arg):
            raise ValueError(f"wrong number of indices for {arg}")


class Einsum(_ArrayExpr):
    """
    The Einsum class represents an operator equivalent to NumPy's
    einsum. Unlike NumPy, it does not evaluate the expression, unless
    .as_explicit() is called.

    Examples
    ========

    >>> from sympy.tensor.array.expressions.einsum_sympy import Einsum
    >>> from sympy.tensor.array.expressions import ArraySymbol
    >>> A = ArraySymbol("A", (4, 4, 4))

    >>> Einsum("abb->a", A)
    Einsum(Str('abb->a'), A)

    >>> Einsum("abb->a", A).as_array_expression()
    ArrayContraction(A, (1, 2))
    >>> Einsum("abc->bca", A).as_array_expression()
    PermuteDims(A, (0 1 2))

    >>> from sympy.abc import x, y, z, t
    >>> from sympy import Array
    >>> B = Array([[x, y], [z, t]])
    >>> B
    [[x, y], [z, t]]

    >>> Einsum("aa->a", B).as_explicit()
    [x, t]
    >>> Einsum("aa->", B).as_explicit()
    t + x
    >>> Einsum("ab->ba", B).as_explicit()
    [[x, z], [y, t]]
    >>> Einsum("ab->a", B).as_explicit()
    [x + y, t + z]
    >>> Einsum("ab->b", B).as_explicit()
    [x + z, t + y]
    """

    def __new__(cls, *args):
        if len(args) == 0:
            raise ValueError("invalid number of arguments")
        path_str: Str
        array_args = []
        if isinstance(args[0], (str, Str)):
            path_src, path_dst = _parse_path_from_string(str(args[0]))
            if isinstance(args[0], str):
                path_str = Str(args[0])
            else:
                path_str = args[0]
            array_args.extend(args[1:])
            _check_path_compatibility_with_array_args(path_src, path_dst, array_args)
        else:
            raise ValueError("first argument must be a string")
        array_args = [_sympify(arg) for arg in array_args]
        obj = _ArrayExpr.__new__(cls, path_str, *array_args)
        path_src = tuple(tuple(i) for i in path_src)
        path_dst = tuple(path_dst)
        obj._path_src = path_src
        obj._path_dst = path_dst
        obj._array_args = tuple(array_args)
        obj._shape = Einsum._get_shape_from_args(path_src, path_dst, array_args)
        return obj

    @staticmethod
    def _get_shape_from_args(path_src: Tuple[tuple], path_dst: tuple, array_args: List[Expr]):
        dims = {}
        for i, arg in zip(path_src, array_args):
            for j, ind2 in enumerate(i):
                dims[ind2] = get_shape(arg)[j]
        return tuple(dims[i] for i in path_dst)

    @property
    def path_string(self) -> str:
        """
        Get the path as a string:

        >>> from sympy.tensor.array.expressions.einsum_sympy import Einsum
        >>> from sympy.tensor.array.expressions import ArraySymbol
        >>> A = ArraySymbol("A", (4, 4, 4))

        >>> Einsum("abc->ac", A).path_string
        'abc->ac'
        """
        return str(self._args[0])

    @property
    def path_src(self):
        return self._path_src

    @property
    def path_dst(self):
        return self._path_dst

    @property
    def array_args(self):
        """
        Return the arrays involved in the Einsum operations.
        """
        return self._array_args

    def as_explicit(self):
        """
        Compute the component-explicit array of the Einsum operations, if
        such computation is possible.
        """
        ret = _convert_einsum_to_sympy_array(
            self.path_string,
            *[arg.as_explicit() if hasattr(arg, "as_explicit") else arg for arg in self.array_args])
        if hasattr(ret, "as_explicit"):
            ret = ret.as_explicit()
        return ret

    def as_array_expression(self):
        """
        Transform the Einsum object into separate contraction, diagonalization
        and dimensional permutation objects.
        """
        return _convert_einsum_to_sympy_array(self.path_string, *self.array_args)

    @property
    def shape(self):
        return self._shape


class _EinsumBuilder:
    def __init__(self, path_src: List[List[int]], path_dst: List[int], *args):
        self.path_src: List[List[int]] = path_src
        self.path_dst: List[int] = path_dst
        self.args: List[Expr] = list(args)

    def reindex(self):
        indices: List[int] = sorted({j for i in self.path_src for j in i})
        self.path_src = [[indices.index(j) for j in i] for i in self.path_src]
        self.path_dst = [indices.index(i) for i in self.path_dst]

    def build(self):
        if len(self.args) == 1 and self.path_src[0] == self.path_dst:
            return self.args[0]
        path = ",".join(["".join(ascii_letters[j] for j in i) for i in self.path_src])
        path += "->"
        path += "".join(ascii_letters[i] for i in self.path_dst)
        return Einsum(path, *self.args)


@singledispatch
def _convert_array_to_einsum(expr) -> _EinsumBuilder:
    raise NotImplementedError()


@_convert_array_to_einsum.register(_ArrayExpr)
def _(expr: _ArrayExpr) -> _EinsumBuilder:
    indices = list(range(get_rank(expr)))
    # Make sure "indices" gets duplicated in memory with [:] !
    return _EinsumBuilder([indices[:]], indices, expr)


@_convert_array_to_einsum.register(PermuteDims)
def _(expr: PermuteDims) -> _EinsumBuilder:
    einsum_builder: _EinsumBuilder = _convert_array_to_einsum(expr.expr)
    if len(set(einsum_builder.path_dst)) != len(einsum_builder.path_dst):
        raise NotImplementedError("repeated destination indices are not supported")
    permutation = expr.permutation
    inv_permutation = permutation**(-1)
    path_dst = [-1 for i in einsum_builder.path_dst]
    for i, e in enumerate(einsum_builder.path_dst):
        new_i = inv_permutation(i)
        path_dst[new_i] = e
    einsum_builder.path_dst = path_dst
    return einsum_builder


@_convert_array_to_einsum.register(ArrayTensorProduct)
def _(expr: ArrayTensorProduct) -> _EinsumBuilder:
    einsum_builder = _EinsumBuilder([], [])
    cumul_dim = 0
    for arg in expr.args:
        arg_ei: _EinsumBuilder = _convert_array_to_einsum(arg)
        einsum_builder.path_src.extend([[j + cumul_dim for j in i] for i in arg_ei.path_src])
        einsum_builder.path_dst.extend([i + cumul_dim for i in arg_ei.path_dst])
        einsum_builder.args.extend(arg_ei.args)
        cumul_dim += get_rank(arg)
    return einsum_builder


@_convert_array_to_einsum.register(ArrayContraction)
def _(expr: ArrayContraction) -> _EinsumBuilder:
    einsum_builder = _convert_array_to_einsum(expr.expr)
    for contraction_tuple in expr.contraction_indices:
        lowest = contraction_tuple[0]
        einsum_builder.path_src = [[lowest if j in contraction_tuple else j for j in i] for i in
                                   einsum_builder.path_src]
        einsum_builder.path_dst = [i for i in einsum_builder.path_dst if i not in contraction_tuple]
    einsum_builder.reindex()
    return einsum_builder


@_convert_array_to_einsum.register(ArrayDiagonal)
def _(expr: ArrayDiagonal) -> _EinsumBuilder:
    einsum_builder = _convert_array_to_einsum(expr.expr)
    for diag_tuple in expr.diagonal_indices:
        lowest = diag_tuple[0]
        einsum_builder.path_src = [[lowest if j in diag_tuple else j for j in i] for i in einsum_builder.path_src]
        # ArrayDiagonal resorts the diagonal indices at the end:
        einsum_builder.path_dst = [i for i in einsum_builder.path_dst if i not in diag_tuple] + [lowest]
    einsum_builder.reindex()
    return einsum_builder


@_convert_array_to_einsum.register(ArrayAdd)
def _(expr: ArrayAdd):
    new_expr = ArrayAdd(*[convert_array_to_einsum(arg) for arg in expr.args])
    indices = list(range(get_rank(new_expr)))
    # Make sure "indices" gets duplicated in memory with [:] !
    return _EinsumBuilder([indices[:]], indices, new_expr)


def convert_array_to_einsum(expr):
    """
    Convert array expression of nested tensor products, contractions and
    diagonalizations into an Einsum object.

    Examples
    ========

    >>> from sympy.tensor.array.expressions.einsum_sympy import convert_array_to_einsum
    >>> from sympy.tensor.array.expressions import ArraySymbol

    Create an array symbol with shape (4, 4, 4), i.e. a three-dimensional array:

    >>> A = ArraySymbol("A", (4, 4, 4))

    Contraction of the 1st and 2nd axes returning a 1-dim vector:

    >>> from sympy.tensor.array.expressions import ArrayContraction
    >>> convert_array_to_einsum(ArrayContraction(A, (0, 1)))
    Einsum(Str('aab->b'), A)

    Diagonalize the 1st and 2nd axes (i.e. only take the diagonal elements,
    discarding off diagonal elements of the first two axes). Remember that
    ArrayDiagonal reorders the axes by putting the diagonalized dimensions at
    the end:

    >>> from sympy.tensor.array.expressions import ArrayDiagonal
    >>> convert_array_to_einsum(ArrayDiagonal(A, (0, 1)))
    Einsum(Str('aab->ba'), A)

    A simple permutation of axes, without reducing the dimension:

    >>> from sympy.tensor.array.expressions import PermuteDims
    >>> convert_array_to_einsum(PermuteDims(A, index_order_old="ijk", index_order_new="kji"))
    Einsum(Str('abc->cba'), A)
    """
    einsum_builder: _EinsumBuilder = _convert_array_to_einsum(expr)
    return einsum_builder.build()
