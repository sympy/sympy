"""
Einstein-summation operator for array expressions.

This module defines :class:`Einsum`, an unevaluated equivalent of NumPy's
``einsum`` operator acting on array expressions, together with the conversion
functions :func:`convert_einsum_to_array` (from a subscript string and
operands to nested array-expression operators) and
:func:`convert_array_to_einsum` (from nested array-expression operators to an
``Einsum`` object).

The supported subscript syntax follows NumPy's ``einsum``:

- explicit mode: ``"ij,jk->ik"``;
- implicit mode (no ``->``): the output indices are the indices appearing
  exactly once, sorted alphabetically, e.g. ``"ij,jk"`` is equivalent to
  ``"ij,jk->ik"``;
- ellipsis: ``"..."`` stands for the dimensions not covered by explicit
  subscript letters, e.g. ``"...ij,...jk->...ik"`` expresses a
  "batched" matrix product.  Unlike NumPy, broadcasting of ellipsis
  dimensions is not supported: all operands using ``"..."`` must have the
  same number (and sizes) of ellipsis dimensions.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from functools import singledispatch
from string import ascii_letters

from sympy.core.basic import Basic
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from sympy.tensor.array.dense_ndim_array import ImmutableDenseNDimArray
from sympy.tensor.array.expressions.array_expressions import (
    _ArrayExpr, ArrayAdd, ArrayContraction, ArrayDiagonal, ArrayTensorProduct,
    PermuteDims, get_ndim, get_shape)


def _parse_einsum_path(path, shapes):
    """Parse and normalize an einsum subscript string.

    Given the subscript string *path* and the list of operand *shapes*,
    return the pair ``(src_specs, dst_spec)`` where ``src_specs`` is a list
    of index strings (one for each operand) and ``dst_spec`` is the output
    index string.  The returned specifications are fully explicit: the
    implicit-mode output is resolved and ellipses are replaced by fresh
    index letters.

    Raises ``ValueError`` for invalid subscript strings and
    ``NotImplementedError`` for ellipsis broadcasting, which is not
    supported.
    """
    path = str(path).replace(" ", "").replace("\t", "")

    if path.count("->") > 1:
        raise ValueError("subscript string contains more than one '->'")

    if "->" in path:
        lhs, dst_spec = path.split("->")
        explicit = True
    else:
        lhs, dst_spec = path, None
        explicit = False

    src_specs = lhs.split(",")
    if len(src_specs) != len(shapes):
        raise ValueError(
            "mismatch between subscripts and number of operands: "
            f"{len(src_specs)} operand subscripts, {len(shapes)} operands")

    def check_spec(spec):
        # A specification may contain letters and at most one ellipsis,
        # written exactly as "...":
        no_ell = spec.replace("...", "", 1)
        if not all(c in ascii_letters for c in no_ell):
            raise ValueError(
                f"invalid characters in einsum subscripts {spec!r}")

    for spec in src_specs:
        check_spec(spec)
    if explicit:
        check_spec(dst_spec)

    # Number of dimensions covered by the ellipsis of each operand:
    ell_ndims = []
    for spec, shape in zip(src_specs, shapes):
        n_letters = len(spec.replace("...", "", 1))
        if "..." in spec:
            n_ell = len(shape) - n_letters
            if n_ell < 0:
                raise ValueError(
                    f"wrong number of indices in {spec!r} for operand of "
                    f"shape {shape}")
            ell_ndims.append(n_ell)
        else:
            if n_letters != len(shape):
                raise ValueError(
                    f"wrong number of indices in {spec!r} for operand of "
                    f"shape {shape}")
            ell_ndims.append(None)

    n_ell_max = max([i for i in ell_ndims if i is not None], default=0)
    for n_ell in ell_ndims:
        if n_ell is not None and n_ell != n_ell_max:
            raise NotImplementedError(
                "broadcasting of ellipsis dimensions is not supported: all "
                "operands using '...' must have the same number of "
                "ellipsis dimensions")

    # Assign fresh letters to the ellipsis dimensions (shared by all
    # operands, as ellipsis dimensions are required to match):
    used_letters = set("".join(src_specs).replace(".", ""))
    if explicit:
        used_letters.update(dst_spec.replace(".", ""))
    free_letters = [c for c in ascii_letters if c not in used_letters]
    if len(free_letters) < n_ell_max:
        raise ValueError("too many einsum indices: only letters a-z and "
                         "A-Z are available")
    ell_letters = "".join(free_letters[:n_ell_max])

    src_specs = [spec.replace("...", ell_letters, 1) for spec in src_specs]

    counted_src = Counter("".join(src_specs))

    if explicit:
        if "..." in dst_spec:
            dst_spec = dst_spec.replace("...", ell_letters, 1)
        elif n_ell_max > 0:
            raise ValueError(
                "operands have more dimensions than subscripts given, but "
                "no '...' ellipsis was provided in the output to cover the "
                "extra dimensions")
        counted_dst = Counter(dst_spec)
        for ind, count in counted_dst.items():
            if count > 1:
                raise ValueError(
                    f"output subscripts include index {ind!r} multiple "
                    "times")
            if ind not in counted_src:
                raise ValueError(
                    f"output subscripts include index {ind!r} which does "
                    "not appear in the input subscripts")
    else:
        # Implicit mode: indices appearing exactly once, sorted
        # alphabetically, preceded by the ellipsis dimensions:
        dst_spec = ell_letters + "".join(sorted(
            ind for ind, count in counted_src.items()
            if count == 1 and ind not in ell_letters))

    # Consistency of the dimensions associated with each index:
    dims = {}
    for spec, shape in zip(src_specs, shapes):
        for ind, dim in zip(spec, shape):
            if ind in dims and dims[ind] != dim:
                raise ValueError(
                    f"dimension mismatch for index {ind!r}: {dims[ind]} "
                    f"and {dim}")
            dims[ind] = dim

    return src_specs, dst_spec


def convert_einsum_to_array(path, *args):
    """
    Convert an einsum subscript string and its operands into nested array
    expression operators (``ArrayTensorProduct``, ``ArrayContraction``,
    ``ArrayDiagonal``, ``PermuteDims``).

    Examples
    ========

    >>> from sympy.tensor.array.expressions import ArraySymbol, convert_einsum_to_array
    >>> A = ArraySymbol("A", (3, 3))
    >>> B = ArraySymbol("B", (3, 3))

    The matrix product:

    >>> convert_einsum_to_array("ij,jk->ik", A, B)
    ArrayContraction(ArrayTensorProduct(A, B), (1, 2))

    The main diagonal of a matrix:

    >>> convert_einsum_to_array("ii->i", A)
    ArrayDiagonal(A, (0, 1))

    The transposition of a matrix, in implicit mode:

    >>> convert_einsum_to_array("ji", A)
    PermuteDims(A, (0 1))
    """
    args = [_sympify_einsum_arg(arg) for arg in args]
    path_src, path_dst = _parse_einsum_path(path, [get_shape(arg) for arg in args])

    counted_src = Counter([j for i in path_src for j in i])
    counted_dst = Counter(path_dst)

    indices_contraction = {i for i, count in counted_src.items() if (count > 1 or i not in path_dst) and i not in counted_dst}
    indices_diagonalization = [i for i, count in counted_src.items() if count > 1 and i in counted_dst]

    path_src_flat = [j for i in path_src for j in i]

    order_after_contraction = [
        i for i in path_src_flat if i not in indices_contraction]
    order_after_diagonalization = [
        i for i in order_after_contraction if i not in indices_diagonalization]
    order_after_diagonalization_part2 = [
        i for i in order_after_contraction if i in indices_diagonalization]
    for i in order_after_diagonalization_part2:
        if i not in order_after_diagonalization:
            order_after_diagonalization.append(i)

    src_indices_pos = defaultdict(list)
    for i, ind in enumerate(path_src_flat):
        src_indices_pos[ind].append(i)

    contraction_indices = []
    for i in sorted(indices_contraction, key=lambda x: src_indices_pos[x][0]):
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
            contraction_indices, diagonalization_indices)
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
        PermuteDims(array_diagonal, index_order_old=order_after_diagonalization, index_order_new=list(path_dst))
        if order_after_diagonalization != list(path_dst) else
        array_diagonal
    )
    return permu


def _sympify_einsum_arg(arg):
    """Sympify an ``Einsum`` operand, converting non-SymPy array-like
    objects (nested lists/tuples, NumPy arrays) to ``ImmutableDenseNDimArray``.
    """
    if isinstance(arg, (list, tuple)):
        return ImmutableDenseNDimArray(arg)
    if hasattr(arg, "__array__") and not isinstance(arg, Basic):
        # NumPy arrays and array-likes:
        return ImmutableDenseNDimArray(arg.tolist())
    return _sympify(arg)


class Einsum(_ArrayExpr):
    """
    The Einsum class represents an operator equivalent to NumPy's
    einsum. Unlike NumPy, it does not evaluate the expression, unless
    ``.as_explicit()`` or ``.doit()`` is called.

    The subscript string is normalized on construction: the output indices
    of implicit-mode strings are made explicit and ellipses are replaced
    by fresh index letters.

    Examples
    ========

    >>> from sympy.tensor.array.expressions import Einsum, ArraySymbol
    >>> A = ArraySymbol("A", (4, 4, 4))

    >>> Einsum("abb->a", A)
    Einsum('abb->a', A)

    >>> Einsum("abb->a", A).as_array_expression()
    ArrayContraction(A, (1, 2))
    >>> Einsum("abc->bca", A).as_array_expression()
    PermuteDims(A, (0 1 2))

    In implicit mode the output indices are the indices appearing exactly
    once, sorted alphabetically:

    >>> B = ArraySymbol("B", (2, 3))
    >>> C = ArraySymbol("C", (3, 4))
    >>> Einsum("ij,jk", B, C)
    Einsum('ij,jk->ik', B, C)

    An ellipsis covers the dimensions not named by explicit index letters,
    so that the same subscript string may be applied to operands of any
    dimension:

    >>> D = ArraySymbol("D", (5, 2, 3))
    >>> E = ArraySymbol("E", (5, 3, 4))
    >>> Einsum("...ij,...jk->...ik", D, E).shape
    (5, 2, 4)

    Explicit arrays are supported:

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
        if not isinstance(args[0], (str, Str)):
            raise ValueError("first argument must be a string")
        array_args = [_sympify_einsum_arg(arg) for arg in args[1:]]
        path_src, path_dst = _parse_einsum_path(
            str(args[0]), [get_shape(arg) for arg in array_args])
        path_str = Str(",".join(path_src) + "->" + path_dst)
        obj = _ArrayExpr.__new__(cls, path_str, *array_args)
        path_src = tuple(tuple(i) for i in path_src)
        path_dst = tuple(path_dst)
        obj._path_src = path_src
        obj._path_dst = path_dst
        obj._array_args = tuple(array_args)
        obj._shape = Einsum._get_shape_from_args(path_src, path_dst, array_args)
        return obj

    @staticmethod
    def _get_shape_from_args(path_src, path_dst, array_args):
        dims = {}
        for i, arg in zip(path_src, array_args):
            for j, ind2 in enumerate(i):
                dims[ind2] = get_shape(arg)[j]
        return tuple(dims[i] for i in path_dst)

    @property
    def path_string(self) -> str:
        """
        Get the normalized subscript string:

        >>> from sympy.tensor.array.expressions import Einsum, ArraySymbol
        >>> A = ArraySymbol("A", (4, 4, 4))

        >>> Einsum("abc->ac", A).path_string
        'abc->ac'

        Implicit-mode subscript strings are normalized to explicit mode:

        >>> Einsum("abc", A).path_string
        'abc->abc'
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
        ret = convert_einsum_to_array(
            self.path_string,
            *[arg.as_explicit() if hasattr(arg, "as_explicit") else arg for arg in self.array_args])
        if hasattr(ret, "as_explicit"):
            ret = ret.as_explicit()
        return ret

    def as_array_expression(self):
        """
        Transform the Einsum object into separate contraction, diagonalization
        and dimensional permutation objects.

        This transformation is not applied recursively to nested ``Einsum``
        objects appearing among the operands; use ``.doit()`` for a
        recursive transformation.
        """
        return convert_einsum_to_array(self.path_string, *self.array_args)

    def doit(self, **hints):
        """
        Rewrite as nested array expression operators, recursively converting
        nested ``Einsum`` objects among the operands as well.

        Examples
        ========

        >>> from sympy.tensor.array.expressions import Einsum, ArraySymbol
        >>> A = ArraySymbol("A", (3, 3))
        >>> B = ArraySymbol("B", (3, 3))
        >>> Einsum("ij,jk->ik", A, Einsum("ab->ba", B)).doit()
        ArrayContraction(ArrayTensorProduct(A, B), (1, 3))
        """
        deep = hints.get("deep", True)
        if deep:
            args = [arg.doit(**hints) if isinstance(arg, Basic) else arg
                    for arg in self.array_args]
        else:
            args = self.array_args
        ret = convert_einsum_to_array(self.path_string, *args)
        if deep and isinstance(ret, Basic):
            ret = ret.doit(**hints)
        return ret

    @property
    def shape(self):
        return self._shape

    def _sympystr(self, printer):
        return "%s(%s, %s)" % (
            type(self).__name__,
            repr(self.path_string),
            ", ".join(printer._print(arg) for arg in self.array_args))

    def _latex(self, printer):
        # The subscript string is typeset verbatim: printed as a
        # mathematical expression, "->" would come out as an arrow-less
        # minus and greater-than sign sequence.
        return r"\operatorname{Einsum}\left(\texttt{%s}, %s\right)" % (
            self.path_string,
            ", ".join(printer._print(arg) for arg in self.array_args))


class _EinsumBuilder:
    def __init__(self, path_src, path_dst, *args):
        self.path_src = path_src
        self.path_dst = path_dst
        self.args = list(args)

    def reindex(self):
        indices = sorted({j for i in self.path_src for j in i})
        self.path_src = [[indices.index(j) for j in i] for i in self.path_src]
        self.path_dst = [indices.index(i) for i in self.path_dst]

    def build(self):
        if len(self.args) == 1 and self.path_src[0] == self.path_dst:
            return self.args[0]
        indices = {j for i in self.path_src for j in i}
        if len(indices) > len(ascii_letters):
            raise ValueError("too many einsum indices: only letters a-z "
                             "and A-Z are available")
        path = ",".join(["".join(ascii_letters[j] for j in i) for i in self.path_src])
        path += "->"
        path += "".join(ascii_letters[i] for i in self.path_dst)
        return Einsum(path, *self.args)


@singledispatch
def _convert_array_to_einsum(expr) -> _EinsumBuilder:
    # Objects that are not expressible with the einsum syntax (e.g.
    # Reshape, ArrayElementwiseApplyFunc, explicit arrays, matrices,
    # scalars) are kept as opaque operands:
    shape = get_shape(expr)
    if shape is not None:
        indices = list(range(len(shape)))
        # Make sure "indices" gets duplicated in memory with [:] !
        return _EinsumBuilder([indices[:]], indices, expr)
    raise NotImplementedError(
        f"cannot convert {type(expr).__name__} object to Einsum")


@_convert_array_to_einsum.register(_ArrayExpr)
def _(expr: _ArrayExpr) -> _EinsumBuilder:
    indices = list(range(get_ndim(expr)))
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
        cumul_dim += get_ndim(arg)
    return einsum_builder


@_convert_array_to_einsum.register(ArrayContraction)
def _(expr: ArrayContraction) -> _EinsumBuilder:
    einsum_builder = _convert_array_to_einsum(expr.expr)
    # The contraction indices are axis positions of the nested expression,
    # while the builder's src/dst entries are einsum index identifiers:
    # translate positions to identifiers through the nested destination.
    old_dst = einsum_builder.path_dst
    removed_positions = set()
    for contraction_tuple in expr.contraction_indices:
        ids = {old_dst[pos] for pos in contraction_tuple}
        lowest = min(ids)
        einsum_builder.path_src = [[lowest if j in ids else j for j in i] for i in
                                   einsum_builder.path_src]
        removed_positions.update(contraction_tuple)
    einsum_builder.path_dst = [
        e for pos, e in enumerate(old_dst) if pos not in removed_positions]
    einsum_builder.reindex()
    return einsum_builder


@_convert_array_to_einsum.register(ArrayDiagonal)
def _(expr: ArrayDiagonal) -> _EinsumBuilder:
    einsum_builder = _convert_array_to_einsum(expr.expr)
    # As in the ArrayContraction case, translate the axis positions of the
    # diagonal indices to einsum index identifiers.  ArrayDiagonal appends
    # the diagonalized axes at the end, one for each diagonal group:
    old_dst = einsum_builder.path_dst
    removed_positions = set()
    appended = []
    for diag_tuple in expr.diagonal_indices:
        ids = {old_dst[pos] for pos in diag_tuple}
        lowest = min(ids)
        einsum_builder.path_src = [[lowest if j in ids else j for j in i] for i in einsum_builder.path_src]
        removed_positions.update(diag_tuple)
        appended.append(lowest)
    einsum_builder.path_dst = [
        e for pos, e in enumerate(old_dst) if pos not in removed_positions
    ] + appended
    einsum_builder.reindex()
    return einsum_builder


@_convert_array_to_einsum.register(ArrayAdd)
def _(expr: ArrayAdd):
    new_expr = ArrayAdd(*[convert_array_to_einsum(arg) for arg in expr.args])
    indices = list(range(get_ndim(new_expr)))
    # Make sure "indices" gets duplicated in memory with [:] !
    return _EinsumBuilder([indices[:]], indices, new_expr)


@_convert_array_to_einsum.register(Einsum)
def _(expr: Einsum) -> _EinsumBuilder:
    indices = list(range(get_ndim(expr)))
    # Make sure "indices" gets duplicated in memory with [:] !
    return _EinsumBuilder([indices[:]], indices, expr)


def convert_array_to_einsum(expr):
    """
    Convert array expression of nested tensor products, contractions and
    diagonalizations into an Einsum object.

    Examples
    ========

    >>> from sympy.tensor.array.expressions import ArraySymbol, convert_array_to_einsum

    Create an array symbol with shape (4, 4, 4), i.e. a three-dimensional array:

    >>> A = ArraySymbol("A", (4, 4, 4))

    Contraction of the 1st and 2nd axes returning a 1-dim vector:

    >>> from sympy.tensor.array.expressions import ArrayContraction
    >>> convert_array_to_einsum(ArrayContraction(A, (0, 1)))
    Einsum('aab->b', A)

    Diagonalize the 1st and 2nd axes (i.e. only take the diagonal elements,
    discarding off diagonal elements of the first two axes). Remember that
    ArrayDiagonal reorders the axes by putting the diagonalized dimensions at
    the end:

    >>> from sympy.tensor.array.expressions import ArrayDiagonal
    >>> convert_array_to_einsum(ArrayDiagonal(A, (0, 1)))
    Einsum('aab->ba', A)

    A simple permutation of axes, without reducing the dimension:

    >>> from sympy.tensor.array.expressions import PermuteDims
    >>> convert_array_to_einsum(PermuteDims(A, index_order_old="ijk", index_order_new="kji"))
    Einsum('abc->cba', A)
    """
    einsum_builder: _EinsumBuilder = _convert_array_to_einsum(expr)
    return einsum_builder.build()
