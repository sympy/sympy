"""Module for differentiation using CSE."""
from __future__ import annotations

from dataclasses import dataclass, replace

from sympy import S, cse, Matrix, Derivative, MatrixBase, Mul, Pow, numbered_symbols
from sympy.core.function import AppliedUndef
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.utilities.iterables import iterable

def _remove_cse_from_derivative(replacements, reduced_expressions):
    """
    This function is designed to postprocess the output of a common subexpression
    elimination (CSE) operation. Specifically, it removes any CSE replacement
    symbols from the arguments of ``Derivative`` terms in the expression. This
    is necessary to ensure that the forward Jacobian function correctly handles
    derivative terms.

    Parameters
    ==========

    replacements : list of (Symbol, expression) pairs
        Replacement symbols and relative common subexpressions that have been
        replaced during a CSE operation.

    reduced_expressions : list of SymPy expressions
        The reduced expressions with all the replacements from the
        replacements list above.

    Returns
    =======

    processed_replacements : list of (Symbol, expression) pairs
        Processed replacement list, in the same format of the
        ``replacements`` input list.

    processed_reduced : list of SymPy expressions
        Processed reduced list, in the same format of the
        ``reduced_expressions`` input list.
    """

    def traverse(node, repl_dict):
        if isinstance(node, Derivative):
            return replace_all(node, repl_dict)
        if not node.args:
            return node
        new_args = [traverse(arg, repl_dict) for arg in node.args]
        return node.func(*new_args)

    def replace_all(node, repl_dict):
        result = node
        while True:
            free_symbols = result.free_symbols
            symbols_dict = {k: repl_dict[k] for k in free_symbols if k in repl_dict}
            if not symbols_dict:
                break
            result = result.xreplace(symbols_dict)
        return result

    repl_dict = dict(replacements)
    processed_replacements = [
        (rep_sym, traverse(sub_exp, repl_dict))
        for rep_sym, sub_exp in replacements
    ]
    processed_reduced = [
        red_exp.__class__([traverse(exp, repl_dict) for exp in red_exp])
        for red_exp in reduced_expressions
    ]

    return processed_replacements, processed_reduced


class _PropagationCache:
    """Cache expression construction during sparse propagation."""

    def __init__(self):
        self._cache = {}
        self._hits = 0
        self._misses = 0

    @property
    def hits(self):
        return self._hits

    @property
    def misses(self):
        return self._misses

    def make(self, func, *args):
        key = (func, tuple(id(arg) for arg in args))
        cached = self._cache.get(key)
        if cached is not None:
            self._hits += 1
            return cached

        self._misses += 1
        result = func(*args)
        self._cache[key] = result
        return result


class _DiffCache:
    """Cache fallback ``diff`` calls during sparse propagation."""

    def __init__(self):
        self._cache = {}
        self._hits = 0
        self._misses = 0

    @property
    def hits(self):
        return self._hits

    @property
    def misses(self):
        return self._misses

    def diff(self, expr, var):
        key = (id(expr), id(var))
        cached = self._cache.get(key)
        if cached is not None:
            self._hits += 1
            return cached
        self._misses += 1
        result = expr.diff(var)
        self._cache[key] = result
        return result


@dataclass
class SparseJacobianIR:
    """Sparse Jacobian intermediate representation.

    This object is the primary internal result of
    ``_forward_jacobian_sparse_cse``. It stores the Jacobian in a COO-like
    sparse form together with row boundaries, the differentiation variables,
    optional intermediate expressions.

    Attributes
    ==========

    shape : tuple[int, int]
        Jacobian shape ``(nrows, ncols)``.
    rows, cols, vals : list
        Parallel arrays describing nonzero Jacobian entries.
    row_slices : list[tuple[int, int]]
        Half-open ``[start, end)`` slices into ``cols`` and ``vals`` for each
        Jacobian row.
    wrt : list
        Differentiation variables in column order.
    intermediates : list[tuple]
        Intermediate expressions referenced by ``vals``.
    """

    shape: tuple[int, int]
    rows: list[int]
    cols: list[int]
    vals: list
    row_slices: list[tuple[int, int]]
    wrt: list
    intermediates: list[tuple]

    def to_matrix(self, matrix_cls):
        """Rebuild the Jacobian as a SymPy matrix."""
        nrows, ncols = self.shape
        entries = {
            (row, col): value
            for row, col, value in zip(self.rows, self.cols, self.vals)
        }
        return matrix_cls.from_dok(nrows, ncols, entries)

    def to_coo_tuple(self):
        """Return a codegen-friendly tuple view."""
        return (
            self.intermediates,
            self.shape,
            self.rows,
            self.cols,
            self.vals,
        )

    def row_vals(self, row_idx):
        """Return ``(col, value)`` pairs for one row."""
        start, end = self.row_slices[row_idx]
        return list(zip(self.cols[start:end], self.vals[start:end]))


def _get_dependencies(expr, wrt):
    """Return column indices that ``expr`` may depend on."""
    wrt_index = {var: i for i, var in enumerate(wrt)}
    wrt_set = set(wrt)

    if all(getattr(var, "is_Symbol", False) for var in wrt_set):
        return {wrt_index[s] for s in expr.free_symbols & wrt_set}

    applied = expr.atoms(AppliedUndef)
    derivatives = expr.atoms(Derivative)
    deps = set()
    for var in wrt:
        if getattr(var, "is_Symbol", False):
            if var in expr.free_symbols:
                deps.add(wrt_index[var])
        elif isinstance(var, AppliedUndef):
            if var in applied:
                deps.add(wrt_index[var])
        elif isinstance(var, Derivative):
            if var in derivatives:
                deps.add(wrt_index[var])
    return deps


def _accumulate_sparse_row(result, col, deriv):
    """Add a derivative contribution to a sparse row map."""
    if deriv is S.Zero:
        return

    current = result.get(col)
    if current is None:
        result[col] = deriv
        return

    new_value = current + deriv
    if new_value is S.Zero:
        result.pop(col, None)
    else:
        result[col] = new_value


def _propagate_add(child_derivs):
    """Propagate a sparse derivative map through an ``Add`` node."""
    result = {}
    for child_deriv in child_derivs:
        for col, deriv in child_deriv.items():
            _accumulate_sparse_row(result, col, deriv)
    return result


def _other_factor_product(node_expr, children, skip_index, cache):
    """Return the product of all factors except ``children[skip_index]``."""
    other_args = tuple(children[k] for k in range(len(children)) if k != skip_index)
    if len(other_args) == 0:
        return S.One
    if len(other_args) == 1:
        return other_args[0]
    return cache.make(Mul, *other_args)


def _propagate_mul(node_expr, children, child_derivs, cache):
    """Propagate a sparse derivative map through a ``Mul`` node."""
    active = [(i, child_deriv) for i, child_deriv in enumerate(child_derivs) if child_deriv]
    if not active:
        return {}

    if len(active) == 1:
        idx, child_deriv = active[0]
        coeff = _other_factor_product(node_expr, children, idx, cache)
        if coeff is S.One:
            return child_deriv
        return {
            col: cache.make(Mul, coeff, deriv)
            for col, deriv in child_deriv.items()
        }

    result = {}
    for idx, child_deriv in active:
        coeff = _other_factor_product(node_expr, children, idx, cache)
        for col, deriv in child_deriv.items():
            contrib = cache.make(Mul, coeff, deriv)
            _accumulate_sparse_row(result, col, contrib)
    return result


def _propagate_unary(node_expr, child, child_deriv, cache):
    """Propagate through common single-argument functions."""
    if not child_deriv:
        return {}

    func = node_expr.func
    if func == exp:
        phi_prime = node_expr
    elif func == sin:
        phi_prime = cache.make(cos, child)
    elif func == cos:
        phi_prime = cache.make(Mul, S.NegativeOne, cache.make(sin, child))
    elif func == log:
        phi_prime = cache.make(Pow, child, S.NegativeOne)
    else:
        phi_prime = node_expr.diff(child)

    if phi_prime is S.Zero:
        return {}
    if phi_prime is S.One:
        return child_deriv

    return {
        col: cache.make(Mul, phi_prime, deriv)
        for col, deriv in child_deriv.items()
    }


def _propagate_pow(node_expr, children, child_derivs, cache):
    """Propagate a sparse derivative map through a ``Pow`` node."""
    base, exp = children
    dbase, dexp = child_derivs

    if not dexp:
        if not dbase:
            return {}
        coeff = cache.make(Mul, exp, cache.make(Pow, base, exp - 1))
        if coeff is S.One:
            return dbase
        return {col: cache.make(Mul, coeff, deriv)
                for col, deriv in dbase.items()}

    if not dbase:
        log_base = cache.make(log, base)
        return {
            col: cache.make(Mul, node_expr, deriv, log_base)
            for col, deriv in dexp.items()
        }

    result = {}
    all_cols = set(dbase) | set(dexp)
    log_base = cache.make(log, base)
    base_inverse = cache.make(Pow, base, S.NegativeOne)
    for col in all_cols:
        value = S.Zero
        if col in dexp:
            value += cache.make(Mul, node_expr, dexp[col], log_base)
        if col in dbase:
            value += cache.make(Mul, node_expr, exp, dbase[col], base_inverse)
        if value is not S.Zero:
            result[col] = value
    return result


def _fallback_local_diff(node_expr, wrt, dcache, columns=None):
    """Fallback local differentiation for unsupported node shapes."""
    if columns is None:
        columns = _get_dependencies(node_expr, wrt)

    result = {}
    for col in columns:
        deriv = dcache.diff(node_expr, wrt[col])
        _accumulate_sparse_row(result, col, deriv)
    return result


def _get_leaf_derivative(expr, wrt_index, replacement_derivatives):
    """"Return sparse derivative for replacement symbols, wrt variables,
    and zero-derivative leaf-like expressions; otherwise return None."""""
    if expr in replacement_derivatives:
        return replacement_derivatives[expr]
    if expr in wrt_index:
        return {wrt_index[expr]: S.One}
    if expr.is_Atom or isinstance(expr, Derivative):
        return {}
    return None


def _get_leaf_support(expr, wrt_index, replacement_support):
    """Return support for leaf-like expressions, or ``None`` to recurse."""
    if expr in replacement_support:
        return replacement_support[expr]
    if expr in wrt_index:
        return frozenset((wrt_index[expr],))
    if expr.is_Atom or isinstance(expr, Derivative):
        return frozenset()
    return None


def _propagate_support(expr, wrt_index, replacement_support, cache=None):
    """Return an over-approximation of the Jacobian support of ``expr``.

    The returned frozenset contains column indices that may be nonzero in the
    derivative of ``expr`` with respect to ``wrt``. This prepass performs only
    set propagation; it does not build derivative expressions.
    """
    leaf_support = _get_leaf_support(expr, wrt_index, replacement_support)
    if leaf_support is not None:
        return leaf_support

    if cache is not None:
        cached = cache.get(expr)
        if cached is not None:
            return cached

    support = frozenset().union(*(
        _propagate_support(child, wrt_index, replacement_support, cache=cache)
        for child in expr.args
    ))
    if cache is not None:
        cache[expr] = support
    return support


def _sparse_derivative(
    expr, wrt, wrt_index, replacement_derivatives, cache, dcache,
    support_map=None, derivative_cache=None,
):
    """Differentiate ``expr`` to a sparse row map.

    The result is a dictionary mapping Jacobian column indices to derivative
    expressions. Known replacement derivatives, sparse propagation rules, local
    fallback differentiation, and an optional per-call memo cache are used to
    keep the traversal sparse.
    """
    leaf_derivative = _get_leaf_derivative(
        expr, wrt_index, replacement_derivatives)
    if leaf_derivative is not None:
        return leaf_derivative

    if derivative_cache is not None:
        cached = derivative_cache.get(expr)
        if cached is not None:
            return cached

    children = expr.args
    child_derivs = []
    for child in children:
        child_deriv = _get_leaf_derivative(
            child, wrt_index, replacement_derivatives)
        if child_deriv is None:
            child_deriv = _sparse_derivative(
                child, wrt, wrt_index, replacement_derivatives, cache, dcache,
                support_map=support_map,
                derivative_cache=derivative_cache)
        child_derivs.append(child_deriv)

    if expr.is_Add:
        result = _propagate_add(child_derivs)
    elif expr.is_Mul:
        result = _propagate_mul(expr, children, child_derivs, cache)
    elif expr.is_Pow:
        result = _propagate_pow(expr, children, child_derivs, cache)
    elif len(children) == 1:
        result = _propagate_unary(expr, children[0], child_derivs[0], cache)
    else:
        columns = support_map.get(expr) if support_map is not None else None
        if columns is None:
            columns = _get_dependencies(expr, wrt)
        result = _fallback_local_diff(
            expr, wrt, dcache, columns=columns)

    if derivative_cache is not None:
        derivative_cache[expr] = result
    return result


def _build_ir(row_maps, wrt_list, intermediates):
    """Build ``SparseJacobianIR`` from sparse row maps."""
    rows = []
    cols = []
    vals = []
    row_slices = []

    for row_idx, row_map in enumerate(row_maps):
        start = len(rows)
        for col in sorted(row_map):
            rows.append(row_idx)
            cols.append(col)
            vals.append(row_map[col])
        row_slices.append((start, len(rows)))

    return SparseJacobianIR(
        shape=(len(row_maps), len(wrt_list)),
        rows=rows,
        cols=cols,
        vals=vals,
        row_slices=row_slices,
        wrt=wrt_list,
        intermediates=list(intermediates),
    )


def _finalize_sparse_jacobian_result(ir, matrix_cls, return_mode="matrix"):
    """Convert the internal IR into the requested external view."""
    if return_mode == "ir":
        return ir

    jacobian = ir.to_matrix(matrix_cls)
    return ir.intermediates, [jacobian], []


def _count_symbol_uses(expr, symbols):
    """Count replacement symbol uses in one expression."""
    counts = {}
    if not expr.free_symbols:
        return counts
    for arg in expr.args:
        child_counts = _count_symbol_uses(arg, symbols)
        for sym, count in child_counts.items():
            counts[sym] = counts.get(sym, 0) + count
    if expr in symbols:
        counts[expr] = counts.get(expr, 0) + 1
    return counts


def _count_all_symbol_uses(replacements, reduced):
    """Count uses of every replacement symbol in one pass."""
    symbols = {sym for sym, _ in replacements}
    counts = {}
    for _, expr in replacements:
        for sym, count in _count_symbol_uses(expr, symbols).items():
            counts[sym] = counts.get(sym, 0) + count
    for expr in reduced:
        for sym, count in _count_symbol_uses(expr, symbols).items():
            counts[sym] = counts.get(sym, 0) + count
    return counts


def _should_extract(expr, reuse_count):
    """Return ``True`` when a shared expression is worth extracting."""
    if expr.func in (sin, cos, exp, log):
        return reuse_count >= 2
    if expr.is_Pow and expr.exp.is_negative:
        return reuse_count >= 2
    if expr.is_Mul and len(expr.args) >= 4:
        return reuse_count >= 2

    ops = expr.count_ops()
    if ops >= 6 and reuse_count >= 2:
        return True
    if ops >= 3 and reuse_count >= 3:
        return True
    return False


def _filter_replacements(replacements, reduced, should_extract_fn, max_rounds=5):
    """Inline low-value CSE replacements back into the outputs."""
    keep = list(replacements)
    reduced = list(reduced)

    for _ in range(max_rounds):
        next_keep = []
        revert = {}
        counts = _count_all_symbol_uses(keep, reduced)

        for sym, expr in keep:
            if should_extract_fn(expr, counts.get(sym, 0)):
                next_keep.append((sym, expr))
            else:
                revert[sym] = expr

        if not revert:
            return keep, reduced

        expanded_revert = {}
        # ``cse()`` replacements are topologically ordered, so expanding in
        # dictionary insertion order is safe here.
        for sym, expr in revert.items():
            expanded_revert[sym] = expr.xreplace(expanded_revert)

        keep = [(sym, expr.xreplace(expanded_revert)) for sym, expr in next_keep]
        reduced = [expr.xreplace(expanded_revert) for expr in reduced]
    return keep, reduced


def _replace_ir_values(ir, vals, intermediates):
    """Return ``ir`` with updated values and intermediates."""
    return replace(ir, vals=vals, intermediates=intermediates)


def _run_value_cse(values, symbols, should_extract_fn=None):
    """Run one value-level CSE pass and filter low-value replacements.

    This helper is shared by row-local and global Jacobian value CSE modes.
    """
    replacements, reduced = cse(
        values,
        order='canonical',
        symbols=symbols,
    )
    if should_extract_fn is not None and replacements:
        replacements, reduced = _filter_replacements(
            replacements, reduced, should_extract_fn)
    return replacements, reduced


def _row_level_cse(ir, should_extract_fn=None):
    """Apply CSE independently to each row of the sparse IR."""
    row_intermediates = []
    vals = []

    for row_idx, (start, end) in enumerate(ir.row_slices):
        row_vals = ir.vals[start:end]
        if len(row_vals) <= 1:
            vals.extend(row_vals)
            continue

        replacements, reduced = _run_value_cse(
            row_vals,
            numbered_symbols(prefix=f'_d{row_idx}_'),
            should_extract_fn=should_extract_fn,
        )
        row_intermediates.extend(replacements)
        vals.extend(reduced)

    return _replace_ir_values(ir, vals, ir.intermediates + row_intermediates)


def _global_cse(ir, should_extract_fn=None):
    """Apply one CSE pass to all sparse Jacobian values."""
    replacements, reduced = _run_value_cse(
        ir.vals,
        numbered_symbols(prefix='_dg_'),
        should_extract_fn=should_extract_fn,
    )

    return _replace_ir_values(ir, reduced, ir.intermediates + replacements)


def _validate_forward_jacobian_input(reduced_expr, wrt):
    """Normalize and validate Jacobian inputs shared by multiple paths."""
    if not isinstance(reduced_expr[0], MatrixBase):
        raise TypeError("``expr`` must be of matrix type")

    if not (reduced_expr[0].shape[0] == 1 or reduced_expr[0].shape[1] == 1):
        raise TypeError("``expr`` must be a row or a column matrix")

    if not iterable(wrt):
        raise TypeError("``wrt`` must be an iterable of variables")
    if not isinstance(wrt, MatrixBase):
        wrt = Matrix(wrt)

    if not (wrt.shape[0] == 1 or wrt.shape[1] == 1):
        raise TypeError("``wrt`` must be a row or a column matrix")

    return wrt


def _forward_jacobian_sparse_cse(
    replacements,
    reduced_expr,
    wrt,
    return_mode="matrix",
    structure_prepass=False,
    value_cse="none",
):
    """
    Compute a sparse forward Jacobian from CSE output.

    This function consumes the ``(replacements, reduced_expr)`` output of
    ``sympy.cse`` and differentiates the reduced expressions with respect to
    ``wrt`` using sparse forward propagation. Internally it builds a
    ``SparseJacobianIR`` and optionally applies value-level CSE on the emitted
    Jacobian entries.

    Parameters
    ==========

    replacements : list[tuple]
        CSE replacement pairs ``(symbol, expr)``.
    reduced_expr : list
        Reduced expressions returned by ``cse``. The first element must be a
        row or column matrix.
    wrt : iterable
        Differentiation variables. May be a matrix or any iterable accepted by
        ``Matrix``.
    return_mode : {"matrix", "ir"}, optional
        ``"matrix"`` returns the legacy-compatible tuple
        ``(replacements, [jacobian_matrix], [])``.
        ``"ir"`` returns a ``SparseJacobianIR`` instance.
    structure_prepass : bool, optional
        Run a support prepass before value propagation to prune obviously zero
        columns for unsupported node shapes.
    value_cse : {"none", "row", "global"}, optional
        Optional value-level CSE applied after sparse propagation.
        ``"row"`` applies CSE independently per Jacobian row.
        ``"global"`` applies one pass across all Jacobian values.

    Returns
    =======

    tuple or SparseJacobianIR
        Return type depends on ``return_mode``.

    """
    if return_mode not in ("matrix", "ir"):
        raise ValueError("``return_mode`` must be 'matrix' or 'ir'")
    if value_cse not in ("none", "row", "global"):
        raise ValueError("``value_cse`` must be 'none', 'row', or 'global'")

    wrt = _validate_forward_jacobian_input(reduced_expr, wrt)
    matrix_cls = reduced_expr[0].__class__

    replacements, reduced_expr = _remove_cse_from_derivative(replacements, reduced_expr)

    wrt_list = list(wrt)
    wrt_index = {var: i for i, var in enumerate(wrt_list)}
    cache = _PropagationCache()
    dcache = _DiffCache()
    replacement_derivatives = {}
    replacement_support = {}
    support_cache = {}
    derivative_cache = {}

    if structure_prepass:
        for sym, expr in replacements:
            replacement_support[sym] = _propagate_support(
                expr, wrt_index, replacement_support, cache=support_cache)

    for sym, expr in replacements:
        if structure_prepass and not replacement_support[sym]:
            replacement_derivatives[sym] = {}
            continue
        replacement_derivatives[sym] = _sparse_derivative(
            expr, wrt_list, wrt_index, replacement_derivatives, cache, dcache,
            support_map=support_cache if structure_prepass else None,
            derivative_cache=derivative_cache,
        )

    row_maps = []
    for expr in reduced_expr[0]:
        if structure_prepass:
            expr_support = _propagate_support(
                expr, wrt_index, replacement_support, cache=support_cache)
            if not expr_support:
                row_maps.append({})
                continue
        row_maps.append(_sparse_derivative(
            expr, wrt_list, wrt_index, replacement_derivatives, cache, dcache,
            support_map=support_cache if structure_prepass else None,
            derivative_cache=derivative_cache,
        ))
    ir = _build_ir(row_maps, wrt_list, replacements)
    if value_cse == "row":
        ir = _row_level_cse(ir, should_extract_fn=_should_extract)
    elif value_cse == "global":
        ir = _global_cse(ir, should_extract_fn=_should_extract)
    return _finalize_sparse_jacobian_result(ir, matrix_cls, return_mode=return_mode)


def _forward_jacobian_cse(replacements, reduced_expr, wrt):
    """
    Core function to compute the Jacobian of an input Matrix of expressions
    through forward accumulation. Takes directly the output of a CSE operation
    (replacements and reduced_expr), and an iterable of variables (wrt) with
    respect to which to differentiate the reduced expression and returns the
    reduced Jacobian matrix and the ``replacements`` list.

    The function also returns a list of precomputed free symbols for each
    subexpression, which are useful in the substitution process.

    Parameters
    ==========

    replacements : list of (Symbol, expression) pairs
        Replacement symbols and relative common subexpressions that have been
        replaced during a CSE operation.

    reduced_expr : list of SymPy expressions
        The reduced expressions with all the replacements from the
        replacements list above.

    wrt : iterable
        Iterable of expressions with respect to which to compute the
        Jacobian matrix.

    Returns
    =======

    replacements : list of (Symbol, expression) pairs
        Replacement symbols and relative common subexpressions that have been
        replaced during a CSE operation. Compared to the input replacement list,
        the output one doesn't contain replacement symbols inside
        ``Derivative``'s arguments.

    jacobian : list of SymPy expressions
        The list only contains one element, which is the Jacobian matrix with
        elements in reduced form (replacement symbols are present).

    precomputed_fs: list
        List of sets, which store the free symbols present in each sub-expression.
        Useful in the substitution process.
    """

    if not isinstance(reduced_expr[0], MatrixBase):
        raise TypeError("``expr`` must be of matrix type")

    if not (reduced_expr[0].shape[0] == 1 or reduced_expr[0].shape[1] == 1):
        raise TypeError("``expr`` must be a row or a column matrix")

    if not iterable(wrt):
        raise TypeError("``wrt`` must be an iterable of variables")

    elif not isinstance(wrt, MatrixBase):
        wrt = Matrix(wrt)

    if not (wrt.shape[0] == 1 or wrt.shape[1] == 1):
        raise TypeError("``wrt`` must be a row or a column matrix")

    replacements, reduced_expr = _remove_cse_from_derivative(replacements, reduced_expr)

    if replacements:
        rep_sym, sub_expr = map(Matrix, zip(*replacements))
    else:
        rep_sym, sub_expr = Matrix([]), Matrix([])

    l_sub, l_wrt, l_red = len(sub_expr), len(wrt), len(reduced_expr[0])

    f1 = reduced_expr[0].__class__.from_dok(l_red, l_wrt,
        {
            (i, j): diff_value
            for i, r in enumerate(reduced_expr[0])
            for j, w in enumerate(wrt)
            if (diff_value := r.diff(w)) != 0
        },
    )

    if not replacements:
        return [], [f1], []

    f2 = Matrix.from_dok(l_red, l_sub,
        {
            (i, j): diff_value
            for i, (r, fs) in enumerate([(r, r.free_symbols) for r in reduced_expr[0]])
            for j, s in enumerate(rep_sym)
            if s in fs and (diff_value := r.diff(s)) != 0
        },
    )

    rep_sym_set = set(rep_sym)
    precomputed_fs = [s.free_symbols & rep_sym_set for s in sub_expr ]

    c_matrix = Matrix.from_dok(1, l_wrt,
                               {(0, j): diff_value for j, w in enumerate(wrt)
                                if (diff_value := sub_expr[0].diff(w)) != 0})

    for i in range(1, l_sub):

        bi_matrix = Matrix.from_dok(1, i,
                                    {(0, j): diff_value for j in range(i + 1)
                                     if rep_sym[j] in precomputed_fs[i]
                                     and (diff_value := sub_expr[i].diff(rep_sym[j])) != 0})

        ai_matrix = Matrix.from_dok(1, l_wrt,
                                    {(0, j): diff_value for j, w in enumerate(wrt)
                                     if (diff_value := sub_expr[i].diff(w)) != 0})

        if bi_matrix._rep.nnz():
            ci_matrix = bi_matrix.multiply(c_matrix).add(ai_matrix)
            c_matrix = Matrix.vstack(c_matrix, ci_matrix)
        else:
            c_matrix = Matrix.vstack(c_matrix, ai_matrix)

    jacobian = f2.multiply(c_matrix).add(f1)
    jacobian = [reduced_expr[0].__class__(jacobian)]

    return replacements, jacobian, precomputed_fs


def _forward_jacobian_norm_in_cse_out(expr, wrt):
    """
    Function to compute the Jacobian of an input Matrix of expressions through
    forward accumulation. Takes a sympy Matrix of expressions (expr) as input
    and an iterable of variables (wrt) with respect to which to compute the
    Jacobian matrix. The matrix is returned in reduced form (containing
    replacement symbols) along with the ``replacements`` list.

    The function also returns a list of precomputed free symbols for each
    subexpression, which are useful in the substitution process.

    Parameters
    ==========

    expr : Matrix
        The vector to be differentiated.

    wrt : iterable
        The vector with respect to which to perform the differentiation.
        Can be a matrix or an iterable of variables.

    Returns
    =======

    replacements : list of (Symbol, expression) pairs
        Replacement symbols and relative common subexpressions that have been
        replaced during a CSE operation. The output replacement list doesn't
        contain replacement symbols inside ``Derivative``'s arguments.

    jacobian : list of SymPy expressions
        The list only contains one element, which is the Jacobian matrix with
        elements in reduced form (replacement symbols are present).

    precomputed_fs: list
        List of sets, which store the free symbols present in each
        sub-expression. Useful in the substitution process.
    """

    replacements, reduced_expr = cse(expr)
    replacements, jacobian, precomputed_fs = _forward_jacobian_cse(replacements, reduced_expr, wrt)

    return replacements, jacobian, precomputed_fs


def _forward_jacobian(expr, wrt):
    """
    Function to compute the Jacobian of an input Matrix of expressions through
    forward accumulation. Takes a sympy Matrix of expressions (expr) as input
    and an iterable of variables (wrt) with respect to which to compute the
    Jacobian matrix.

    Explanation
    ===========

    Expressions often contain repeated subexpressions. Using a tree structure,
    these subexpressions are duplicated and differentiated multiple times,
    leading to inefficiency.

    Instead, if a data structure called a directed acyclic graph (DAG) is used
    then each of these repeated subexpressions will only exist a single time.
    This function uses a combination of representing the expression as a DAG and
    a forward accumulation algorithm (repeated application of the chain rule
    symbolically) to more efficiently calculate the Jacobian matrix of a target
    expression ``expr`` with respect to an expression or set of expressions
    ``wrt``.

    Note that this function is intended to improve performance when
    differentiating large expressions that contain many common subexpressions.
    For small and simple expressions it is likely less performant than using
    SymPy's standard differentiation functions and methods.

    Parameters
    ==========

    expr : Matrix
        The vector to be differentiated.

    wrt : iterable
        The vector with respect to which to do the differentiation.
        Can be a matrix or an iterable of variables.

    See Also
    ========

    Direct Acyclic Graph : https://en.wikipedia.org/wiki/Directed_acyclic_graph
    """

    replacements, reduced_expr = cse(expr)

    if replacements:
        rep_sym, _ = map(Matrix, zip(*replacements))
    else:
        rep_sym = Matrix([])

    replacements, jacobian, precomputed_fs = _forward_jacobian_cse(replacements, reduced_expr, wrt)

    if not replacements: return jacobian[0]

    sub_rep = dict(replacements)
    for i, ik in enumerate(precomputed_fs):
        sub_dict = {j: sub_rep[j] for j in ik}
        sub_rep[rep_sym[i]] = sub_rep[rep_sym[i]].xreplace(sub_dict)

    return jacobian[0].xreplace(sub_rep)
