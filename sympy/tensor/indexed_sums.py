from __future__ import print_function, division

import collections
import copy
import itertools

from .indexed import (IndexedBase, Indexed, Idx, IndexException,
                      DeltaIndexedBase)

from sympy.core import (Number, sympify, S, Add, Mul, Pow, Wild, Integer,
                        preorder_traversal, Symbol, Tuple)
from sympy.core.cache import cacheit
from sympy.core.compatibility import range, default_sort_key
from sympy.core.function import Function, expand
from sympy.simplify import simplify
from sympy.utilities.decorator import doctest_depends_on


# TODO: Implement IndexedSum, allowing explicit specification of indices that
# should/shouldn't be summed?


_MAX_NUM_SIMPLIFY_ITERATIONS = 3
_MAX_INDEX_COUNT_TO_PERMUTE = 5


class IndexConformanceException(Exception):
    pass


@doctest_depends_on(modules=('numpy',))
class EinsteinSum(Function):
    """Class to enable implicit summation of repeated indices.

    Simply enclose any expression involving ``Indexed`` objects in
    ``EinsteinSum`` to invoke a form of the `Einstein summation convention`_,
    where repeating an index implies a sum over that index. By default the
    "outer" indices, i.e. non-summation indices, of the resulting expression are
    ordered lexicographically. To override this, specify the desired ordering
    using the ``outer`` keyword argument to the constructor.

    Once an ``EinsteinSum`` expression is created, its index contraction
    structure can be determined using the ``index_structure`` property.

    Note that co/contravariant indices are not supported using ``EinsteinSum``;
    see ``sympy.tensor.tensor``.

    .. _Einstein summation convention:
        http://en.wikipedia.org/wiki/Einstein_notation

    Examples
    ========

    Create an ``EinsteinSum`` that multiplies a vector by a matrix:

    >>> from sympy import symbols, IndexedBase, EinsteinSum
    >>> A, x = symbols('A x', cls=IndexedBase)
    >>> i, j = symbols('i j')
    >>> ein_sum = EinsteinSum(A[i, j] * x[j]); ein_sum
    EinsteinSum(x[j]*A[i, j])

    ``EinsteinSum`` objects can be simplified:

    >>> from sympy import simplify
    >>> k = symbols('k')
    >>> ein_sum = EinsteinSum(A[i, j] * x[i] * x[j] + A[j, k] * x[k] * x[j])
    >>> simplify(ein_sum)
    EinsteinSum(2*x[i]*x[j]*A[i, j])

    Code Generation
    ===============

    ``EinsteinSum`` objects support NumPy code generation using
    ``sympy.lambdify``. For example, a matrix multiplication function can be
    created from symbolic expressions:

    >>> from sympy import symbols, IndexedBase, EinsteinSum, lambdify
    >>> import numpy as np
    >>> i, j, k = symbols('i j k')
    >>> A, B = symbols('A B', cls=IndexedBase)
    >>> m, n = 3, 4
    >>> A_array = np.random.rand(m*n).reshape(m, n)
    >>> B_array = np.random.rand(n*m).reshape(n, m)
    >>> expr = EinsteinSum(A[i, j]*B[j, k])
    >>> func = lambdify([A, B], expr)
    >>> np.allclose(func(A_array, B_array), A_array.dot(B_array))
    True

    To obtain the transpose instead, change the outer index ordering (defaults
    to lexicographic):

    >>> expr = EinsteinSum(A[i, j]*B[j, k], outer=(k, i))
    >>> func_trans = lambdify([A, B], expr)
    >>> np.allclose(func_trans(A_array, B_array), (A_array.dot(B_array)).T)
    True

    Competitive speedwise with default matrix multiply for small arrays (default
    is faster due to optimizations specific to matrix multiplication, e.g.
    `Strassen algorithm`_):

    >>> from timeit import timeit
    >>> m, n = 30, 40
    >>> A_array = np.random.rand(m*n).reshape(m, n)
    >>> B_array = np.random.rand(n*m).reshape(n, m)
    >>> timeit(lambda: func(A_array, B_array), number=10000)  # doctest: +SKIP
    0.1936895000108052
    >>> timeit(lambda: A_array.dot(B_array), number=10000)  # doctest: +SKIP
    0.058242397994035855

    If the ``EinsteinSum`` involves one or more Kronecker delta ``Indexed``
    objects (made via ``DeltaIndexedBase()``), they should *not* be included in
    the ``args`` parameter of ``lambdify``. Instead, they should have ``Idx``
    instances as indices, each possessing numerical shape information. This is
    how ``lambdify`` knows how to automatically construct corresponding identity
    matrices.

    >>> from sympy import Idx, DeltaIndexedBase
    >>> N = 3
    >>> i, j, k, l = symbols('i j k l', cls=Idx, range=N)
    >>> delta = DeltaIndexedBase()
    >>> expr = EinsteinSum(delta[i, j] * delta[k, l])
    >>> func = lambdify([], expr)
    >>> np.allclose(func(), np.einsum('kl,ij->ijkl', np.eye(N), np.eye(N)))
    True

    .. _Strassen algorithm:
        http://en.wikipedia.org/wiki/Strassen_algorithm

    """
    def __new__(cls, arg, outer=None):
        if isinstance(arg, EinsteinSum):
            return EinsteinSum(arg.expr, outer=outer)

        index_structure = cls._get_index_structure(arg)
        computed_outer = index_structure['outer']
        if outer and tuple(outer) != tuple(computed_outer):
            if (
                sorted(list(outer), key=default_sort_key)
                != sorted(list(computed_outer), key=default_sort_key)
            ):
                msg = ("Supplied outer indices, {0!s}, should be a permutation "
                       "of computed outer indices, {1!s}, but they aren't")
                raise ValueError(msg.format(outer, computed_outer))
            obj = Function.__new__(cls, sympify(arg), Tuple(*outer))
        else:
            obj = Function.__new__(cls, sympify(arg))
            outer = computed_outer

        obj._index_structure = index_structure
        obj._index_structure['outer'] = outer

        return obj

    @property
    def expr(self):
        return self.args[0]

    @cacheit
    def _lambda_str(self):
        """Print out a NumPy version of an ``EinsteinSum``

        Helper method for ``lambdify()``.

        The constructed string uses ``numpy.einsum`` to evaluate tensor
        expressions numerically. Using ``Idx`` instances with numerical upper
        bounds as outer indices is required if any ``DeltaIndexedBase``
        instances are present (so that corresponding identity matrices can be
        constructed automatically).

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy import IndexedBase, EinsteinSum, DeltaIndexedBase
        >>> import numpy as np
        >>> i, j, k = symbols('i j k')
        >>> Q, A, x = symbols('Q A x', cls=IndexedBase)
        >>> delta = DeltaIndexedBase()
        >>> expr = EinsteinSum(Q[i, j, k]*A[i, j]*x[k] + delta[i, j]*x[i]*x[j])
        >>> expr._lambda_str()
        "einsum('i,i->', x, x) + einsum('i,jk,jki->', x, A, Q)"

        """
        outer_indices = self.index_structure['outer']

        ein_sum = self._eval_simplify()

        args = []
        for arg in preorder_traversal(ein_sum):
            if (
                isinstance(arg, IndexedBase)
                and not isinstance(arg, DeltaIndexedBase)
            ):
                args.append(arg)

        # Produce a component of a string that will be passed to numpy.einsum.
        # This component specifies the indices of the resultant tensor.
        outer_string = '->{0}'.format(
            ''.join([str(index) for index in outer_indices]))

        # Arrays args will be organized below as [arrays from args, identity
        # matrices]. array_pos corresponds to an array's position in this list.
        # E.g. if args = [a, b, c] and there is one Kronecker delta tensor, then
        # the array_pos of b is 1, the delta tensor's array pos is 3, etc.

        # Make a map from each argument (IndexedBase instance) to its
        # corresponding array_pos. E.g. if args = [a, b, c], this would be
        # {a: 0, b: 1, c: 2}.
        arg_to_array_pos = dict([(arg, pos) for pos, arg in enumerate(args)])

        # Find and organize (by dimension) all DeltaIndexedBase instances.
        delta_dims = set()
        for arg in preorder_traversal(ein_sum):
            if (
                isinstance(arg, Indexed)
                and isinstance(arg.base, DeltaIndexedBase)
            ):
                if (
                    not arg.shape or not isinstance(arg.shape[0], Number)
                    or not isinstance(arg.shape[1], Number)
                ):
                    msg = ("Delta tensors appearing in EinsteinSum must have "
                           "numerical shapes: {0!s}")
                    raise ValueError(msg.format(arg))
                delta_dims.add(arg.shape[0])
        delta_dims = sorted(list(delta_dims))
        delta_args = ["eye({dim!s})".format(dim=dim) for dim in delta_dims]
        num_args = len(args)
        # Make a map from each delta tensor's dimension to its corresponding
        # array_pos. E.g. if args = [a, b, c] and there is just one delta tensor
        # with shape = (2, 2), this would be {2: 3}.
        delta_dim_to_array_pos = dict([
            (dim, num_args + pos) for pos, dim in enumerate(delta_dims)])

        prefactors, decomp = ein_sum._decompose(ein_sum.expr)

        # Loop through EinsteinSum's monomial terms. For each, produce 1) an
        # ordering of arrays to later pass to np.einsum (e.g. [1, 1, 0]
        # meaning pass array 1, then array 1 again, then array 0), and 2) a
        # contraction specification string for np.einsum (e.g. "i,k,ij->jk").
        einsum_term_strings = []
        term_array_orders = []
        for term_decomp in decomp:
            # Loop through multiplicative factors in each monomial.
            einsum_factor_strings = []
            term_array_order = []
            for factor in term_decomp:
                assert isinstance(factor, Indexed)
                einsum_factor_strings.append(
                    "".join([str(index) for index in factor.indices]))
                if isinstance(factor.base, DeltaIndexedBase):
                    term_array_order.append(
                        delta_dim_to_array_pos[factor.shape[0]])
                else:
                    term_array_order.append(arg_to_array_pos[factor.base])
            einsum_term_strings.append(
                ",".join(einsum_factor_strings) + outer_string)
            term_array_orders.append(term_array_order)

        # Function to evaluate each monomial.
        def make_term_string(prefactor, args, order, einsum_string):
            args_string = ", ".join([str(args[pos]) for pos in order])
            prefactor_string = ""
            if prefactor != 1:
                prefactor_string = "{0!s} * ".format(prefactor)
            term_string = ("{pref}einsum('{es}', {args})")
            return term_string.format(pref=prefactor_string,
                                      es=einsum_string,
                                      args=args_string)

        term_strings = []
        # Append identity matrices to array arguments.
        args = list(args) + delta_args
        iterable = zip(prefactors, term_array_orders, einsum_term_strings)
        for prefactor, arg_order, einsum_string in iterable:
            term_strings.append(
                make_term_string(prefactor, args, arg_order, einsum_string))
        return " + ".join(term_strings)

    @property
    def index_structure(self):
        """Determine an ``EinsteinSum``'s index structure.

        Returns
        =======

        dict : with the following elements:
            ``'monomial_list'``
                A list of the EinsteinSum's monomial terms. Sum of these should
                be equivalent to ``self.expr``.
            ``'outer'``
                A list of the EinsteinSum's outer indices.
            ``'inner_list'``
                A list of lists, each enumerating a corresponding monomial's
                inner indices.

        Raises
        ======

        IndexConformanceException
            If two monomials in a sum have differing outer indices, or if an
            index is repeated more than once in a monomial.
        IndexException
            If an ``IndexedBase`` appears without indices.

        Notes
        =====

        Inner indices are those that are implicitly summed according to the
        summation convention; outer indices are those that are not summed and
        therefore constitute the indices of the resultant tensor.

        Examples
        ========

        >>> from sympy import symbols, IndexedBase, EinsteinSum
        >>> i, j, k = symbols('i j k')
        >>> A, B = symbols('A B', cls=IndexedBase)
        >>> ein_sum = EinsteinSum(A[i, j]*B[j, k] + 2*A[i, k])
        >>> structure = ein_sum.index_structure
        >>> structure['monomial_list']
        [A[i, j]*B[j, k], 2*A[i, k]]
        >>> structure['outer']
        [i, k]
        >>> structure['inner_list']
        [[j], []]

        """
        # Make a copy so that EinsteinSum is immutable. (Otherwise user could
        # mutate returned dictionary.)
        return copy.deepcopy(self._index_structure)

    @classmethod
    def _get_index_structure(cls, expr):
        """Compute index structure (see index_structure property)"""
        term_list = []
        inner_list = []
        for term, outer, inner in cls._get_index_structure_recurse(expr):
            term_list.append(term)
            inner_list.append(inner)

        term_list, inner_list = zip(*sorted(zip(term_list, inner_list),
                                            key=default_sort_key))
        term_list = list(term_list)
        inner_list = list(inner_list)

        return {'monomial_list': term_list, 'outer': outer,
                'inner_list': inner_list}

    @classmethod
    def _get_monomial_indices(cls, monomial):
        """List all indices appearing in a monomial with repetition.

        Helper for _get_index_structure_recurse.

        """
        if isinstance(monomial, Mul):
            # Build up list of indices of each term.
            indices_list = []
            for arg in monomial.args:
                indices_list.extend(cls._get_monomial_indices(arg))
        elif isinstance(monomial, Pow):
            # Get indices of base and duplicate each p times, where p is the
            # power.
            base, power = monomial.args
            indices_list = cls._get_monomial_indices(base)
            if indices_list and not isinstance(power, Integer):
                msg = "Only integral powers are allowed, not: {0!s}"
                raise ValueError(msg.format(power))
            if indices_list:
                indices_list *= power
        elif isinstance(monomial, Indexed):
            indices_list = list(monomial.indices)
        elif isinstance(monomial, Add):
            if monomial.has(Indexed):
                raise ValueError("Not a monomial: {0!s}".format(monomial))
            indices_list = []
        else:
            # Assume everything else has no indices.
            indices_list = []

        indices_list.sort(key=default_sort_key)
        return indices_list

    @classmethod
    def _classify_indices(cls, index_list):
        """Classify indices as inner or outer.

        Helper for _get_index_structure_recurse.

        Parameters
        ==========

        ``index_list`` : ``list``
            A list of indices appearing in a monomial with repetition.

        Returns
        =======

        tuple : ``set`` of outer indices, ``set`` of inner indices

        """
        outer = set()
        inner = set()
        index_counts = {}
        for index in index_list:
            if index in index_counts:
                index_counts[index] += 1
            else:
                index_counts[index] = 1

        for index in index_counts:
            if index_counts[index] == 1:
                outer.add(index)
            elif index_counts[index] == 2:
                inner.add(index)
            else:
                msg = "> 2 occurrences of index {0!s}: {1!s}"
                raise IndexConformanceException(msg.format(index, index_list))
        return outer, inner

    @classmethod
    def _get_index_structure_recurse(cls, expr):
        """Return a list of tuples describing a tensor expression.

        Helper for _get_index_structure.

        """
        expr = expand(expr)
        if isinstance(expr, Add):
            index_structure = []
            for arg in expr.args:
                index_structure.extend(cls._get_index_structure_recurse(arg))

            # Get outer indices of each term and ensure they're all the same.
            prev_outer = None
            for _, outer, inner in index_structure:
                if prev_outer is not None and outer != prev_outer:
                    msg = "Inconsistent index structure across sum: {0!s}"
                    raise IndexConformanceException(msg.format(expr))
                prev_outer = outer

            return index_structure
        else:
            index_list = cls._get_monomial_indices(expr)
            outer, inner = cls._classify_indices(index_list)
            return [(expr,
                     sorted(list(outer), key=default_sort_key),
                     sorted(list(inner), key=default_sort_key))]

    @classmethod
    def _decompose(cls, expr):
        """Decompose an ``EinsteinSum`` into constituent terms and factors.

        Returns
        =======

        tuple : prefactors (``list``), decomp (``list`` of ``list``s)
            Each sublist in ``decomp`` corresponds to a monomial, and enumerates
            constituent ``Indexed`` object factors contained in each monomial.
            ``prefactors`` is a list of prefactors (everything other than
            Indexed objects), one for each monomial.

        Helper method for _lambda_str.

        """
        expr = expand(expr)
        if isinstance(expr, Add):
            args = expr.args
        else:
            args = [expr]

        decomp = []
        prefactors = []
        for arg in args:
            prefactor, monomial_decomp = cls._decompose_monomial(arg)
            decomp.append(monomial_decomp)
            prefactors.append(prefactor)

        return prefactors, decomp

    @classmethod
    def _decompose_monomial(cls, monomial):
        """Helper for ``_decompose()``: decompose a single monomial term."""
        prefactor = S.One
        monomial_decomp = []
        if isinstance(monomial, Mul):
            # Tabulate terms in Mul.
            for arg in monomial.args:
                arg_prefactor, arg_decomp = cls._decompose_monomial(arg)
                prefactor *= arg_prefactor
                monomial_decomp.extend(arg_decomp)
        elif isinstance(monomial, Pow):
            # Get decomp of base, then modify accounting for power.
            base, power = monomial.args
            base_prefactor, base_decomp = cls._decompose_monomial(base)
            if base_decomp and not isinstance(power, Integer):
                msg = "Only integral powers are allowed, not: {0!s}"
                raise ValueError(msg.format(power))
            prefactor *= base_prefactor ** power
            monomial_decomp = base_decomp * power
        elif isinstance(monomial, Indexed):
            monomial_decomp.append(monomial)
        elif isinstance(monomial, Add):
            if monomial.has(Indexed):
                raise ValueError("Not a monomial: {0!s}".format(monomial))
            prefactor *= monomial
        elif isinstance(monomial, IndexedBase):
                msg = ("No IndexedBase objects should be present w/o indices: "
                       "{0!s}")
                raise IndexException(msg.format(monomial))
        else:
            # Assume everything else has no Indexed objects.
            prefactor *= monomial

        return prefactor, monomial_decomp

    def _eval_derivative(self, wrt):
        return type(self)(self.expr.diff(wrt)).simplify_deltas()

    def _eval_simplify(self, **kwargs):
        old = self
        for iteration in range(_MAX_NUM_SIMPLIFY_ITERATIONS):
            new = old.simplify_deltas()
            new = new.canonicalize_inner()
            new = type(self)(simplify(new.expr, **kwargs))
            if new == old:
                break
            old = new
        return new

    def canonicalize_inner(self):
        """Rename and reorder inner (summation) indices into a canonical form.

        Since inner indices are "dummies", this only cosmetically changes the
        expression. Renaming only uses existing inner indices -- no new indices
        are created. Note that all inner indices must be interchangeable in
        order to be rearranged: they must all be of the same type and must all
        have the same range information if of type ``Idx``.

        Examples
        ========

        >>> from sympy import EinsteinSum, IndexedBase, symbols
        >>> i, j, l, m = symbols('i j l m')
        >>> L, h, Q = symbols('L h Q', cls=IndexedBase)
        >>> expr = EinsteinSum(L[i, j] * h[i] * h[j] + L[j, i] * h[i] * h[j])
        >>> expr.canonicalize_inner()
        EinsteinSum(2*h[i]*h[j]*L[i, j])
        >>> expr = EinsteinSum(L[i, j] * h[j] + Q[i, l, m] * h[l] * h[m])
        >>> expr.canonicalize_inner()
        EinsteinSum(h[j]*h[l]*Q[i, j, l] + h[j]*L[i, j])

        """
        structure = self.index_structure
        inner_list = structure["inner_list"]
        monomial_list = structure["monomial_list"]
        new_monomial_list = monomial_list[:]

        # Construct a canonical list of inner indices.
        canonical_inner = set().union(*inner_list)
        canonical_inner = sorted(list(canonical_inner), key=default_sort_key)
        if not canonical_inner:
            return type(self)(self.expr)

        # Ensure all inner indices are interchangeable. Otherwise don't attempt
        # to simplify.
        inner_type = type(canonical_inner[0])
        for inner in canonical_inner:
            if not isinstance(inner, inner_type):
                return type(self)(self.expr)
            if inner_type == Idx and not (
                inner.lower == canonical_inner[0].lower and
                inner.upper == canonical_inner[0].upper
            ):
                return type(self)(self.expr)

        # Construct a list of monomials using a minimal set of inner indices.
        for pos in range(len(inner_list)):
            inner = sorted(list(inner_list[pos]), key=default_sort_key)
            sub = zip(inner, canonical_inner[:len(inner)])
            new_monomial_list[pos] = monomial_list[pos].subs(sub,
                                                             simultaneous=True)

        # Iteratively reorder each monomial's inner indices to minimize hash.
        for pos, monomial in enumerate(new_monomial_list):
            num_inner = len(inner_list[pos])
            if num_inner > _MAX_INDEX_COUNT_TO_PERMUTE:
                continue
            inner = canonical_inner[:num_inner]

            monomial_permutations = [
                monomial.subs(zip(inner, new_inner), simultaneous=True)
                for new_inner in itertools.permutations(inner)]
            new_monomial_list[pos] = min(monomial_permutations,
                                         key=default_sort_key)

        return type(self)(sum(new_monomial_list))

    def simplify_deltas(self):
        """Simplify trivial contractions involving Kronecker deltas.

        Returns
        =======

        EinsteinSum : An equivalent, simplified ``EinsteinSum`` object.

        Examples
        ========

        >>> from sympy import EinsteinSum, IndexedBase, DeltaIndexedBase
        >>> from sympy import symbols, S
        >>> i, j = symbols('i j')
        >>> v = IndexedBase('v')
        >>> delta = DeltaIndexedBase()
        >>> expr = EinsteinSum(v[i] * delta[i, j])
        >>> expr.simplify_deltas()
        EinsteinSum(v[j])

        """
        expr = self._perform_delta_contractions(self.expr)
        return type(self)(expr)

    @classmethod
    def _perform_delta_contractions(cls, expr):
        """Evaluate trivial Kronecker delta contractions.

        Helper method for ``simplify_deltas()``.

        """
        expr = expand(expr)
        if isinstance(expr, Add):
            new_args = [cls._perform_delta_contractions(arg)
                        for arg in expr.args]
            return Add(*new_args)
        elif isinstance(expr, Mul):
            # First determine which factors have which indices: form a map from
            # index to position in args where an Indexed object resides
            # possessing that index.
            args = list(expr.args)
            index_to_arg_ids = collections.defaultdict(list)
            for arg_id, arg in enumerate(args):
                if not isinstance(arg, Indexed):
                    # Might be a Pow (handled below). Pow objects can be handled
                    # separately (without e.g. incorporating into the map we're
                    # building).
                    args[arg_id] = cls._perform_delta_contractions(arg)
                    continue
                for index in arg.indices:
                    index_to_arg_ids[index].append(arg_id)

            # Build up a list of argument positions ("ids") where deltas we have
            # simplified reside, so we can get rid of them later.
            arg_ids_to_delete = set()

            def eliminate_delta(dummy_index, delta_id, other_id):
                if delta_id in arg_ids_to_delete:
                    # Already handled this delta.
                    return
                if delta_id == other_id:
                    # Internal delta contraction. Ignore.
                    return
                delta = args[delta_id]
                other_tensor = args[other_id]
                other_index = (delta.indices[0]
                               if delta.indices[1] == dummy_index
                               else delta.indices[1])
                args[other_id] = other_tensor.subs(dummy_index, other_index)
                arg_ids_to_delete.add(delta_id)

            # Loop through indices we have tabulated and handle each one.
            for index in sorted(index_to_arg_ids.keys(),
                                key=default_sort_key):
                tensor_ids = index_to_arg_ids[index]
                if len(tensor_ids) == 1:
                    # No repeated dummy index. Ignore.
                    continue
                elif len(tensor_ids) > 2:
                    # Inconsistent with Einstein convention.
                    msg = "Index {0!s} repeated > 2 times in {1!s}"
                    raise IndexConformanceException(msg.format(index, expr))

                # Handling a dummy index. See if either tensor possessing the
                # dummy index is a Kronecker delta; if so, perform contraction
                # and eliminate it.
                if isinstance(args[tensor_ids[0]].base, DeltaIndexedBase):
                    eliminate_delta(index, tensor_ids[0], tensor_ids[1])
                    continue
                elif isinstance(args[tensor_ids[1]].base, DeltaIndexedBase):
                    eliminate_delta(index, tensor_ids[1], tensor_ids[0])
                    continue

            new_args = [args[arg_id] for arg_id in range(len(args))
                        if arg_id not in arg_ids_to_delete]
            return Mul(*new_args)
        elif (isinstance(expr, Pow) and isinstance(expr.args[0], Indexed)
              and isinstance(expr.args[0].base, DeltaIndexedBase)):
            delta_base = expr.args[0].base
            wild1 = Wild("wild1")
            wild2 = Wild("wild2")
            return expr.replace(delta_base[wild1, wild2] ** 2,
                                delta_base[wild1, wild1])
        else:
            return expr


def get_indices(expr):
    """Determine the non-summation indices present in a tensor expression.

    Note that unless passed an ``EinsteinSum`` object, this function does NOT
    assume implicit sums over repeated indices.

    Returns
    =======

    set : Set of resultant "outer" indices appearing in ``expr``.

    Examples
    ========

    >>> from sympy import IndexedBase, EinsteinSum, symbols, get_indices
    >>> A, x = symbols('A x', cls=IndexedBase)
    >>> i, j = symbols('i j')
    >>> get_indices(A[i, j] * x[j])
    set([i, j])
    >>> get_indices(A[i, j] + x[j])
    set([i, j])
    >>> get_indices(2*EinsteinSum(A[i, j] * x[j]))
    set([i])

    """
    if isinstance(expr, EinsteinSum):
        return set(expr.index_structure['outer'])
    elif isinstance(expr, Indexed):
        # Get all indices.
        indices = set()
        for index_expr in expr.indices:
            for arg in preorder_traversal(index_expr):
                if isinstance(arg, Symbol) or isinstance(arg, Idx):
                    indices.add(arg)

        # However, Idx objects each contain constituent Symbols. Ignore these.
        idx_symbols = set()
        for index in indices:
            if isinstance(index, Idx):
                for arg in preorder_traversal(index):
                    if isinstance(arg, Symbol):
                        idx_symbols.add(arg)

        return indices - idx_symbols
    else:
        return set().union(*[get_indices(arg) for arg in expr.args])


def _index_has_shape(index):
    return isinstance(index, Idx) and isinstance(index.upper, Number)
