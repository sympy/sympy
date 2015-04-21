from __future__ import print_function, division

import collections
import itertools

from .indexed import (IndexedBase, Indexed, Idx, IndexException,
                      DeltaIndexedBase)

from sympy.core import (Symbol, Number, sympify, S, Add, Mul, Pow,
                        Wild, Integer, preorder_traversal)
from sympy.core.cache import cacheit
from sympy.core.compatibility import range, default_sort_key
from sympy.core.function import Function, expand
from sympy.printing.repr import srepr
from sympy.simplify import simplify
from sympy.utilities.decorator import doctest_depends_on


# TODO: Implement IndexedSum, allowing explicit specification of indices that
# should/shouldn't be summed?


_MAX_NUM_SIMPLIFY_ITERATIONS = 3
_MAX_INDEX_COUNT_TO_PERMUTE = 5


class IndexConformanceException(Exception):
    pass


class EinsteinSum(Function):
    """Class to enable implicit summation of repeated indices.

    Simply enclose any expression involving ``Indexed`` objects in
    ``EinsteinSum`` to invoke a form of the `Einstein summation convention`_,
    where repeating an index implies a sum over that index.

    Once an ``EinsteinSum`` expression is created, its index contraction
    structure can be determined using the ``index_structure`` property.

    An ``EinsteinSum`` can be converted into a function accepting NumPy arrays
    using the ``numpify`` method.

    Note that co/contravariant indices are not supported using ``EinsteinSum``;
    see ``sympy.tensor.tensor``.

    .. _Einstein summation convention:
        http://en.wikipedia.org/wiki/Einstein_notation

    Examples
    ========

    Create an ``EinsteinSum`` that effects matrix multiplication:

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

    """
    def __new__(cls, arg):
        if isinstance(arg, EinsteinSum):
            return arg
        obj = Function.__new__(cls, sympify(arg))
        # Compute this now to ensure valid index structure. (It's decorated with
        # @cacheit.)
        obj.index_structure
        return obj

    @property
    def expr(self):
        return self.args[0]

    @doctest_depends_on(modules=('numpy',))
    def numpify(self, args, outer_indices, dtype=float):
        """Make a function accepting NumPy arrays out of an ``EinsteinSum``.

        The constructed function is then simply a wrapper for ``numpy.einsum``.
        Using ``Idx`` instances with numerical upper bounds as outer indices is
        required if any ``DeltaIndexedBase`` instances are present (so that
        corresponding identity matrices can be constructed automatically).

        Note that there can be no non-numerical elements of ``EinsteinSum``
        objects that call ``numpify`` other than ``Indexed`` objects -- making
        functions that also take other symbolic objects as parameters is not yet
        supported.

        Parameters
        ==========

        args : ``list`` (or iterable)
            List of all constituent ``IndexedBase`` instances specifying the
            argument NumPy arrays the constructed function will take, in order.
            This does not include ``DeltaIndexedBase`` instances, which are
            handled automatically.
        outer_indices : ``list`` (or iterable)
            List of outer (i.e. non-summation) indices present, which specifies
            the index order of the output of the constructed function.
        dtype : a number type (either ``float`` or ``complex``)
            Number type constructed function should use for computation.

        Examples
        ========

        Matrix multiplication
        ---------------------

        >>> from sympy import symbols, IndexedBase, EinsteinSum
        >>> import numpy as np
        >>> i, j, k = symbols('i j k')
        >>> A, B = symbols('A B', cls=IndexedBase)
        >>> m, n = 3, 4
        >>> A_array = np.random.rand(m*n).reshape(m, n)
        >>> B_array = np.random.rand(n*m).reshape(n, m)
        >>> expr = EinsteinSum(A[i, j]*B[j, k])
        >>> func = expr.numpify([A, B], [i, k])
        >>> np.allclose(func(A_array, B_array), A_array.dot(B_array))
        True

        Competitive speedwise with default matrix multiply for large arrays:

        >>> from timeit import timeit
        >>> m, n = 3, 4
        >>> A_array = np.random.rand(m*n).reshape(m, n)
        >>> B_array = np.random.rand(n*m).reshape(n, m)
        >>> timeit(lambda: func(A_array, B_array), number=10)  # doctest: +SKIP
        0.21692680500564165
        >>> timeit(lambda: A_array.dot(B_array), number=10)  # doctest: +SKIP
        0.3139198549906723

        Kronecker deltas
        ----------------

        If the ``EinsteinSum`` involves one or more Kronecker delta ``Indexed``
        object (made via ``DeltaIndexedBase()``), they should *not* be included
        in ``args``. Instead, they should have ``Idx`` instances as indices,
        each possessing numerical shape information. This is how ``numpify``
        knows how to construct corresponding identity matrices.

        >>> from sympy import Idx, DeltaIndexedBase
        >>> N = 3
        >>> i, j, k, l = symbols('i j k l', cls=Idx, range=N)
        >>> delta = DeltaIndexedBase()
        >>> expr = EinsteinSum(delta[i, j] * delta[k, l])
        >>> func = expr.numpify([], [i, j, k, l])
        >>> np.allclose(func(), np.einsum('kl,ij->ijkl', np.eye(N), np.eye(N)))
        True

        See Also
        ========

        sympy.utilities.autowrap

        """
        try:
            import numpy as np
        except ImportError as e:
            msg = ("NumPy is required for EinsteinSum.numpify but it couldn't "
                   "be imported: {0!s}".format(e))
            raise ImportError(msg)

        self = self._eval_simplify()
        self._check_numpify_args(args, outer_indices)

        # Produce a string component that we will pass to numpy.einsum below.
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
        for arg in preorder_traversal(self):
            if (
                isinstance(arg, Indexed)
                and isinstance(arg.base, DeltaIndexedBase)
            ):
                delta_dims.add(arg.shape[0])
        delta_dims = sorted(list(delta_dims))
        delta_arrays = [np.eye(dim) for dim in delta_dims]
        num_args = len(args)
        # Make a map from each delta tensor's dimension to its corresponding
        # array_pos. E.g. if args = [a, b, c] and there is just one delta tensor
        # with shape = (2, 2), this would be {2: 3}.
        delta_dim_to_array_pos = dict([
            (dim, num_args + pos) for pos, dim in enumerate(delta_dims)])

        prefactors, decomp = self._decompose(self.expr)
        prefactors = [dtype(prefactor) for prefactor in prefactors]

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
        def term_func(prefactor, arrays, order, einsum_string):
            return prefactor * np.einsum(einsum_string,
                                         *[arrays[pos] for pos in order])

        # Create function that evaluates entire EinsteinSum = sum of monomials.
        def func(*arrays):
            # Append identity matrices to array arguments.
            arrays = list(arrays) + delta_arrays
            result = 0.
            iterable = zip(prefactors, term_array_orders, einsum_term_strings)
            for prefactor, array_order, einsum_string in iterable:
                result += term_func(
                    prefactor, arrays, array_order, einsum_string)
            return result
        func.__doc__ = """Tensor function created with EinsteinSum.numpify.

        Parameters
        ==========

        {arrays} : each a ``numpy.ndarray`` or object that can be converted into
            one via ``numpy.array``.

        Returns
        ========

        ``np.ndarray``
            Result of:
                {expr!s}
            With index ordering:
                {outer!s}

        """.format(
            expr=self.expr, arrays=", ".join([str(arg) for arg in args]),
            outer=outer_indices)
        return func

    def _check_numpify_args(self, args, outer_indices):
        """Various checks to ensure arguments to numpify are valid."""
        expr = self.expr

        # Check that args are all IndexedBase objects. Ensure none are
        # DeltaIndexedBase objects.
        for arg in args:
            if not isinstance(arg, IndexedBase):
                raise TypeError("args should only have IndexedBase objects")
            if isinstance(arg, DeltaIndexedBase):
                raise TypeError(
                    "No need to specify Kronecker deltas as arguments; just "
                    "ensure relevant indices are instances of Idx and have the "
                    "right numerical upper bounds")

        # Check supplied outer indices.
        index_structure = self.index_structure
        actual_outer_indices = set(index_structure['outer'])
        if set(outer_indices) != actual_outer_indices:
            msg = ("Supplied outer_indices = {0} should be a permutation of "
                   "true outer indices = {1} but isn't")
            raise IndexConformanceException(msg.format(outer_indices,
                                                       actual_outer_indices))

        # Gather various objects.
        indexed_base_objects = set()
        all_indices = set()
        indexed_objects = set()
        for arg in preorder_traversal(expr):
            if isinstance(arg, Indexed):
                if not isinstance(arg.base, DeltaIndexedBase):
                    indexed_objects.add(arg)
                    indexed_base_objects.add(arg.base)
                    all_indices |= set(arg.indices)
                elif not all(
                    [_index_has_shape(index) for index in arg.indices]
                ):
                    msg = ("Kronecker delta tensors specified via "
                           "DeltaIndexedBase should have Idx indices with "
                           "numerical ranges: {0}")
                    raise TypeError(msg.format(arg))

        # Check that all IndexedBase objects appearing in the EinsteinSum have
        # been mentioned in args. (Extras in args are okay though.)
        if not set(args) >= indexed_base_objects:
            msg = ("Supplied args list, {0}, should list all IndexedBase "
                   "objects appearing in EinsteinSum, {1}, but it doesn't")
            raise ValueError(msg.format(args, indexed_base_objects))

        # Ensure there are no free symbols that aren't indices.
        replaced_expr = expr.replace(Indexed, lambda *_: S.One)
        if not isinstance(replaced_expr, Number):
            msg = ("Only numbers and Indexed objects are allowed when calling "
                   "EinsteinSum.numpify(): {0}")
            raise ValueError(msg.format(expr))

        # Ensure all indices are symbolic (not e.g. numbers).
        for index in all_indices:
            if not (isinstance(index, Symbol) or isinstance(index, Idx)):
                msg = "Index is neither a Symbol nor an Idx: {0}"
                raise TypeError(msg.format(srepr(index)))

    @property
    @cacheit
    def index_structure(self):
        """Determine an ``EinsteinSum``'s index structure.

        Returns
        =======

        dict : with the following elements:
            ``'monomial_list'``
                A list of the EinsteinSum's monomial terms. Sum of these should
                be equivalent to ``self.expr``.
            ``'outer'``
                A set of the EinsteinSum's outer indices.
            ``'inner_list'``
                A list of sets, each listing a corresponding monomial's inner
                indices.

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
        set([i, k])
        >>> structure['inner_list']
        [set([j]), set()]

        """
        term_list = []
        inner_list = []
        for term, outer, inner in _get_index_structure(self.expr,
                                                       einstein_notation=True):
            term_list.append(term)
            inner_list.append(inner)

        term_list, inner_list = zip(*sorted(zip(term_list, inner_list),
                                            key=default_sort_key))
        term_list = list(term_list)
        inner_list = list(inner_list)

        return {'monomial_list': term_list, 'outer': outer,
                'inner_list': inner_list}

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

        Helper method for numpify.

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


def _get_monomial_indices(monomial):
    """List all indices appearing in a monomial with repetition.

    Helper for _get_index_structure."""
    if isinstance(monomial, Mul):
        # Build up list of indices of each term.
        indices_list = []
        for arg in monomial.args:
            indices_list.extend(_get_monomial_indices(arg))
    elif isinstance(monomial, Pow):
        # Get indices of base and duplicate each p times, where p is the
        # power.
        base, power = monomial.args
        indices_list = _get_monomial_indices(base)
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


def _classify_indices(index_list):
    """Classify indices as inner or outer.

    Helper for _get_index_structure.

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


def _get_index_structure(expr, einstein_notation):
    """Return a list of tuples describing a tensor expression."""
    expr = expand(expr)
    if isinstance(expr, Add):
        index_structure = []
        for arg in expr.args:
            index_structure.extend(_get_index_structure(arg, einstein_notation))

        # Get outer indices of each term and ensure they're all the same.
        prev_outer = None
        for _, outer, inner in index_structure:
            if prev_outer is not None and outer != prev_outer:
                msg = "Inconsistent index structure across sum: {0!s}"
                raise IndexConformanceException(msg.format(expr))
            prev_outer = outer

        return index_structure
    else:
        index_list = _get_monomial_indices(expr)
        if einstein_notation:
            outer, inner = _classify_indices(index_list)
        else:
            outer = set(index_list)
            inner = set()
        return [(expr, outer, inner)]


def get_indices(expr):
    """Determine the outer indices present in a tensor expression.

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
        return expr.index_structure['outer']
    elif isinstance(expr, Indexed):
        return set(expr.indices)
    else:
        return set().union(*[get_indices(arg) for arg in expr.args])


def _index_has_shape(index):
    return isinstance(index, Idx) and isinstance(index.upper, Number)
