from __future__ import print_function, division

import collections

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

    >>> from sympy import symbols, IndexedBase, EinsteinSum
    >>> A, x = symbols('A x', cls=IndexedBase)
    >>> i, j = symbols('i j')
    >>> ein_sum = EinsteinSum(A[i, j] * x[j]); ein_sum
    EinsteinSum(x[j]*A[i, j])

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
        Any shape information in the component ``Indexed`` objects is ignored
        by this method.

        Note that there can be no non-numerical elements of ``EinsteinSum``
        objects that call ``numpify`` other than ``Indexed`` objects -- making
        functions that also take other symbolic objects as parameters is not yet
        supported.

        Parameters
        ==========

        args : ``list`` (or iterable)
            List of all constituent ``IndexedBase`` instances specifying the
            argument NumPy arrays the constructed function will take, in order.
        outer_indices : ``list`` (or iterable)
            Index order of output of constructed function.
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

        Competitive speedwise with default matrix multiply:

        >>> from timeit import timeit
        >>> m, n = 300, 400
        >>> A_array = np.random.rand(m*n).reshape(m, n)
        >>> B_array = np.random.rand(n*m).reshape(n, m)
        >>> timeit(func(A_array, B_array), number=10)  # doctest: +SKIP
        0.21692680500564165
        >>> timeit(A_array.dot(B_array), number=10)  # doctest: +SKIP
        0.3139198549906723

        Kronecker deltas
        ----------------

        Currently if the ``EinsteinSum`` involves a Kronecker delta ``Indexed``
        object (made via ``DeltaIndexedBase()``), then it must be referenced
        in ``args``, and an appropriately shaped NumPy identity matrix must be
        supplied explicitly as an argument to the constructed function. If
        multiple different Kronecker deltas are required, they should be given
        different explicit labels, otherwise they will be determined to be the
        same. In the future this should be fixed so that identity matrices are
        supplied automatically.

        >>> from sympy import DeltaIndexedBase
        >>> v = symbols('v', cls=IndexedBase)
        >>> delta = DeltaIndexedBase()
        >>> m = 3
        >>> v_array = np.random.rand(m)
        >>> delta_array = np.eye(m)
        >>> expr = EinsteinSum(delta[i, j]*v[j])
        >>> func = expr.numpify([delta, v], [i])
        >>> np.allclose(func(delta_array, v_array), v_array)
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

        self._check_numpify_args(args, outer_indices)

        # Produce a string component that we will pass to numpy.einsum below.
        # This component specifies the indices of the resultant tensor.
        outer_string = '->{0}'.format(
            ''.join([str(index) for index in outer_indices]))

        # Make a map from argument (IndexedBase instance) to argument's
        # position in args. E.g. if args = [a, b, c], would be
        # {a: 0, b: 1, c: 2}.
        base_to_arg_position = dict([(arg, arg_id)
                                     for arg_id, arg in enumerate(args)])

        def arrange_args(args, order):
            """Permute args according to order."""
            return [args[index] for index in order]

        prefactors, decomp = self._decompose(self.expr)
        prefactors = [dtype(prefactor) for prefactor in prefactors]

        # Loop through EinsteinSum's monomial terms. For each, produce 1) an
        # ordering of arguments to later pass to np.einsum (e.g. [1, 1, 0]
        # meaning pass arg 1, then arg 1 again, then arg 0), and 2) a
        # contraction specification string for np.einsum (e.g. "i,k,ij->jk").
        einsum_term_strings = []
        term_arg_orders = []
        for term_decomp in decomp:
            # Loop through multiplicative factors in each monomial.
            term_arg_orders.append([base_to_arg_position[factor.base]
                                    for factor in term_decomp])
            einsum_factor_strings = []
            for factor in term_decomp:
                assert isinstance(factor, Indexed)
                einsum_factor_strings.append(
                    "".join([str(index) for index in factor.indices]))
            einsum_term_strings.append(
                ",".join(einsum_factor_strings) + outer_string)

        # Create functions to evaluate each monomial. Need to have this factory
        # function to avoid using loop variable in closure (classic Python
        # mistake).
        def make_term_func(einsum_string):
            def term_func(*args):
                return np.einsum(einsum_string, *args, dtype=dtype)
            term_func.einsum_string = einsum_string
            return term_func
        term_funcs = []
        for einsum_string in einsum_term_strings:
            term_funcs.append(make_term_func(einsum_string))

        # Create function that evaluates entire EinsteinSum = sum of monomials.
        def func(*args):
            result = 0.
            iterable = zip(term_funcs, term_arg_orders, prefactors)
            for term_func, arg_order, prefactor in iterable:
                arranged_args = arrange_args(args, arg_order)
                result += prefactor * term_func(*arranged_args)
            return result
        func.__doc__ = """Tensor function created with EinsteinSum.numpify.

        Parameters
        ==========

        {args} : each a ``numpy.ndarray`` or object that can be converted into
            one via ``numpy.array``.

        Returns
        ========

        ``np.ndarray``
            Result of:
                {expr!s}
            With index ordering:
                {outer!s}

        """.format(
            expr=self.expr, args=", ".join([str(arg) for arg in args]),
            outer=outer_indices)
        func.term_funcs = term_funcs
        return func

    def _check_numpify_args(self, args, outer_indices):
        """Various checks to ensure arguments to numpify are valid."""
        expr = self.expr

        # Check that args are all IndexedBase objects.
        for arg in args:
            if not isinstance(arg, IndexedBase):
                raise TypeError("args should only have IndexedBase objects")

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
                indexed_objects.add(arg)
                indexed_base_objects.add(arg.base)
                all_indices |= set(arg.indices)

        # Check that all IndexedBase objects appearing in the EinsteinSum have
        # been mentioned in args.
        if not set(args) == indexed_base_objects:
            msg = ("Supplied args list, {0}, should list all IndexedBase "
                   "objects appearing in EinsteinSum, {1}, but it doesn't")
            raise ValueError(msg.format(args, indexed_base_objects))

        # Ensure there are no free symbols that aren't indices.
        substitutions = zip(indexed_objects,
                            [S.One for _ in range(len(indexed_objects))])
        if not isinstance(expr.subs(substitutions), Number):
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
        ein_sum = self.simplify_deltas()
        expr = simplify(ein_sum.expr, **kwargs)
        return type(self)(expr)

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
