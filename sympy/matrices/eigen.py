from __future__ import division, print_function

from types import FunctionType

from mpmath.libmp.libmpf import prec_to_dps

from sympy.core.compatibility import default_sort_key
from sympy.core.logic import fuzzy_and, fuzzy_or
from sympy.core.numbers import Float
from sympy.core.symbol import Dummy
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.polys import roots
from sympy.simplify import nsimplify, simplify as _simplify
from sympy.utilities.exceptions import SymPyDeprecationWarning

from .common import (MatrixError, NonSquareMatrixError,
    NonPositiveDefiniteMatrixError)

from .utilities import _iszero


# This functions is a candidate for caching if it gets implemented for matrices.
def _eigenvals(M, error_when_incomplete=True, dotprodsimp=None, **flags):
    r"""Return eigenvalues using the Berkowitz agorithm to compute
    the characteristic polynomial.

    Parameters
    ==========

    error_when_incomplete : bool, optional
        If it is set to ``True``, it will raise an error if not all
        eigenvalues are computed. This is caused by ``roots`` not returning
        a full list of eigenvalues.

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    simplify : bool or function, optional
        If it is set to ``True``, it attempts to return the most
        simplified form of expressions returned by applying default
        simplification method in every routine.

        If it is set to ``False``, it will skip simplification in this
        particular routine to save computation resources.

        If a function is passed to, it will attempt to apply
        the particular function as simplification method.

    rational : bool, optional
        If it is set to ``True``, every floating point numbers would be
        replaced with rationals before computation. It can solve some
        issues of ``roots`` routine not working well with floats.

    multiple : bool, optional
        If it is set to ``True``, the result will be in the form of a
        list.

        If it is set to ``False``, the result will be in the form of a
        dictionary.

    Returns
    =======

    eigs : list or dict
        Eigenvalues of a matrix. The return format would be specified by
        the key ``multiple``.

    Raises
    ======

    MatrixError
        If not enough roots had got computed.

    NonSquareMatrixError
        If attempted to compute eigenvalues from a non-square matrix.

    Examples
    ========

    >>> from sympy.matrices import Matrix
    >>> M = Matrix(3, 3, [0, 1, 1, 1, 0, 0, 1, 1, 1])
    >>> M.eigenvals()
    {-1: 1, 0: 1, 2: 1}

    See Also
    ========

    MatrixDeterminant.charpoly
    eigenvects

    Notes
    =====

    Eigenvalues of a matrix `A` can be computed by solving a matrix
    equation `\det(A - \lambda I) = 0`
    """

    simplify = flags.get('simplify', False) # Collect simplify flag before popped up, to reuse later in the routine.
    multiple = flags.get('multiple', False) # Collect multiple flag to decide whether return as a dict or list.
    rational = flags.pop('rational', True)

    if not M:
        return {}

    if rational:
        M = M.applyfunc(
            lambda x: nsimplify(x, rational=True) if x.has(Float) else x)

    if M.is_upper or M.is_lower:
        if not M.is_square:
            raise NonSquareMatrixError()

        diagonal_entries = [M[i, i] for i in range(M.rows)]

        if multiple:
            eigs = diagonal_entries

        else:
            eigs = {}

            for diagonal_entry in diagonal_entries:
                if diagonal_entry not in eigs:
                    eigs[diagonal_entry] = 0

                eigs[diagonal_entry] += 1

    else:
        flags.pop('simplify', None)  # pop unsupported flag

        if isinstance(simplify, FunctionType):
            eigs = roots(M.charpoly(x=Dummy('x'), simplify=simplify,
                    dotprodsimp=dotprodsimp), **flags)
        else:
            eigs = roots(M.charpoly(x=Dummy('x'), dotprodsimp=dotprodsimp), **flags)

    # make sure the algebraic multiplicity sums to the
    # size of the matrix
    if error_when_incomplete and (sum(eigs.values()) if
            isinstance(eigs, dict) else len(eigs)) != M.cols:
        raise MatrixError("Could not compute eigenvalues for {}".format(M))

    # Since 'simplify' flag is unsupported in roots()
    # simplify() function will be applied once at the end of the routine.
    if not simplify:
        return eigs
    if not isinstance(simplify, FunctionType):
        simplify = _simplify

    # With 'multiple' flag set true, simplify() will be mapped for the list
    # Otherwise, simplify() will be mapped for the keys of the dictionary
    if not multiple:
        return {simplify(key): value for key, value in eigs.items()}
    else:
        return [simplify(value) for value in eigs]


# This functions is a candidate for caching if it gets implemented for matrices.
def _eigenvects(M, error_when_incomplete=True, iszerofunc=_iszero,
        dotprodsimp=None, **flags):
    """Return list of triples (eigenval, multiplicity, eigenspace).

    Parameters
    ==========

    error_when_incomplete : bool, optional
        Raise an error when not all eigenvalues are computed. This is
        caused by ``roots`` not returning a full list of eigenvalues.

    iszerofunc : function, optional
        Specifies a zero testing function to be used in ``rref``.

        Default value is ``_iszero``, which uses SymPy's naive and fast
        default assumption handler.

        It can also accept any user-specified zero testing function, if it
        is formatted as a function which accepts a single symbolic argument
        and returns ``True`` if it is tested as zero and ``False`` if it
        is tested as non-zero, and ``None`` if it is undecidable.

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    simplify : bool or function, optional
        If ``True``, ``as_content_primitive()`` will be used to tidy up
        normalization artifacts.

        It will also be used by the ``nullspace`` routine.

    chop : bool or positive number, optional
        If the matrix contains any Floats, they will be changed to Rationals
        for computation purposes, but the answers will be returned after
        being evaluated with evalf. The ``chop`` flag is passed to ``evalf``.
        When ``chop=True`` a default precision will be used; a number will
        be interpreted as the desired level of precision.

    Returns
    =======
    ret : [(eigenval, multiplicity, eigenspace), ...]
        A ragged list containing tuples of data obtained by ``eigenvals``
        and ``nullspace``.

        ``eigenspace`` is a list containing the ``eigenvector`` for each
        eigenvalue.

        ``eigenvector`` is a vector in the form of a ``Matrix``. e.g.
        a vector of length 3 is returned as ``Matrix([a_1, a_2, a_3])``.

    Raises
    ======

    NotImplementedError
        If failed to compute nullspace.

    Examples
    ========

    >>> from sympy.matrices import Matrix
    >>> M = Matrix(3, 3, [0, 1, 1, 1, 0, 0, 1, 1, 1])
    >>> M.eigenvects()
    [(-1, 1, [Matrix([
    [-1],
    [ 1],
    [ 0]])]), (0, 1, [Matrix([
    [ 0],
    [-1],
    [ 1]])]), (2, 1, [Matrix([
    [2/3],
    [1/3],
    [  1]])])]

    See Also
    ========

    eigenvals
    MatrixSubspaces.nullspace
    """

    def eigenspace(eigenval):
        """Get a basis for the eigenspace for a particular eigenvalue"""

        m   = M - M.eye(M.rows) * eigenval
        ret = m.nullspace(iszerofunc=iszerofunc, dotprodsimp=dotprodsimp)

        # the nullspace for a real eigenvalue should be
        # non-trivial.  If we didn't find an eigenvector, try once
        # more a little harder
        if len(ret) == 0 and simplify:
            ret = m.nullspace(iszerofunc=iszerofunc, simplify=True, dotprodsimp=dotprodsimp)
        if len(ret) == 0:
            raise NotImplementedError(
                    "Can't evaluate eigenvector for eigenvalue %s" % eigenval)

        return ret

    simplify = flags.get('simplify', True)

    if not isinstance(simplify, FunctionType):
        simpfunc = _simplify if simplify else lambda x: x

    primitive = flags.get('simplify', False)
    chop      = flags.pop('chop', False)

    flags.pop('multiple', None)  # remove this if it's there

    has_floats = M.has(Float) # roots doesn't like Floats, so replace them with Rationals

    if has_floats:
        M = M.applyfunc(lambda x: nsimplify(x, rational=True))

    eigenvals = M.eigenvals(rational=False,
            error_when_incomplete=error_when_incomplete,
            dotprodsimp=dotprodsimp, **flags)

    ret = [(val, mult, eigenspace(val)) for val, mult in
                sorted(eigenvals.items(), key=default_sort_key)]

    if primitive:
        # if the primitive flag is set, get rid of any common
        # integer denominators
        def denom_clean(l):
            from sympy import gcd
            return [(v / gcd(list(v))).applyfunc(simpfunc) for v in l]

        ret = [(val, mult, denom_clean(es)) for val, mult, es in ret]

    if has_floats:
        # if we had floats to start with, turn the eigenvectors to floats
        ret = [(val.evalf(chop=chop), mult, [v.evalf(chop=chop) for v in es])
                for val, mult, es in ret]

    return ret


def _is_diagonalizable_with_eigen(M, reals_only=False, dotprodsimp=None):
    """See _is_diagonalizable. This function returns the bool along with the
    eigenvectors to avoid calculating them again in functions like
    ``diagonalize``."""

    if not M.is_square:
        return False, []

    eigenvecs = M.eigenvects(simplify=True, dotprodsimp=dotprodsimp)

    for val, mult, basis in eigenvecs:
        if reals_only and not val.is_real: # if we have a complex eigenvalue
            return False, eigenvecs

        if mult != len(basis): # if the geometric multiplicity doesn't equal the algebraic
            return False, eigenvecs

    return True, eigenvecs

def _is_diagonalizable(M, reals_only=False, dotprodsimp=None, **kwargs):
    """Returns ``True`` if a matrix is diagonalizable.

    Parameters
    ==========

    reals_only : bool, optional
        If ``True``, it tests whether the matrix can be diagonalized
        without complex numbers.

        If ``False``, it tests whether the matrix can be diagonalized
        with complex numbers.

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification
        is used during matrix multiplications to control expression
        blowup and thus speed up calculation.

    Examples
    ========

    Example of a diagonalizable matrix:

    >>> from sympy import Matrix
    >>> M = Matrix([[1, 2, 0], [0, 3, 0], [2, -4, 2]])
    >>> M.is_diagonalizable()
    True

    Example of a non-diagonalizable matrix:

    >>> M = Matrix([[0, 1], [0, 0]])
    >>> M.is_diagonalizable()
    False

    Example of a matrix which is diagonalizable with a diagonal matrix
    with complex entries, but not with real entries:

    >>> M = Matrix([[0, 1], [-1, 0]])
    >>> M.is_diagonalizable(reals_only=False)
    True
    >>> M.is_diagonalizable(reals_only=True)
    False

    See Also
    ========

    is_diagonal
    diagonalize
    """

    if 'clear_cache' in kwargs:
        SymPyDeprecationWarning(
            feature='clear_cache',
            deprecated_since_version=1.4,
            issue=15887
        ).warn()

    if 'clear_subproducts' in kwargs:
        SymPyDeprecationWarning(
            feature='clear_subproducts',
            deprecated_since_version=1.4,
            issue=15887
        ).warn()

    if not M.is_square:
        return False

    if all(e.is_real for e in M) and M.is_symmetric():
        return True

    if all(e.is_complex for e in M) and M.is_hermitian:
        return True

    return _is_diagonalizable_with_eigen(M, reals_only=reals_only,
            dotprodsimp=dotprodsimp)[0]


def _diagonalize(M, reals_only=False, sort=False, normalize=False,
        dotprodsimp=None):
    """
    Return (P, D), where D is diagonal and

        D = P^-1 * M * P

    where M is current matrix.

    Parameters
    ==========

    reals_only : bool. Whether to throw an error if complex numbers are need
                    to diagonalize. (Default: False)

    sort : bool. Sort the eigenvalues along the diagonal. (Default: False)

    normalize : bool. If True, normalize the columns of P. (Default: False)

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    Examples
    ========

    >>> from sympy.matrices import Matrix
    >>> M = Matrix(3, 3, [1, 2, 0, 0, 3, 0, 2, -4, 2])
    >>> M
    Matrix([
    [1,  2, 0],
    [0,  3, 0],
    [2, -4, 2]])
    >>> (P, D) = M.diagonalize()
    >>> D
    Matrix([
    [1, 0, 0],
    [0, 2, 0],
    [0, 0, 3]])
    >>> P
    Matrix([
    [-1, 0, -1],
    [ 0, 0, -1],
    [ 2, 1,  2]])
    >>> P.inv() * M * P
    Matrix([
    [1, 0, 0],
    [0, 2, 0],
    [0, 0, 3]])

    See Also
    ========

    is_diagonal
    is_diagonalizable
    """

    if not M.is_square:
        raise NonSquareMatrixError()

    is_diagonalizable, eigenvecs = _is_diagonalizable_with_eigen(M,
                reals_only=reals_only, dotprodsimp=dotprodsimp)

    if not is_diagonalizable:
        raise MatrixError("Matrix is not diagonalizable")

    if sort:
        eigenvecs = sorted(eigenvecs, key=default_sort_key)

    p_cols, diag = [], []

    for val, mult, basis in eigenvecs:
        diag   += [val] * mult
        p_cols += basis

    if normalize:
        p_cols = [v / v.norm() for v in p_cols]

    return M.hstack(*p_cols), M.diag(*diag)


def _eval_is_positive_definite(M, method="eigen", dotprodsimp=None):
    """Algorithm dump for computing positive-definiteness of a
    matrix.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    method : str, optional
        Specifies the method for computing positive-definiteness of
        a matrix.

        If ``'eigen'``, it computes the full eigenvalues and decides
        if the matrix is positive-definite.

        If ``'CH'``, it attempts computing the Cholesky
        decomposition to detect the definitiveness.

        If ``'LDL'``, it attempts computing the LDL
        decomposition to detect the definitiveness.
    """

    if M.is_hermitian:
        if method == 'eigen':
            eigen = M.eigenvals(dotprodsimp=dotprodsimp)
            args  = [x.is_positive for x in eigen.keys()]

            return fuzzy_and(args)

        elif method == 'CH':
            try:
                M.cholesky(hermitian=True)
            except NonPositiveDefiniteMatrixError:
                return False

            return True

        elif method == 'LDL':
            try:
                M.LDLdecomposition(hermitian=True)
            except NonPositiveDefiniteMatrixError:
                return False

            return True

        else:
            raise NotImplementedError()

    elif M.is_square:
        M_H = (M + M.H) / 2

        return M_H._eval_is_positive_definite(method=method,
                dotprodsimp=dotprodsimp)

def _is_positive_definite(M):
    return M._eval_is_positive_definite()

def _is_positive_semidefinite(M):
    if M.is_hermitian:
        eigen = M.eigenvals()
        args  = [x.is_nonnegative for x in eigen.keys()]

        return fuzzy_and(args)

    elif M.is_square:
        return ((M + M.H) / 2).is_positive_semidefinite

    return None

def _is_negative_definite(M):
    if M.is_hermitian:
        eigen = M.eigenvals()
        args  = [x.is_negative for x in eigen.keys()]

        return fuzzy_and(args)

    elif M.is_square:
        return ((M + M.H) / 2).is_negative_definite

    return None

def _is_negative_semidefinite(M):
    if M.is_hermitian:
        eigen = M.eigenvals()
        args  = [x.is_nonpositive for x in eigen.keys()]

        return fuzzy_and(args)

    elif M.is_square:
        return ((M + M.H) / 2).is_negative_semidefinite

    return None

def _is_indefinite(M):
    if M.is_hermitian:
        eigen        = M.eigenvals()
        args1        = [x.is_positive for x in eigen.keys()]
        any_positive = fuzzy_or(args1)
        args2        = [x.is_negative for x in eigen.keys()]
        any_negative = fuzzy_or(args2)

        return fuzzy_and([any_positive, any_negative])

    elif M.is_square:
        return ((M + M.H) / 2).is_indefinite

    return None

_doc_positive_definite = \
    r"""Finds out the definiteness of a matrix.

    Examples
    ========

    An example of numeric positive definite matrix:

    >>> from sympy import Matrix
    >>> A = Matrix([[1, -2], [-2, 6]])
    >>> A.is_positive_definite
    True
    >>> A.is_positive_semidefinite
    True
    >>> A.is_negative_definite
    False
    >>> A.is_negative_semidefinite
    False
    >>> A.is_indefinite
    False

    An example of numeric negative definite matrix:

    >>> A = Matrix([[-1, 2], [2, -6]])
    >>> A.is_positive_definite
    False
    >>> A.is_positive_semidefinite
    False
    >>> A.is_negative_definite
    True
    >>> A.is_negative_semidefinite
    True
    >>> A.is_indefinite
    False

    An example of numeric indefinite matrix:

    >>> A = Matrix([[1, 2], [2, 1]])
    >>> A.is_positive_definite
    False
    >>> A.is_positive_semidefinite
    False
    >>> A.is_negative_definite
    False
    >>> A.is_negative_semidefinite
    False
    >>> A.is_indefinite
    True

    Notes
    =====

    Definitiveness is not very commonly discussed for non-hermitian
    matrices.

    However, computing the definitiveness of a matrix can be
    generalized over any real matrix by taking the symmetric part:

    `A_S = 1/2 (A + A^{T})`

    Or over any complex matrix by taking the hermitian part:

    `A_H = 1/2 (A + A^{H})`

    And computing the eigenvalues.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Definiteness_of_a_matrix#Eigenvalues

    .. [2] http://mathworld.wolfram.com/PositiveDefiniteMatrix.html

    .. [3] Johnson, C. R. "Positive Definite Matrices." Amer.
        Math. Monthly 77, 259-264 1970.
    """

_is_positive_definite.__doc__     = _doc_positive_definite
_is_positive_semidefinite.__doc__ = _doc_positive_definite
_is_negative_definite.__doc__     = _doc_positive_definite
_is_negative_semidefinite.__doc__ = _doc_positive_definite
_is_indefinite.__doc__            = _doc_positive_definite


def _jordan_form(M, calc_transform=True, dotprodsimp=None, **kwargs):
    """Return ``(P, J)`` where `J` is a Jordan block
    matrix and `P` is a matrix such that

        ``M == P*J*P**-1``

    Parameters
    ==========

    calc_transform : bool
        If ``False``, then only `J` is returned.

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    chop : bool
        All matrices are converted to exact types when computing
        eigenvalues and eigenvectors.  As a result, there may be
        approximation errors.  If ``chop==True``, these errors
        will be truncated.

    Examples
    ========

    >>> from sympy.matrices import Matrix
    >>> M = Matrix([[ 6,  5, -2, -3], [-3, -1,  3,  3], [ 2,  1, -2, -3], [-1,  1,  5,  5]])
    >>> P, J = M.jordan_form()
    >>> J
    Matrix([
    [2, 1, 0, 0],
    [0, 2, 0, 0],
    [0, 0, 2, 1],
    [0, 0, 0, 2]])

    See Also
    ========

    jordan_block
    """

    if not M.is_square:
        raise NonSquareMatrixError("Only square matrices have Jordan forms")

    chop       = kwargs.pop('chop', False)
    mat        = M
    has_floats = M.has(Float)

    if has_floats:
        try:
            max_prec = max(term._prec for term in M._mat if isinstance(term, Float))
        except ValueError:
            # if no term in the matrix is explicitly a Float calling max()
            # will throw a error so setting max_prec to default value of 53
            max_prec = 53

        # setting minimum max_dps to 15 to prevent loss of precision in
        # matrix containing non evaluated expressions
        max_dps = max(prec_to_dps(max_prec), 15)

    def restore_floats(*args):
        """If ``has_floats`` is `True`, cast all ``args`` as
        matrices of floats."""

        if has_floats:
            args = [m.evalf(n=max_dps, chop=chop) for m in args]
        if len(args) == 1:
            return args[0]

        return args

    # cache calculations for some speedup
    mat_cache = {}

    def eig_mat(val, pow):
        """Cache computations of ``(M - val*I)**pow`` for quick
        retrieval"""

        if (val, pow) in mat_cache:
            return mat_cache[(val, pow)]

        if (val, pow - 1) in mat_cache:
            mat_cache[(val, pow)] = mat_cache[(val, pow - 1)].multiply(
                    mat_cache[(val, 1)], dotprodsimp=dotprodsimp)
        else:
            mat_cache[(val, pow)] = (mat - val*M.eye(M.rows)).pow(pow,
                    dotprodsimp=dotprodsimp)

        return mat_cache[(val, pow)]

    # helper functions
    def nullity_chain(val, algebraic_multiplicity):
        """Calculate the sequence  [0, nullity(E), nullity(E**2), ...]
        until it is constant where ``E = M - val*I``"""

        # mat.rank() is faster than computing the null space,
        # so use the rank-nullity theorem
        cols    = M.cols
        ret     = [0]
        nullity = cols - eig_mat(val, 1).rank(dotprodsimp=dotprodsimp)
        i       = 2

        while nullity != ret[-1]:
            ret.append(nullity)

            if nullity == algebraic_multiplicity:
                break

            nullity  = cols - eig_mat(val, i).rank(dotprodsimp=dotprodsimp)
            i       += 1

            # Due to issues like #7146 and #15872, SymPy sometimes
            # gives the wrong rank. In this case, raise an error
            # instead of returning an incorrect matrix
            if nullity < ret[-1] or nullity > algebraic_multiplicity:
                raise MatrixError(
                    "SymPy had encountered an inconsistent "
                    "result while computing Jordan block: "
                    "{}".format(M))

        return ret

    def blocks_from_nullity_chain(d):
        """Return a list of the size of each Jordan block.
        If d_n is the nullity of E**n, then the number
        of Jordan blocks of size n is

            2*d_n - d_(n-1) - d_(n+1)"""

        # d[0] is always the number of columns, so skip past it
        mid = [2*d[n] - d[n - 1] - d[n + 1] for n in range(1, len(d) - 1)]
        # d is assumed to plateau with "d[ len(d) ] == d[-1]", so
        # 2*d_n - d_(n-1) - d_(n+1) == d_n - d_(n-1)
        end = [d[-1] - d[-2]] if len(d) > 1 else [d[0]]

        return mid + end

    def pick_vec(small_basis, big_basis):
        """Picks a vector from big_basis that isn't in
        the subspace spanned by small_basis"""

        if len(small_basis) == 0:
            return big_basis[0]

        for v in big_basis:
            _, pivots = M.hstack(*(small_basis + [v])).echelon_form(
                    with_pivots=True, dotprodsimp=dotprodsimp)

            if pivots[-1] == len(small_basis):
                return v

    # roots doesn't like Floats, so replace them with Rationals
    if has_floats:
        mat = mat.applyfunc(lambda x: nsimplify(x, rational=True))

    # first calculate the jordan block structure
    eigs = mat.eigenvals(dotprodsimp=dotprodsimp)

    # make sure that we found all the roots by counting
    # the algebraic multiplicity
    if sum(m for m in eigs.values()) != mat.cols:
        raise MatrixError("Could not compute eigenvalues for {}".format(mat))

    # most matrices have distinct eigenvalues
    # and so are diagonalizable.  In this case, don't
    # do extra work!
    if len(eigs.keys()) == mat.cols:
        blocks     = list(sorted(eigs.keys(), key=default_sort_key))
        jordan_mat = mat.diag(*blocks)

        if not calc_transform:
            return restore_floats(jordan_mat)

        jordan_basis = [eig_mat(eig, 1).nullspace(dotprodsimp=dotprodsimp)[0]
                for eig in blocks]
        basis_mat    = mat.hstack(*jordan_basis)

        return restore_floats(basis_mat, jordan_mat)

    block_structure = []

    for eig in sorted(eigs.keys(), key=default_sort_key):
        algebraic_multiplicity = eigs[eig]
        chain = nullity_chain(eig, algebraic_multiplicity)
        block_sizes = blocks_from_nullity_chain(chain)

        # if block_sizes =       = [a, b, c, ...], then the number of
        # Jordan blocks of size 1 is a, of size 2 is b, etc.
        # create an array that has (eig, block_size) with one
        # entry for each block
        size_nums = [(i+1, num) for i, num in enumerate(block_sizes)]

        # we expect larger Jordan blocks to come earlier
        size_nums.reverse()

        block_structure.extend(
            (eig, size) for size, num in size_nums for _ in range(num))

    jordan_form_size = sum(size for eig, size in block_structure)

    if jordan_form_size != M.rows:
        raise MatrixError(
            "SymPy had encountered an inconsistent result while "
            "computing Jordan block. : {}".format(M))

    blocks     = (mat.jordan_block(size=size, eigenvalue=eig) for eig, size in block_structure)
    jordan_mat = mat.diag(*blocks)

    if not calc_transform:
        return restore_floats(jordan_mat)

    # For each generalized eigenspace, calculate a basis.
    # We start by looking for a vector in null( (A - eig*I)**n )
    # which isn't in null( (A - eig*I)**(n-1) ) where n is
    # the size of the Jordan block
    #
    # Ideally we'd just loop through block_structure and
    # compute each generalized eigenspace.  However, this
    # causes a lot of unneeded computation.  Instead, we
    # go through the eigenvalues separately, since we know
    # their generalized eigenspaces must have bases that
    # are linearly independent.
    jordan_basis = []

    for eig in sorted(eigs.keys(), key=default_sort_key):
        eig_basis = []

        for block_eig, size in block_structure:
            if block_eig != eig:
                continue

            null_big   = (eig_mat(eig, size)).nullspace(dotprodsimp=dotprodsimp)
            null_small = (eig_mat(eig, size - 1)).nullspace(dotprodsimp=dotprodsimp)

            # we want to pick something that is in the big basis
            # and not the small, but also something that is independent
            # of any other generalized eigenvectors from a different
            # generalized eigenspace sharing the same eigenvalue.
            vec      = pick_vec(null_small + eig_basis, null_big)
            new_vecs = [eig_mat(eig, i).multiply(vec, dotprodsimp=dotprodsimp)
                    for i in range(size)]

            eig_basis.extend(new_vecs)
            jordan_basis.extend(reversed(new_vecs))

    basis_mat = mat.hstack(*jordan_basis)

    return restore_floats(basis_mat, jordan_mat)


def _left_eigenvects(M, **flags):
    """Returns left eigenvectors and eigenvalues.

    This function returns the list of triples (eigenval, multiplicity,
    basis) for the left eigenvectors. Options are the same as for
    eigenvects(), i.e. the ``**flags`` arguments gets passed directly to
    eigenvects().

    Examples
    ========

    >>> from sympy.matrices import Matrix
    >>> M = Matrix([[0, 1, 1], [1, 0, 0], [1, 1, 1]])
    >>> M.eigenvects()
    [(-1, 1, [Matrix([
    [-1],
    [ 1],
    [ 0]])]), (0, 1, [Matrix([
    [ 0],
    [-1],
    [ 1]])]), (2, 1, [Matrix([
    [2/3],
    [1/3],
    [  1]])])]
    >>> M.left_eigenvects()
    [(-1, 1, [Matrix([[-2, 1, 1]])]), (0, 1, [Matrix([[-1, -1, 1]])]), (2,
    1, [Matrix([[1, 1, 1]])])]

    """

    eigs = M.transpose().eigenvects(**flags)

    return [(val, mult, [l.transpose() for l in basis]) for val, mult, basis in eigs]


def _singular_values(M, dotprodsimp=None):
    """Compute the singular values of a Matrix

    Examples
    ========

    >>> from sympy import Matrix, Symbol
    >>> x = Symbol('x', real=True)
    >>> M = Matrix([[0, 1, 0], [0, x, 0], [-1, 0, 0]])
    >>> M.singular_values()
    [sqrt(x**2 + 1), 1, 0]

    See Also
    ========

    condition_number
    """

    if M.rows >= M.cols:
        valmultpairs = M.H.multiply(M, dotprodsimp=dotprodsimp) \
                .eigenvals(dotprodsimp=dotprodsimp)
    else:
        valmultpairs = M.multiply(M.H, dotprodsimp=dotprodsimp) \
                .eigenvals(dotprodsimp=dotprodsimp)

    # Expands result from eigenvals into a simple list
    vals = []

    for k, v in valmultpairs.items():
        vals += [sqrt(k)] * v  # dangerous! same k in several spots!

    # Pad with zeros if singular values are computed in reverse way,
    # to give consistent format.
    if len(vals) < M.cols:
        vals += [M.zero] * (M.cols - len(vals))

    # sort them in descending order
    vals.sort(reverse=True, key=default_sort_key)

    return vals
