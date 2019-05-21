from sympy import Basic, Symbol, Tuple, Expr
from sympy.matrices.immutable import ImmutableMatrix
from sympy.stats.rv import RandomSymbol
from sympy.stats.index_rv import IndexedRandomSymbol
from sympy.core.compatibility import string_types

class RandomMatrix(ImmutableMatrix):

    def __new__(cls, *args, **kwargs):
        rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
        flat_list = list(flat_list)
        if not any(elem.atoms(RandomSymbol) for elem in flat_list):
            raise ValueError("Random matrix should have at least one RandomSymbol.")
        if not isinstance(flat_list, Tuple):
            flat_list = Tuple(*flat_list)
        return Basic.__new__(cls, rows, cols, flat_list)

    def rv(self, sym):
        if isinstance(sym, string_types):
            sym = Symbol(sym)
        return IndexedRandomSymbol(sym, self)

    def joint_eig_dist(self, *args, **kwargs):
        raise NotImplementedError()

    def compute_expectation(self, *args, **kwargs):
        raise NotImplementedError()

def matrix_rv(sym, *args, **kwargs):
    """
    >>> from sympy.stats.random_matrix import matrix_rv
    >>> from sympy.stats import Normal
    >>> from sympy import Matrix
    >>> X = Normal('X', 0, 1)
    >>> R = matrix_rv('R', [[X, X, X], [X, X, X], [X, X, X]])
    >>> M = Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
    >>> R.shape
    (2, 3)
    >>> R[0, 0]
    X
    >>> type(R[0, 0])
    <class 'sympy.stats.rv.RandomSymbol'>
    >>> N = matrix_rv('N', R + M) # use matrix_rv for expressions having random matrix variables
    >>> N.type
    Matrix([[X + 1, X + 2, X + 3], [X + 1, X + 2, X + 3], [X + 1, X + 2, X + 3]])
    """
    if len(args) > 1:
        raise ValueError("Too many inputs, please pass an expr, list or list of lists.")
    val = args[0]
    if isinstance(val, Expr):
        subsd = {at: at.type for at in val.atoms(IndexedRandomSymbol)}
        return RandomMatrix(val.subs(subsd), **kwargs).rv(sym)
    return RandomMatrix(*args, **kwargs)
