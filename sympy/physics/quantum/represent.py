"""Logic for representing operators in state in various bases."""

from sympy import S, Add, Mul, Matrix, Pow

from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.tensorproduct import TensorProduct

__all__ = [
    'represent'
]

#-----------------------------------------------------------------------------
# Represent
#-----------------------------------------------------------------------------


def represent(expr, basis, **options):
    """Represent the quantum expression in the given basis.

    In quantum mechanics abstract states and operators can be represented in
    various basis sets. Under this operator, states become vectors or
    functions and operators become matrices or differential operators. This
    function is the top-level interface for this action.

    This function walks the sympy expression tree looking for ``QExpr``
    instances that have a ``_represent`` method. This method is then called
    and the object is replaced by the representation returned by this method.
    By default, the ``_represent`` method will dispatch to other methods
    that handle the representation logic for a particular basis set. The
    naming convention for these methods is the following::

        def _represent_FooBasis(self, e, basis, **options)

    This function will have the logic for representing instances of its class
    in the basis set having a class named ``FooBasis``.

    Parameters
    ==========
    expr  : Expr
        The expression to represent.
    basis : Operator, basis set
        An object that contains the information about the basis set. If an
        operator is used, the basis is assumed to be the orthonormal
        eigenvectors of that operator. In general though, the basis argument
        can be any object that contains the basis set information.
    options : dict
        Key/value pairs of options that are passed to the underlying method
        that does finds the representation. These options can be used to
        control how the representation is done. For example, this is where
        the size of the basis set would be set.

    Returns
    =======
    e : Expr
        The sympy expression of the represented quantum expression.

    Examples
    ========

    Here we subclass ``Operator`` and ``Ket`` to create the z-spin operator
    and its spin 1/2 up eigenstate. By definining the ``_represent_SzOp``
    method, the ket can be represented in the z-spin basis.

        >>> from sympy.physics.quantum import Operator, represent, Ket
        >>> from sympy import Matrix

        >>> class SzUpKet(Ket):
        ...     def _represent_SzOp(self, basis, **options):
        ...         return Matrix([1,0])
        ...     
        >>> class SzOp(Operator):
        ...     pass
        ... 
        >>> sz = SzOp('Sz')
        >>> up = SzUpKet('up')
        >>> represent(up, sz)
        [1]
        [0]
    """
    if isinstance(expr, QExpr):
        return expr._represent(basis, **options)
    elif isinstance(expr, Add):
        result = S.Zero
        for args in expr.args:
            if not result:
                result = represent(args, basis, **options)
            else:
                result += represent(args, basis, **options)
        return result
    elif isinstance(expr, Pow):
        return represent(expr.base, basis, **options)**expr.exp
    elif isinstance(expr, TensorProduct):
        new_args = [represent(arg, basis, **options) for arg in expr.args]
        return TensorProduct(*new_args)
    elif not isinstance(expr, Mul):
        return expr

    if not isinstance(expr, Mul):
        raise TypeError('Mul expected, got: %r' % expr)

    result = S.One
    for arg in reversed(expr.args):
        result = represent(arg, basis, **options)*result
    if isinstance(result, Matrix):
        if result.shape == (1,1):
            result = result[0]
    return result
