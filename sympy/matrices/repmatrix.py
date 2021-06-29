from sympy.core.sympify import _sympify, SympifyError

from .matrices import MatrixBase


def _compare_sequence(a, b):
    """Compares the elements of a list/tuple `a`
    and a list/tuple `b`.  `_compare_sequence((1,2), [1, 2])`
    is True, whereas `(1,2) == [1, 2]` is False"""
    if type(a) is type(b):
        # if they are the same type, compare directly
        return a == b
    # there is no overhead for calling `tuple` on a
    # tuple
    return tuple(a) == tuple(b)


class RepMatrix(MatrixBase):

    def __eq__(self, other):
        # Skip sympify for mutable matrices...
        if not isinstance(other, RepMatrix):
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented

        if isinstance(other, RepMatrix):
            return self._rep.unify_eq(other._rep)

        self_shape = getattr(self, 'shape', None)
        other_shape = getattr(other, 'shape', None)
        if self_shape == None:
            1/0
        if None in (self_shape, other_shape):
            return False
        if self_shape != other_shape:
            return False

        if isinstance(other, MatrixBase):
            1/0
            from sympy.matrices.dense import Matrix
            return _compare_sequence(self._flat, Matrix(other)._flat)

    @property
    def _flat(self):
        return self._rep.to_sympy().to_list_flat()
