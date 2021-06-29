from sympy.core.sympify import _sympify, SympifyError

from .matrices import MatrixBase


class RepMatrix(MatrixBase):

    def __eq__(self, other):
        # Skip sympify for mutable matrices...
        if not isinstance(other, RepMatrix):
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
            if not isinstance(other, RepMatrix):
                return NotImplemented

        return self._rep.unify_eq(other._rep)

    @property
    def _flat(self):
        return self._rep.to_sympy().to_list_flat()
