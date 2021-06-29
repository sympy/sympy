from sympy.core.kind import NumberKind, UndefinedKind
from sympy.core.sympify import _sympify, SympifyError
from sympy.polys.domains import ZZ, QQ, EXRAW

from .matrices import MatrixBase, MatrixKind


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

    def _eval_todok(self):
        return self._rep.to_sympy().to_dok()

    def _eval_values(self):
        return list(self.todok().values())

    @property
    def kind(self):
        domain = self._rep.domain
        if domain in (ZZ, QQ):
            element_kind = NumberKind
        elif domain == EXRAW:
            kinds = set(e.kind for e in self.values())
            if len(kinds) == 1:
                [element_kind] = kinds
            else:
                element_kind = UndefinedKind
        else: # pragma: no cover
            raise RuntimeError("Domain should only be ZZ, QQ or EXRAW")
        return MatrixKind(element_kind)
