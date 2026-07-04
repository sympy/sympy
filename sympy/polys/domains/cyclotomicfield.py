"""Implementation of :class:`CyclotomicField` class. """
from __future__ import annotations

from typing import TYPE_CHECKING
from functools import cached_property
from sympy.core.numbers import AlgebraicNumber
from sympy.external.gmpy import lcm
from sympy.polys.domains.algebraicfield import AlgebraicField
from sympy.polys.domains.rationalfield import QQ
from sympy.utilities import public

if TYPE_CHECKING:
    from sympy.core.expr import Expr

@public
class CyclotomicField(
    AlgebraicField
):
    r"""Cyclotomic field :ref:`QQ<zeta>`

    Parameters
    ==========

    n : int
        Construct the nth cyclotomic field.
    ss : boolean, optional (default=False)
        If True, append *n* as a subscript on the alias string.
    alias : str, optional (default="zeta")
        Symbol name for the generator.
    gen : :py:class:`~.Symbol`, optional (default=None)
        Desired variable for the cyclotomic polynomial that defines the
        field. If ``None``, a dummy variable will be used.
    root_index : int, optional (default=-1)
        Specifies which root of the polynomial is desired. The ordering is
        as defined by the :py:class:`~.ComplexRootOf` class. The default of
        ``-1`` selects the root $\mathrm{e}^{2\pi i/n}$.

    Examples
    ========

    >>> from sympy import CyclotomicField
    >>> K = CyclotomicField(5)
    >>> K.to_sympy(K([-1, 1]))
    1 - zeta
    >>> L = CyclotomicField(7, ss=True)
    >>> a = L.to_sympy(L([-1, 1]))
    >>> print(a)
    1 - zeta7
    """
    _conductor: int

    is_CyclotomicField = is_Cyclotomic = True

    def __init__(self, n: int, ss: bool = False, alias: str = "zeta",
                 gen: Expr | None = None, root_index: int = -1
                 ):
        from sympy.polys.specialpolys import cyclotomic_poly
        from sympy.polys.rootoftools import CRootOf

        if ss:
            alias += str(n)

        root = CRootOf(cyclotomic_poly(n, gen), root_index)
        alpha = AlgebraicNumber(root, alias=alias)

        super().__init__(QQ, alpha, alias=alias)

        self._conductor = n

    @property
    def conductor(self) -> int:
        """Return the conductor of the cyclotomic field.

        Examples
        ========

        >>> from sympy import CyclotomicField
        >>> CyclotomicField(12).conductor
        12
        """
        return self._conductor

    @cached_property
    def conjugate(self):
        # ext.conjugate() == ext**(conductor - 1)
        dom = self.dom
        rep = [dom.one] + [dom.zero] * (self.conductor - 1)
        conj = self.dtype(rep, self.mod.to_list(), dom)

        def conjugate(a):
            v = self.zero
            for c in a.rep:
                v = v * conj + c
            return v

        return conjugate

    def cyclotomic_field(self, n: int, ss: bool = False, alias: str = "zeta",
                         gen: Expr | None = None, root_index: int = -1
                         ) -> CyclotomicField:
        return CyclotomicField(int(lcm(n, self.conductor)), ss=ss, alias=alias,
                               gen=gen, root_index=root_index)
