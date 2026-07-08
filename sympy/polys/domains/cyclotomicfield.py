"""Implementation of :class:`CyclotomicField` class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from functools import cached_property
from sympy.core.numbers import AlgebraicNumber
from sympy.external.gmpy import lcm
from sympy.polys.densetools import dup_eval
from sympy.polys.domains.algebraicfield import AlgebraicField
from sympy.polys.domains.rationalfield import QQ
from sympy.utilities import public

if TYPE_CHECKING:
    from sympy.core.expr import Expr


@public
class CyclotomicField(AlgebraicField):
    r"""Cyclotomic field :ref:`QQ(zeta_n)`

    Parameters
    ==========

    n : int
        Construct the nth cyclotomic field.
    ss : boolean, optional (default=True)
        If True, append *n* as a subscript on the alias string.
    alias : str, optional (default="zeta")
        Symbol name for the generator.
    root_index : int, optional (default=-1)
        Specifies which root of the polynomial is desired. The ordering is
        as defined by the :py:class:`~.ComplexRootOf` class. The default of
        ``-1`` selects the root $\mathrm{e}^{2\pi i/n}$.

    Examples
    ========

    >>> from sympy import CyclotomicField
    >>> K = CyclotomicField(5)
    >>> K.to_sympy(K([-1, 1]))
    1 - zeta5
    >>> K.zeta_order
    5
    """

    _zeta_order: int

    is_CyclotomicField = is_Cyclotomic = True

    def __init__(
        self,
        n: int,
        ss: bool = True,
        alias: str = "zeta",
        root_index: int = -1,
    ):
        from sympy.polys.specialpolys import cyclotomic_poly
        from sympy.polys.rootoftools import CRootOf

        if ss:
            alias += str(n)

        root = CRootOf(cyclotomic_poly(n), root_index)
        alpha = AlgebraicNumber(root, alias=alias)

        super().__init__(QQ, alpha, alias=alias)

        self._zeta_order = n

    @property
    def zeta_order(self) -> int:
        """Return the order of the unity root of the cyclotomic field.

        Examples
        ========

        >>> from sympy import CyclotomicField
        >>> CyclotomicField(12).zeta_order
        12
        """
        return self._zeta_order

    @cached_property
    def conjugate(self):
        # ext.conjugate() == ext**(order - 1)
        dom = self.dom
        rep = [dom.one] + [dom.zero] * (self.zeta_order - 1)
        conj = self.dtype(rep, self.mod.to_list(), dom)

        def conjugate(a):
            return dup_eval(a.rep, conj, dom)

        return conjugate

    def cyclotomic_field(
        self,
        n: int,
        ss: bool = True,
        alias: str = "zeta",
        gen: Expr | None = None,
        root_index: int = -1,
    ) -> CyclotomicField:
        return CyclotomicField(
            int(lcm(n, self.zeta_order)),
            ss=ss,
            alias=alias,
            root_index=root_index,
        )
