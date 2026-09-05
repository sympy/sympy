"""Lifting phase of cylindrical algebraic decomposition.

Given the projection factor sets `P_1, \\ldots, P_n` of a set of
polynomials (see :func:`~.projection_sets`), the lifting phase builds the
decomposition level by level. The real line is split at the real roots of
`P_1` into *sections* (the roots) and *sectors* (the open intervals between
them). Over every cell `c` of `\\mathbb{R}^{k-1}` the polynomials of `P_k`
are evaluated at the sample point of `c`; their real roots split the
cylinder `c \\times \\mathbb{R}` into sections and sectors again, each with
its own sample point. The polynomials are sign-invariant on every cell, so
the sign of every input polynomial on a cell is read off its sample point.
"""
from __future__ import annotations

from sympy.core.numbers import Rational
from sympy.polys.polyerrors import PolynomialError
from sympy.polys.polytools import Poly

from .projection import projection_sets, _to_polys
from .samplepoints import (SamplePoint, compare_real, rational_between,
    rational_below, rational_above)


class NotWellOriented(PolynomialError):
    """Raised when a projection factor vanishes identically on a cell of
    positive dimension, so that McCallum's projection is not guaranteed to
    give a sign-invariant decomposition."""


class CADCell:
    """A cell of a cylindrical algebraic decomposition.

    Attributes
    ==========

    index : tuple of int
        The position of the cell in the decomposition, one integer per
        level: the `k`-th entry is odd for a sector and even for a section
        of the cylinder over the parent cell, counted from 1 upwards.
    point : tuple of Expr
        The coordinates of an exact sample point of the cell: rational
        numbers for the sectors and :class:`~.ComplexRootOf` roots (possibly
        times a rational number) for the sections.
    sample : SamplePoint
        The same sample point with its coordinates in a common algebraic
        field, for exact sign evaluations.
    parent : CADCell or None
        The cell of the previous level over which this cell lies.
    """

    __slots__ = ('index', 'point', 'sample', 'parent', '_signs', 'signs')

    def __init__(self, index, point, sample, parent, signs):
        self.index = index
        self.point = point
        self.sample = sample
        self.parent = parent
        # signs of the projection factors of this level at the sample point
        self._signs = signs
        # signs of the input polynomials, filled by CAD for the top level
        self.signs = None

    @property
    def level(self):
        return len(self.index)

    @property
    def dimension(self):
        """The dimension of the cell: the number of sectors in its index."""
        return sum(i % 2 for i in self.index)

    @property
    def is_section(self):
        """Whether the cell is a section (a root) over its parent."""
        return self.index[-1] % 2 == 0

    def __repr__(self):
        return "CADCell(%s, %s)" % (self.index, self.point)

    def _factor_sign(self, level, i):
        """Sign of the ``i``-th projection factor of level ``level`` at the
        sample point of this cell."""
        cell = self
        while cell.level > level:
            cell = cell.parent
        return cell._signs[i]


class CAD:
    """A cylindrical algebraic decomposition of `\\mathbb{R}^n`, sign-invariant
    for a set of polynomials.

    Instances are built by :func:`cylindrical_algebraic_decomposition`.

    Attributes
    ==========

    gens : tuple of Symbol
        The variables, in the order of the levels.
    polys : list of Poly
        The input polynomials.
    projection : list of list of Poly
        The projection factor sets ``P_1, ..., P_n``.
    method : str
        The projection operator used, ``'mccallum'`` or ``'hong'``.
    cells : list of CADCell
        The cells of `\\mathbb{R}^n`, in lexicographic order of their indices.
        The ``signs`` attribute of each cell is the tuple of signs of the
        input polynomials on the cell.
    """

    def __init__(self, gens, polys, projection, method, levels):
        self.gens = tuple(gens)
        self.polys = polys
        self.projection = projection
        self.method = method
        self._levels = levels
        self.cells = levels[-1]

    def __len__(self):
        return len(self.cells)

    def __iter__(self):
        return iter(self.cells)

    def __repr__(self):
        return "CAD(%d cells, %s)" % (len(self.cells), ", ".join(map(str, self.gens)))

    def cells_at(self, level):
        """The cells of `\\mathbb{R}^k` for ``level = k``, ``1 <= k <= n``."""
        if not 1 <= level <= len(self.gens):
            raise ValueError("level must be between 1 and %d" % len(self.gens))
        return self._levels[level - 1]

    def children(self, cell):
        """The cells of the next level lying over ``cell``."""
        if cell.level >= len(self.gens):
            return []
        return [c for c in self._levels[cell.level] if c.parent is cell]


def _merge_roots(roots, new):
    """Merge the sorted lists of distinct real roots ``roots`` and ``new``."""
    result = []
    i = j = 0
    while i < len(roots) and j < len(new):
        c = compare_real(roots[i], new[j])
        if c < 0:
            result.append(roots[i])
            i += 1
        elif c > 0:
            result.append(new[j])
            j += 1
        else:
            result.append(roots[i])
            i += 1
            j += 1
    result.extend(roots[i:])
    result.extend(new[j:])
    return result


def _lift(projection, gens, method):
    """Build the cells of all levels from the projection factor sets."""
    root = CADCell((), (), SamplePoint(), None, ())
    levels = []
    current = [root]
    for k, polys in enumerate(projection, start=1):
        level_gens = gens[:k]
        cells = []
        for parent in current:
            roots = []
            nullified = set()
            for i, f in enumerate(polys):
                r = parent.sample.real_roots(f, level_gens)
                if r is None:
                    if method == 'mccallum' and parent.dimension > 0:
                        raise NotWellOriented(
                            "%s vanishes identically on a cell of dimension %d"
                            % (f.as_expr(), parent.dimension))
                    nullified.add(i)
                else:
                    roots = _merge_roots(roots, r)
            values = []
            if not roots:
                values.append((Rational(0), None))
            else:
                values.append((rational_below(roots[0]), None))
                for j, r in enumerate(roots):
                    values.append((r, r))
                    if j + 1 < len(roots):
                        values.append((rational_between(r, roots[j + 1]), None))
                values.append((rational_above(roots[-1]), None))
            for j, (value, section_root) in enumerate(values, start=1):
                sample = parent.sample.extend(value)
                signs = []
                for i, f in enumerate(polys):
                    if i in nullified:
                        signs.append(0)
                    else:
                        signs.append(sample.sign(f, level_gens))
                cells.append(CADCell(parent.index + (j,), parent.point + (value,),
                                     sample, parent, tuple(signs)))
        levels.append(cells)
        current = cells
    return levels


def _level_of(f, gens):
    for k in range(len(gens), 0, -1):
        if f.degree(gens[k - 1]) > 0:
            return k
    return 0


def _factor_signs(polys, projection, gens):
    """For every input polynomial, the sign of its constant factor and the
    positions ``(level, index, exponent)`` of its irreducible factors in
    the projection sets."""
    positions = {}
    for k, level in enumerate(projection, start=1):
        for i, g in enumerate(level):
            positions[Poly(g.as_expr(), *gens)] = (k, i)
    result = []
    for f in polys:
        coeff, factors = f.factor_list()
        c = 0 if f.is_zero else (1 if coeff > 0 else -1)
        items = []
        for g, e in factors:
            if g.LC() < 0:
                g = -g
                if e % 2:
                    c = -c
            items.append(positions[g] + (e,))
        result.append((c, items))
    return result


def cylindrical_algebraic_decomposition(polys, gens, method=None):
    """Cylindrical algebraic decomposition sign-invariant for ``polys``.

    Parameters
    ==========

    polys : list of Expr or Poly
        Polynomials with rational coefficients in ``gens``.
    gens : list of Symbol
        The variables. The decomposition is cylindrical with respect to this
        order: the cells of `\\mathbb{R}^n` project onto cells of
        `\\mathbb{R}^{n-1}` in the first ``n - 1`` variables, and so on.
    method : ``'mccallum'``, ``'hong'`` or None
        The projection operator, see :mod:`sympy.polys.cad.projection`. By
        default McCallum's projection is tried first and Hong's is used if
        the polynomials are not well-oriented for it. With
        ``method='mccallum'`` :class:`NotWellOriented` is raised instead.

    Returns
    =======

    A :class:`CAD` whose ``cells`` are sign-invariant for every polynomial:
    each cell is a connected set on which every input polynomial is
    constantly positive, negative or zero, and every point of
    `\\mathbb{R}^n` belongs to exactly one cell.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.polys.cad import cylindrical_algebraic_decomposition
    >>> cad = cylindrical_algebraic_decomposition([x**2 + y**2 - 1], [x, y])
    >>> cad
    CAD(13 cells, x, y)
    >>> for cell in cad:
    ...     print(cell.index, cell.point, cell.signs)
    (1, 1) (-2, 0) (1,)
    (2, 1) (-1, -1) (1,)
    (2, 2) (-1, 0) (0,)
    (2, 3) (-1, 1) (1,)
    (3, 1) (0, -2) (1,)
    (3, 2) (0, -1) (0,)
    (3, 3) (0, 0) (-1,)
    (3, 4) (0, 1) (0,)
    (3, 5) (0, 2) (1,)
    (4, 1) (1, -1) (1,)
    (4, 2) (1, 0) (0,)
    (4, 3) (1, 1) (1,)
    (5, 1) (2, 0) (1,)

    The cells of index ``(2, 2)``, ``(3, 2)``, ``(3, 4)`` and ``(4, 2)``
    are the two points and the two arcs of the circle, ``(3, 3)`` is the
    open disc and the others cover its complement.
    """
    gens = list(gens)
    if not gens:
        raise ValueError("at least one generator is needed")
    polys = _to_polys(polys, gens)
    if method is None:
        methods = ['mccallum', 'hong']
    elif method in ('mccallum', 'hong'):
        methods = [method]
    else:
        raise ValueError("unknown projection method %r" % (method,))

    for m in methods:
        projection = projection_sets(polys, gens, method=m)
        try:
            levels = _lift(projection, gens, m)
        except NotWellOriented:
            if m == methods[-1]:
                raise
            continue
        break

    factor_signs = _factor_signs(polys, projection, gens)
    for cell in levels[-1]:
        signs = []
        for c, items in factor_signs:
            s = c
            for level, i, e in items:
                s *= cell._factor_sign(level, i)**e
            signs.append(s)
        cell.signs = tuple(signs)
    return CAD(gens, polys, projection, m, levels)
