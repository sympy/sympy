from __future__ import annotations

from sympy.core.numbers import Rational
from sympy.polys.polytools import Poly
from sympy.polys.rootoftools import CRootOf
from sympy.polys.cad.lifting import (cylindrical_algebraic_decomposition, CAD,
    CADCell, NotWellOriented, _merge_roots)
from sympy.polys.cad.samplepoints import SamplePoint, compare_real
from sympy.testing.pytest import raises
from sympy.abc import x, y, z, w


def _check_signs(cad):
    """The signs of the input polynomials on each cell, which are derived
    from the signs of the projection factors, must agree with a direct
    exact evaluation at the sample point."""
    for cell in cad:
        direct = tuple(cell.sample.sign(f, cad.gens) for f in cad.polys)
        assert direct == cell.signs, (cell, direct)
        assert cell.sample.as_exprs() != () or cell.point == ()
        for value, coord in zip(cell.point, cell.sample.as_exprs()):
            assert abs(value.evalf(30) - coord.evalf(30)) < 1e-25


def _check_structure(cad):
    """Indices are lexicographic, parents are consistent and the sample
    points of the children lie over the parent's sample point."""
    n = len(cad.gens)
    for k in range(1, n + 1):
        cells = cad.cells_at(k)
        assert [c.index for c in cells] == sorted(c.index for c in cells)
        assert all(c.level == k for c in cells)
        if k > 1:
            parents = cad.cells_at(k - 1)
            for c in cells:
                assert c.parent in parents
                assert c.index[:-1] == c.parent.index
                assert c.point[:-1] == c.parent.point
            for p in parents:
                children = cad.children(p)
                assert children and [c.index[-1] for c in children] == list(range(1, len(children) + 1))
                # sections alternate with sectors
                assert len(children) % 2 == 1
                for c in children:
                    assert c.is_section == (c.index[-1] % 2 == 0)
                # sections are ordered
                sections = [c.point[-1] for c in children if c.is_section]
                for a, b in zip(sections, sections[1:]):
                    assert compare_real(a, b) < 0
        for c in cells:
            assert c.dimension == sum(i % 2 for i in c.index)
    assert cad.cells is cad.cells_at(n)


def test_merge_roots():
    r = CRootOf(x**2 - 2, 1)
    assert _merge_roots([], []) == []
    assert _merge_roots([1], []) == [1]
    assert _merge_roots([Rational(-1), r], [Rational(1), r, 2]) == [-1, 1, r, 2]
    assert _merge_roots([0, 3], [1, 2]) == [0, 1, 2, 3]


def test_univariate():
    cad = cylindrical_algebraic_decomposition([x**2 - 2, x - 1], [x])
    assert len(cad) == 7
    assert cad.gens == (x,)
    r0, r1 = CRootOf(x**2 - 2, 0), CRootOf(x**2 - 2, 1)
    assert [c.point for c in cad] == [(-2,), (r0,), (0,), (1,), (Rational(5, 4),), (r1,), (2,)]
    assert [c.signs for c in cad] == [(1, -1), (0, -1), (-1, -1), (-1, 0),
                                       (-1, 1), (0, 1), (1, 1)]
    assert [c.index for c in cad] == [(1,), (2,), (3,), (4,), (5,), (6,), (7,)]
    assert [c.dimension for c in cad] == [1, 0, 1, 0, 1, 0, 1]
    assert cad.projection == [[Poly(x**2 - 2, x), Poly(x - 1, x)]]
    assert cad.method == 'mccallum'
    _check_signs(cad)
    _check_structure(cad)

    cad = cylindrical_algebraic_decomposition([x**2 + 1], [x])
    assert len(cad) == 1 and cad.cells[0].point == (0,) and cad.cells[0].signs == (1,)
    cad = cylindrical_algebraic_decomposition([], [x])
    assert len(cad) == 1 and cad.cells[0].point == (0,) and cad.cells[0].signs == ()
    cad = cylindrical_algebraic_decomposition([3, -x**2], [x])
    assert [c.signs for c in cad] == [(1, -1), (1, 0), (1, -1)]
    cad = cylindrical_algebraic_decomposition([0], [x])
    assert [c.signs for c in cad] == [(0,)]
    cad = cylindrical_algebraic_decomposition([(x - 1)**2*(x + 1)], [x])
    assert [c.point for c in cad] == [(-2,), (-1,), (0,), (1,), (2,)]
    assert [c.signs for c in cad] == [(-1,), (0,), (1,), (0,), (1,)]
    cad = cylindrical_algebraic_decomposition([x/2 - Rational(1, 3)], [x])
    assert [c.point for c in cad] == [(0,), (Rational(2, 3),), (1,)]
    assert [c.signs for c in cad] == [(-1,), (0,), (1,)]


def test_circle():
    cad = cylindrical_algebraic_decomposition([x**2 + y**2 - 1], [x, y])
    assert repr(cad) == "CAD(13 cells, x, y)"
    assert len(cad) == 13
    assert [c.index for c in cad] == [
        (1, 1), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5),
        (4, 1), (4, 2), (4, 3), (5, 1)]
    assert [c.point for c in cad] == [
        (-2, 0), (-1, -1), (-1, 0), (-1, 1), (0, -2), (0, -1), (0, 0), (0, 1),
        (0, 2), (1, -1), (1, 0), (1, 1), (2, 0)]
    assert [c.signs[0] for c in cad] == [1, 1, 0, 1, 1, 0, -1, 0, 1, 1, 0, 1, 1]
    assert [c.dimension for c in cad] == [2, 1, 0, 1, 2, 1, 2, 1, 2, 1, 0, 1, 2]
    assert len(cad.cells_at(1)) == 5
    assert [c.point for c in cad.cells_at(1)] == [(-2,), (-1,), (0,), (1,), (2,)]
    assert cad.children(cad.cells_at(1)[2]) == cad.cells[4:9]
    assert cad.children(cad.cells[0]) == []
    raises(ValueError, lambda: cad.cells_at(0))
    raises(ValueError, lambda: cad.cells_at(3))
    assert repr(cad.cells[6]) == "CADCell((3, 3), (0, 0))"
    assert isinstance(cad.cells[6].sample, SamplePoint)
    assert isinstance(cad.cells[6], CADCell)
    assert isinstance(cad, CAD)
    assert list(cad) == cad.cells
    _check_signs(cad)
    _check_structure(cad)


def test_circle_and_line():
    circle, line = x**2 + y**2 - 1, x - y
    cad = cylindrical_algebraic_decomposition([circle, line], [x, y])
    assert len(cad) == 47
    assert [len(cad.children(c)) for c in cad.cells_at(1)] == [3, 5, 7, 5, 7, 5, 7, 5, 3]
    assert {c.signs for c in cad} == {
        (1, -1), (1, 1), (1, 0), (0, -1), (0, 1), (0, 0), (-1, -1), (-1, 1), (-1, 0)}
    # the two intersection points
    both = [c for c in cad if c.signs == (0, 0)]
    r0, r1 = CRootOf(2*x**2 - 1, 0), CRootOf(2*x**2 - 1, 1)
    assert [c.point for c in both] == [(r0, r0), (r1, r1)]
    assert [c.dimension for c in both] == [0, 0]
    # the open regions
    regions = [c for c in cad if c.dimension == 2]
    assert len(regions) == 16
    assert {c.signs for c in regions} == {(1, -1), (1, 1), (-1, -1), (-1, 1)}
    _check_signs(cad)
    _check_structure(cad)
    # the order of the polynomials does not change the cells
    cad2 = cylindrical_algebraic_decomposition([line, circle], [x, y])
    assert [c.point for c in cad2] == [c.point for c in cad]
    assert [c.signs for c in cad2] == [c.signs[::-1] for c in cad]
    # the other variable order gives the symmetric decomposition
    cad3 = cylindrical_algebraic_decomposition([circle, line], [y, x])
    assert [c.point for c in cad3] == [c.point for c in cad]
    assert [c.signs for c in cad3] == [(s, -t) for s, t in [c.signs for c in cad]]


def test_parabola_and_circle():
    cad = cylindrical_algebraic_decomposition([y - x**2, x**2 + y**2 - 1], [x, y])
    assert len(cad) == 47
    assert [c.point[0] for c in cad.cells_at(1)] == [
        -2, -1, Rational(-7, 8), CRootOf(x**4 + x**2 - 1, 0), 0,
        CRootOf(x**4 + x**2 - 1, 1), Rational(7, 8), 1, 2]
    both = [c for c in cad if c.signs == (0, 0)]
    assert [c.index for c in both] == [(4, 4), (6, 4)]
    assert [c.point[1] for c in both] == [CRootOf(x**2 + x - 1, 1)]*2
    # the sample point of the intersection lies in a field of degree 4
    assert both[0].sample.degree == 4
    assert both[0].sample.sign(y - x**2, [x, y]) == 0
    _check_signs(cad)
    _check_structure(cad)


def test_polynomials_in_fewer_variables():
    cad = cylindrical_algebraic_decomposition([x - 1], [x, y])
    assert [c.point for c in cad] == [(0, 0), (1, 0), (2, 0)]
    assert [c.signs for c in cad] == [(-1,), (0,), (1,)]
    cad = cylindrical_algebraic_decomposition([y - 1], [x, y])
    assert [c.point for c in cad] == [(0, 0), (0, 1), (0, 2)]
    assert [c.signs for c in cad] == [(-1,), (0,), (1,)]
    cad = cylindrical_algebraic_decomposition([x*y], [x, y])
    assert [c.point for c in cad] == [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    assert [c.signs for c in cad] == [(a*b,) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    cad = cylindrical_algebraic_decomposition([x*y - 1], [x, y])
    assert [c.point for c in cad] == [
        (-1, -2), (-1, -1), (-1, 0), (0, 0), (1, 0), (1, 1), (1, 2)]
    assert [c.signs for c in cad] == [(1,), (0,), (-1,), (-1,), (-1,), (0,), (1,)]
    _check_signs(cad)
    _check_structure(cad)
    cad = cylindrical_algebraic_decomposition([], [x, y, z])
    assert [c.point for c in cad] == [(0, 0, 0)]
    assert cad.cells[0].index == (1, 1, 1)


def test_three_variables():
    sphere, saddle = x**2 + y**2 + z**2 - 1, z - x*y
    cad = cylindrical_algebraic_decomposition([sphere, saddle], [x, y, z])
    assert len(cad) == 137
    assert cad.method == 'mccallum'
    assert {c.signs for c in cad} == {(s, t) for s in (-1, 0, 1) for t in (-1, 0, 1)}
    _check_signs(cad)
    _check_structure(cad)
    # Whitney umbrella
    cad = cylindrical_algebraic_decomposition([x**2 - y**2*z], [x, y, z])
    assert len(cad) == 21
    _check_signs(cad)
    _check_structure(cad)


def test_not_well_oriented():
    # x*w + y vanishes identically on the line x = y = 0 in (x, y, z)
    raises(NotWellOriented, lambda: cylindrical_algebraic_decomposition(
        [x*w + y], [x, y, z, w], method='mccallum'))
    cad = cylindrical_algebraic_decomposition([x*w + y], [x, y, z, w])
    assert cad.method == 'hong'
    assert len(cad) == 21
    assert {c.signs for c in cad} == {(-1,), (0,), (1,)}
    zero = [c for c in cad if c.signs == (0,)]
    assert [c.point for c in zero] == [
        (-1, -1, 0, -1), (-1, 0, 0, 0), (-1, 1, 0, 1), (0, 0, 0, 0),
        (1, -1, 0, 1), (1, 0, 0, 0), (1, 1, 0, -1)]
    assert [c.dimension for c in zero] == [3, 2, 3, 2, 3, 2, 3]
    _check_signs(cad)
    _check_structure(cad)
    cad2 = cylindrical_algebraic_decomposition([x*w + y], [x, y, z, w], method='hong')
    assert [c.point for c in cad2] == [c.point for c in cad]
    # the same polynomial in three variables is well-oriented
    cad = cylindrical_algebraic_decomposition([x*z + y], [x, y, z], method='mccallum')
    assert cad.method == 'mccallum' and len(cad) == 21
    _check_signs(cad)


def test_hong_method():
    circle, line = x**2 + y**2 - 1, x - y
    cad = cylindrical_algebraic_decomposition([circle, line], [x, y], method='hong')
    assert cad.method == 'hong'
    # Hong's projection adds the factor x, so there are more cells
    assert len(cad) == 61
    assert {c.signs for c in cad if c.dimension == 2} == {
        (1, -1), (1, 1), (-1, -1), (-1, 1)}
    _check_signs(cad)
    _check_structure(cad)


def test_errors():
    raises(ValueError, lambda: cylindrical_algebraic_decomposition([x], []))
    raises(ValueError, lambda: cylindrical_algebraic_decomposition([x], [x], method='collins'))
    raises(Exception, lambda: cylindrical_algebraic_decomposition([x + z], [x, y]))
    raises(Exception, lambda: cylindrical_algebraic_decomposition([x**Rational(1, 2)], [x]))
