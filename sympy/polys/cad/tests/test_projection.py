from __future__ import annotations
from sympy.core.numbers import Rational
from sympy.polys.polytools import Poly
from sympy.polys.polyerrors import PolynomialError
from sympy.polys.cad.projection import (squarefree_basis, mccallum_projection,
    hong_projection, projection_sets, _reducta, _psc, _coefficients)
from sympy.testing.pytest import raises
from sympy.abc import x, y, z, w


def _exprs(polys):
    return {p.as_expr() for p in polys}


def test_squarefree_basis():
    assert squarefree_basis([], [x, y]) == [[], []]
    assert squarefree_basis([3, -7], [x, y]) == [[], []]
    assert squarefree_basis([(x**2 - 1)*(x + y)**2, 6], [x, y]) == [
        [Poly(x - 1, x, y), Poly(x + 1, x, y)], [Poly(x + y, x, y)]]
    # sign normalization and deduplication
    assert squarefree_basis([-x + 1, 2*x - 2, x - 1], [x]) == [[Poly(x - 1, x)]]
    # rational coefficients are cleared
    assert squarefree_basis([x**2/4 + y**2/9 - 1], [x, y]) == [
        [], [Poly(9*x**2 + 4*y**2 - 36, x, y)]]
    # levels are determined by the last generator present
    assert squarefree_basis([y**2 - 2, x*z - 1, w], [x, y, z, w]) == [
        [], [Poly(y**2 - 2, x, y, z, w)], [Poly(x*z - 1, x, y, z, w)],
        [Poly(w, x, y, z, w)]]
    assert squarefree_basis([Poly(x**2 - y, x, y)], [x, y]) == [
        [], [Poly(x**2 - y, x, y)]]
    raises(PolynomialError, lambda: squarefree_basis([x + z], [x, y]))


def test_coefficients_and_reducta():
    f = Poly(x*y**2 + (x - 1)*y + 3, y, x)
    assert _coefficients(f) == [Poly(x, x), Poly(x - 1, x), Poly(3, x)]
    assert _reducta(f) == [f, Poly((x - 1)*y + 3, y, x), Poly(3, y, x)]
    g = Poly(x*y**2 + 3, y, x)
    assert _reducta(g) == [g, Poly(3, y, x)]
    h = Poly(x**3 - 2*x, x)
    assert _coefficients(h) == [1, 0, -2, 0]
    assert _reducta(h) == [h, Poly(-2*x, x)]


def test_psc():
    f = Poly(x**2 + y**2 - 1, y, x)
    assert _psc(f, f.diff(y)) == [Poly(4*x**2 - 4, x), Poly(2, x)]
    g = Poly(y - x, y, x)
    assert _psc(f, g) == [Poly(2*x**2 - 1, x), Poly(1, x)]
    u = Poly(x**2 + 1, x)
    assert _psc(u, Poly(x**2 - 1, x)) == [4, 0, 1]


def test_mccallum_projection():
    circle = Poly(x**2 + y**2 - 1, x, y)
    assert mccallum_projection([circle], y) == [Poly(x - 1, x), Poly(x + 1, x)]
    assert mccallum_projection([circle], x) == [Poly(y - 1, y), Poly(y + 1, y)]
    line = Poly(x - y, x, y)
    assert _exprs(mccallum_projection([circle, line], y)) == {x - 1, x + 1, 2*x**2 - 1}
    assert mccallum_projection([line], y) == []
    # coefficients are taken down to the first nonzero constant one
    f = Poly((x**2 - 1)*y**2 + x*y + 5, x, y)
    assert _exprs(mccallum_projection([f], y)) == {x - 1, x + 1, x, 19*x**2 - 20}
    g = Poly(x*y**3 + (x - 2)*y + x**2, x, y)
    proj = _exprs(mccallum_projection([g], y))
    assert x in proj and x - 2 in proj
    assert g.reorder(y, x).discriminant().as_expr() in proj or any(
        p != x and p != x - 2 for p in proj)
    # a vanishing coefficient does not stop the search
    h = Poly(x*y**2 + 3, x, y)
    assert _exprs(mccallum_projection([h], y)) == {x}
    # three variables
    sphere = Poly(x**2 + y**2 + z**2 - 1, x, y, z)
    saddle = Poly(z - x*y, x, y, z)
    assert _exprs(mccallum_projection([sphere, saddle], z)) == {
        x**2 + y**2 - 1, x**2*y**2 + x**2 + y**2 - 1}


def test_hong_projection():
    circle = Poly(x**2 + y**2 - 1, x, y)
    assert hong_projection([circle], y) == [Poly(x - 1, x), Poly(x + 1, x)]
    # the constant term of a linear polynomial is a reductum, so it is added
    line = Poly(x - y, x, y)
    assert _exprs(hong_projection([circle, line], y)) == {x - 1, x + 1, x, 2*x**2 - 1}
    # Hong's projection contains McCallum's
    sphere = Poly(x**2 + y**2 + z**2 - 1, x, y, z)
    saddle = Poly(z - x*y, x, y, z)
    cubic = Poly(x*z**3 + y*z + 1, x, y, z)
    for polys in [[sphere, saddle], [sphere, cubic], [saddle, cubic],
                  [sphere, saddle, cubic]]:
        assert _exprs(mccallum_projection(polys, z)) <= _exprs(hong_projection(polys, z))
    assert _exprs(hong_projection([cubic], z)) == {x, y, 27*x + 4*y**3}
    assert _exprs(mccallum_projection([cubic], z)) == {x, y, 27*x + 4*y**3}
    f = Poly((x**2 - 1)*y**2 + x*y + 5, x, y)
    assert _exprs(hong_projection([f], y)) == {x - 1, x + 1, x, 19*x**2 - 20}


def test_projection_sets():
    circle, line = x**2 + y**2 - 1, x - y
    P1, P2 = projection_sets([circle, line], [x, y])
    assert _exprs(P1) == {x - 1, x + 1, 2*x**2 - 1}
    assert P1[0].gens == (x,)
    assert _exprs(P2) == {circle, line}
    assert P2[0].gens == (x, y)
    H1, H2 = projection_sets([circle, line], [x, y], method='hong')
    assert _exprs(H1) == {x - 1, x + 1, x, 2*x**2 - 1} and H2 == P2
    raises(ValueError, lambda: projection_sets([circle], [x, y], method='collins'))
    raises(ValueError, lambda: projection_sets([circle], []))

    assert projection_sets([3], [x, y]) == [[], []]
    assert projection_sets([x**2 - 2, y], [x, y]) == [
        [Poly(x**2 - 2, x)], [Poly(y, x, y)]]
    # polynomials in fewer variables are placed at their own level
    assert projection_sets([x*y*z + 1], [x, y, z]) == [
        [Poly(x, x)], [Poly(y, x, y)], [Poly(x*y*z + 1, x, y, z)]]
    P1, P2, P3, P4 = projection_sets([x*w + y], [x, y, z, w])
    assert (P1, P2, P3) == ([Poly(x, x)], [Poly(y, x, y)], [])
    assert P4 == [Poly(x*w + y, x, y, z, w)]

    sphere, saddle = x**2 + y**2 + z**2 - 1, z - x*y
    P1, P2, P3 = projection_sets([sphere, saddle], [x, y, z])
    assert _exprs(P3) == {sphere, x*y - z}
    assert _exprs(P2) == {x**2 + y**2 - 1, x**2*y**2 + x**2 + y**2 - 1}
    assert _exprs(P1) == {x - 1, x + 1, x**2 + 1, x}
    H1, H2, H3 = projection_sets([sphere, saddle], [x, y, z], method='hong')
    assert _exprs(P1) <= _exprs(H1) and _exprs(P2) <= _exprs(H2) and H3 == P3
    assert len(H1) == len(set(H1)) and len(H2) == len(set(H2))

    # rational coefficients
    assert projection_sets([x**2/4 + y**2/9 - 1], [x, y]) == [
        [Poly(x - 2, x), Poly(x + 2, x)], [Poly(9*x**2 + 4*y**2 - 36, x, y)]]
    assert projection_sets([Rational(1, 2)*x - y], [x, y]) == [
        [], [Poly(x - 2*y, x, y)]]
