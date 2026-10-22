from __future__ import annotations

from sympy import symbols, And, Lt, Gt, Eq, S, Ge, Le
from sympy.polys.cad import projection, cad, solve_cad, Exists, ForAll, solve_qe

def test_projection():
    x, y = symbols('x y')
    f = x**2 + y**2 - 1
    # Projection of f with respect to y should include coeffs and discriminant with respect to y
    # f = y**2 + (x**2 - 1). coeffs: [1, 0, x**2 - 1]. disc: -4*(x**2 - 1)
    proj = projection([f], y)
    proj_exprs = [p.as_expr() for p in proj]
    # Should contain a polynomial that is a multiple of x**2 - 1
    assert any(p.has(x**2 - 1) or p.has(1 - x**2) for p in proj_exprs)


def test_cad_1d():
    x = symbols('x')
    # 1D CAD on x**2 - 1
    # Roots: -1, 1.
    # Should have 5 cells: sample points should be below -1, at -1, between -1 and 1, at 1, above 1.
    cells = cad([x**2 - 1], [x])
    assert len(cells) == 5
    sample_points = [float(cell.sample_point[0]) for cell in cells]
    assert sample_points[0] < -1
    assert sample_points[1] == -1
    assert -1 < sample_points[2] < 1
    assert sample_points[3] == 1
    assert sample_points[4] > 1


def test_cad_2d():
    x, y = symbols('x y')
    # 2D CAD on circle x**2 + y**2 - 1
    cells = cad([x**2 + y**2 - 1], [x, y])
    # The number of cells should be finite and non-empty
    assert len(cells) > 0
    # Let's check some sample points
    # If x = 0, y**2 - 1 has roots -1, 1.
    # So there should be a cell at (0, 0)
    assert any(cell.sample_point == (0, 0) for cell in cells)


def test_solve_cad():
    x, y = symbols('x y')
    # Solve unit disk: x**2 + y**2 < 1
    formula = Lt(x**2 + y**2, 1)
    matching_cells = solve_cad(formula, [x, y])

    # Verify that the cell with sample point (0, 0) is matching
    assert any(cell.sample_point == (0, 0) for cell in matching_cells)
    # Verify that cells far away do not match
    for cell in matching_cells:
        sx, sy = cell.sample_point
        assert sx**2 + sy**2 < 1


def test_solve_cad_triangle():
    x, y = symbols('x y')
    # Solve a standard triangle: x > 0 and y > 0 and x + y < 1
    formula = And(Gt(x, 0), Gt(y, 0), Lt(x + y, 1))
    matching_cells = solve_cad(formula, [x, y])

    assert len(matching_cells) > 0
    # The point (1/2, 1/4) should be in the solution
    assert any(cell.sample_point == (S(1)/2, S(1)/4) for cell in matching_cells)

    # Check that all matching cells are indeed in the triangle
    for cell in matching_cells:
        sx, sy = cell.sample_point
        assert sx > 0
        assert sy > 0
        assert sx + sy < 1


def test_solve_cad_algebraic_lifting():
    x, y = symbols('x y')
    # Solve system with algebraic roots: x**2 - 2 == 0 and y**2 - x == 0
    # Roots: (sqrt(2), 2**(1/4)) and (sqrt(2), -2**(1/4))
    formula = And(Eq(x**2 - 2, 0), Eq(y**2 - x, 0))
    matching_cells = solve_cad(formula, [x, y])

    assert len(matching_cells) == 2
    # Verify that the coordinates are algebraic roots
    coords = [cell.sample_point for cell in matching_cells]
    assert any(c[0]**2 == 2 and c[1]**2 == c[0] for c in coords)


def test_solve_qe():
    x, y = symbols('x y')
    # 1. Existential quantifier: Exists(y, x**2 + y**2 < 1)
    # This should eliminate y and return x**2 - 1 < 0 (equivalent to -1 < x < 1)
    formula1 = Exists(y, Lt(x**2 + y**2, 1))
    res1 = solve_qe(formula1)
    assert res1 == Lt(x**2 - 1, 0)

    # 2. Universal quantifier: ForAll(y, y**2 >= x)
    # This should eliminate y and return x <= 0
    formula2 = ForAll(y, Ge(y**2 - x, 0))
    res2 = solve_qe(formula2)
    assert res2 == Le(x, 0)
