from sympy.categories.diagram_drawing import _GrowableGrid
from sympy.categories import DiagramGrid, Object, NamedMorphism, Diagram

def test_GrowableGrid():
    grid = _GrowableGrid(1, 2)

    # Check dimensions.
    assert grid.width == 1
    assert grid.height == 2

    # Check initialisation of elements.
    assert grid[0, 0] == None
    assert grid[1, 0] == None

    # Check assignment to elements.
    grid[0, 0] = 1
    grid[1, 0] = "two"

    assert grid[0, 0] == 1
    assert grid[1, 0] == "two"

    # Check appending a row.
    grid.append_row()

    assert grid.width == 1
    assert grid.height == 3

    assert grid[0, 0] == 1
    assert grid[1, 0] == "two"
    assert grid[2, 0] == None

    # Check appending a column.
    grid.append_column()
    assert grid.width == 2
    assert grid.height == 3

    assert grid[0, 0] == 1
    assert grid[1, 0] == "two"
    assert grid[2, 0] == None

    assert grid[0, 1] == None
    assert grid[1, 1] == None
    assert grid[2, 1] == None

    grid = _GrowableGrid(1, 2)
    grid[0, 0] = 1
    grid[1, 0] = "two"

    # Check prepending a row.
    grid.prepend_row()
    assert grid.width == 1
    assert grid.height == 3

    assert grid[0, 0] == None
    assert grid[1, 0] == 1
    assert grid[2, 0] == "two"

    # Check prepending a column.
    grid.prepend_column()
    assert grid.width == 2
    assert grid.height == 3

    assert grid[0, 0] == None
    assert grid[1, 0] == None
    assert grid[2, 0] == None

    assert grid[0, 1] == None
    assert grid[1, 1] == 1
    assert grid[2, 1] == "two"

def test_DiagramGrid():
    # Set up some objects and morphisms.
    A = Object("A")
    B = Object("B")
    C = Object("C")
    D = Object("D")
    E = Object("E")

    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    h = NamedMorphism(D, A, "h")
    k = NamedMorphism(D, B, "k")

    # A triangle.
    d = Diagram([f, g], {g * f:"unique"})
    grid = DiagramGrid(d)

    assert grid.width == 2
    assert grid.height == 2
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[1, 0] == C
    assert grid[1, 1] is None

    # A simple diagram.
    d = Diagram([f, g, h, k])
    grid = DiagramGrid(d)

    assert grid.width == 3
    assert grid.height == 2
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] == D
    assert grid[1, 0] is None
    assert grid[1, 1] == C
    assert grid[1, 2] is None

    # A chain of morphisms.
    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    h = NamedMorphism(C, D, "h")
    k = NamedMorphism(D, E, "k")
    d = Diagram([f, g, h, k])
    grid = DiagramGrid(d)

    assert grid.width == 3
    assert grid.height == 3
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] is None
    assert grid[1, 0] is None
    assert grid[1, 1] == C
    assert grid[1, 2] == D
    assert grid[2, 0] is None
    assert grid[2, 1] is None
    assert grid[2, 2] == E

    # A square.
    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, D, "g")
    h = NamedMorphism(A, C, "h")
    k = NamedMorphism(C, D, "k")
    d = Diagram([f, g, h, k])
    grid = DiagramGrid(d)

    assert grid.width == 2
    assert grid.height == 2
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[1, 0] == C
    assert grid[1, 1] == D
