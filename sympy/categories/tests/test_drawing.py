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

    # A strange diagram which resulted from a typo when creating a
    # test for five lemma, but which allowed to stop one extra problem
    # in the algorithm.
    A = Object("A")
    B = Object("B")
    C = Object("C")
    D = Object("D")
    E = Object("E")
    A_ = Object("A'")
    B_ = Object("B'")
    C_ = Object("C'")
    D_ = Object("D'")
    E_ = Object("E'")

    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    h = NamedMorphism(C, D, "h")
    i = NamedMorphism(D, E, "i")

    # These 4 morphisms should be between primed objects.
    j = NamedMorphism(A, B, "j")
    k = NamedMorphism(B, C, "k")
    l = NamedMorphism(C, D, "l")
    m = NamedMorphism(D, E, "m")

    o = NamedMorphism(A, A_, "o")
    p = NamedMorphism(B, B_, "p")
    q = NamedMorphism(C, C_, "q")
    r = NamedMorphism(D, D_, "r")
    s = NamedMorphism(E, E_, "s")

    d = Diagram([f, g, h, i, j, k, l, m, o, p, q, r, s])
    grid = DiagramGrid(d)

    assert grid.width == 3
    assert grid.height == 4
    assert grid[0, 0] is None
    assert grid[0, 1] == A
    assert grid[0, 2] == A_
    assert grid[1, 0] == C
    assert grid[1, 1] == B
    assert grid[1, 2] == B_
    assert grid[2, 0] == C_
    assert grid[2, 1] == D
    assert grid[2, 2] == E
    assert grid[3, 0] is None
    assert grid[3, 1] == D_
    assert grid[3, 2] == E_

    # A cube.
    A1 = Object("A1")
    A2 = Object("A2")
    A3 = Object("A3")
    A4 = Object("A4")
    A5 = Object("A5")
    A6 = Object("A6")
    A7 = Object("A7")
    A8 = Object("A8")

    # The top face of the cube.
    f1 = NamedMorphism(A1, A2, "f1")
    f2 = NamedMorphism(A1, A3, "f2")
    f3 = NamedMorphism(A2, A4, "f3")
    f4 = NamedMorphism(A3, A4, "f3")

    # The bottom face of the cube.
    f5 = NamedMorphism(A5, A6, "f5")
    f6 = NamedMorphism(A5, A7, "f6")
    f7 = NamedMorphism(A6, A8, "f7")
    f8 = NamedMorphism(A7, A8, "f8")

    # The remaining morphisms.
    f9 = NamedMorphism(A1, A5, "f9")
    f10 = NamedMorphism(A2, A6, "f10")
    f11 = NamedMorphism(A3, A7, "f11")
    f12 = NamedMorphism(A4, A8, "f11")

    d = Diagram([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12])
    grid = DiagramGrid(d)

    assert grid.width == 4
    assert grid.height == 3
    assert grid[0, 0] is None
    assert grid[0, 1] == A5
    assert grid[0, 2] == A6
    assert grid[0, 3] is None
    assert grid[1, 0] is None
    assert grid[1, 1] == A1
    assert grid[1, 2] == A2
    assert grid[1, 3] is None
    assert grid[2, 0] == A7
    assert grid[2, 1] == A3
    assert grid[2, 2] == A4
    assert grid[2, 3] == A8
