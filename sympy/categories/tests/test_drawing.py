from sympy.categories.diagram_drawing import _GrowableGrid
from sympy.categories import DiagramGrid, Object, NamedMorphism, Diagram
from sympy import FiniteSet

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

    # A one-morphism diagram.
    d = Diagram([f])
    grid = DiagramGrid(d)

    assert grid.width == 2
    assert grid.height == 1
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid.morphisms == {f:FiniteSet()}

    # A triangle.
    d = Diagram([f, g], {g * f:"unique"})
    grid = DiagramGrid(d)

    assert grid.width == 2
    assert grid.height == 2
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[1, 0] == C
    assert grid[1, 1] is None
    assert grid.morphisms == {f:FiniteSet(), g:FiniteSet(),
                              g * f:FiniteSet("unique")}

    # A triangle with a "loop" morphism.
    l_A = NamedMorphism(A, A, "l_A")
    d = Diagram([f, g, l_A])
    grid = DiagramGrid(d)

    assert grid.width == 2
    assert grid.height == 2
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[1, 0] is None
    assert grid[1, 1] == C
    assert grid.morphisms == {f:FiniteSet(), g:FiniteSet(), l_A:FiniteSet()}

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
    assert grid.morphisms == {f:FiniteSet(), g:FiniteSet(), h:FiniteSet(),
                              k:FiniteSet()}

    assert str(grid) == '[[Object("A"), Object("B"), Object("D")], ' \
           '[None, Object("C"), None]]'

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
    assert grid.morphisms == {f:FiniteSet(), g:FiniteSet(), h:FiniteSet(),
                              k:FiniteSet()}

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
    assert grid.morphisms == {f:FiniteSet(), g:FiniteSet(), h:FiniteSet(),
                              k:FiniteSet()}

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
    assert grid[2, 2] == D_
    assert grid[3, 0] is None
    assert grid[3, 1] == E
    assert grid[3, 2] == E_

    morphisms = {}
    for m in [f, g, h, i, j, k, l, m, o, p, q, r, s]:
        morphisms[m] = FiniteSet()
    assert grid.morphisms == morphisms

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

    morphisms = {}
    for m in [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12]:
        morphisms[m] = FiniteSet()
    assert grid.morphisms == morphisms

    # A line diagram.
    A = Object("A")
    B = Object("B")
    C = Object("C")
    D = Object("D")
    E = Object("E")

    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    h = NamedMorphism(C, D, "h")
    i = NamedMorphism(D, E, "i")
    d = Diagram([f, g, h, i])
    grid = DiagramGrid(d, layout="sequential")

    assert grid.width == 5
    assert grid.height == 1
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] == C
    assert grid[0, 3] == D
    assert grid[0, 4] == E
    assert grid.morphisms == {f:FiniteSet(), g:FiniteSet(), h:FiniteSet(),
                              i:FiniteSet()}

    # Test the transposed version.
    grid = DiagramGrid(d, layout="sequential", transpose=True)

    assert grid.width == 1
    assert grid.height == 5
    assert grid[0, 0] == A
    assert grid[1, 0] == B
    assert grid[2, 0] == C
    assert grid[3, 0] == D
    assert grid[4, 0] == E
    assert grid.morphisms == {f:FiniteSet(), g:FiniteSet(), h:FiniteSet(),
                              i:FiniteSet()}

    # A pullback.
    m1 = NamedMorphism(A, B, "m1")
    m2 = NamedMorphism(A, C, "m2")
    s1 = NamedMorphism(B, D, "s1")
    s2 = NamedMorphism(C, D, "s2")
    f1 = NamedMorphism(E, B, "f1")
    f2 = NamedMorphism(E, C, "f2")
    g = NamedMorphism(E, A, "g")

    d = Diagram([m1, m2, s1, s2, f1, f2], {g: "unique"})
    grid = DiagramGrid(d)

    assert grid.width == 3
    assert grid.height == 2
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] == E
    assert grid[1, 0] == C
    assert grid[1, 1] == D
    assert grid[1, 2] is None

    morphisms = {g:FiniteSet("unique")}
    for m in [m1, m2, s1, s2, f1, f2]:
        morphisms[m] = FiniteSet()
    assert grid.morphisms == morphisms

    # Test the pullback with sequential layout, just for stress
    # testing.
    grid = DiagramGrid(d, layout="sequential")

    assert grid.width == 5
    assert grid.height == 1
    assert grid[0, 0] == D
    assert grid[0, 1] == B
    assert grid[0, 2] == A
    assert grid[0, 3] == C
    assert grid[0, 4] == E
    assert grid.morphisms == morphisms

    # Test a pullback with object grouping.
    grid = DiagramGrid(d, groups=FiniteSet(E, FiniteSet(A, B, C, D)))

    assert grid.width == 3
    assert grid.height == 2
    assert grid[0, 0] == E
    assert grid[0, 1] == A
    assert grid[0, 2] == B
    assert grid[1, 0] is None
    assert grid[1, 1] == C
    assert grid[1, 2] == D
    assert grid.morphisms == morphisms

    # Five lemma, actually.
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

    j = NamedMorphism(A_, B_, "j")
    k = NamedMorphism(B_, C_, "k")
    l = NamedMorphism(C_, D_, "l")
    m = NamedMorphism(D_, E_, "m")

    o = NamedMorphism(A, A_, "o")
    p = NamedMorphism(B, B_, "p")
    q = NamedMorphism(C, C_, "q")
    r = NamedMorphism(D, D_, "r")
    s = NamedMorphism(E, E_, "s")

    d = Diagram([f, g, h, i, j, k, l, m, o, p, q, r, s])
    grid = DiagramGrid(d)

    assert grid.width == 5
    assert grid.height == 3
    assert grid[0, 0] is None
    assert grid[0, 1] == A
    assert grid[0, 2] == A_
    assert grid[0, 3] is None
    assert grid[0, 4] is None
    assert grid[1, 0] == C
    assert grid[1, 1] == B
    assert grid[1, 2] == B_
    assert grid[1, 3] == C_
    assert grid[1, 4] is None
    assert grid[2, 0] == D
    assert grid[2, 1] == E
    assert grid[2, 2] is None
    assert grid[2, 3] == D_
    assert grid[2, 4] == E_

    morphisms = {}
    for m in [f, g, h, i, j, k, l, m, o, p, q, r, s]:
        morphisms[m] = FiniteSet()
    assert grid.morphisms == morphisms

    # Test the five lemma with object grouping.
    grid = DiagramGrid(d, FiniteSet(
        FiniteSet(A, B, C, D, E), FiniteSet(A_, B_, C_, D_, E_)))

    assert grid.width == 6
    assert grid.height == 3
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] is None
    assert grid[0, 3] == A_
    assert grid[0, 4] == B_
    assert grid[0, 5] is None
    assert grid[1, 0] is None
    assert grid[1, 1] == C
    assert grid[1, 2] == D
    assert grid[1, 3] is None
    assert grid[1, 4] == C_
    assert grid[1, 5] == D_
    assert grid[2, 0] is None
    assert grid[2, 1] is None
    assert grid[2, 2] == E
    assert grid[2, 3] is None
    assert grid[2, 4] is None
    assert grid[2, 5] == E_
    assert grid.morphisms == morphisms

    # Test the five lemma with object grouping, but mixing containers
    # to represent groups.
    grid = DiagramGrid(d, [(A, B, C, D, E), set([A_, B_, C_, D_, E_])])

    assert grid.width == 6
    assert grid.height == 3
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] is None
    assert grid[0, 3] == A_
    assert grid[0, 4] == B_
    assert grid[0, 5] is None
    assert grid[1, 0] is None
    assert grid[1, 1] == C
    assert grid[1, 2] == D
    assert grid[1, 3] is None
    assert grid[1, 4] == C_
    assert grid[1, 5] == D_
    assert grid[2, 0] is None
    assert grid[2, 1] is None
    assert grid[2, 2] == E
    assert grid[2, 3] is None
    assert grid[2, 4] is None
    assert grid[2, 5] == E_
    assert grid.morphisms == morphisms

    # Test the five lemma with object grouping and hints.
    grid = DiagramGrid(d, {
        FiniteSet(A, B, C, D, E): {"layout": "sequential",
                                   "transpose": True},
        FiniteSet(A_, B_, C_, D_, E_): {"layout": "sequential",
                                        "transpose": True}},
                       transpose=True)

    assert grid.width == 5
    assert grid.height == 2
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] == C
    assert grid[0, 3] == D
    assert grid[0, 4] == E
    assert grid[1, 0] == A_
    assert grid[1, 1] == B_
    assert grid[1, 2] == C_
    assert grid[1, 3] == D_
    assert grid[1, 4] == E_
    assert grid.morphisms == morphisms

    # A two-triangle disconnected diagram.
    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    f_ = NamedMorphism(A_, B_, "f")
    g_ = NamedMorphism(B_, C_, "g")
    d = Diagram([f, g, f_, g_], {g * f: "unique", g_ * f_: "unique"})
    grid = DiagramGrid(d)

    assert grid.width == 4
    assert grid.height == 2
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] == A_
    assert grid[0, 3] == B_
    assert grid[1, 0] == C
    assert grid[1, 1] is None
    assert grid[1, 2] == C_
    assert grid[1, 3] is None
    assert grid.morphisms == {f: FiniteSet(), g: FiniteSet(), f_: FiniteSet(),
                              g_: FiniteSet(), g * f: FiniteSet("unique"),
                              g_ * f_: FiniteSet("unique")}

    # A two-morphism disconnected diagram.
    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(C, D, "g")
    d = Diagram([f, g])
    grid = DiagramGrid(d)

    assert grid.width == 4
    assert grid.height == 1
    assert grid[0, 0] == A
    assert grid[0, 1] == B
    assert grid[0, 2] == C
    assert grid[0, 3] == D
    assert grid.morphisms == {f: FiniteSet(), g: FiniteSet()}

    # Test a one-object diagram.
    f = NamedMorphism(A, A, "f")
    d = Diagram([f])
    grid = DiagramGrid(d)

    assert grid.width == 1
    assert grid.height == 1
    assert grid[0, 0] == A

    # Test a two-object disconnected diagram.
    g = NamedMorphism(B, B, "g")
    d = Diagram([f, g])
    grid = DiagramGrid(d)

    assert grid.width == 2
    assert grid.height == 1
    assert grid[0, 0] == A
    assert grid[0, 1] == B

    # Test a diagram in which even growing a pseudopod does not
    # eventually help.
    F = Object("F")
    f1 = NamedMorphism(A, B, "f1")
    f2 = NamedMorphism(A, C, "f2")
    f3 = NamedMorphism(A, D, "f3")
    f4 = NamedMorphism(A, E, "f4")
    f5 = NamedMorphism(A, A_, "f5")
    f6 = NamedMorphism(A, B_, "f6")
    f7 = NamedMorphism(A, C_, "f7")
    f8 = NamedMorphism(A, D_, "f8")
    f9 = NamedMorphism(A, E_, "f9")
    f10 = NamedMorphism(A, F, "f10")
    d = Diagram([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10])
    grid = DiagramGrid(d)

    assert grid.width == 5
    assert grid.height == 3
    assert grid[0, 0] == E
    assert grid[0, 1] == C
    assert grid[0, 2] == C_
    assert grid[0, 3] == E_
    assert grid[0, 4] == F
    assert grid[1, 0] == D
    assert grid[1, 1] == A
    assert grid[1, 2] == A_
    assert grid[1, 3] is None
    assert grid[1, 4] is None
    assert grid[2, 0] == D_
    assert grid[2, 1] == B
    assert grid[2, 2] == B_
    assert grid[2, 3] is None
    assert grid[2, 4] is None

    morphisms = {}
    for f in [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]:
        morphisms[f] = FiniteSet()
    assert grid.morphisms == morphisms
