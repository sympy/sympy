from sympy.categories import Object, Morphism

def test_object():
    A = Object("A")

    assert A.name == "A"

def test_morphism():
    A = Object("A")
    B = Object("B")
    C = Object("C")
    D = Object("D")

    f = Morphism(A, B, "f")
    g = Morphism(B, C, "g")
    h = Morphism(C, D, "h")

    assert f.name == "f"
    assert f.domain == A
    assert f.codomain == B
    assert f.components == [f]

    assert f * g == None
    assert f * f == None

    k = g.compose(f, "k")

    assert k.domain == A
    assert k.codomain == C
    assert k.name == "k"
    assert k.components == [f, g]

    k = g * f
    p = h * g
    u = h * g * f

    assert k.domain == A
    assert k.codomain == C
    assert k.name == ""
    assert k.components == [f, g]

    assert h * k == u
    assert p * f == u

    assert u.domain == A
    assert u.codomain == D
    assert u.name == ""
    assert u.components == [f, g, h]

    u1 = u.flatten()

    assert u1.domain == A
    assert u1.codomain == D
    assert u1.name == ""
    assert u1.components == [u1]

    u1 = u.flatten("u")

    assert u1.domain == A
    assert u1.codomain == D
    assert u1.name == "u"
    assert u1.components == [u1]
