from sympy.categories import Object, Morphism, Diagram
from sympy.utilities.pytest import XFAIL, raises
from sympy import FiniteSet, EmptySet

def test_object():
    A = Object("A")

    assert A.name == "A"

    assert A == Object("A")
    assert A != Object("A1")
    assert Object("") != A
    assert Object("") != Object("")

    assert hash(A) == hash(Object("A"))

    assert A != None

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

    assert f != None

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

    assert f == Morphism(A, B, "f")
    assert f != g
    assert f != Morphism(A, B, "")
    assert Morphism(A, B, "") != Morphism(A, B, "")

    assert hash(f) == hash(Morphism(A, B, "f"))

    id_A = Morphism(A, A, identity=True)
    id_B = Morphism(B, B, identity=True)

    assert id_A.identity == True
    assert id_A == Morphism(A, A, name="f", identity=True)
    assert id_A != Morphism(A, A, name="f")

    assert id_A * id_A == id_A
    assert f * id_A == f
    assert id_B * f == f

    raises(ValueError, lambda: Morphism(A, B, identity=True))

def test_diagram():
    A = Object("A")
    B = Object("B")
    C = Object("C")

    f = Morphism(A, B, "f")
    g = Morphism(B, C, "g")
    id_A = Morphism(A, A, "1_A")
    id_B = Morphism(B, B, "1_B")

    empty = EmptySet()

    # Test the addition of identities.
    d1 = Diagram()
    d1.add_premise(f)

    assert d1.list_objects() == FiniteSet(A, B)
    assert d1.hom(A, B) == (FiniteSet(f), empty)
    assert d1.hom(A, A) == (FiniteSet(Morphism(A, A, identity=True)), empty)
    assert d1.hom(B, B) == (FiniteSet(Morphism(B, B, identity=True)), empty)

    # Test the addition of composites.
    d2 = Diagram()
    d2.add_premise(f)
    d2.add_premise(g)
    homAC = d2.hom(A, C)[0]

    assert d2.list_objects() == FiniteSet(A, B, C)
    assert g * f in d2.premises.keys()

    # Test equality, inequality and hash.
    d11 = Diagram()
    d11.add_premise(f)

    assert d1 == d11
    assert d1 != d2
    assert hash(d1) == hash(d11)

    d11 = Diagram()
    d11.add_premise(f, "unique")

    assert d1 != d11

    # Make sure that (re-)adding composites (with new properties)
    # works as expected.
    d = Diagram()
    d.add_premise(f)
    d.add_premise(g)
    d.add_conclusion(g * f, "unique")

    assert d.conclusions[g * f] == FiniteSet("unique")

    # Check how the properties of composite morphisms are computed.
    d = Diagram()
    d.add_premise(f, "unique", "isomorphism")
    d.add_premise(g, "unique")

    assert d.premises[g * f] == FiniteSet("unique")

    # Check that conclusion morphisms with new objects are not allowed.
    d = Diagram()
    d.add_premise(f)
    d.add_conclusion(g)

    assert d.conclusions == {}
