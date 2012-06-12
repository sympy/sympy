from sympy.categories import Object, Morphism, IdentityMorphism, Diagram, Category
from sympy.utilities.pytest import XFAIL, raises
from sympy import FiniteSet, EmptySet, Dict, Tuple

def test_object():
    A = Object("A")

    assert A.name == "A"

    assert A == Object("A")
    assert A != Object("A1")

    assert hash(A) == hash(Object("A"))

    assert A != None

    raises(ValueError, lambda: Object(""))

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
    assert f.components == Tuple(f)

    assert f != None
    assert f == f
    assert f == Morphism(A, B, "f")
    assert f != g

    raises(ValueError, lambda: f * g)
    raises(ValueError, lambda: f * f)

    k = g.compose(f, "k")

    assert k.domain == A
    assert k.codomain == C
    assert k.name == "k"
    assert k.components == Tuple(f, g)

    k = g * f
    p = h * g
    u = h * g * f

    assert k.domain == A
    assert k.codomain == C
    assert k.name == ""
    assert k.components == Tuple(f, g)

    assert h * k == u
    assert p * f == u

    assert u.domain == A
    assert u.codomain == D
    assert u.name == ""
    assert u.components == Tuple(f, g, h)

    u1 = u.flatten()

    assert u1.domain == A
    assert u1.codomain == D
    assert u1.name == ""
    assert u1.components == Tuple(u1)

    u1 = u.flatten("u")

    assert u1.domain == A
    assert u1.codomain == D
    assert u1.name == "u"
    assert u1.components == Tuple(u1)

    assert f == Morphism(A, B, "f")
    assert f != g
    assert f != Morphism(A, B, "")
    assert Morphism(A, B, "") != Morphism(A, B, "")

    assert hash(f) == hash(Morphism(A, B, "f"))

    id_A = Morphism(A, A, identity=True)
    id_B = Morphism(B, B, identity=True)

    assert type(id_A) == IdentityMorphism
    assert id_A == IdentityMorphism(A, "id_A")

    assert id_A.is_identity == True
    assert id_A.components == Tuple(id_A)
    assert id_A == Morphism(A, A, name="f", identity=True)
    assert hash(id_A) == hash(Morphism(A, A, name="f", identity=True))
    assert id_A != Morphism(A, A, name="f")
    assert id_A != id_B

    assert id_A * id_A == id_A
    assert f * id_A == f
    assert id_B * f == f

    raises(ValueError, lambda: Morphism(A, B, identity=True))

    f = Morphism(A, B)
    assert f != Morphism(A, B)
    assert f == f

    raises(TypeError, lambda: f.compose(None))
    raises(TypeError, lambda: id_A.compose(None))
    raises(TypeError, lambda: f * None)
    raises(TypeError, lambda: id_A * None)

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
    d1 = Diagram([f])

    assert d1.objects == FiniteSet(A, B)
    assert d1.hom(A, B) == (FiniteSet(f), empty)
    assert d1.hom(A, A) == (FiniteSet(Morphism(A, A, identity=True)), empty)
    assert d1.hom(B, B) == (FiniteSet(Morphism(B, B, identity=True)), empty)

    # Test the addition of composites.
    d2 = Diagram([f, g])
    homAC = d2.hom(A, C)[0]

    assert d2.objects == FiniteSet(A, B, C)
    assert g * f in d2.premises.keys()

    # Test equality, inequality and hash.
    d11 = Diagram([f])

    assert d1 == d11
    assert d1 != d2
    assert hash(d1) == hash(d11)

    d11 = Diagram({f:"unique"})
    assert d1 != d11

    # Make sure that (re-)adding composites (with new properties)
    # works as expected.
    d = Diagram([f, g], {g * f:"unique"})
    assert d.conclusions[g * f] == FiniteSet("unique")

    # Check how the properties of composite morphisms are computed.
    d = Diagram({f:["unique", "isomorphism"], g:"unique"})
    assert d.premises[g * f] == FiniteSet("unique")

    # Check that conclusion morphisms with new objects are not allowed.
    d = Diagram([f], [g])
    assert d.conclusions == Dict({})

    # Test an empty diagram.
    d = Diagram()
    assert d.premises == Dict({})
    assert d.conclusions == Dict({})
    assert d.objects == empty

    # Check a SymPy Dict object.
    d = Diagram(Dict({f:FiniteSet("unique", "isomorphism"), g:"unique"}))
    assert d.premises[g * f] == FiniteSet("unique")

def test_category():
    A = Object("A")
    B = Object("B")
    C = Object("C")

    f = Morphism(A, B, "f")
    g = Morphism(B, C, "g")

    d1 = Diagram([f, g])
    d2 = Diagram([f])

    K = Category("K", commutative=[d1, d2])

    assert K.commutative == FiniteSet(d1, d2)
