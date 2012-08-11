from sympy.categories import (Object, Morphism, IdentityMorphism,
                              NamedMorphism, CompositeMorphism,
                              Diagram, Category, Implication)
from sympy.categories.baseclasses import Class
from sympy.utilities.pytest import XFAIL, raises
from sympy import FiniteSet, EmptySet, Dict, Tuple

def test_morphisms():
    A = Object("A")
    B = Object("B")
    C = Object("C")
    D = Object("D")

    # Test the base morphism.
    f = NamedMorphism(A, B, "f")
    assert f.domain == A
    assert f.codomain == B
    assert f == NamedMorphism(A, B, "f")

    # Test identities.
    id_A = IdentityMorphism(A)
    id_B = IdentityMorphism(B)
    assert id_A.domain == A
    assert id_A.codomain == A
    assert id_A == IdentityMorphism(A)
    assert id_A != id_B

    # Test named morphisms.
    g = NamedMorphism(B, C, "g")
    assert g.name == "g"
    assert g != f
    assert g == NamedMorphism(B, C, "g")
    assert g != NamedMorphism(B, C, "f")

    # Test composite morphisms.
    assert f == CompositeMorphism(f)

    k = g.compose(f)
    assert k.domain == A
    assert k.codomain == C
    assert k.components == Tuple(f, g)
    assert g * f == k
    assert CompositeMorphism(f, g) == k

    assert CompositeMorphism(g * f) == g * f

    # Test the associativity of composition.
    h = NamedMorphism(C, D, "h")

    p = h * g
    u = h * g * f

    assert h * k == u
    assert p * f == u
    assert CompositeMorphism(f, g, h) == u

    # Test flattening.
    u2 = u.flatten("u")
    assert isinstance(u2, NamedMorphism)
    assert u2.name == "u"
    assert u2.domain == A
    assert u2.codomain == D

    # Test identities.
    assert f * id_A == f
    assert id_B * f == f
    assert id_A * id_A == id_A
    assert CompositeMorphism(id_A) == id_A

    # Test bad compositions.
    raises(ValueError, lambda: f * g)

    raises(TypeError, lambda: f.compose(None))
    raises(TypeError, lambda: id_A.compose(None))
    raises(TypeError, lambda: f * None)
    raises(TypeError, lambda: id_A * None)

    raises(TypeError, lambda: CompositeMorphism(f, None, 1))

    raises(ValueError, lambda: NamedMorphism(A, B, ""))
    raises(NotImplementedError, lambda: Morphism(A, B))

def test_diagram():
    A = Object("A")
    B = Object("B")
    C = Object("C")

    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    id_A = IdentityMorphism(A)
    id_B = IdentityMorphism(B)
    id_C = IdentityMorphism(C)

    empty = EmptySet()

    # Test identities.
    d1 = Diagram(f)

    assert set(d1.generators) == set([f, id_A, id_B])
    assert d1.generators_properties == Dict({f: empty, id_A: empty, id_B: empty})
    assert d1.is_finite == True
    assert set(d1.morphisms) == set([f, id_A, id_B])
    assert d1.objects == FiniteSet(A, B)
    assert set(d1.hom(A, B)) == set([f])
    assert set(d1.hom(A, A)) == set([id_A])
    assert set(d1.hom(B, B)) == set([id_B])

    # Test construction from an iterable.
    assert d1 == Diagram([f])

    # Test some basic generator simplifications.
    assert d1 == Diagram(id_A, f)
    assert d1 == Diagram(f, f)

    # Test composites.
    d2 = Diagram([f, g])

    assert d2.is_finite == True
    assert d2.objects == FiniteSet(A, B, C)
    assert g * f in d2
    assert set(d2.hom(A, C)) == set([g * f])
    assert d2 == Diagram(f, g)

    d = Diagram(g * f)
    assert set(d.generators) == set([id_A, id_C, g * f])
    assert d.generators_properties == Dict({g * f: empty, id_A: empty,
                                            id_C: empty})
    assert d.is_finite == True
    assert set(d.morphisms) == set([id_A, id_C, g * f])
    assert d.objects == FiniteSet(A, C)

    # Test equality, inequality and hash.
    d11 = Diagram([f])

    assert d1 == d11
    assert d1 != d2
    assert hash(d1) == hash(d11)

    d11 = Diagram({f: "unique"})
    assert d1 != d11

    # Make sure that composites with properties work as expected.
    d = Diagram({f: empty, g: empty, g * f: "unique"})
    assert d[g * f] == FiniteSet("unique")
    assert g * f in d

    # Check how the properties of composite morphisms are computed.
    d = Diagram({f: ["unique", "isomorphism"], g: "unique"})
    assert d[g * f] == FiniteSet("unique")
    assert g * f in d

    # Test an empty diagram.
    d = Diagram()
    assert set(d.generators) == set([])
    assert set(d.morphisms) == set([])
    assert d.objects == empty
    assert d.is_finite == True

    # Check a SymPy Dict object.
    d = Diagram(Dict({f: FiniteSet("unique", "isomorphism"), g: "unique"}))
    assert d[g * f] == FiniteSet("unique")

    # Check subdiagrams.
    d = Diagram({f: empty, g: empty, g * f:"unique"})

    d1 = Diagram([f])
    assert d1 <= d
    assert not (d1 >= d)

    d = Diagram([NamedMorphism(B, A, "f'")])
    assert not (d1 <= d)
    assert not (d1 >= d)

    d1 = Diagram({f: empty, g: empty, g * f: ["unique", "something"]})
    assert not (d1 <= d)
    assert not (d1 >= d)

    d = Diagram({f: "blooh"})
    d1 = Diagram({f: "bleeh"})
    assert not (d1 <= d)
    assert not (d1 >= d)

    d = Diagram({f: "unique", g: empty, g * f: "veryunique"})
    d1 = d.subdiagram_from_objects(FiniteSet(A, B))
    assert d1 == Diagram({f: "unique"})
    raises(ValueError, lambda: d.subdiagram_from_objects(FiniteSet(A, Object("D"))))

    # Test how identities with properties work.
    d = Diagram({id_A: "unique", f: []})
    assert set(d.generators) == set([id_A, id_B, f])
    assert d.generators_properties == Dict(
        {id_A: FiniteSet("unique"), id_B: FiniteSet(), f: FiniteSet()})
    assert d.is_finite == True
    assert set(d.morphisms) == set([id_A, id_B, f])
    assert d[id_A] == FiniteSet("unique")
    assert d[id_B] == FiniteSet()
    assert d[f] == FiniteSet()
    assert id_A in d
    assert id_C not in d

    d = Diagram(f, id_C)
    assert id_C in d

    # Test the dictionary-like interface of ``Diagram``.
    d = Diagram({f: [], g: [], g * f: "unique"})
    assert set(m for m in d) == set(m for m in d.morphisms)
    assert len(d) == 6
    assert f in d
    assert d[g * f] == FiniteSet("unique")
    assert d[g] == FiniteSet()

    # Test how property-less composites are simplified out of
    # generators.
    d = Diagram(f, g, g * f)
    assert set(d.generators) == set([f, g, id_A, id_B, id_C])
    assert set(d.morphisms) == set([id_A, id_B, id_C, f, g, g * f])
    assert d == Diagram(f, g)

def test_category():
    A = Object("A")
    B = Object("B")
    C = Object("C")

    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")

    d1 = Diagram([f, g])
    d2 = Diagram([f])

    objects = d1.objects | d2.objects

    K = Category("K", objects, commutative_diagrams=[d1, d2])

    assert K.name == "K"
    assert K.objects == Class(objects)
    assert K.commutative_diagrams == FiniteSet(d1, d2)

    raises(ValueError, lambda: Category(""))

def test_implication():
    A = Object("A")
    B = Object("B")
    C = Object("C")

    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")

    premise = Diagram(g, f)
    conclusion = Diagram({f: [], g: [], g * f: "unique"})

    imp = Implication(premise, conclusion)

    # Generic tests.
    assert imp.premise == premise
    assert imp.conclusion == conclusion
    assert imp.to_diagram() == conclusion

    # Flattening and adding an attribute to conclusions.
    flattened_implication = Diagram({g: "conclusion", f: "conclusion",
                                     g * f: ["unique", "conclusion"]})
    assert imp.to_diagram("conclusion") == flattened_implication

    # Extra objects in conclusion.
    raises(ValueError, lambda: Implication(Diagram(f), Diagram(g * f)))

    # Test diff.
    assert imp.diff() == FiniteSet(g * f)

    h = NamedMorphism(C, A, "h")
    imp = Implication(premise, Diagram(h))
    assert imp.diff() == FiniteSet(h)
