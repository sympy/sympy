from sympy.categories import (Object, Morphism, IdentityMorphism,
                              NamedMorphism, CompositeMorphism,
                              DerivedMorphism, Diagram, Category,
                              Implication)
from sympy.categories.baseclasses import Class
from sympy.utilities.pytest import XFAIL, raises
from sympy import FiniteSet, EmptySet, Dict, Tuple
from itertools import islice

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
    assert len(k) == 2
    assert tuple(k) == k.components
    assert f in k
    assert g in k

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
    assert isinstance(u2, DerivedMorphism)
    assert u2.name == "u"
    assert u2.domain == A
    assert u2.codomain == D
    assert u2.original_morphism == u

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

    # Test derived morphisms.
    assert DerivedMorphism(A, B, "f", f) != f

    g = DerivedMorphism(A, B, "g", f)
    assert g.original_morphism == f

    raises(ValueError, lambda: DerivedMorphism(A, B, "", f))
    raises(ValueError, lambda: DerivedMorphism(A, B, "f", None))

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
    assert d1.is_finite
    assert d1.is_hom_set_finite(A, B)
    assert d1.is_hom_set_finite(A, A)
    assert d1.is_hom_set_finite(B, B)

    assert not d1.is_hom_set_empty(A, B)
    assert d1.is_hom_set_empty(B, A)

    assert set(d1.morphisms) == set([f, id_A, id_B])
    assert set(d1.expanded_generators) == set(d1.morphisms)
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

    assert d2.is_finite
    assert d2.objects == FiniteSet(A, B, C)
    assert g * f in d2
    assert set(d2.hom(A, C)) == set([g * f])
    assert d2 == Diagram(f, g)
    assert set(d2.expanded_generators) == set(d2.morphisms)

    raises(ValueError, lambda: Diagram([g * f]))
    raises(ValueError, lambda: Diagram(g * f))

    # Test equality, inequality and hash.
    d11 = Diagram([f])

    assert d1 == d11
    assert d1 != d2
    assert hash(d1) == hash(d11)

    d11 = Diagram({f: "unique"})
    assert d1 != d11
    assert d11.get(f) == FiniteSet("unique")
    assert d11.get(g) is None
    raises(ValueError, lambda: d11[g])

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
    assert set(d.expanded_generators) == set([])
    assert set(d.morphisms) == set([])
    assert d.objects == empty
    assert d.is_finite

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
    d1 = d.subdiagram_from_objects([A, B])
    assert d1 == Diagram({f: "unique"})
    raises(ValueError, lambda: d.subdiagram_from_objects([A, Object("D")]))

    # Test how identities with properties work.
    d = Diagram({id_A: "unique", f: []})
    assert set(d.generators) == set([id_A, id_B, f])
    assert d.generators_properties == Dict(
        {id_A: FiniteSet("unique"), id_B: FiniteSet(), f: FiniteSet()})
    assert d.is_finite
    assert set(d.morphisms) == set([id_A, id_B, f])
    assert set(d.morphisms) == set(d.expanded_generators)
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

    # Further tests for cycle detection.
    h = NamedMorphism(C, A, "h")
    d = Diagram(f, g, h)
    assert not d.is_finite
    assert not d.is_hom_set_finite(A, B)
    assert not d.is_hom_set_finite(B, C)
    assert not d.is_hom_set_finite(C, A)
    assert not d.is_hom_set_finite(A, A)
    assert not d.is_hom_set_finite(B, B)
    assert not d.is_hom_set_finite(C, C)
    assert set(islice(d, 6)) == set(d.generators)
    assert set(d.expanded_generators) == set(islice(d, 18))
    raises(TypeError, lambda: len(d))

    assert not d.is_hom_set_empty(A, B)

    # The same diagram, but with two extra morphisms sticking outward.
    D = Object("D")
    E = Object("E")
    k = NamedMorphism(D, A, "k")
    m = NamedMorphism(E, D, "m")
    d = Diagram(f, g, h, k, m)
    assert not d.is_finite
    assert not d.is_hom_set_finite(A, B)
    assert d.is_hom_set_finite(E, D)
    assert d.is_hom_set_finite(E, E)
    assert not d.is_hom_set_finite(D, A)
    assert not d.is_hom_set_finite(E, A)
    assert d.is_hom_set_empty(D, E)

    # Test condensation.
    d = Diagram(f, g, h)
    c = d.condensation
    assert c.objects == FiniteSet((d,))
    assert c.generators_properties == Dict({
        IdentityMorphism(d): FiniteSet()})
    assert set(c.morphisms) == set(c.generators)
    assert c.is_finite

    # Further test condensation for a diagram with two strongly
    # connected components (`\{A, B, C\}` and `\{A', B', C'\}`).
    A_ = Object("A'")
    B_ = Object("B'")
    C_ = Object("C'")
    f_ = NamedMorphism(A_, B_, "f'")
    g_ = NamedMorphism(B_, C_, "g'")
    h_ = NamedMorphism(C_, A_, "h'")

    k1 = NamedMorphism(A, A_, "k1")
    k2 = NamedMorphism(B, C_, "k2")

    d = Diagram({f: [], g: [], h: [],
                 f_: [], g_: [], h_: [],
                 k1: [], k2: "blaster"})
    d1 = Diagram(f, g, h)
    d2 = Diagram(f_, g_, h_)
    assert not d.is_finite

    c = d.condensation
    assert c.objects == FiniteSet(d1, d2)
    assert c.generators_properties == Dict({
        IdentityMorphism(d1): FiniteSet(),
        IdentityMorphism(d2): FiniteSet(),
        DerivedMorphism(d1, d2, "xi_1", k1): FiniteSet(),
        DerivedMorphism(d1, d2, "xi_2", k2): FiniteSet("blaster")
        })
    assert set(c.generators) == set(c.morphisms)
    assert c.is_finite

    # Test loop morphisms.
    h = NamedMorphism(A, A, "h")
    d = Diagram(f, g, h)

    assert not d.is_finite
    assert h * h in d
    assert f * h in d
    assert f * h * h in d

    s = set(islice(d, 9))
    assert f * h in s
    assert h * h in s

    assert not d.is_hom_set_finite(A, B)
    assert not d.is_hom_set_finite(A, A)
    assert d.is_hom_set_finite(B, B)
    assert set(d.expanded_generators) == set([
        g, f, g * f, h, f * h, g * f * h, id_A, id_B, id_C])

    raises(ValueError, lambda: Diagram(f, g, f * h * h))
    raises(ValueError, lambda: Diagram({f: [], g: [], f * h * h: []}))
    raises(ValueError, lambda: Diagram(f, g, h, f * h * h))

    # Test expanded generators with composites among generators.
    d = Diagram({f: [], g: [], h: [], g * f: "big"})
    assert set(d.expanded_generators) == set([
        g, f, g * f, h, f * h, g * f * h, id_A, id_B, id_C])

    # Test expanded generators with two loops starting at the same
    # object.
    h_ = NamedMorphism(A, A, "h'")
    d = Diagram(f, g, h, h_)

    egen = set(d.expanded_generators)

    assert g * f in egen
    assert g * f * h in egen
    assert g * f * h_ in egen
    assert g * f * h * h_ in egen
    assert g * f * h_ * h in egen

    assert f in egen
    assert f * h in egen
    assert f * h_ in egen
    assert f * h * h_ in egen
    assert f * h_ * h in egen

    assert f * h * h not in egen
    assert f * h * h_ * h not in egen
    assert f * h_ * h_ * h not in egen

    # To test expanded generators with some complex cycles, we will
    # need more objects.
    A1 = Object("A1")
    A2 = Object("A2")
    A3 = Object("A3")
    A4 = Object("A4")
    A5 = Object("A5")
    A6 = Object("A6")
    A7 = Object("A7")

    # Test two cycles starting at the same object.
    f1 = NamedMorphism(A1, A2, "f1")
    f2 = NamedMorphism(A2, A3, "f2")
    f3 = NamedMorphism(A3, A2, "f3")
    f4 = NamedMorphism(A2, A4, "f4")
    f5 = NamedMorphism(A4, A2, "f5")
    f6 = NamedMorphism(A2, A5, "f6")
    d = Diagram(f1, f2, f3, f4, f5, f6)

    egen = set(d.expanded_generators)

    assert f2 * f1 in egen
    assert f2 * f3 * f2 * f1 in egen
    assert f2 * f5 * f4 * f1 in egen
    assert f2 * f5 * f4 * f3 * f2 * f1 in egen
    assert f2 * f3 * f2 * f5 * f4 * f1 in egen

    assert f2 * f3 * f2 * f3 * f2 * f1 not in egen
    assert f2 * f5 * f4 * f5 * f4 * f1 not in egen
    assert f2 * f5 * f4 * f3 * f2 * f5 * f4 * f3 * f2 * f1 not in egen

    # Test a cycle which starts in the middle of another cycle.
    f1 = NamedMorphism(A1, A2, "f1")
    f2 = NamedMorphism(A2, A3, "f2")
    f3 = NamedMorphism(A3, A4, "f3")
    f4 = NamedMorphism(A4, A5, "f4")
    f5 = NamedMorphism(A5, A6, "f5")
    f6 = NamedMorphism(A6, A3, "f6")
    f7 = NamedMorphism(A5, A2, "f7")
    f8 = NamedMorphism(A2, A7, "f8")
    d = Diagram(f1, f2, f3, f4, f5, f6, f7, f8)

    egen = set(d.expanded_generators)
    assert f8 * f1 in egen
    assert f8 * f7 * f4 * f3 * f2 * f1 in egen
    assert f8 * f7 * f4 * f3 * f6 * f5 * f4 * f3 * f2 * f1 in egen

    assert f8 * f7 * f4 * f3 * f2 * f7 * f4 * f3 * f2 * f1 not in \
           egen
    assert f8 * f7 * f4 * f3 * f6 * f5 * f4 * f3 * f6 * f5 * f4 * \
           f3 * f2 * f1 not in egen
    assert f8 * f7 * f4 * f3 * f6 * f5 * f4 * f3 * f2 * f7 * f4 * \
           f3 * f6 * f5 * f4 * f3 * f2 * f1 not in egen

    # Test two cycles containing the same edge.
    f1 = NamedMorphism(A1, A2, "f1")
    f2 = NamedMorphism(A2, A3, "f2")
    f3 = NamedMorphism(A3, A4, "f3")
    f4 = NamedMorphism(A3, A5, "f4")
    f5 = NamedMorphism(A5, A2, "f5")
    f6 = NamedMorphism(A3, A6, "f6")
    f7 = NamedMorphism(A6, A2, "f7")
    d = Diagram(f1, f2, f3, f4, f5, f6, f7)

    egen = set(d.expanded_generators)
    assert f3 * f2 * f1 in egen
    assert f3 * f2 * f5 * f4 * f2 * f1 in egen
    assert f3 * f2 * f7 * f6 * f2 * f1 in egen
    assert f3 * f2 * f5 * f4 * f2 * f7 * f6 * f2 * f1 in egen

    assert f3 * f2 * f5 * f4 * f2 * f5 * f4 * f2 * f1 not in egen
    assert f3 * f2 * f7 * f6 * f2 * f7 * f6 * f2 * f1 not in egen
    assert f3 * f2 * f5 * f4 * f2 * f7 * f6 * f2 * f5 * f4 * f2 * \
           f7 * f6 * f2 * f1 not in egen

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
    raises(ValueError, lambda: Implication(Diagram(f), Diagram(g, f, g * f)))

    # Test diff.
    assert imp.diff() == FiniteSet([g * f])

    h = NamedMorphism(C, A, "h")
    imp = Implication(premise, Diagram(h))
    assert imp.diff() == FiniteSet(h)
