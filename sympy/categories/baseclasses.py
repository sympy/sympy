from sympy.core import (Set, Basic, FiniteSet, EmptySet, Dict, Symbol,
                        Tuple)
from sympy.core.compatibility import iterable
from itertools import product
from sympy import cacheit

class Class(Set):
    r"""
    The base class for any kind of class in the set-theoretic sense.

    In axiomatic set theories, everything is a class.  A class which
    can be a member of another class is a set.  A class which is not a
    member of another class is a proper class.  The class `\{1, 2\}`
    is a set; the class of all sets is a proper class.

    This class is essentially a synonym for :class:`sympy.core.Set`.
    The goal of this class is to assure easier migration to the
    eventual proper implementation of set theory.
    """
    is_proper = False

class Object(Symbol):
    """
    The base class for any kind of object in an abstract category.

    While technically any instance of :class:`Basic` will do, this
    class is the recommended way to create abstract objects in
    abstract categories.
    """

class Morphism(Basic):
    """
    The base class for any morphism in an abstract category.

    In abstract categories, a morphism is an arrow between two
    category objects.  The object where the arrow starts is called the
    domain, while the object where the arrow ends is called the
    codomain.

    Two morphisms between the same pair of objects are considered to
    be the same morphisms.  To distinguish between morphisms between
    the same objects use :class:`NamedMorphism`.

    It is prohibited to instantiate this class.  Use one of the
    derived classes instead.

    See Also
    ========

    IdentityMorphism, NamedMorphism, CompositeMorphism
    """
    def __new__(cls, domain, codomain):
        raise(NotImplementedError(
            "Cannot instantiate Morphism.  Use derived classes instead."))

    @property
    def domain(self):
        """
        Returns the domain of the morphism.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = NamedMorphism(A, B, "f")
        >>> f.domain
        Object("A")

        """
        return self.args[0]

    @property
    def codomain(self):
        """
        Returns the codomain of the morphism.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = NamedMorphism(A, B, "f")
        >>> f.codomain
        Object("B")

        """
        return self.args[1]

    def compose(self, other):
        r"""
        Composes self with the supplied morphism.

        The order of elements in the composition is the usual order,
        i.e., to construct `g\circ f` use ``g.compose(f)``.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> g * f
        CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"),
        NamedMorphism(Object("B"), Object("C"), "g")))
        >>> (g * f).domain
        Object("A")
        >>> (g * f).codomain
        Object("C")

        """
        return CompositeMorphism(other, self)

    def __mul__(self, other):
        r"""
        Composes self with the supplied morphism.

        The semantics of this operation is given by the following
        equation: ``g * f == g.compose(f)`` for composable morphisms
        ``g`` and ``f``.

        See Also
        ========

        compose
        """
        return self.compose(other)

class IdentityMorphism(Morphism):
    """
    Represents an identity morphism.

    An identity morphism is a morphism with equal domain and codomain,
    which acts as an identity with respect to composition.

    Examples
    ========

    >>> from sympy.categories import Object, NamedMorphism, IdentityMorphism
    >>> A = Object("A")
    >>> B = Object("B")
    >>> f = NamedMorphism(A, B, "f")
    >>> id_A = IdentityMorphism(A)
    >>> id_B = IdentityMorphism(B)
    >>> f * id_A == f
    True
    >>> id_B * f == f
    True

    See Also
    ========

    Morphism
    """
    def __new__(cls, domain):
        return Basic.__new__(cls, domain, domain)

class NamedMorphism(Morphism):
    """
    Represents a morphism which has a name.

    Names are used to distinguish between morphisms which have the
    same domain and codomain: two named morphisms are equal if they
    have the same domains, codomains, and names.

    Examples
    ========

    >>> from sympy.categories import Object, NamedMorphism
    >>> A = Object("A")
    >>> B = Object("B")
    >>> f = NamedMorphism(A, B, "f")
    >>> f
    NamedMorphism(Object("A"), Object("B"), "f")
    >>> f.name
    'f'

    See Also
    ========

    Morphism
    """
    def __new__(cls, domain, codomain, name):
        if not name:
            raise ValueError("Empty morphism names not allowed.")

        return Basic.__new__(cls, domain, codomain, Symbol(name))

    @property
    def name(self):
        """
        Returns the name of the morphism.

        Examples
        ========
        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = NamedMorphism(A, B, "f")
        >>> f.name
        'f'

        """
        return self.args[2].name

class CompositeMorphism(Morphism):
    r"""
    Represents a morphism which is a composition of other morphisms.

    Two composite morphisms are equal if the morphisms they were
    obtained from (components) are the same and were listed in the
    same order.

    The arguments to the constructor for this class should be listed
    in diagram order: to obtain the composition `g\circ f` from the
    instances of :class:`Morphism` ``g`` and ``f`` use
    ``CompositeMorphism(f, g)``.

    Examples
    ========

    >>> from sympy.categories import Object, NamedMorphism, CompositeMorphism
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = NamedMorphism(A, B, "f")
    >>> g = NamedMorphism(B, C, "g")
    >>> g * f
    CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"),
    NamedMorphism(Object("B"), Object("C"), "g")))
    >>> CompositeMorphism(f, g) == g * f
    True

    """
    @staticmethod
    def _add_morphism(t, morphism):
        """
        Intelligently adds ``morphism`` to tuple ``t``.

        If ``morphism`` is a composite morphism, its components are
        added to the tuple.  If ``morphism`` is an identity, nothing
        is added to the tuple.

        No composability checks are performed.
        """
        if isinstance(morphism, CompositeMorphism):
            # ``morphism`` is a composite morphism; we have to
            # denest its components.
            return t + morphism.components
        elif isinstance(morphism, IdentityMorphism):
            # ``morphism`` is an identity.  Nothing happens.
            return t
        else:
            return t + Tuple(morphism)

    def __new__(cls, *components):
        if components and not isinstance(components[0], Morphism):
            # Maybe the user has explicitly supplied a list of
            # morphisms.
            return CompositeMorphism.__new__(cls, *components[0])

        normalised_components = Tuple()

        # TODO: Fix the unpythonicity.
        for i in xrange(len(components) - 1):
            current = components[i]
            following = components[i + 1]

            if not isinstance(current, Morphism) or \
                   not isinstance(following, Morphism):
                raise TypeError("All components must be morphisms.")

            if current.codomain != following.domain:
                raise ValueError("Uncomposable morphisms.")

            normalised_components = CompositeMorphism._add_morphism(
                normalised_components, current)

        # We haven't added the last morphism to the list of normalised
        # components.  Add it now.
        normalised_components = CompositeMorphism._add_morphism(
            normalised_components, components[-1])

        if not normalised_components:
            # If ``normalised_components`` is empty, only identities
            # were supplied.  Since they all were composable, they are
            # all the same identities.
            return components[0]
        elif len(normalised_components) == 1:
            # No sense to construct a whole CompositeMorphism.
            return normalised_components[0]

        return Basic.__new__(cls, normalised_components)

    @property
    def components(self):
        """
        Returns the components of this composite morphism.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> (g * f).components
        (NamedMorphism(Object("A"), Object("B"), "f"),
        NamedMorphism(Object("B"), Object("C"), "g"))

        """
        return self.args[0]

    @property
    def domain(self):
        """
        Returns the domain of this composite morphism.

        The domain of the composite morphism is the domain of its
        first component.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> (g * f).domain
        Object("A")

        """
        return self.components[0].domain

    @property
    def codomain(self):
        """
        Returns the codomain of this composite morphism.

        The codomain of the composite morphism is the codomain of its
        last component.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> (g * f).codomain
        Object("C")

        """
        return self.components[-1].codomain

    def flatten(self, new_name):
        """
        Forgets the composite structure of this morphism.

        If ``new_name`` is not empty, returns a :class:`NamedMorphism`
        with the supplied name, otherwise returns a :class:`Morphism`.
        In both cases the domain of the new morphism is the domain of
        this composite morphism and the codomain of the new morphism
        is the codomain of this composite morphism.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> (g * f).flatten("h")
        NamedMorphism(Object("A"), Object("C"), "h")

        """
        return NamedMorphism(self.domain, self.codomain, new_name)

    def __len__(self):
        """
        Returns the number of components of this
        :class:`CompositeMorphism`.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> len(g * f)
        2

        """
        return len(self.components)

    def __contains__(self, other):
        """
        Checks if ``other`` is a component of this
        :class:`CompositeMorphism`.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> g in (g * f)
        True

        """
        return other in self.components

    def __iter__(self):
        """
        Returns an iterator over the components of this
        :class:`CompositeMorphism`.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> from sympy import pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> pprint([m for m in g * f])
        [f:A-->B, g:B-->C]

        """
        return iter(self.components)

class Category(Basic):
    r"""
    An (abstract) category.

    A category [JoyOfCats] is a quadruple `\mbox{K} = (O, \hom, id,
    \circ)` consisting of

    * a (set-theoretical) class `O`, whose members are called
      `K`-objects,

    * for each pair `(A, B)` of `K`-objects, a set `\hom(A, B)` whose
      members are called `K`-morphisms from `A` to `B`,

    * for a each `K`-object `A`, a morphism `id:A\rightarrow A`,
      called the `K`-identity of `A`,

    * a composition law `\circ` associating with every `K`-morphisms
      `f:A\rightarrow B` and `g:B\rightarrow C` a `K`-morphism `g\circ
      f:A\rightarrow C`, called the composite of `f` and `g`.

    Composition is associative, `K`-identities are identities with
    respect to composition, and the sets `\hom(A, B)` are pairwise
    disjoint.

    This class knows nothing about its objects and morphisms.
    Concrete cases of (abstract) categories should be implemented as
    classes derived from this one.

    Certain instances of :class:`Diagram` can be asserted to be
    commutative in a :class:`Category` by supplying the argument
    ``commutative_diagrams`` in the constructor.

    Examples
    ========

    >>> from sympy.categories import Object, NamedMorphism, Diagram, Category
    >>> from sympy import FiniteSet
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = NamedMorphism(A, B, "f")
    >>> g = NamedMorphism(B, C, "g")
    >>> d = Diagram(f, g)
    >>> K = Category("K", commutative_diagrams=[d])
    >>> K.commutative_diagrams == FiniteSet((d,))
    True

    See Also
    ========
    Diagram
    """
    def __new__(cls, name, objects=EmptySet(), commutative_diagrams=EmptySet()):
        if not name:
            raise ValueError("A Category cannot have an empty name.")

        new_category = Basic.__new__(cls, Symbol(name), Class(objects),
                                     FiniteSet(commutative_diagrams))
        return new_category

    @property
    def name(self):
        """
        Returns the name of this category.

        Examples
        ========

        >>> from sympy.categories import Category
        >>> K = Category("K")
        >>> K.name
        'K'

        """
        return self.args[0].name

    @property
    def objects(self):
        """
        Returns the class of objects of this category.

        Examples
        ========

        >>> from sympy.categories import Object, Category
        >>> from sympy import FiniteSet
        >>> A = Object("A")
        >>> B = Object("B")
        >>> K = Category("K", FiniteSet(A, B))
        >>> K.objects
        Class({Object("A"), Object("B")})

        """
        return self.args[1]

    @property
    def commutative_diagrams(self):
        """
        Returns the :class:`FiniteSet` of diagrams which are known to
        be commutative in this category.

        >>> from sympy.categories import Object, NamedMorphism, Diagram, Category
        >>> from sympy import FiniteSet
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> K = Category("K", commutative_diagrams=[d])
        >>> K.commutative_diagrams == FiniteSet((d,))
        True

        """
        return self.args[2]

    def hom(self, A, B):
        raise NotImplementedError(
            "hom-sets are not implemented in Category.")

    def all_morphisms(self):
        raise NotImplementedError(
            "Obtaining the class of morphisms is not implemented in Category.")

class Diagram(Basic):
    """
    This class represents a diagram in a certain category.

    Overview
    ========

    Informally, a diagram is a collection of objects of a category and
    certain morphisms between them.  Identity morphisms, as well as
    all composites of morphisms included in the diagram, belong to the
    diagram.  For a more formal approach to this notion see
    [Pare1970].

    A :class:`Diagram` represents a mapping between morphisms and the
    :class:`FiniteSet`'s of their properties.

    >>> from sympy.categories import Object, NamedMorphism, Diagram
    >>> from sympy import FiniteSet, pprint, default_sort_key, Dict
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = NamedMorphism(A, B, "f")
    >>> g = NamedMorphism(B, C, "g")
    >>> d = Diagram(f, g, g * f)
    >>> pprint(d)
    {id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: EmptySet(), f:A-->B: Em
    ptySet(), g:B-->C: EmptySet()}
    >>> pprint(list(d.hom(A, B)))
    [f:A-->B]

    .. note::

       No checks are carried out of whether the supplied objects and
       morphisms do belong to the same category.

    All possible morphism compositions are considered to belong to the
    diagram.

    >>> d == Diagram(f, g)
    True
    >>> g * f in Diagram(f, g)
    True
    >>> pprint(list(d.hom(A, C)))
    [g*f:A-->C]

    Morphisms can be assigned some properties which can later be
    retrieved using the :class:`dict`-like interface of the
    :class:`Diagram`:

    >>> d = Diagram({f: [], g: "unique"})
    >>> pprint(d)
    {id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: EmptySet(), f:A-->B: Em
    ptySet(), g:B-->C: {unique}}
    >>> d[g]
    {unique}
    >>> d[f]
    EmptySet()

    The set of properties of such a composite morphism is the
    intersection of the sets of properties of its components:

    >>> d = Diagram({f: ["unique", "subobject"], g: ["unique"]})
    >>> d[g * f]
    {unique}

    One can construct subdiagrams from sets of objects.

    >>> subd = d.subdiagram_from_objects([A, B])
    >>> pprint(subd)
    {id:A-->A: EmptySet(), id:B-->B: EmptySet(), f:A-->B: {subobject, unique}}
    >>> d.is_subdiagram(subd)
    True
    >>> subd <= d
    True

    More Details
    ============

    From the morphisms supplied to its constructor, :class:`Diagram`
    constructs a set of morphisms which is referred to as
    _generators_.  :class:`Diagram` adds the identities of all objects
    referred to by the supplied morphism.

    >>> pprint(Diagram(f))
    {id:A-->A: EmptySet(), id:B-->B: EmptySet(), f:A-->B: EmptySet()}
    >>> pprint(sorted(Diagram(f).generators, key=default_sort_key))
    [id:A-->A, id:B-->B, f:A-->B]
    >>> pprint(Diagram(f).generators_properties)
    {id:A-->A: EmptySet(), id:B-->B: EmptySet(), f:A-->B: EmptySet()}

    On the other hand, property-less composites of other generators
    are removed::

    >>> pprint(Diagram(f, g, g * f))
    {id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: EmptySet(), f:A-->B: Em
    ptySet(), g:B-->C: EmptySet()}

    To get all morphisms of a :class:`Diagram`, one can either use
    ``Diagram.morphisms``, or directly enumerate the :class:`Diagram`.

    >>> d = Diagram(g, f)
    >>> g * f in d.generators
    False
    >>> g * f in d.morphisms
    True
    >>> g * f in d
    True
    >>> pprint(sorted(d.morphisms, key=default_sort_key))
    [g*f:A-->C, id:A-->A, id:B-->B, id:C-->C, f:A-->B, g:B-->C]
    >>> pprint(sorted(d, key=default_sort_key))
    [g*f:A-->C, id:A-->A, id:B-->B, id:C-->C, f:A-->B, g:B-->C]
    >>> len(d)
    6

    :class:`Diagram` can cope with infinite diagrams, i.e., diagrams
    which contain an infinite number of morphisms.  Such diagrams can
    be obtained by creating cycles.

    >>> d = Diagram(g, f)
    >>> d.is_finite
    True
    >>> f_ = NamedMorphism(B, A, "f'")
    >>> d = Diagram(f, f_)
    >>> d.is_finite
    False

    It is also possible to check the finiteness of a certain
    `\hom`-set::

    >>> d = Diagram(g, f)
    >>> d.is_hom_set_finite(A, C)
    True
    >>> d = Diagram(f, f_)
    >>> d.is_hom_set_finite(A, A)
    False

    Further, it is possible to check the emptiness of a certain
    `\hom`-set::

    >>> d = Diagram(f)
    >>> d.is_hom_set_empty(A, B)
    False
    >>> d.is_hom_set_empty(B, A)
    True
    >>> d = Diagram(f, f_)
    >>> d.is_hom_set_empty(B, A)
    False

    Computing the length of an infinite diagram is an error.  However,
    one can still enumerate the morphisms of the diagram.  This
    enumeration will never stop, so one can explicitly retrieve a
    certain number of morphisms.

    >>> from itertools import islice
    >>> slice = islice(d, 8)
    >>> pprint(sorted(slice, key=default_sort_key))
    [f'*f:A-->A, f*f':B-->B, f*f'*f:A-->B, f'*f*f':B-->A, id:A-->A, id:B-->B, f:A-
    ->B, f':B-->A]

    .. note::

       No particular ordering is guaranteed for any of the lists
       produced by :class:`Diagram`.

    Even More Details.  Graph Theory
    ================================

    Obviously, a diagram is a graph-like structure: it has objects
    connected with arrows.  A diagram is indeed very similar to a
    directed multigraph, because there can be more than one morphism
    between the same two objects.  Note however that, formally
    speaking, only the generators of the diagram form a multigraph;
    that's because the total number of morphisms belonging to the
    diagram may be infinite.

    :class:`Diagram` offers a couple graph-theoretical tools which may
    prove useful in handling infinite cases.

    Observe that there are only two cases in which a diagram can be
    infinite:

    * the multigraph defined by the generators has cycles, or

    * there are non-identity loop morphisms among the generators
      (i.e., non-identity morphisms with the same domain and
      codomain).

    When the multigraph defined by the generators has cycles, it may
    be useful to know the strongly connected components of that
    multigraph.

    Consider the following infinite diagram::

    >>> h = NamedMorphism(C, A, "h")
    >>> D = Object("D")
    >>> k = NamedMorphism(D, A, "k")
    >>> d = Diagram(f, g, h, k)
    >>> pprint(d)
    {id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: EmptySet(), id:D-->D: E
    mptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet(), h:C-->A: EmptySet(), k:D-
    ->A: EmptySet()}

    The multigraph of the generators has two strongly connected
    components in this case: `\{A, B, C\}` and `\{D\}`.

    >>> component1 = d.subdiagram_from_objects([A, B, C])
    >>> d.objects_components == Dict({A: component1, B: component1, C: component1,
    ...                               D: D})
    True

    Remark that one-vertex strongly-connected components are not
    represented as :class:`Diagram`'s, but just as objects.

    A step further is computing the condensation of the diagram.  The
    condensation of a digraph is the graph obtained by contracting
    every strongly connecting component of the digraph to a vertex and
    by connecting two such new vertices with an edge if there is an
    edge between two original vertices located in the corresponding
    strongly connected components [WCon].  This condensation of a
    digraph is a directed acyclic graph [WCon].

    If in a :class:`Diagram` there is a morphism `f:A\rightarrow B`
    with properties ``props``, the condensation will contain a
    morphism between the connected components containing `A` and `B`;
    this morphism will have the same name and the same properties as
    the original morphism.

    Consider the following slight modification of the previous
    example::

    >>> d = Diagram({f: [], g: [], h: [], k: "blaster"})
    >>> component1 = d.subdiagram_from_objects([A, B, C])
    >>> condensation_k = NamedMorphism(D, component1, "k")
    >>> d.condensation == Diagram({condensation_k: "blaster"})
    True

    Another powerful, graph-theory based, infinite diagram analysis
    tool is the notion of an _expanded generator_.  Consider the
    following diagram with a loop morphism::

    >>> h = NamedMorphism(A, A, "h")
    >>> d = Diagram(f, g, h)

    It is obviously infinite, because it includes the loop morphism
    `h`.  Note, however, that the set of morphisms of this diagram has
    a very simple structure::

    >>> pprint(set(islice(d, 18)))
    set([h*h:A-->A, f*h:A-->B, g*f:A-->C, h*h*h:A-->A, f*h*h:A-->B, g*f*h:A-->C, h
    *h*h*h:A-->A, f*h*h*h:A-->B, g*f*h*h:A-->C, h*h*h*h*h:A-->A, f*h*h*h*h:A-->B,
    g*f*h*h*h:A-->C, id:A-->A, id:B-->B, id:C-->C, h:A-->A, f:A-->B, g:B-->C])

    That is, we have the morphisms `f`, `g`, and `g\circ f` and then
    we obtain the other morphisms by "appending" arbitrary long "runs"
    of `h` to the end of `f` or `g\circ f`.  Remark that all morphisms
    in the family `[f\circ h] = \{f\circ h^n\mid n\in\mathbb{N}\}`
    (where `\mathbb{N}` does _not_ include zero) are very tightly
    related.  Firstly, either the whole family `[f\circ h]` is
    included in a diagram, or neither of these morphisms belongs to a
    diagram.  Secondly, since the properties of composite morphisms
    are computed as intersections of the properties of the components,
    all morphisms in `[f\circ h]` will have the _same_ set of
    properties.  That is, when analysing infinite diagrams, we are
    generally more interested in considering such families of
    morphisms rather than attempting to handle the infinite multitude
    as a whole.

    An _expanded generator_ is a the shortest composite morphism in a
    family of morphisms of the form described in the previous
    paragraph.

    Let's check that the expanded generators of the diagram ``d`` are
    those which we expected::

    >>> pprint(set(d.expanded_generators))
    set([f*h:A-->B, g*f:A-->C, g*f*h:A-->C, id:A-->A, id:B-->B, id:C-->C, h:A-->A,
    f:A-->B, g:B-->C])

    .. note::

       While we have only discussed loops, the same reasoning applies
       to the situations when there are cycles in the multigraph of
       generators.

    Expanded generators have a number of properties which are very
    useful in practise.  Notice firstly that an expanded generator
    does not pass through any loop or cycle in the multigraph of
    generators more than once.  An immediate consequence of this
    observation and the fact that we only allow finite sets of
    generators is that the set of expanded generators is always
    _finite_.

    .. note::

       :class:`Diagram` will accept as generators only those morphism
       which can be expanded generators::

       >>> Diagram(f, g, f * h * h)
       Traceback (most recent call last):
          ...
       ValueError: All generators must be expanded generators.

    References
    ==========

    [Pare1970] B. Pareigis: Categories and functors.  Academic Press,
    1970.

    [WCon] http://en.wikipedia.org/wiki/Condensation_%28graph_theory%29
    """
    def  __new__(cls, *args):
        """
        Creates a new instance of :class:`Diagram`.

        If no arguments are supplied, an empty diagram is created.

        If the first argument is a dictionary, it is interpreted as a
        mapping between the morphisms which form the diagram and their
        properties.  If the first argument is an iterable, but not a
        dictionary, it is interpreted as a collection of morphisms,
        each of which has no properties.  Otherwise, all arguments are
        interpreted as morphisms.

        Raises a :class:`ValueError` if it encounters a morphism which
        is not an expanded generator.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> from sympy.categories import IdentityMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> IdentityMorphism(A) in d
        True
        >>> g * f in d
        True
        >>> d = Diagram({f: [], g: [], g * f: "unique"})
        >>> d[g * f]
        {unique}

        """
        generators = {}
        objects = set([])

        if args:
            first_arg = args[0]

            if isinstance(first_arg, dict) or isinstance(first_arg, Dict):
                # The user has supplied a dictionary of morphisms and
                # their properties.  Assure that the properties are
                # converted to ``FiniteSet``'s.
                generators = {}
                for morphism, props in first_arg.items():
                    if not Diagram._is_expanded_generator(morphism):
                        raise ValueError(
                            "All generators must be expanded generators.")

                    generators[morphism] = FiniteSet(props)
            elif iterable(first_arg) and not isinstance(
                first_arg, CompositeMorphism):
                # The user has supplied a list of morphisms, none of
                # which have any properties.
                empty = EmptySet()
                for morphism in first_arg:
                    if not Diagram._is_expanded_generator(morphism):
                        raise ValueError(
                            "All generators must be expanded generators.")

                    generators[morphism] = empty
            else:
                # Attempt to interpret ``args`` as a list of
                # morphisms.
                return Diagram(args)

        for morphism, props in generators.items():
            # Drop property-less composites of other morphisms in the
            # generators.
            if isinstance(morphism, CompositeMorphism) and not props:
                if all(component in generators for component
                       in morphism.components):
                    del generators[morphism]

                    continue

            # This morphism is all right, register the objects it
            # involves.
            objects.update([morphism.domain, morphism.codomain])

        # Add identity morphisms for those objects for which they have
        # not been already added.  We need this to be able to have
        # isolated objects in diagrams.
        #
        # We need one-object diagrams for some technical reasons, for
        # example, in diagram drawing.  This does not violate the
        # actual notion of a ``Diagram``, however.  While it is indeed
        # somewhat unusual to have one-object diagrams, it is
        # perfectly fine from the formal perspective.
        for obj in objects:
            id_obj = IdentityMorphism(obj)
            if id_obj not in generators:
                generators[id_obj] = FiniteSet()

        return Basic.__new__(cls, Dict(generators), FiniteSet(objects))

    @property
    def generators(self):
        """
        Returns the list of generators of this diagram.

        .. note::

           No particular ordering of generators is guaranteed.

        Returns the :class:`Dict` mapping the morphisms included in
        this :class:`Diagram` to their properties.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import pprint, default_sort_key
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = NamedMorphism(A, B, "f")
        >>> d = Diagram(f)
        >>> sorted_gens = sorted(d.generators, key=default_sort_key)
        >>> pprint(sorted_gens, use_unicode=False)
        [id:A-->A, id:B-->B, f:A-->B]

        """
        return self.args[0].keys()

    @property
    def generators_properties(self):
        """
        Returns the dictionary mapping the generators of this
        :class:`Diagram` to the :class:`FiniteSet`'s of their
        properties.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = NamedMorphism(A, B, "f")
        >>> d = Diagram(f)
        >>> pprint(d.generators_properties, use_unicode=False)
        {id:A-->A: EmptySet(), id:B-->B: EmptySet(), f:A-->B: EmptySet()}

        """
        return self.args[0]

    @cacheit
    def _build_generator_adj_lists(self):
        """
        Returns the adjacency lists of the graph defined by the
        generators.

        This will skip identities and loop morphisms.
        """
        adj_lists = {}

        for obj in self.objects:
            adj_lists[obj] = []

        for m in self.generators:
            if m.domain == m.codomain:
                # The generators may contain some identity morphisms,
                # which were explicitly supplied with properties.
                # Yet, since they are identities, they do not add new
                # morphisms to the diagram, so they should be skipped
                # here.
                #
                # We will also skip loop morphisms here.  Loop
                # morphisms are always handled in a special way,
                # without looking at the adjacency lists.
                continue
            adj_lists[m.domain].append(m.codomain)

        return adj_lists

    @staticmethod
    def _dfs_markup_postorder(obj, adj_lists, component_idx, obj_indices,
                              postorder_indices=None, counter_ref=None):
        """
        Traverses the graph given by the adjacency lists ``adj_lists``
        starting from the object ``obj``.  On traversal it marks all
        visited objects by setting to ``component_idx`` in
        ``object_indices``.  In postorder, it marks the vertices with
        the values of the counter ``counter_ref`` in
        ``postorder_indices``.

        ``counter_ref`` is a one-element list, containing the value of
        the counter.  This trick allows to have one global counter,
        unique for all recursive call of this function, without
        polluting ``self`` with additional instance variables.
        """
        obj_indices[obj] = component_idx

        for adj_obj in adj_lists[obj]:
            if obj_indices[adj_obj] is None:
                Diagram._dfs_markup_postorder(
                    adj_obj, adj_lists, component_idx, obj_indices,
                    postorder_indices, counter_ref)

        if postorder_indices:
            postorder_indices[obj] = counter_ref[0]
            counter_ref[0] += 1

    @cacheit
    def _build_strongly_connected_components(self):
        """
        Finds the strongly connected components of the directed
        multigraph defined by the generators of this diagram.  Returns
        a :class:`Dict` mapping each object to the number of the
        connected component it makes part of.

        This method returns a :class:`Dict` instead of the built-in
        :class:`dict` because we need to return an immutable object to
        make use of SymPy's cache.

        This implements Kosaraju's algorithm [Wikipedia].  [Sedg2002]
        contains another explanation and some theoretical background.
        I have followed [Sedg2002] in writing this method.  I chose
        Kosaraju's algorithm instead of Tarjan's algorithm [WTar]
        because we have the edges of the graph (the generators), and
        not its adjacency matrix.

        References
        ==========

        [WKos] http://en.wikipedia.org/wiki/Kosaraju%27s_algorithm

        [Sedg2002] Robert Sedgewick, Algorithms in C++ Graph
        Algorithms, Part 5.  Third Edition.  Addison-Wesley, 2002.
        Section Digraphs and DAGs, Listing 19.10 and the corresponding
        explanations.

        [WTar] http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm

        """
        # This dictionary will eventually contain the result, i.e.,
        # will map objects to numbers of the corresponding component
        # indices.
        obj_indices = {}

        # These dictionary will store the adjacency lists of the graph
        # of generators with reversed arrows.
        rev_adj_lists = {}

        # This dictionary will store the numbers of vertices, as they
        # are traversed in postorder by DFS.
        postorder_indices = {}

        for obj in self.objects:
            obj_indices[obj] = None
            postorder_indices[obj] = None
            rev_adj_lists[obj] = []

        # Build the adjacency lists.
        adj_lists = self._build_generator_adj_lists()

        for v in adj_lists:
            for w in adj_lists[v]:
                rev_adj_lists[w].append(v)

        # This one-element list holds a global counter for the DFS
        # traversal (cf. the docstring of
        # ``Diagram._dfs_markup_postorder``).
        counter_ref = [0]

        # At first, DFS the reverse graph.
        #
        # Since the original graph and the graph with reversed arrows
        # have the same strongly connected components, it doesn't
        # matter which one of them we traverse first.  However, let's
        # follow [Sedg2002] as closely as we can.
        #
        # Note that we don't need to start the traversal from sources
        # only.  Even if start the traversal from a vertex `w`, which
        # is not a source, and thus don't reach the vertex `v` for
        # which the edge `(v, w)` exists, we will still arrive at `v`
        # at a later time, and will mark it with a _greater_ postorder
        # index, which is right anyway.
        for obj in self.objects:
            if obj_indices[obj] is None:
                # In this run, we only use ``obj_indices`` to assure
                # that DFS doesn't visit the same vertices twice; that
                # is, the value of ``component_idx`` is of no
                # importance.
                Diagram._dfs_markup_postorder(
                    obj, rev_adj_lists, -1, obj_indices,
                    postorder_indices, counter_ref)

        for obj in self.objects:
            obj_indices[obj] = None

        # We would like now to traverse the objects in decreasing
        # direction of the postorder markup value in
        # ``postorder_indices``.
        nobjects = len(self.objects)
        ordered_objects = [None] * nobjects
        for obj in self.objects:
            idx = nobjects - 1 - postorder_indices[obj]
            ordered_objects[idx] = obj

        # We will use this to count through the detected strongly
        # connected components.
        component_idx = 0

        # Now, DFS the original graph (now the direction of edges will
        # be _reversed_ as compared to the previous run).
        #
        # Note that direction of traversal here is a bit different
        # from that in the Listing 19.10 in [Sedg2002]; however, it
        # does correspond to what is described in the text and on
        # [WKos].
        for obj in ordered_objects:
            if obj_indices[obj] is None:
                Diagram._dfs_markup_postorder(
                    obj, adj_lists, component_idx, obj_indices)

                component_idx += 1

        return Dict(obj_indices)

    @cacheit
    def _get_loop_morphisms(self):
        """
        Returns a tuple containing all loop morphisms, excluding
        identities.

        The purpose of having this method is to benefit from SymPy's
        cache in order to avoid extracting loop morphisms anew,
        whenever it should be necessary.
        """
        return tuple(m for m in self.generators
                     if not isinstance(m, IdentityMorphism) and
                     (m.domain == m.codomain))

    @property
    @cacheit
    def is_finite(self):
        """
        Checks whether this diagram contains a finite number of
        morphisms.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = NamedMorphism(A, B, "f")
        >>> d = Diagram(f)
        >>> d.is_finite
        True
        >>> f_ = NamedMorphism(B, A, "f'")
        >>> d = Diagram(f, f_)
        >>> d.is_finite
        False

        """
        if self._get_loop_morphisms():
            # A diagram with loop morphisms is always infinite.
            return False

        obj_components = self._build_strongly_connected_components()

        # If the graph defined by the generators of this diagram has a
        # strongly connected component which includes more that one
        # vertex, this graph has cycles.  Let's see if there are such
        # components.

        # ``used_index[i] == True`` if there is a vertex which
        # belongs to the connected component with index ``i``.
        used_index = [False] * len(self.objects)
        for obj in obj_components:
            component_idx = obj_components[obj]
            if used_index[component_idx]:
                # ``component_idx`` is a connected component with at
                # least two vertices.
                return False
            else:
                used_index[component_idx] = True

        return True

    @property
    @cacheit
    def objects_components(self):
        """
        Returns a :class:`Dict` mapping objects to the strongly
        connected component they belong to.  Strongly connected
        components are represented as :class:`Diagram`'s.

        One-vertex strongly connected components are represented using
        the corresponding objects instead of one-object diagrams.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import Dict
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> D = Object("D")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> h = NamedMorphism(C, A, "h")
        >>> k = NamedMorphism(D, A, "k")
        >>> d = Diagram(f, g, h, k)
        >>> component1 = d.subdiagram_from_objects([A, B, C])
        >>> d.objects_components == Dict({A: component1, B: component1,
        ...                               C: component1, D: D})
        True

        """
        # Build the actual sets of objects corresponding to each
        # strongly connected component.
        component_indices = self._build_strongly_connected_components()

        # ``component_objects[i]`` will hold the list of objects
        # contained in the strongly connected component with index
        # ``i``.
        component_objects = {}
        for obj, idx in component_indices.items():
            if idx not in component_objects:
                component_objects[idx] = [obj]
            else:
                component_objects[idx].append(obj)

        # For each nontrivial component, set up a diagram from the
        # objects.  Also set up the mapping between objects and the
        # corresponding diagram.
        object_diagrams = {}
        for idx, objects in component_objects.items():
            if len(objects) > 1:
                diagram = self.subdiagram_from_objects(objects)
                for obj in objects:
                    object_diagrams[obj] = diagram
            else:
                # Every object belongs to a strongly connected
                # component, so if ``objects`` does not have more than
                # one element, it has exactly one element.
                obj = objects[0]
                object_diagrams[obj] = obj

        return Dict(object_diagrams)

    @property
    @cacheit
    def condensation(self):
        """
        Returns the :class:`Diagram` which represents the condensation
        of this :class:`Diagram`.

        The condensation of a digraph is the graph obtained by
        contracting every strongly connecting component of the digraph
        to a vertex and by connecting two such new vertices with an
        edge if there is an edge between two original vertices located
        in the corresponding strongly connected components [WCon].
        This condensation of a digraph is a directed acyclic graph
        [WCon].

        This method generalises the notion for directed multigraphs.

        If in the current diagram there is a morphism `f:A\rightarrow
        B` with properties ``props``, the condensation will contain a
        morphism between the connected components containing `A` and
        `B`, this morphism will have the name ``f`` and it will have
        the same properties as the original morphism.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> D = Object("D")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> h = NamedMorphism(C, A, "h")
        >>> k = NamedMorphism(D, A, "k")
        >>> d = Diagram({f: [], g: [], h: [], k: "blaster"})
        >>> component1 = d.subdiagram_from_objects([A, B, C])
        >>> condensation_k = NamedMorphism(D, component1, "k")
        >>> d.condensation == Diagram({condensation_k: "blaster"})
        True

        References
        ==========

        [WCon] http://en.wikipedia.org/wiki/Condensation_%28graph_theory%29
        """
        objects_components = self.objects_components

        new_generators = {}

        for generator, props in self.generators_properties.items():
            dom_component = objects_components[generator.domain]
            cod_component = objects_components[generator.codomain]
            if dom_component == cod_component:
                # ``generator`` belongs to a connected component.
                continue
            else:
                new_generator = NamedMorphism(dom_component, cod_component,
                                              generator.name)
                new_generators[new_generator] = props

        # Add an identity morphism for every component, so that all
        # components make it into the diagram.
        for component in objects_components.values():
            new_generators[IdentityMorphism(component)] = []

        return Diagram(new_generators)

    def _morphisms_in_rounds(self):
        """
        Returns a generator that enumerates all morphisms belonging to
        the diagram, and yields the number of the generation round as
        well.
        """
        for m in self.generators:
            yield (m, 1)

        # Next we will generate morphisms in rounds; at each round we
        # will generate morphisms from combinations of a fixed number
        # of generators.

        # If this diagram is finite, we want to only go on until there
        # are composites.  Otherwise, we just go on and on.
        #
        # The point of having a lambda is to avoid checking
        # ``is_finite`` on every iteration.
        if self.is_finite:
            halt_condition = lambda ncomposites: ncomposites == 0
        else:
            halt_condition = lambda x: False

        # Here we will store the previous round of morphisms.
        # Remember, we have just yielded all generators, so that's the
        # previous round here.
        previous_morphisms = [gen for gen in self.generators
                              if not isinstance(gen, IdentityMorphism)]

        round = 2

        ncomposites = 1
        while not halt_condition(ncomposites):
            ncomposites = 0

            # Here we will store the current round of morphisms.
            current_morphisms = []

            # We need to yield all possible combinations between
            # non-trivial generators and ``previous_morphisms``.
            for (previous, current) in product(previous_morphisms,
                                               self.generators):
                if not isinstance(current, IdentityMorphism) and \
                       (current.domain == previous.codomain):
                    composite = current * previous

                    # The generators may already include this
                    # composite, if has been given some properties.
                    if composite not in self.generators_properties:
                        yield (composite, round)
                        ncomposites += 1
                        current_morphisms.append(composite)

            previous_morphisms = current_morphisms
            round += 1

    @property
    def morphisms(self):
        """
        Returns a generator that provides all morphisms belonging to
        this diagram.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import pprint, default_sort_key
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> pprint(sorted(d.morphisms, key=default_sort_key))
        [g*f:A-->C, id:A-->A, id:B-->B, id:C-->C, f:A-->B, g:B-->C]

        """
        for (morphism, round) in self._morphisms_in_rounds():
            yield morphism

    @staticmethod
    def _is_expanded_generator(morphism):
        """
        Checks whether the ``morphism`` is an expanded generator.  Any
        non-composite morphism is an expanded generator.

        For examples of composite expanded generators, consider the
        following situation:

                +----> C ----+
                |            |
                |           \|/
        A ----> B            D ----> F
               /|\           |
                |            |
                +----- E ----+

        A->B->C->D->F is an expanded generator, so is
        A->B->C->D->E->B->C->D->F.  However, the morphisms
        A->B->C->D->E->B and A->B->C->D->E->B->C->D->E->B->C->D->F are
        not expanded generators.

        See the docstring of this class for more details.
        """
        if not isinstance(morphism, CompositeMorphism):
            return True

        # To see if ``composite`` is an expanded generator, we will we
        # count how much it passes through the objects it involves.
        # By taking a look at the examples in the docstring it is
        # rather easy to see the following:
        #
        #   * a non-loop expanded generator will pass at most twice
        #     through each object it involves, except for endpoints,
        #     thought which it passes only once;
        #
        #   * a loop expanded generator will pass at most twice
        #     through each object it involves.
        #
        # Thus, any expanded generator will not pass more than twice
        # through each object it involves.

        objects_involved = Diagram._get_involved_objects(morphism)

        object_pass_count = dict.fromkeys(objects_involved, 0)
        object_pass_count[morphism.domain] = 1
        for component in morphism:
            object_pass_count[component.codomain] += 1

        return all(object_pass_count[obj] <= 2 for obj in objects_involved)

    @property
    def expanded_generators(self):
        """
        Enumerates the expanded generators of this :class:`Diagram`.

        Expanded generators include all morphisms which are obtained
        by composing generators, without repeating any one of them.
        The set of expanded generators is always finite.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> h = NamedMorphism(A, A, "h")
        >>> d = Diagram(f, g, h)
        >>> pprint(set(d.expanded_generators))
        set([f*h:A-->B, g*f:A-->C, g*f*h:A-->C, id:A-->A, id:B-->B, id:C-->C, h:A-->A,
        f:A-->B, g:B-->C])

        """
        # Any composite expanded generator passes through every object
        # it involves at most twice.  The same applies to the
        # generators it involves: it passes through each of them at
        # most twice.  Therefore, any expanded morphism is of a length
        # less or equal to twice the number of generators.
        #
        # We know that ``self.morphisms`` produces morphisms in rounds
        # of the same length in terms of generators, i.e., at first,
        # it produces the generators, then all composites of pairs of
        # generators, then all composites of triples of generators,
        # etc.  This means that once we have arrived at a non-expanded
        # morphism of length greater than the limit length, we have
        # exhausted all expanded generator morphisms.

        # Keep in mind that ``CompositeMorphism`` always stores the
        # flattened list of components, so we cannot just do ``2 *
        # len(self.generators)``.
        limit_length = 0
        for gen in self.generators:
            if isinstance(gen, CompositeMorphism):
                limit_length += len(gen)
            else:
                limit_length += 1
        limit_length *= 2

        for morphism in self.morphisms:
            if self._is_expanded_generator(morphism):
                yield morphism
            elif len(morphism) > limit_length:
                return

    @property
    def objects(self):
        """
        Returns the :class:`FiniteSet` of objects that appear in this
        diagram.

        Examples
        ========
        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> d.objects
        {Object("A"), Object("B"), Object("C")}

        """
        return self.args[1]

    def hom(self, A, B):
        """
        Returns a generator which enumerates all morphisms between the
        objects ``A`` and ``B``.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> pprint(list(d.hom(A, C)), use_unicode=False)
        [g*f:A-->C]

        See Also
        ========
        Object, Morphism
        """
        # For infinite diagrams, we rely on the fact that
        # ``morphisms`` assures a sufficiently BFS-like traversal.
        return (m for m in self.morphisms if (m.domain == A) and
                (m.codomain == B))

    @staticmethod
    def _get_involved_objects(morphism):
        """
        Given a morphism, returns the :class:`set` of objects it
        involves.

        A simple (not composite) morphism involves its domain and
        codomain.  A composite morphism involves the domains and the
        codomains of each of its components.
        """
        objects = set([morphism.domain, morphism.codomain])

        if isinstance(morphism, CompositeMorphism):
            for component in morphism.components:
                objects.add(component.domain)

        return objects

    def is_hom_set_finite(self, A, B):
        """
        Returns ``True`` if the `\hom`-set `\hom(A, B)` contains a
        finite number of morphisms and ``False`` otherwise.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = NamedMorphism(A, B, "f")
        >>> d = Diagram(f)
        >>> d.is_hom_set_finite(A, B)
        True
        >>> f_ = NamedMorphism(B, A, "f'")
        >>> d = Diagram(f, f_)
        >>> d.is_hom_set_finite(A, B)
        False

        """
        objects_components = self.objects_components
        component_A = objects_components[A]
        component_B = objects_components[B]
        condensation = self.condensation

        condensation_loop_objects = set([objects_components[m.domain] for m
                                         in self._get_loop_morphisms()])

        # If there is a path from ``component_A`` to ``component_B``
        # in the condensation that passes through a non-trivial
        # strongly connected component, `\hom(A, B)` is infinite.
        # This is so because the a nontrivial strongly connected
        # component contains at least one cycle.
        #
        # Further, if there is a path from ``component_A`` to
        # ``component_B`` in the condensation that passes through an
        # object with a loop, `\hom(A, B)` is infinite, because one
        # can use that loop morphism to produce as many
        # (theoretically) different morphisms as one wishes.
        #
        # Note that the condensation is surely acyclic, so it's
        # finite.  We will therefore abuse the property ``morphisms``
        # a bit, for the sake of clarity (hopefully).

        for m in condensation.hom(component_A, component_B):
            components_involved = self._get_involved_objects(m)

            if components_involved & condensation_loop_objects:
                # ``m`` passes through a loop object.
                return False

            for component in components_involved:
                if isinstance(component, Diagram):
                    # ``m`` passes through a nontrivial strongly
                    # connected component.
                    return False

        return True

    def is_hom_set_empty(self, A, B):
        """
        Returns ``True`` if the `\hom`-set `\hom(A, B)` contains at
        least one morphism and ``False`` otherwise.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> d.is_hom_set_empty(A, C)
        False
        >>> d.is_hom_set_empty(C, A)
        True

        """
        # We will start a DFS at ``A`` and see if we can get to ``B``.

        adj_lists = self._build_generator_adj_lists()

        # ``visited[obj] == True`` means that ``obj`` was visited
        # during the DFS traversal of the graph.
        visited = {}
        for obj in self.objects:
           visited[obj] = None

        Diagram._dfs_markup_postorder(A, adj_lists, True, visited)

        # If ``B`` was visited (``visited[B] == True``), the
        # `\hom`-set `\hom(A, B)` is not empty, so we have to return
        # ``False``.
        return visited[B] is not True

    def is_subdiagram(self, other):
        """
        Checks whether ``other`` is a subdiagram of ``self``.  Diagram
        `D'` is a subdiagram of `D` if all morphisms of `D'` are
        contained in the morphisms of `D`.  The morphisms contained
        both in `D'` and `D` should have the same properties for `D'`
        to be a subdiagram of `D`.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> d1 = Diagram(f)
        >>> d.is_subdiagram(d1)
        True
        >>> d1.is_subdiagram(d)
        False

        """
        # Check the inclusion of generators.  This suffices since all
        # other morphisms in diagrams are composites of the
        # generators.
        for m in other.generators:
            if (m not in self.generators_properties) or \
               not self.generators_properties[m].subset(
                other.generators_properties[m]):
                return False

        return True

    def __le__(self, other):
        """
        Checks if this :class:`Diagram` is a subdiagram of ``other``.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> d <= d
        True
        >>> d_ = Diagram(f)
        >>> d_ <= d
        True

        See Also
        ========

        is_subdiagram

        """
        return other.is_subdiagram(self)

    def __ge__(self, other):
        """
        Checks if ``other`` is a subdiagram of this :class:`Diagram`.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> d >= d
        True
        >>> d_ = Diagram(f)
        >>> d >= d_
        True

        See Also
        ========

        is_subdiagram

        """
        return self.is_subdiagram(other)

    def subdiagram_from_objects(self, objects):
        """
        If ``objects`` is a subset of the objects of ``self``, returns
        a :class:`Diagram` which has as morphisms all those morphisms
        of ``self`` which have a domains and codomains in ``objects``.
        Properties are preserved.

        Examples
        ========
        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram({f: "unique", g * f: "veryunique"})
        >>> d1 = d.subdiagram_from_objects([A, B])
        >>> d1 == Diagram({f: "unique"})
        True

        """
        objects = FiniteSet(objects)

        if not self.objects.subset(objects):
            raise ValueError("Supplied objects should all belong to the diagram.")

        # It suffices to filter the generators, because all other
        # morphisms are obtained as compositions of generators.

        new_generators = {}
        for morphism, props in self.generators_properties.items():
            if (morphism.domain in objects) and (morphism.codomain in objects):
                new_generators[morphism] = props

        return Diagram(new_generators)

    def __iter__(self):
        """
        Produces an iterator over all morphisms in this
        :class:`Diagram`.

        This returns the same generator as ``Diagram.morphisms``.

        Example
        =======

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import FiniteSet, pprint, default_sort_key
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> pprint(sorted(d, key=default_sort_key))
        [g*f:A-->C, id:A-->A, id:B-->B, id:C-->C, f:A-->B, g:B-->C]

        See Also
        ========

        morphisms

        """
        return iter(self.morphisms)

    def __contains__(self, morphism):
        """
        Checks whether ``morphism`` is contained in the diagram.

        Example
        =======

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> g * f in d
        True

        """
        if isinstance(morphism, CompositeMorphism) and \
               morphism not in self.generators_properties:
            # Composite morphisms always store the flattened-down
            # version of components, so it suffices to check the
            # inclusion of each of the component.
            return all(component in self.generators for component
                       in morphism.components)
        else:
            return morphism in self.generators

    def __len__(self):
        """
        If this :class:`Diagram` is finite, returns the number of
        morphisms it includes.

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> len(d)
        6

        """
        if self.is_finite:
            return len(list(self.morphisms))
        else:
            raise TypeError("This diagram is infinite.")

    def _morphism_properties(self, morphism):
        """
        Returns the properties of ``morphism`` in this diagram.  Does
        not check whether ``morphism`` indeed belongs to ``self``.
        """
        if morphism in self.generators_properties:
            return self.generators_properties[morphism]
        elif isinstance(morphism, CompositeMorphism):
            # This is a composite morphisms which was not explicitly
            # included in the generators.
            #
            # Remember that the properties of the composite morphism are
            # the intersection of the properties of its components.
            resulting_props =  self.generators_properties[morphism.components[0]]
            for component in morphism.components:
                resulting_props &= self.generators_properties[component]
            return resulting_props

    def __getitem__(self, morphism):
        r"""
        Retrieves the properties of the supplied ``morphism``, if it
        belongs to the diagram.  Throws :class:`ValueError` if
        ``morphism`` does not belong to this diagram.

        Example
        =======

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram({f: "unique", g: []})
        >>> d[f]
        {unique}

        """
        # People say that ``all`` is implemented directly in C, so
        # this check should not be too resource-intensive.  The code
        # is however clearer in this form.
        if morphism not in self:
            raise ValueError("%s does not belong to this diagram." %
                             str(morphism))

        return self._morphism_properties(morphism)

    def get(self, morphism, default=None):
        """
        Retrieves the properties of the supplied ``morphism``, if it
        belongs to the diagram.  Returns the value of ``default`` if
        ``morphism`` does not belong to this diagram.

        Example
        =======

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram({f: "unique", g: []})
        >>> d.get(f)
        {unique}
        >>> d.get(NamedMorphism(B, A, "f'")) is None
        True

        """
        # People say that ``all`` is implemented directly in C, so
        # this check should not be too resource-intensive.  The code
        # is however clearer in this form.
        if morphism not in self:
            return default

        return self._morphism_properties(morphism)

class Implication(Basic):
    r"""
    Represents an category theoretic implication in terms of diagrams.

    In category theoretic reasoning, a commutative diagram is often
    accompanied by a statement of the following kind: "if such
    morphisms with such properties exist, then such morphisms which
    such properties exist and the diagram is commutative".  To
    represent this, an :class:`Implication` includes two
    :class:`Diagram`'s: one to represent the premise of the
    implication, and the other to represent the conclusion.  The
    objects of the conclusion :class:`Diagram` should be a subset of
    the premise :class:`Diagram`.

    For example, the trivial fact that for every two morphisms
    `f:A\rightarrow B` and `g:B\rightarrow C` there exists the unique
    composite `g\circ f:A\rightarrow C` can be expressed by
    establishing the following implication:

    >>> from sympy.categories import Object, NamedMorphism, Diagram, Implication
    >>> from sympy import FiniteSet, pprint
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = NamedMorphism(A, B, "f")
    >>> g = NamedMorphism(B, C, "g")
    >>> premise = Diagram(f, g)
    >>> conclusion = Diagram({g * f: "unique"})
    >>> imp = Implication(premise, conclusion)
    >>> pprint(imp)
    {id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: EmptySet(), f:A-->B: Em
    ptySet(), g:B-->C: EmptySet()} ==> {g*f:A-->C: {unique}}

    Thus, the diagrams included in an :class:`Implication` should be
    interpreted in the following way: in the situation described by
    the premise diagram, there exists the morphisms with the
    properties described in the conclusion.  Notice however that a
    :class:`Diagram` includes all possible compositions between the
    morphisms supplied at creation:

    >>> pprint(imp.conclusion)
    {g*f:A-->C: {unique}, id:A-->A: EmptySet(), id:C-->C: EmptySet()}

    To get the morphisms that present some interest, use the method
    ``diff``.  This method returns the set of those morphisms which
    either appear only in the conclusion, or have non-empty properties
    in the conclusion, different from their properties in the premise

    >>> pprint(imp.diff())
    {g*f:A-->C}

    In some situations, it is useful to flatten the implication by
    squashing the premises and the conclusions in a single
    :class:`Diagram`.  This can be achieved via the method
    ``to_diagram``:

    >>> pprint(imp.to_diagram())
    {g*f:A-->C: {unique}, id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: Em
    ptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet()}
    >>> imp.to_diagram() == Diagram({g: [], f: [], g * f: "unique"})
    True

    See Also
    ========
    Diagram
    """
    def __new__(cls, premise, conclusion):
        r"""
        Constructs a new instance of :class:`Implication` from two
        :class:`Diagram`'s which represent the premise and the conclusion
        of the implication.

        The objects of the conclusion :class:`Diagram` must be a subset of
        the objects of the premise :class:`Diagram`.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram, Implication
        >>> from sympy import FiniteSet, pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> premise = Diagram(f, g)
        >>> conclusion = Diagram({g * f: "unique"})
        >>> imp = Implication(premise, conclusion)
        >>> pprint(imp)
        {id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: EmptySet(), f:A-->B: Em
        ptySet(), g:B-->C: EmptySet()} ==> {g*f:A-->C: {unique}}

        """
        if not premise.objects.subset(conclusion.objects):
            raise ValueError(
                "The conclusion must include the same objects as the premise.")

        return Basic.__new__(cls, premise, conclusion)

    @property
    def premise(self):
        """
        Returns the premise of this :class:`Implication`.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram, Implication
        >>> from sympy import FiniteSet, pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> premise = Diagram(f, g)
        >>> conclusion = Diagram({g * f: "unique"})
        >>> imp = Implication(premise, conclusion)
        >>> pprint(imp.premise)
        {id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: EmptySet(), f:A-->B: Em
        ptySet(), g:B-->C: EmptySet()}


        """
        return self.args[0]

    @property
    def conclusion(self):
        """
        Returns the conclusion of this :class:`Implication`.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram, Implication
        >>> from sympy import FiniteSet, pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> premise = Diagram(f, g)
        >>> conclusion = Diagram({g * f: "unique"})
        >>> imp = Implication(premise, conclusion)
        >>> pprint(imp.conclusion)
        {g*f:A-->C: {unique}, id:A-->A: EmptySet(), id:C-->C: EmptySet()}


        """
        return self.args[1]

    def to_diagram(self, conclusion_property=None):
        """
        Merges the premise and the conclusion of this implication into
        a single :class:`Diagram`.  If ``conclusion_property`` is supplied,
        every morphism from the conclusion will have ``conclusion_property``
        appended to its properties.

        If a morphism occurs both in the premise and in the conclusion, in
        the resulting :class:`Diagram` it will have the properties it has
        in the conclusion.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram, Implication
        >>> from sympy import FiniteSet, pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> premise = Diagram(f, g)
        >>> conclusion = Diagram({g * f: "unique"})
        >>> imp = Implication(premise, conclusion)
        >>> pprint(imp.to_diagram())
        {g*f:A-->C: {unique}, id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: Em
        ptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet()}

        """
        new_generators = dict(self.premise.generators_properties)
        for morphism, props in self.conclusion.generators_properties.items():
            new_props = props
            if conclusion_property and not isinstance(morphism, IdentityMorphism):
                new_props |= FiniteSet(conclusion_property)
            new_generators[morphism] = new_props

        return Diagram(new_generators)

    def diff(self):
        """
        Returns the :class:`FiniteSet` of morphisms which appear in the
        conclusion but do not appear in the premise, or which have non-empty
        properties in conclusion, different from the properties in the premise.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram, Implication
        >>> from sympy import FiniteSet, pprint
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> premise = Diagram(f, g)
        >>> conclusion = Diagram({g * f: "unique"})
        >>> imp = Implication(premise, conclusion)
        >>> pprint(imp.diff())
        {g*f:A-->C}

        """
        diff_morphisms = set([])
        for morphism, props in self.conclusion.generators_properties.items():
            if morphism not in self.premise:
                diff_morphisms.add(morphism)
            elif props and (props != self.premise[morphism]):
                diff_morphisms.add(morphism)
        return FiniteSet(diff_morphisms)
