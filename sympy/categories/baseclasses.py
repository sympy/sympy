from sympy.core import Set, Basic, EmptySet

class Class(Set):
    """
    The base class for any kind of class in set-theoretic sense.

    In axiomatic set theories, everything is a class.  A class which
    can is a member of another class is a set.  A class which is not a
    member of another class is a proper class.  The class {1, 2} is a
    set; the class of all sets is a proper class.

    This class is essentially a synonym for :class:`sympy.core.Set`.
    The goal of this class is to assure easier migration to the
    eventual proper implementation of set theory.
    """
    is_proper = False

class Object(Basic):
    """
    The base class for any kind of object in an abstract category.

    While concrete categories may have some concrete SymPy classes as
    object types, in abstract categories only the name of an object is
    known.

    Two objects with the same name are the same object.  An unnamed
    object is not equal to any other object.

    Examples
    ========

    >>> from sympy.categories import Object
    >>> Object("A") == Object("A")
    True

    >>> Object("A") == Object("")
    False

    """
    def __init__(self, name=""):
        self.name = name

    def __eq__(self, obj):
        if (not obj.name) or (not self.name):
            return False

        return self.name == obj.name

    def __ne__(self, obj):
        if (not obj.name) or (not self.name):
            return True

        return self.name != obj.name

    def __hash__(self):
        return hash(self.name)

class Morphism(Basic):
    """
    The base class for any kind of morphism in an abstract category.

    In abstract categories, a morphism is an arrow between two
    category objects.  The object where the arrow starts is called the
    domain, while the object where the arrow ends is called the
    codomain.

    Two simple (not composed) morphisms with the same name, domain,
    and codomain are the same morphisms.  A simple unnamed morphism is
    not equal to any other morphism.

    Two composed morphisms are equal if they have the same components,
    in the same order (which guarantees the equality of domains and
    codomains).  The names of such composed morphisms are not taken in
    consideration at comparison.

    Morphisms with the same domain and codomain can be defined to be
    identity morphisms.  Identity morphisms with the same (co)domains
    are equal.  Identity morphisms are identities with respect to
    composition.

    Examples
    ========

    >>> from sympy.categories import Object, Morphism
    >>> A = Object("A"); B = Object("B"); C = Object("C")
    >>> f = Morphism(A, B, "f")
    >>> g = Morphism(B, C, "g")

    >>> f == Morphism(A, B, "f")
    True
    >>> f == Morphism(A, B, "")
    False

    >>> f * g is None
    True
    >>> g * f
    Morphism(Object("B"), Object("C"), "g") *
    Morphism(Object("A"), Object("B"), "f")

    >>> id_A = Morphism(A, A, identity=True)
    >>> id_A == Morphism(A, A, identity=True)
    True

    >>> f * id_A == f
    True

    """
    def __new__(cls, domain, codomain, name="", identity=False):
        new_morphism = Basic.__new__(cls, domain, codomain, name, identity)

        new_morphism.domain = domain
        new_morphism.codomain = codomain
        new_morphism.name = name

        new_morphism.components = [new_morphism]

        new_morphism.identity = identity

        if identity and (domain != codomain):
            raise ValueError(
                "identity morphisms must have the same domain and codomain")

        return new_morphism

    def compose(self, g, new_name=""):
        """
        If ``self`` is a morphism from B to C and ``g`` is a morphism
        from A to B, returns the morphism from A to C which results
        from the composition of these morphisms.  Otherwise, returns
        ``None``.

        If either ``self`` or ``g`` are morphisms resulted from some
        previous composition, components in the resulting morphism
        will be the concatenation of ``g.components`` and
        ``self.components``, in this order.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> A = Object("A"); B = Object("B"); C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")

        >>> f.compose(g) is None
        True

        >>> g.compose(f)
        Morphism(Object("B"), Object("C"), "g") *
        Morphism(Object("A"), Object("B"), "f")

        >>> (g.compose(f, "h")).name
        'h'

        """
        if g.codomain != self.domain:
            return None

        if self.identity:
            return g
        if g.identity:
            return self

        composite = Morphism(g.domain, self.codomain, new_name)
        composite.components = g.components + self.components

        return composite

    def __mul__(self, g):
        """
        Returns the result of the composition of ``self`` with the
        argument, if this composition is defined.

        The semantics of multiplication is as follows: ``f * g =
        f.compose(g)``.

        See Also
        =======
        compose
        """
        return self.compose(g)

    def flatten(self, new_name=""):
        """
        If ``self`` resulted from composition of other morphisms,
        returns a new morphism without any information about the
        morphisms it resulted from.

        Note that comparing ``self`` with the new morphism need NOT
        return ``True``.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> A = Object("A"); B = Object("B"); C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")

        >>> (g * f).flatten("h")
        Morphism(Object("A"), Object("C"), "h")

        See Also
        ========
        compose
        """
        return Morphism(self.domain, self.codomain, new_name)

    def __eq__(self, g):
        if self.identity and g.identity:
            # All identities are equal.
            return self.domain == g.domain
        elif self.identity or g.identity:
            # One of the morphisms is an identity, but not both.
            return False

        if (len(self.components) == 1) and (len(g.components) == 1):
            # We are comparing two simple morphisms.
            if (not self.name) or (not g.name):
                return False

            return (self.name == g.name) and \
                   (self.domain == g.domain) and \
                   (self.codomain == g.codomain)
        else:
            # One of the morphisms is composed.  Compare the
            # components.
            for (self_component, g_component) in zip(self.components, g.components):
                if self_component != g_component:
                    return False
            return True

    def __ne__(self, g):
        return not (self == g)

    def __hash__(self):
        return hash((self.name, self.domain, self.codomain))

class Category(Basic):
    """
    An (abstract) category.

    A category [JoyOfCats] is a quadruple K = (O, hom, id, *)
    consisting of

    * a (set-theoretical) class O, whose members are called K-objects,

    * for each pair (A, B) of K-objects, a set hom(A, B) whose members
      are called K-morphisms from A to B,

    * for a each K-object A, a morphism id:A->A, called the K-identity
      on A,

    * a composition law associating with every K-morphisms f:A->B and
      g:B->C a K-morphism g*f:A->C, called the composite of f and g.

    Composition is associative, K-identities are identities with
    respect to composition, and the sets hom(A, B) are pairwise
    disjoint.

    This class nothing about its objects and morphisms.  Concrete
    cases of (abstract) categories should be implemented as classes
    derived from this one.
    """
    objects = None
    name = ""

    def __init__(self, name, objects=EmptySet()):
        self.name = name
        self.objects = objects

    def hom(self, A, B):
        raise NotImplementedError(
            "hom-sets are not implemented in Category.")

    def all_morphisms(self):
        raise NotImplementedError(
            "Obtaining the class of morphisms is not implemented in Category.")

class Diagram(Basic):
    r"""
    Represents a diagram in a certain category.

    Informally, a diagram is a collection of objects of a category and
    certain morphisms between them.  A diagram is still a monoid with
    respect to morphism composition; i.e., identity morphisms, as well
    as all composites of morphisms included in the diagram belong to
    the diagram.  For a more formal approach to this notion see
    [Pare1970].

    A commutative diagram is often accompanied by a statement of the
    following kind: "if such morphisms with such properties exist,
    then such morphisms which such properties exist and the diagram is
    commutative".  To represent this, an instance of :class:`Diagram`
    includes a list of morphisms which are the premises and another
    list of conclusions.  ``premises`` and ``conclusions`` associate
    morphisms belonging to the corresponding categories with the
    :class:`FiniteSet`'s of their properties.

    While ``premises`` and ``conclusions`` are public attributes of
    this class, it is highly discouraged to directly modify them.
    This may result in inconsistent behaviour.  Use ``add_premise``
    and ``add_conclusion`` instead.

    No checks are carried out of whether the supplied object and
    morphisms do belong to one and the same category.

    Examples
    ========
    TODO: Add examples.

    References
    ==========
    [Pare1970] B. Pareigis: Categories and functors.  Academic Press,
    1970.
    """
    def __init__(self):
        self.premises = {}
        self.conclusions = {}

    def add_premise(self, morphism, *props):
        """
        Adds a morphism and its properties to the premises of this
        diagram.

        ``morphism`` should be compatible with :class:`Morphism`.  A
        property is a string.

        When another morphisms is added to the diagram, all necessary
        morphisms are added to keep the diagram a monoid with respect
        to composition: the identity morphisms of the domain and the
        codomain, as well as all possible compositions with the
        morphisms already included in the diagram.  The set of
        properties of a composite morphism is the intersection of the
        sets of properties of its components.

        Examples
        ========
        TODO: Add examples.

        See Also
        ========
        Morphism
        """
        pass

    def add_conclusion(self, morphism, *props):
        """
        Adds a morphism and its properties to the conclusions of this
        diagram.

        ``morphism`` should be compatible with :class:`Morphism`.  A
        property is a string.  The domain and codomain of the morphism
        should be among the domains and codomains of the morphisms
        listed as the premises of this diagram.

        When another morphisms is added to the diagram, all necessary
        morphisms are added to keep the diagram a monoid with respect
        to composition: the identity morphisms of the domain and the
        codomain, as well as all possible compositions with the
        morphisms already included in the diagram.  The set of
        properties of a composite morphism is the intersection of the
        sets of properties of its components.

        Examples
        ========
        TODO: Add examples.

        See Also
        ========
        Morphism, add_premise
        """
        pass

    def list_objects(self):
        """
        Returns the set of objects which appear in this diagram.

        A :class:`Diagram` does not explicitly store a list of objects
        and constructs it from the domains and codomains of the
        included morphisms.

        Examples
        ========
        TODO: Add examples.

        See Also
        ========
        Object
        """
        pass

    def hom(self, A, B):
        """
        Returns a 2-tuple of sets of morphisms between objects A and
        B: one set of morphisms listed as premises, and the other set
        of morphisms listed as conclusions.

        Examples
        ========
        TODO: Add examples.

        See Also
        ========
        Object, Morphism
        """
        pass

    def __eq__(self, other):
        pass

    def __eq__(self, other):
        pass

    def hash(self):
        pass
