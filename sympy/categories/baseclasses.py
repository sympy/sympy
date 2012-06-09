from sympy.core import Set, Basic, FiniteSet, EmptySet

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
        if not isinstance(obj, Object):
            return False

        if (not obj.name) or (not self.name):
            return False

        return self.name == obj.name

    def __ne__(self, obj):
        if not isinstance(obj, Object):
            return True

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
        if not isinstance(g, Morphism):
            return False

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

    Certain :class:`Diagram`s can be asserted to be commutative in a
    :class:`Category`, using the method :method:`assert_commutative`.

    Examples
    ========

    >>> from sympy.categories import Object, Morphism, Diagram, Category
    >>> from sympy import FiniteSet

    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")

    >>> f = Morphism(A, B, "f")
    >>> g = Morphism(B, C, "g")

    >>> d = Diagram()
    >>> d.add_premise(f)
    >>> d.add_premise(g)

    >>> K = Category("K")
    >>> K.assert_commutative(d)

    >>> K.commutative == FiniteSet(d)
    True

    See Also
    ========
    Diagram
    """
    objects = None
    name = ""

    def __init__(self, name, objects=EmptySet()):
        self.name = name
        self.objects = objects

        self.commutative = EmptySet()

    def hom(self, A, B):
        raise NotImplementedError(
            "hom-sets are not implemented in Category.")

    def all_morphisms(self):
        raise NotImplementedError(
            "Obtaining the class of morphisms is not implemented in Category.")

    def assert_commutative(self, *diagrams):
        """
        Asserts that certain diagrams are commutative in this
        category.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism, Diagram, Category
        >>> from sympy import FiniteSet

        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")

        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")

        >>> d = Diagram()
        >>> d.add_premise(f)
        >>> d.add_premise(g)

        >>> K = Category("K")
        >>> K.assert_commutative(d)

        >>> K.commutative == FiniteSet(d)
        True
        """
        self.commutative |= FiniteSet(diagrams)

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

    >>> from sympy.categories import Object, Morphism, Diagram
    >>> from sympy import FiniteSet
    >>> A = Object("A"); B = Object("B"); C = Object("C")
    >>> f = Morphism(A, B, "f"); g = Morphism(B, C, "g")
    >>> d = Diagram()

    >>> d.add_premise(f)
    >>> d.add_premise(g)

    >>> Morphism(A, A, identity=True) in d.premises.keys()
    True

    >>> g * f in d.premises.keys()
    True

    >>> d.add_conclusion(g * f, "unique")
    >>> d.conclusions[g * f] == FiniteSet("unique")
    True

    References
    ==========
    [Pare1970] B. Pareigis: Categories and functors.  Academic Press,
    1970.
    """
    def __init__(self):
        self.premises = {}
        self.conclusions = {}

    @staticmethod
    def _set_dict_union(dictionary, key, value):
        """
        If ``key`` is in ``dictionary``, set the new value of ``key``
        to be the union between the old value and ``value``.
        Otherwise, set the value of ``key`` to ``value.

        Returns ``True`` if the key already was in the dictionary and
        ``False`` otherwise.

        If ``key`` is ``None``, returns True and does nothing.
        """
        if not key:
            return True

        if key in dictionary:
            dictionary[key] = dictionary[key] | value
            return True
        else:
            dictionary[key] = value
            return False

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

        >>> from sympy.categories import Object, Morphism, Diagram
        >>> from sympy import FiniteSet
        >>> A = Object("A"); B = Object("B")
        >>> f = Morphism(A, B, "f")
        >>> d = Diagram()
        >>> d.add_premise(f)

        >>> f in d.premises.keys()
        True

        >>> Morphism(A, A, identity=True) in d.premises.keys()
        True

        See Also
        ========
        Morphism
        """
        props = FiniteSet(props)
        if self._set_dict_union(self.premises, morphism, props) == False:
            # We have just added a new morphism.

            if morphism.identity:
                return

            empty = EmptySet()

            # Add the identity morphisms for the domain and the
            # codomain.
            id_dom = Morphism(morphism.domain, morphism.domain, identity=True)
            id_cod = Morphism(morphism.codomain, morphism.codomain, identity=True)
            self._set_dict_union(self.premises, id_dom, empty)
            self._set_dict_union(self.premises, id_cod, empty)

            # Add all possible compositions
            for existing_morphism in self.premises.keys():
                left = morphism * existing_morphism
                right = existing_morphism * morphism

                self._set_dict_union(self.premises, left, props)
                self._set_dict_union(self.premises, right, props)

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

        >>> from sympy.categories import Object, Morphism, Diagram
        >>> from sympy import FiniteSet
        >>> A = Object("A"); B = Object("B"); C = Object("C")
        >>> f = Morphism(A, B, "f"); g = Morphism(B, C, "g")
        >>> d = Diagram()
        >>> d.add_premise(f)
        >>> d.add_premise(g)

        >>> d.add_conclusion(g * f, "unique")
        >>> d.conclusions[g * f] == FiniteSet("unique")
        True

        See Also
        ========
        Morphism, add_premise
        """
        objects = self.list_objects()
        if (morphism.domain not in objects) or \
           (morphism.codomain not in objects):
            return

        props = FiniteSet(props)
        if self._set_dict_union(self.conclusions, morphism, props) == False:
            # We have just added a new morphism.

            if morphism.identity:
                return

            empty = EmptySet()

            # We don't need to add more identities, because the domain
            # and the codomain of ``morphism`` are already in the
            # premises.

            # Add all possible compositions
            for existing_morphism in self.conclusions.keys():
                left = morphism * existing_morphism
                right = existing_morphism * morphism

                self._set_dict_union(self.conclusions, left, props)
                self._set_dict_union(self.conclusions, right, props)

    def list_objects(self):
        """
        Returns the set of objects which appear in this diagram.

        A :class:`Diagram` does not explicitly store a list of objects
        and constructs it from the domains and codomains of the
        included morphisms.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism, Diagram
        >>> from sympy import FiniteSet
        >>> A = Object("A"); B = Object("B"); C = Object("C")
        >>> f = Morphism(A, B, "f"); g = Morphism(B, C, "g")
        >>> d = Diagram()

        >>> d.add_premise(f)
        >>> d.add_premise(g)

        >>> d.list_objects() == FiniteSet(A, B, C)
        True

        See Also
        ========
        Object
        """
        objects = EmptySet()

        # Note that morphisms in the conclusions list cannot introduce
        # new objects.
        for morphism in self.premises.keys():
            objects |= FiniteSet(morphism.domain, morphism.codomain)

        return objects

    def hom(self, A, B):
        """
        Returns a 2-tuple of sets of morphisms between objects A and
        B: one set of morphisms listed as premises, and the other set
        of morphisms listed as conclusions.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism, Diagram
        >>> from sympy import FiniteSet
        >>> A = Object("A"); B = Object("B"); C = Object("C")
        >>> f = Morphism(A, B, "f"); g = Morphism(B, C, "g")
        >>> d = Diagram()
        >>> d.add_premise(f)
        >>> d.add_premise(g)
        >>> d.add_conclusion(g * f, "unique")

        >>> d.hom(A, C) == (FiniteSet(g * f), FiniteSet(g * f))
        True

        See Also
        ========
        Object, Morphism
        """
        premises = EmptySet()
        conclusions = EmptySet()

        for morphism in self.premises.keys():
            if (morphism.domain == A) and (morphism.codomain == B):
                premises |= FiniteSet(morphism)
        for morphism in self.conclusions.keys():
            if (morphism.domain == A) and (morphism.codomain == B):
                conclusions |= FiniteSet(morphism)

        return (premises, conclusions)

    def __eq__(self, other):
        if not isinstance(other, Diagram):
            return False

        return (self.premises == other.premises) and \
               (self.conclusions == other.conclusions)

    def __ne__(self, other):
        return not (self == other)
