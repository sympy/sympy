from sympy.core import (Set, Basic, FiniteSet, EmptySet, Dict, Symbol,
                        Dummy, Tuple)

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

def _make_symbol(name):
    """
    If ``name`` is not empty, creates a :class:`Symbol` with this
    name.  Otherwise creates a :class:`Dummy`.
    """
    if name:
        return Symbol(name)
    else:
        return Dummy("")

class Object(Basic):
    """
    The base class for any kind of object in an abstract category.

    While concrete categories may have some concrete SymPy classes as
    object types, in abstract categories only the name of an object is
    known.  Anonymous objects are not allowed.

    Two objects with the same name are the same object.

    Examples
    ========

    >>> from sympy.categories import Object
    >>> Object("A") == Object("A")
    True

    """
    def __new__(cls, name):
        if not name:
            raise ValueError("Anonymous Objects are not allowed.")

        return Basic.__new__(cls, Symbol(name))

    @property
    def name(self):
        """
        Returns the name of this object.

        Examples
        ========
        >>> from sympy.categories import Object
        >>> A = Object("A")
        >>> A.name
        'A'

        """
        return self._args[0].name

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
    composition.  Identity morphisms are instances of
    :class:`IdentityMorphism`.

    Examples
    ========

    >>> from sympy.categories import Object, Morphism
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = Morphism(A, B, "f")
    >>> g = Morphism(B, C, "g")
    >>> f == Morphism(A, B, "f")
    True
    >>> f == Morphism(A, B, "")
    False
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
        if identity and (domain != codomain):
            raise ValueError(
                "identity morphisms must have the same domain and codomain")

        # The last component of self.args represents the components of
        # this morphism.
        if identity:
            return IdentityMorphism(domain, name)
        else:
            return Basic.__new__(cls, domain, codomain,
                                 _make_symbol(name), Tuple())

    @property
    def domain(self):
        """
        Returns the domain of this morphism.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = Morphism(A, B)
        >>> f.domain
        Object("A")

        """
        return self.args[0]

    @property
    def codomain(self):
        """
        Returns the codomain of this morphism.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = Morphism(A, B)
        >>> f.codomain
        Object("B")

        """
        return self.args[1]

    @property
    def name(self):
        """
        Returns the name of this morphism.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = Morphism(A, B, "f")
        >>> f.name
        'f'

        """
        return self.args[2].name

    @property
    def is_identity(self):
        """
        Is ``True`` if this morphism is known to be an identity
        morphism.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = Morphism(A, B, "f")
        >>> f.is_identity
        False
        >>> id_A = Morphism(A, A, identity=True)
        >>> id_A.is_identity
        True

        """
        return False

    @property
    def components(self):
        """
        Returns the components of this morphisms.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> from sympy import Tuple
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")
        >>> (g * f).components == Tuple(f, g)
        True

        """
        components = self.args[3]
        if not components:
            return Tuple(self)
        else:
            return components

    def compose(self, g, new_name=""):
        """
        If ``self`` is a morphism from `B` to `C` and ``g`` is a
        morphism from `A` to `B`, returns the morphism from `A` to `C`
        which results from the composition of these morphisms.

        If either ``self`` or ``g`` are morphisms resulted from some
        previous composition, components in the resulting morphism
        will be the concatenation of ``g.components`` and
        ``self.components``, in this order.

        Instead of ``f.compose(g)`` it is possible to write ``f * g``.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")
        >>> g.compose(f)
        Morphism(Object("B"), Object("C"), "g") *
        Morphism(Object("A"), Object("B"), "f")
        >>> (g.compose(f, "h")).name
        'h'

        """
        if not isinstance(g, Morphism):
            raise TypeError("Morphisms can only be composed with morphisms.")

        if g.codomain != self.domain:
            raise ValueError("Uncomponsable morphisms.")

        if g.is_identity:
            return self

        # We don't really know whether the new morphism is an identity
        # (even if g.domain == self.codomain), so let's suppose it's
        # not an identity.
        return Basic.__new__(Morphism, g.domain, self.codomain,
                             _make_symbol(new_name), g.components +
                             self.components)

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
        try:
            return self.compose(g)
        except TypeError:
            return NotImplemented

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
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")
        >>> (g * f).flatten("h")
        Morphism(Object("A"), Object("C"), "h")

        See Also
        ========
        compose
        """

        return Morphism(self.domain, self.codomain, new_name)

    def __eq__(self, other):
        if other is self:
            return True

        if not isinstance(other, Morphism):
            return False

        if self.is_identity and other.is_identity:
            # All identities are equal.
            return self.domain == other.domain
        elif self.is_identity or other.is_identity:
            # One of the morphisms is an identity, but not both.
            return False

        if (len(self.components) == 1) and (len(other.components) == 1):
            # We are comparing two simple morphisms.
            if (not self.name) or (not other.name):
                return False

            return (self.name == other.name) and \
                   (self.domain == other.domain) and \
                   (self.codomain == other.codomain)
        else:
            # One of the morphisms is composed.  Compare the
            # components.
            for (self_component, other_component) in \
                    zip(self.components, other.components):
                if self_component != other_component:
                    return False
            return True

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.name, self.domain, self.codomain))

class IdentityMorphism(Morphism):
    """
    An identity morphism.

    An identity morphism is a morphism with equal domain and codomain,
    which acts as an identity with respect to composition.

    Examples
    ========

    >>> from sympy.categories import Object, Morphism, IdentityMorphism
    >>> A = Object("A")
    >>> B = Object("B")
    >>> f = Morphism(A, B, "f")
    >>> id_A = IdentityMorphism(A)
    >>> f * id_A == f
    True

    """
    def __new__(cls, domain, name=""):
        return Basic.__new__(cls, domain, _make_symbol(name))

    @property
    def domain(self):
        """
        Returns the domain of this identity morphism.

        Examples
        ========

        >>> from sympy.categories import Object, IdentityMorphism
        >>> A = Object("A")
        >>> id_A = IdentityMorphism(A)
        >>> id_A.domain
        Object("A")

        """
        return self.args[0]

    @property
    def codomain(self):
        """
        Returns the codomain of this identity morphism.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism
        >>> A = Object("A")
        >>> id_A = Morphism(A, A, identity=True)
        >>> id_A.codomain
        Object("A")

        """
        return self.args[0]

    @property
    def name(self):
        """
        Returns the name of this identity morphism.

        Examples
        ========

        >>> from sympy.categories import Object, IdentityMorphism
        >>> A = Object("A")
        >>> id_A = IdentityMorphism(A, "id_A")
        >>> id_A.name
        'id_A'

        """
        return self.args[1].name

    @property
    def is_identity(self):
        """
        Is ``True`` if this morphism is known to be an identity
        morphism.

        Examples
        ========

        >>> from sympy.categories import Object, IdentityMorphism
        >>> A = Object("A")
        >>> id_A = IdentityMorphism(A)
        >>> id_A.is_identity
        True

        """
        return True

    @property
    def components(self):
        r"""
        Returns the components of this morphisms.

        Since this is an identity morphisms, it always has itself as
        the only component.

        Examples
        ========

        >>> from sympy.categories import Object, IdentityMorphism
        >>> A = Object("A")
        >>> id_A = IdentityMorphism(A, 'id_A')
        >>> id_A.components
        (IdentityMorphism(Object("A"), "id_A"),)

        """
        return Tuple(self)

    def compose(self, g, new_name=""):
        """
        If ``g`` is a morphism with codomain equal to the domain of
        this identity morphisms, returns ``g``.

        The argument ``new_name`` is not used.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism, IdentityMorphism
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = Morphism(A, B, "f")
        >>> id_A = IdentityMorphism(A)
        >>> f.compose(id_A) == f
        True

        """
        if not isinstance(g, Morphism):
            raise TypeError("Morphisms can only be composed with morphisms.")

        if g.codomain != self.domain:
            raise ValueError("Uncomponsable morphisms.")

        return g

    def __hash__(self):
        return hash(self.domain)

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
    ``commutative`` in the constructor.

    Examples
    ========

    >>> from sympy.categories import Object, Morphism, Diagram, Category
    >>> from sympy import FiniteSet
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = Morphism(A, B, "f")
    >>> g = Morphism(B, C, "g")
    >>> d = Diagram([f, g])
    >>> K = Category("K", commutative=[d])
    >>> K.commutative == FiniteSet(d)
    True

    See Also
    ========
    Diagram
    """
    def __new__(cls, name, objects=EmptySet(), commutative=EmptySet()):
        if not name:
            raise ValueError("A Category cannot have an empty name.")

        new_category = Basic.__new__(cls, Symbol(name), objects,
                                     FiniteSet(commutative))
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
        {Object("B"), Object("A")}

        """
        return self.args[1]

    @property
    def commutative(self):
        """
        Returns the :class:`FiniteSet` of diagrams which are known to
        be commutative in this category.

        >>> from sympy.categories import Object, Morphism, Diagram, Category
        >>> from sympy import FiniteSet
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")
        >>> d = Diagram([f, g])
        >>> K = Category("K", commutative=[d])
        >>> K.commutative == FiniteSet(d)
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
    includes a collection of morphisms which are the premises and
    another collection of conclusions.  ``premises`` and
    ``conclusions`` associate morphisms belonging to the corresponding
    categories with the :class:`FiniteSet`'s of their properties.

    The set of properties of a composite morphism is the intersection
    of the sets of properties of its components.  The domain and
    codomain of a conclusion morphism should be among the domains and
    codomains of the morphisms listed as the premises of a diagram.

    No checks are carried out of whether the supplied object and
    morphisms do belong to one and the same category.

    Examples
    ========

    >>> from sympy.categories import Object, Morphism, Diagram
    >>> from sympy import FiniteSet
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = Morphism(A, B, "f")
    >>> g = Morphism(B, C, "g")
    >>> d = Diagram([f, g])
    >>> Morphism(A, A, identity=True) in d.premises.keys()
    True
    >>> g * f in d.premises.keys()
    True
    >>> d = Diagram([f, g], {g * f:"unique"})
    >>> d.conclusions[g * f] == FiniteSet("unique")
    True

    References
    ==========
    [Pare1970] B. Pareigis: Categories and functors.  Academic Press,
    1970.
    """
    @staticmethod
    def _set_dict_union(dictionary, key, value):
        """
        If ``key`` is in ``dictionary``, set the new value of ``key``
        to be the union between the old value and ``value``.
        Otherwise, set the value of ``key`` to ``value.

        Returns ``True`` if the key already was in the dictionary and
        ``False`` otherwise.
        """
        if key in dictionary:
            dictionary[key] = dictionary[key] | value
            return True
        else:
            dictionary[key] = value
            return False

    @staticmethod
    def _add_morphism(morphisms, morphism, props, add_identities=True):
        """
        Adds a morphism and its attributes to the supplied dictionary
        ``morphisms``.  If ``add_identities`` is True, also adds the
        identity morphisms for the domain and the codomain of
        ``morphism``.
        """
        if Diagram._set_dict_union(morphisms, morphism, props) == False:
            # We have just added a new morphism.

            if morphism.is_identity:
                return

            if add_identities:
                empty = EmptySet()

                id_dom = Morphism(morphism.domain, morphism.domain, identity=True)
                id_cod = Morphism(morphism.codomain, morphism.codomain, identity=True)
                Diagram._set_dict_union(morphisms, id_dom, empty)
                Diagram._set_dict_union(morphisms, id_cod, empty)

            for existing_morphism, existing_props in morphisms.items():
                new_props = existing_props & props
                if morphism.domain == existing_morphism.codomain:
                    left = morphism * existing_morphism
                    Diagram._set_dict_union(morphisms, left, new_props)
                if morphism.codomain == existing_morphism.domain:
                    right = existing_morphism * morphism
                    Diagram._set_dict_union(morphisms, right, new_props)

    def __new__(cls, *args):
        premises = {}
        conclusions = {}

        # Here we will keep track of the objects which appear in the
        # premises.
        objects = EmptySet()

        if len(args) >= 1:
            # We've got some premises in the arguments.
            premises_arg = args[0]

            if isinstance(premises_arg, list):
                # The user has supplied a list of morphisms, none of
                # which have any attributes.
                empty = EmptySet()

                for morphism in premises_arg:
                    objects |= FiniteSet(morphism.domain, morphism.codomain)
                    Diagram._add_morphism(premises, morphism, empty)
            elif isinstance(premises_arg, dict) or isinstance(premises_arg, Dict):
                # The user has supplied a dictionary of morphisms and
                # their properties.
                for morphism, props in premises_arg.items():
                    objects |= FiniteSet(morphism.domain, morphism.codomain)
                    Diagram._add_morphism(premises, morphism, FiniteSet(props))

        if len(args) >= 2:
            # We also have some conclusions.
            conclusions_arg = args[1]

            if isinstance(conclusions_arg, list):
                # The user has supplied a list of morphisms, none of
                # which have any attributes.
                empty = EmptySet()

                for morphism in conclusions_arg:
                    # Check that no new objects appear in conclusions.
                    if (morphism.domain in objects) and \
                       (morphism.codomain in objects):
                        # No need to add identities this time.
                        Diagram._add_morphism(conclusions, morphism, empty, False)
            elif isinstance(conclusions_arg, dict) or \
                     isinstance(conclusions_arg, Dict):
                # The user has supplied a dictionary of morphisms and
                # their properties.
                for morphism, props in conclusions_arg.items():
                    # Check that no new objects appear in conclusions.
                    if (morphism.domain in objects) and \
                       (morphism.codomain in objects):
                        # No need to add identities this time.
                        Diagram._add_morphism(conclusions, morphism,
                                           FiniteSet(props), False)

        return Basic.__new__(cls, Dict(premises), Dict(conclusions), objects)

    @property
    def premises(self):
        """
        Returns the premises of this diagram.

        Examples
        ========
        >>> from sympy.categories import Object, Morphism, Diagram
        >>> from sympy import EmptySet, Dict
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = Morphism(A, B, "f")
        >>> id_A = Morphism(A, A, identity=True)
        >>> id_B = Morphism(B, B, identity=True)
        >>> d = Diagram([f])
        >>> d.premises == Dict({f:EmptySet(), id_A:EmptySet(), id_B:EmptySet()})
        True

        """
        return self.args[0]

    @property
    def conclusions(self):
        """
        Returns the conclusions of this diagram.

        Examples
        ========
        >>> from sympy.categories import Object, Morphism, Diagram
        >>> from sympy import FiniteSet
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")
        >>> d = Diagram([f, g])
        >>> Morphism(A, A, identity=True) in d.premises.keys()
        True
        >>> g * f in d.premises.keys()
        True
        >>> d = Diagram([f, g], {g * f:"unique"})
        >>> d.conclusions[g * f] == FiniteSet("unique")
        True

        """
        return self.args[1]

    @property
    def objects(self):
        """
        Returns the :class:`FiniteSet` of objects that appear in this
        diagram.

        Examples
        ========
        >>> from sympy.categories import Object, Morphism, Diagram
        >>> from sympy import FiniteSet, Dict
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")
        >>> d = Diagram([f, g])
        >>> d.objects == FiniteSet(A, B, C)
        True

        """
        return self.args[2]

    def hom(self, A, B):
        """
        Returns a 2-tuple of sets of morphisms between objects A and
        B: one set of morphisms listed as premises, and the other set
        of morphisms listed as conclusions.

        Examples
        ========

        >>> from sympy.categories import Object, Morphism, Diagram
        >>> from sympy import FiniteSet
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = Morphism(A, B, "f")
        >>> g = Morphism(B, C, "g")
        >>> d = Diagram([f, g], {g * f: "unique"})
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
