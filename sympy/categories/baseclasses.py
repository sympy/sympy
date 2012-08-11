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
    TODO: Rewrite all the docstrings in this class.

    Represents a diagram in a certain category.

    Informally, a diagram is a collection of objects of a category and
    certain morphisms between them.  Identity morphisms, as well as
    all composites of morphisms included in the diagram belong to the
    diagram.  For a more formal approach to this notion see
    [Pare1970].

    A :class:`Diagram` stores a mapping between morphisms and the
    :class:`FiniteSet`'s of their properties.  All possible morphism
    compositions, as well as the components of composite morphisms are
    also added to the diagram.  No properties are assigned to such
    morphisms by default.  The set of properties of a composite
    morphism is the intersection of the sets of properties of its
    components.

    No checks are carried out of whether the supplied objects and
    morphisms do belong to one and the same category.

    Examples
    ========

    >>> from sympy.categories import Object, NamedMorphism, Diagram
    >>> from sympy import FiniteSet, pprint, default_sort_key
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = NamedMorphism(A, B, "f")
    >>> g = NamedMorphism(B, C, "g")
    >>> d = Diagram(f, g)
    >>> morphisms = sorted(d, key=default_sort_key)
    >>> pprint(morphisms, use_unicode=False)
    [g*f:A-->C, id:A-->A, id:B-->B, id:C-->C, f:A-->B, g:B-->C]
    >>> pprint(d.morphisms, use_unicode=False)
    {g*f:A-->C: EmptySet(), id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C:
    EmptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet()}

    References
    ==========
    [Pare1970] B. Pareigis: Categories and functors.  Academic Press,
    1970.
    """
    def  __new__(cls, *args):
        """
        TODO: Document.

        TODO: Mention that identities are dropped.
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
                    generators[morphism] = FiniteSet(props)
            elif iterable(first_arg):
                # The user has supplied a list of morphisms, none of
                # which have any properties.
                empty = EmptySet()
                for morphism in first_arg:
                    generators[morphism] = empty
            else:
                # Attempt to interpret ``args`` as a list of
                # morphisms.
                return Diagram(args)

        for morphism in generators:
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
        TODO: Rewrite

        Returns the :class:`Dict` mapping the morphisms included in
        this :class:`Diagram` to their properties.

        Examples
        ========
        >>> from sympy.categories import Object, NamedMorphism
        >>> from sympy.categories import IdentityMorphism, Diagram
        >>> from sympy import pretty
        >>> A = Object("A")
        >>> B = Object("B")
        >>> f = NamedMorphism(A, B, "f")
        >>> id_A = IdentityMorphism(A)
        >>> id_B = IdentityMorphism(B)
        >>> d = Diagram(f)
        >>> print pretty(d.morphisms, use_unicode=False)
        {id:A-->A: EmptySet(), id:B-->B: EmptySet(), f:A-->B: EmptySet()}

        """
        return self.args[0].keys()

    @property
    def generators_properties(self):
        """
        TODO: Describe.
        """
        return self.args[0]

    @staticmethod
    def _check_cycles_in_generators(obj, adj_lists, pre_visited, post_visited):
        """
        Traverses the graph supplied in ``adj_lists`` starting from
        ``object``, marking the visited objects in preorder in
        ``pre_visited`` and in postorder in ``post_visited``.  Reports
        ``False`` if it never finds a back edge and ``True``
        otherwise.

        References
        ==========

        [Sedg2002] Robert Sedgewick, Algorithms in C++ Graph
        Algorithms, Part 5.  Third Edition.  Addison-Wesley, 2002.
        Section Digraphs and DAGs.

        [Wikipedia] http://en.wikipedia.org/wiki/Depth-first_search#Output_of_a_depth-first_search

        """
        pre_visited[obj] = True

        for adj_obj in adj_lists[obj]:
            if not pre_visited[adj_obj]:
                # We haven't yet been there.
                if Diagram._check_cycles_in_generators(
                    adj_obj, adj_lists, pre_visited, post_visited):
                    return True
            elif not post_visited[adj_obj]:
                # ``(obj, adj_obj)`` is a back edge ([Wikipedia],
                # [Sedg2002], Listing 19.2), and we have detected a
                # cycle.
                #
                # Essentially, this means that ``adj_obj`` actually
                # precedes ``obj`` in the DFS tree.  This means that
                # one can go from ``adj_obj`` to ``obj`` (because
                # ``adj_obj`` precedes ``obj``), and then back to
                # ``adj_obj``, thus closing a cycle.
                #
                # The reason this check works is the following.
                # ``post_visited`` is filled in during _postorder_
                # traversal, that is, when we are getting up the DFS
                # tree.  ``post_visited[adj_obj] == False`` means that
                # we haven't yet reached ``adj_obj`` on our way back.
                # Therefore, when we were going down the DFS tree, we
                # got to ``obj`` through ``adj_obj``.  (The best way
                # to understand this is to actually draw some simple
                # graphs.)
                return True

        post_visited[obj] = True
        # No cycles found yet.
        return False

    @property
    @cacheit
    def is_finite(self):
        """
        Checks whether this diagram contains a finite number of
        morphisms.

        TODO: Expand
        """
        # This method essentially checks whether there are cycles in
        # the digraph defined by the generators.  Note that it is
        # sufficient to only check the generators, since the other
        # (derived) morphisms essentially form a transitive closure of
        # the digraph and thus cannot contribute new cycles.
        #
        # For references, see ``_check_cycles_in_generators``.  Note
        # however that the code in that method and in ``is_infinite``
        # is an adaptation of what can be found in [Sedg2002].

        # In this dictionary, vertices will be marked as visited
        # during preorder traversal.
        pre_visited = {}

        # In this dictionary, vertices will be marked as visited
        # during preorder traversal.
        post_visited = {}

        # We will also count the number of morphisms which go into
        # each object.
        indegrees = {}

        # This dictionary will store the adjacency lists of the graph
        # of generators.
        adj_lists = {}

        for obj in self.objects:
            adj_lists[obj] = []
            pre_visited[obj] = False
            post_visited[obj] = False
            indegrees[obj] = 0

        # Build the adjacency lists of the graph of generators.
        for m in self.generators:
            if isinstance(m, IdentityMorphism):
                # The generators may contain some identity morphisms,
                # which were explicitly supplied with properties.
                # Yet, since they are identities, they do not add new
                # morphisms to the diagram, so they should be skipped
                # here.
                continue
            adj_lists[m.domain].append(m.codomain)
            indegrees[m.codomain] += 1

        # Sort the objects by the number of morphisms which go into
        # them.  This will allow us to start the traversal with
        # sources, if there are sources.
        objects_by_indegree = sorted(
            self.objects,key=lambda obj: indegrees[obj])

        for obj in objects_by_indegree:
            if not pre_visited[obj]:
                # Look for cycles in the DFS tree which starts at this
                # object.
                if self._check_cycles_in_generators(
                    obj, adj_lists, pre_visited, post_visited):
                    # There are cycles, this diagram is infinite.
                    return False

        # No cycles, this is a finite diagram.
        return True

    def _composites_by_length(self, length, identities=True):
        """
        Given the desired ``length`` (component count), generates all
        possible compositions of morphisms from ``generators`` of this
        length.  For ``length == 1``, this method will produce exactly
        the list of the generators of this :class:`Diagram`.  For
        ``length == 0`` will produce nothing.

        If ``identities == True``, and ``length == 1``, this function
        will produce _all_ generators, including identities.
        Otherwise it will only produce non-identity morphisms.
        """
        if not length:
            return
        if length == 1:
            for m in self.generators:
                if isinstance(m, IdentityMorphism) and not identities:
                    continue
                yield m
        else:
            first_half_len = length/2
            second_half_len = length - first_half_len

            first_half = self._composites_by_length(first_half_len,
                                                    identities=False)
            second_half = self._composites_by_length(second_half_len,
                                                    identities=False)

            for (first, second) in product(first_half, second_half):
                # We only need to check either ``first * second`` or
                # ``second * first``, and it doesn't matter which of
                # them, because we go through all possible
                # combinations.
                #
                # Note that some composites may have already been
                # supplied as generators; no need to yield them twice.
                if (second.domain == first.codomain):
                    composite = second * first
                    if composite not in self.generators_properties:
                        yield composite

    @property
    def morphisms(self):
        """
        Generates all the morphisms of this diagram.

        TODO: Expand
        """
        length = 1

        # If this diagram is finite, we want to only go on until there
        # are composites.  Otherwise, we just go on and on.
        #
        # The point of having a lambda is to avoid checking
        # ``is_finite`` on every iteration.
        if self.is_finite:
            halt_condition = lambda ncomposites: ncomposites == 0
        else:
            halt_condition = lambda x: False

        ncomposites = 1
        while not halt_condition(ncomposites):
            ncomposites = 0
            for composite in self._composites_by_length(length):
                ncomposites += 1
                yield composite
            length += 1

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
        TODO: Rewrite

        Returns the :class:`FiniteSet` of morphisms between the
        objects ``A`` and ``B``.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import pretty
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> print pretty(d.hom(A, C), use_unicode=False)
        {g*f:A-->C}

        See Also
        ========
        Object, Morphism
        """
        # For infinite diagrams, we rely on the fact that
        # ``morphisms`` assures a sufficiently BFS-like traversal.
        return (m for m in self.morphisms if (m.domain == A) and
                (m.codomain == B))

    def is_subdiagram(self, other):
        """
        TODO: Check

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

    def subdiagram_from_objects(self, objects):
        """
        TODO: Check

        If ``objects`` is a subset of the objects of ``self``, returns
        a diagram which has as morphisms all those morphisms of
        ``self`` which have a domains and codomains in ``objects``.
        Properties are preserved.

        Examples
        ========
        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import FiniteSet
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram({f: "unique", g * f: "veryunique"})
        >>> d1 = d.subdiagram_from_objects(FiniteSet(A, B))
        >>> d1 == Diagram({f: "unique"})
        True
        """
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
        TODO: Check

        Produces an iterator over the underlying dictionary of
        morphisms.

        Example
        =======

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import FiniteSet, pretty, default_sort_key
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> d = Diagram(f, g)
        >>> sorted(d, key=default_sort_key)
        [CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"),
        NamedMorphism(Object("B"), Object("C"), "g"))),
        IdentityMorphism(Object("A")), IdentityMorphism(Object("B")),
        IdentityMorphism(Object("C")), NamedMorphism(Object("A"),
        Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g")]

        """
        return iter(self.morphisms)

    def __contains__(self, morphism):
        """
        Checks whether ``morphism`` is contained in the diagram.

        Example
        =======

        >>> from sympy.categories import Object, NamedMorphism, Diagram
        >>> from sympy import FiniteSet, pretty
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
        TODO: Document.
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
        >>> from sympy import FiniteSet, pretty
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
        >>> from sympy import FiniteSet, pretty
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
    {g*f:A-->C: EmptySet(), id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C:
    EmptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet()} ==> {g*f:A-->C: {unique}
    }

    Thus, the diagrams included in an :class:`Implication` should be
    interpreted in the following way: in the situation described by
    the premise diagram, there exists the morphisms with the
    properties described in the conclusion.  Notice however that a
    :class:`Diagram` includes all possible compositions between the
    morphisms supplied at creation:

    >>> pprint(imp.conclusion)
    {g*f:A-->C: {unique}, id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: Em
    ptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet()}

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
    >>> imp.to_diagram() == Diagram({g * f: "unique"})
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
        {g*f:A-->C: EmptySet(), id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C:
        EmptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet()} ==> {g*f:A-->C: {unique}
        }

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
        {g*f:A-->C: EmptySet(), id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C:
        EmptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet()}

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
        {g*f:A-->C: {unique}, id:A-->A: EmptySet(), id:B-->B: EmptySet(), id:C-->C: Em
        ptySet(), f:A-->B: EmptySet(), g:B-->C: EmptySet()}

        """
        return self.args[1]

    def to_diagram(self, conclusion_property=None):
        """
        TODO: Check

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
