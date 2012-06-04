from sympy.core import Set, Basic

class Class(Set):
    """
    The base class for any kind of class in set-theoretic sense.

    In axiomatic set theories, everything is a class.  A class which
    can is a member of another class is a set.  A class which is not a
    member of another class is a proper class.  The class {1, 2} is a
    set; the class of all sets is a proper class.

    This class is essentially a synonym for sympy.core.Set.  The goal
    of this class is to assure easier migration to the eventual proper
    implementation of set theory.
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
    """
    def __init__(self, domain, codomain, name=""):
        self.domain = domain
        self.codomain = codomain
        self.name = name

        self.components = [self]

    def compose(self, g, new_name=""):
        """
        If self is a morphism from B to C and g is a morphism from A
        to B, returns the morphism from A to C which results from the
        composition of these morphisms.  Otherwise, returns None.

        If either self or g are morphisms resulted from some previous
        composition, components in the resulting morphism will be the
        concatenation of g.components and self.components, in this
        order.

        Examples
        ========
        TODO: Add examples.
        """
        if g.codomain != self.domain:
            return None

        composite = Morphism(g.domain, self.codomain, new_name)
        composite.components = g.components + self.components

        return composite

    def __mul__(self, g):
        """
        Returns the result of the composition of self with the
        argument, if this composition is defined.

        The semantics of multiplication is as follows: f * g =
        f.compose(g).

        See Also
        =======
        compose
        """
        return self.compose(g)

    def flatten(self, new_name=""):
        """
        If self resulted from composition of other morphisms, returns
        a new morphism without any information about the morphisms it
        resulted from.

        Note that comparing self with the new morphism need NOT return
        True.

        Examples
        ========
        TODO: Add examples.

        See Also
        ========
        compose
        """
        return Morphism(self.domain, self.codomain, new_name)

    def __eq__(self, g):
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
