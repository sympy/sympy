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
    """
    def __init__(self, name=""):
        self.name = name

class Morphism(Basic):
    """
    The base class for any kind of morphism in an abstract category.

    In abstract categories, a morphism is an arrow between two
    category objects.  The object where the arrow starts is called the
    domain, while the object where the arrow ends is called the
    codomain.
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

        If new_name is not an empty string, the new morphism will have
        the name new_name.  Otherwise the name of the new morphism
        will the juxtaposition of the names of morphisms in the
        components list, in reversed order, interspersed with '*'.

        Examples
        ========
        TODO: Add examples.
        """
        if g.codomain != self.domain:
            return None

        composite = Morphism(g.domain, self.codomain, new_name)
        composite.components = g.components + self.components

        if not new_name:
            for component in reversed(composite.components):
                composite.name += component.name + "*"
            composite.name = composite.name[:-1]

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

        If new_name is not an empty string, the new morphism will have
        the name new_name.  Otherwise the name of the new morphism
        will be the juxtaposition of the names of morphisms in
        self.components, in reversed order, interspersed with *.

        Examples
        ========
        TODO: Add examples.

        See Also
        ========
        compose
        """
        flattened = Morphism(self.domain, self.codomain, new_name)

        if not new_name:
            for component in reversed(self.components):
                flattened.name += component.name + "*"
            flattened.name = flattened.name[:-1]

        return flattened
