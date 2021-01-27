"""Trait for implementing domain elements. """


from sympy.utilities import public

@public
class DomainElement:
    """
    Represents an element of a domain.

    Mix in this trait into a class whose instances should be recognized as
    elements of a domain. Method ``parent()`` gives that domain.
    """

    def parent(self):
        """Get the domain associated with ``self`` """
        raise NotImplementedError("abstract method")
