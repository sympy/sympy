"""Trait for implementing domain elements. """

class DomainElement(object):
    """
    Represents an element of a domain.

    Mix in this trait into a class which instances should be recognized as
    elements of a domain. Method ``parent()`` gives that domain.

    """

    def parent(self):
        raise NotImplementedError("abstract method")
