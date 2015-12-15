class FreeGroup(Group):
    r"""
    Called  with  a  positive  integer rank, FreeGroup
    returns a free group on rank generators. If optional
    argument name is given then the  generators are printed
    as name1, name2 etc., that is, each name is the
    concatenation of the string name and an integer from 1
    to range. The default for name is the string "f".
    """
    is_group = True

    def __new__(cls, *args, **kwargs):
        if not(len(args) == 1 or isinstance(args[0], int)):
            raise ValueError("")

        obj = Basic.__new__(cls, *args)
        obj._generators = args
        obj._order = S.Infinity
        obj._center = []
        obj._degree = obj._generators
        obj._is_abelian = False

        return obj

    def __init__(self, *args):
        if isinstance(args[0], int):
            self.gens_assign(*args)

    def gens_assign(self, *args):
        f = var('f0:%s' %args)

    def __repr__(self):
        str_form = '<free group on the generators[ '
        for i in range(len(self)):
            str_form += 'f%s, ' %i
        str_form = str_form[:-2]
        str_form += ' ]>'
        return str_form

    def __len__(self):
        if isinstance(self.args[0], int):
            return self.args[0]

    def __getitem__(self):
        return self._generators[i]


class Group(Basic):
    pass


class GroupElem(Basic):
    """
    Group element definition

    """

    def __init__(self, elem, group):
        if not isinstance(group, Group):
            raise TypeError("group is not a Group")
        if not elem in group.Set:
            raise TypeError("elem is not an element of group")
        self.elem = elem
        self.group = group

    def __str__(self):
        return str(self.elem)

    def __eq__(self, other):
        """
        Two GroupElems are equal if they represent the same element,
        regardless of the Groups they belong to
        """

        if not isinstance(other, GroupElem):
            raise TypeError("other is not a GroupElem")
        return self.elem == other.elem

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.elem)

    def __mul__(self, other):
        """
        If other is a group element, returns self * other.
        If other = n is an int, and self is in an abelian group, returns self**n
        """
        if self.group.is_abelian() and isinstance(other, (int, long)):
            return self ** other

        if not isinstance(other, GroupElem):
            raise TypeError("other must be a GroupElem, or an int " \
                            "(if self's group is abelian)")
        try:
            return GroupElem(self.group.bin_op((self.elem, other.elem)), \
                             self.group)
        # This can return a TypeError in Funcion.__call__ if self and other
        # belong to different Groups. So we see if we can make sense of this
        # operation the other way around.
        except TypeError:
            return other.__rmul__(self)

    def __rmul__(self, other):
        """
        If other is a group element, returns other * self.
        If other = n is an int, and self is in an abelian group, returns self**n
        """
        if self.group.is_abelian() and isinstance(other, (int, long)):
            return self ** other

        if not isinstance(other, GroupElem):
            raise TypeError("other must be a GroupElem, or an int " \
                            "(if self's group is abelian)")

        return GroupElem(self.group.bin_op((other.elem, self.elem)), self.group)

    def __add__(self, other):
        """Returns self + other for Abelian groups"""
        if self.group.is_abelian():
            return self * other
        raise TypeError("not an element of an abelian group")

    def __pow__(self, n, modulo=None):
        """
        Returns self**n

        modulo is included as an argument to comply with the API, and ignored
        """
        if not isinstance(n, (int, long)):
            raise TypeError("n must be an int or a long")

        if n == 0:
            return self.group.e
        elif n < 0:
            return self.group.inverse(self) ** -n
        elif n % 2 == 1:
            return self * (self ** (n - 1))
        else:
            return (self * self) ** (n / 2)

    def __neg__(self):
        """Returns self ** -1 if self is in an abelian group"""
        if not self.group.is_abelian():
            raise TypeError("self must be in an abelian group")
        return self ** (-1)

    def __sub__(self, other):
        """Returns self * (other ** -1) if self is in an abelian group"""
        if not self.group.is_abelian():
            raise TypeError("self must be in an abelian group")
        return self * (other ** -1)

    def order(self):
        """Returns the order of self in the Group"""
        return len(self.group.generate([self]))
