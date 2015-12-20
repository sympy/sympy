## ------------- How the API is formed of `FreeGroupElm`? ------------------ ##
# FreeGroup( [wfilt, ]rank[, name] ) one example could be


# I don't know what the `wfilt` argument is all about ???????

# First API form is:  FreeGroup( rank )  here `rank` is any +ve integer

# FreeGroup( rank ) for any positive integer `rank`, returns a `FreeGroup`
# on `rank` generators.

# For example:
# gap> FreeGroup(3)
# <free group on the generators [ f1, f2, f3 ]>
# In Python i think similar functionality could be added
# >>> FreeGroup(3)
# <free group on the generators [ f1, f2, f3 ]>


# Second form of API is:   `FreeGroup(rank, "name")` here `rank` is any +ve integer and the string `name` is
# used for printing seems like only. Like `name1`, `name2` ... `name_rank` so `rank` number of generators.
# If the optional argument `name` is given then the generators
# are printed as `name1`, `name2` ... so on. So that's a concatenation of
# the string `name` provided by user and an integer from `1` to `range`.
# The default for `name` is the string `f`.

# For example:
# gap> FreeGroup(3, "apple");
# <free group on the generators [ apple1, apple2, apple3 ]>

# In this seems like "apple" string serves no other purpose other than representing it's name in the `__str__`
# form as we do it in Python.


# Third form of API is:   `FreeGroup("string1", "string2", "string3", ...)   here we just specify the
# name of string for each of generators.

# Fourth form of API is:   `FreeGroup(infinity, name, init)`

# So i think this has implication that the `FreeGroupElm` class can have the API of
# def __new__(cls, )


class FreeGroup(Basic):
    def __new__(cls, *args, **kwargs):

        obj = Basic.__new__(cls, *args, **kwargs)

        # (1) or (2) the First API form with `rank` and may be
        # provided with a `string`
        if isinstance(args[0], int) and args[0] >= 0:
            # (1) form of the API used here
            if len(args) == 1:
                obj._str = None
            # (2) form of the API used here
            elif len(args) == 2:
                obj._str = list([args[1]])
            # otherwise raise ValueError
            else:
                raise ValueError("Invalid arguments")
            obj._rank = args[0]
            obj._as_str = False

        # (3) API
        elif all([isinstance(i, str) for i in args]):
            obj._str = list(args)
            obj._rank = len(args)
            obj._as_str = True
        else:
            raise ValueError("Invalid arguments")

        # (4) API form with `rank` being `Infinity
        # TODO

        return obj

    def __getitem__(self, i):
        if i >= self.rank:
            raise IndexError("No such generator exists")
        elif self._str is None:
            return FreeGroupElm(self, i, 1, "f")
        elif len(self._str) == 1:
            return FreeGroupElm(self, i, 1, self._str[0])
        else:
            return FreeGroupElm(self, i, 1, self._str[i])

    def __str__(self):
        str_form = "<free group on the generators "
        gens = self.generators
        str_form += str(gens) + ">"
        return str_form

    __repr__ = __str__

    @property
    def as_str(self):
        return self._as_str

    def order(self):
        if self.rank == 0:
            return 1
        else:
            return S.Infinity

    @property
    def elements(self):
        if self.rank == 0:
            # A universal Identity element is returned
            return IdentityElm
        else:
            raise ValueError("Groups contains infinitely many "
                            ", hence can't be represented")

    @property
    def generators(self):
        if self.rank == 0:
            return list()
        else:
            return list([self[i] for i in range(self.rank)])

    @property
    def rank(self):
        """
        In group theory, the rank of a group G, denoted rank(G),
        can refer to the smallest cardinality of a generating set
        for G, that is

        \operatorname{rank}(G)=\min\{ |X|: X\subseteq G, \langle X\rangle =G\}.

        """
        return self._rank

    @property
    def is_abelian(self):
        if self.rank == 0:
            return True
        return False

# Now about the API of `FreeGroupElm` it is going to take in the
# The important things for elements of the same FreeGroup include
# things like the number.
# gap> f:=FreeGroup(4);;
# gap> f.1 < f.2;
# true
#
# gap> f.2 < f.3
# true

# quite strange are the properties

# def __new__(cls, group, number)

class FreeGroupElm(Basic):
    """
    Represents an element of FreeGroup. Other than
    the Identity element, all have Infinite order.

    """
    is_Identity = None

    def __new__(cls, free_group, index, pow_val, str_form="f"):
        # here `index` represents the max value of `generator` in
        # `FreeGroupElm` for example `f0**8*f2**5` has the
        # `index` value of `2`. since `f2` is there
        # while `pow_val` represents the `pow_val` of `index`
        # element in `FreeGroupElm`. in the above expression
        # `pow_val` is 5, since `f2**5` is there in `FreeGroupElm`
        obj = Basic.__new__(cls, free_group, index, pow_val, str_form)
        obj._group = free_group
        obj._index = index
        obj._pow_val = pow_val
        obj._str_form = str_form
        return obj

    @property
    def group(self):
        return self._group

    @property
    def str_form(self):
        return self._str_form

    @property
    def pow_val(self):
        return self._pow_val

    @property
    def group_index(self):
        return self._index

    def __str__(self):
        if self.pow_val != 1:
            return self.str_form + str(self.group_index) + "**" + str(self.pow_val)
        return self.str_form + str(self.group_index)

    __repr__ = __str__

    def __pow__(self, other):
        return FreeGroupElm(self.group, self.group_index, other, self.str_form)

    def __mul__(self, other):
        pass

    def order(self):
        if self.is_Identity:
            return 1
        else:
            return S.Infinity


class Group(Basic):
    pass


class GroupElm(Basic):
    """
    Group element definition

    """

    def __init__(self, elem, group=None):
        if group is None:
            raise ValueError("")
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


class GroupIdentity(Basic):
    """
    Represents the idenity for any Group.
    """
    is_Identity = True
    is_Abelian = True

    def __repr__(self):
        return "<identity ...>"

    def order(self):
        return 1

    def __mul__(self, other):
        return other

    __rmul__ = __mul__

    def __pow__(self, other):
        return self


def One(group):
    """
    Returns the Identity Element of Group `group`.
    """
    pass


class InfiniteGenerators(Basic):
    """
    Returns an InfiniteList of Generators.
    It prints the Infinite list of generators in the form
    of a limited with two elements [ f0, f1 ... ]
    """
    pass
