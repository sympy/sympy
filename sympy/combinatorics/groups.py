from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.core.compatibility import as_int
from sympy.utilities.iterables import flatten

## ------------ How the API is formed of `FreeGroupElm`? --------------- ##
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


# Second form of API is:   `FreeGroup(rank, "name")` here `rank` is any +ve
# integer and the string `name` is used for printing seems like only. Like
# `name1`, `name2` ... `name_rank` so `rank` number of generators. If the optional
# argument `name` is given then the generators are printed as `name1`, `name2` ...
# so on. So that's a concatenation of the string `name` provided by user and an
# integer from `1` to `range`. The default for `name` is the string `f`.

# For example:
# gap> FreeGroup(3, "apple");
# <free group on the generators [ apple1, apple2, apple3 ]>

# In this seems like "apple" string serves no other purpose other than representing
# it's name in the `__str__`form as we do it in Python.


# Third form of API is:   `FreeGroup("string1", "string2", "string3", ...)
# here we just specify the name of string for each of generators.

# Fourth form of API is:   `FreeGroup(infinity, name, init)`

# So i think this has implication that the `FreeGroupElm` class can have the API of
# def __new__(cls, )


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


# Second form of API is:   `FreeGroup(rank, "name")` here `rank` is any +ve
# integer and the string `name` is used for printing seems like only. Like
# `name1`, `name2` ... `name_rank` so `rank` number of generators. If the optional
# argument `name` is given then the generators are printed as `name1`, `name2` ...
# so on. So that's a concatenation of the string `name` provided by user and an integer
# from `1` to `range`. The default for `name` is the string `f`.

# For example:
# gap> FreeGroup(3, "apple");
# <free group on the generators [ apple1, apple2, apple3 ]>

# In this seems like "apple" string serves no other purpose other than representing
# it's name in the `__str__` form as we do it in Python.


# Third form of API is:   `FreeGroup("string1", "string2", "string3", ...)   here we just
# specify the name of string for each of generators.

# Fourth form of API is:   `FreeGroup(infinity, name, init)`

# So i think this has implication that the `FreeGroupElm` class can have the API of
# def __new__(cls, )


class FreeGroup(Basic):
    """

    """
    is_associative = True
    is_group = True
    is_FreeGroup = True
    is_PermutationGroup = False

    def __new__(cls, *args, **kwargs):
        """
        Called with a positive integer rank, FreeGroup returns
        a free group on rank generators. If the optional argument
        `name` is given then the generators are printed as name1,
        name2 etc., that is, each name is the concatenation of the
        string name and an integer from 1 to `range`. The default
        for name is the string "f". Called in the second form,
        FreeGroup returns a free group on as many generators as
        arguments, printed as name1, name2 etc.

        Called in the second form, FreeGroup returns a free group on
        as many generators as arguments, printed as name1, name2 etc.
        """

        obj = Basic.__new__(cls, *args, **kwargs)

        if isinstance(args[0], int) and args[0] >= 0:

            # (1) form of the API used here
            if len(args) == 1:
                obj._str = "f"

            # (2) form of the API used here
            elif len(args) == 2:
                obj._str = args[1]
            else:
                raise ValueError("Invalid arguments")
            obj._rank = args[0]

        # (3) form of API, where all the generators of `FreeGroup`
        # given as strings.
        elif all([isinstance(i, str) for i in args]):
            pass

        elif args[0] is S.Infinity:
            pass

        return obj

    @property
    def generators(self):
        """Returns the generators of the FreeGroup

        Examples
        ========

        >>> f = FreeGroup( 3, "swapnil" )
        >>> f.generators
        [swapnil0, swapnil1, swapnil2]

        """
        gens = []
        for i in range(self.rank):
            gens.append(FreeGroupElm(self, [(i, 1)], self.as_str))
        return gens

    def __getitem__(self, i):
        return self.generators[i]

    def __contains__(self, i):
        """
        Return True if `i` is contained in FreeGroup.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4, "swap" )
        >>> g = FreeGroup( 4, "tul" )

        >>> f[0]**2*f[3] in f
        True
        >>> f[0]**2*f[3] in g
        False

        """
        if not isinstance(i, FreeGroupElm):
            raise TypeError("FreeGroup contains elements of type `FreeGroupElm`")
        return self.contains(i)

    def __len__(self):
        return self.rank

    def __str__(self):
        if self.rank > 30:
            str_form = "<free group with %s generators>" % self.rank
        else:
            str_form = "<free group on the generators "
            gens = self.generators
            str_form += str(gens) + ">"
        return str_form

    __repr__ = __str__

    def __eq__(self, other):
        """No `FreeGroup` if equal to any other `FreeGroup`.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4, "swapnil" )
        >>> g = FreeGroup( 4, "swapnil" )
        >>> f == g
        False
        >>> f == f
        True

        """
        if self is other:
            return True
        else:
            return False

    @property
    def as_str(self):
        return self._str

    def order(self):
        if self.rank == 0:
            return 1
        else:
            return S.Infinity

    @property
    def elements(self):
        if self.rank == 0:
            # A set containing universal Identity element is returned
            return set([IdentityElm])
        else:
            raise ValueError("Group contains infinitely many elements"
                            ", hence can't be represented")

    @property
    def rank(self):
        r"""
        In group theory, the rank of a group G, denoted rank(G),
        can refer to the smallest cardinality of a generating set
        for G, that is

        \operatorname{rank}(G)=\min\{ |X|: X\subseteq G, \langle X\rangle =G\}.

        """
        return self._rank

    @property
    def is_abelian(self):
        """Test if the group is Abelian.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> f.is_abelian
        False

        >>> g = FreeGroup( 0 )
        >>> g.is_abelian
        True

        """
        if self.rank == 0:
            return True
        return False

    def contains(self, g):
        """Test if Free Group element ``g`` belong to self, ``G``.

        In mathematical terms any linear combination of generators
        of a Free Group is contained in it.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> f.contains(f[0]**3*f[1]**2)
        True

        """
        if not isinstance(i, FreeGroupElm):
            return False
        elif self != i.group:
            return False
        else:
            return True

    def assign_variables(self):
        """
        If  self is a group, whose generators are represented by symbols
        (for example a free group, a finitely presented group or a pc group)
        this function assigns these generators to global variables with the
        same names.
        """
        pass


class FreeGroupElm(Basic):
    """
    Represents an element of FreeGroup. Other than
    the Identity element, all have Infinite order.
    """
    is_Identity = None
    is_AssocWord = True

    def __new__(cls, free_group, array_form, str_expr):
        obj = Basic.__new__(cls, free_group, array_form, str_expr)

        # obj._array_form is used internally by the methods
        # of `FreeGroupElm`
        obj._array_form = array_form
        return obj

    @property
    def group(self):
        """
        Returns the `FreeGroup` on which the element itself is
        defined.

        Examples
        ========

        >>> f = FreeGroup(4)
        >>> g = FreeGroup(4, "swapnil")
        >>> (f[0]**2*f[2]).group
        <free group on the generators [f0, f1, f2, f3]>

        >>> (g[1]**2*g[3]*g[1]).group
        <free group on the generators [swapnil0, swapnil1, swapnil2, swapnil3]>

        """
        return self.args[0]

    @property
    def str_expr(self):
        return self.args[2]

    @property
    def is_Identity(self):
        if self.array_form == list():
            return True
        else:
            return False

    @property
    def array_form(self):
        """
        SymPy provides two different internal kinds of representation
        of assosciative words. The first one is called the `array_form`
        which is a list containing `tuples` as its elements, where the
        size of each tuple is two. At the first position the tuple
        contains the `generator-index`, while at the second position
        of tuple contains the exponent of that generator at the position.
        Since elements (aka words) don't commute, the indexing of list
        makes that property to stay.

        The structure in `array_form` of `FreeGroupElm` is shown below,

        [ ( index_of_gen , exponent ), ( , ), ... ( , ) ]

        Note
        ====

        The representations are for internal use only, though for its
        clarity it has not been declared internal in a Python way.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> (f[0]*f[2]).array_form
        [(0, 1), (2, 1)]
        >>> (f[0]**2*f[2]*f[1]*f[0]**2).array_form
        [(0, 2), (2, 1), (1, 1), (0, 2)]

        See Also
        ========

        letter_repr

        """
        return self._array_form

    @property
    def letter_form(self):
        """
        The  letter  representation  of an `FreeGroupElm` is as a
        list of integers, each entry corresponding to a group
        generator. Inverses of the generators are represented by
        negative numbers.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroupElm( 4 )
        >>> (f[0]**3).letter_form
        [1, 1, 1]
        >>> (f[0]**2*f[3]**-2*f[0]*f[1]**-4).letter_form
        [1, 1, -4, -4, 1, -2, -2, -2, -2]
        >>> (f[0]**-2*f[1]**3*f[3]).letter_form
        [-1, -1, 2, 2, 2, 4]

        See Also
        ========

        array_form

        """
        return flatten([[i + 1]*j if j > 0 else [-i - 1]*(-j)
                        for i, j in self.array_form])

    def __str__(self):
        str_form = ""
        array_form = self.array_form
        for i in range(len(array_form)):
            if i == len(array_form) - 1:
                if array_form[i][1] == 1:
                    str_form += self.str_expr + str(array_form[i][0])
                else:
                    str_form += self.str_expr + str(array_form[i][0]) + \
                            "**" + str(array_form[i][1])
            else:
                if array_form[i][1] == 1:
                    str_form += self.str_expr + str(array_form[i][0]) + "*"
                else:
                    str_form += self.str_expr + str(array_form[i][0]) + \
                            "**" + str(array_form[i][1]) + "*"
        return str_form

    __repr__ = __str__

    def __pow__(self, other):
        if not isinstance(other, int):
            raise TypeError("exponent of type: int expected not "
                             "of type: %s" % type(other))
        if other == 0:
            return IdentityElm()

        if other < 0:
            other = -other
            return (self.inverse())**other

        if len(self.array_form) == 1:
            new_array = [(i, other*j) for i, j in self.array_form]
            return FreeGroupElm(self.group, new_array,
                    str_expr=self.str_expr)
        else:
            result = self
            for i in range(other - 1):
                result = result*self
            return result

    def __mul__(self, other):
        """Returns the product of elements belonging to the same `FreeGroup`.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4, "swapnil" )
        >>> f[0]*f[1]**2*f[1]**-4
        swapnil0*swapnil1**-2
        >>> f[2]*f[1]**-2
        swapnil2*swapnil1**-2

        """
        if not isinstance(other, (FreeGroupElm, IdentityElm)) or \
                self.group != other.group:
            raise TypeError("only FreeGroup elements of same FreeGroup can "
                             "be multiplied")
        new_array = self.array_form[:-1]
        if self.array_form[-1][0] == other.array_form[0][0]:
            new_array.append((self.array_form[-1][0], self.array_form[-1][1] + other.array_form[0][1]))
            new_array += other.array_form[1:]
        else:
            new_array.append(self.array_form[-1])
            new_array += other.array_form[:]
        a = new_array
        for i in range(len(a) - 1, -1, -1):
            if a[i][1] == 0:
                a[i:i + 1] = []
        return FreeGroupElm(self.group, new_array, str_expr=self.str_expr)

    def __div__(self, other):
        return self*(other.inverse())

    def __len__(self):
        """
        For a FreeGroup element self, `len()` returns the number of
        letters in self.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup(2); gens = f.generators
        >>> a := gens[0];; b := gens[1];;w := a**5*b*a**2*b**-4*a;;
        >>> w
        a**5*b*a**2*b**-4*a
        >>> len( w )
        13
        >>> len( a**17 )
        17
        >>> len( w**0 )
        0

        """
        return sum([abs(j) for i, j in self.array_form])

    def inverse(self):
        """
        Returns the inverse of a `FreeGroupElm` element

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 2, "swapnil" )
        >>> f[0].inverse()
        swapnil0**-1
        >>> (f[0]*f[1]).inverse()
        swapnil1**-1*swapnil0**-1

        """
        new_array = [(i, -j) for i, j in self.array_form[::-1]]
        return FreeGroupElm(self.group, new_array, self.str_expr)

    def __eq__(self, other):
        """
        Two  associative words are equal if they are words over the
        same alphabet and if they are sequences of the same letters.
        This is equivalent to saying that the external representations
        of the words are equal.
        There is no "universal" empty word, every alphabet has its own
        empty word.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f= FreeGroup( 2, "swapnil" )
        >>> f
        <free group on the generators [swapnil0, swapnil1]>
        >>> f[0] == f[1]
        False
        >>> f[0]*f[1] == f[1]/f[1]*f[0]*f[1]
        True
        >>> f[0]*f[1] == f[1]*f[0]
        False

        """
        if not isinstance(other, FreeGroupElm):
            return False
        if self.group != other.group:
            return False
        if self.array_form != other.array_form:
            return False
        return True

    def __lt__(self, other):
        """
        The  ordering  of  associative  words is defined by length and
        lexicography (this ordering is called short-lex ordering), that
        is, shorter words are smaller than longer words, and words of the
        same length are compared w.r.t. the lexicographical ordering induced
        by the ordering of generators. Generators  are  sorted  according
        to the order in which they were created. If the generators are
        invertible then each generator g is larger than its inverse g**-1,
        and g**-1 is larger than every generator that is smaller than g.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> a = FreeGroup( 4 )
        >>> a[1] < a[0]
        False
        >>> a[0] < a[0].inverse()
        False
        """
        if not isinstance(other, FreeGroupElm) or self.group != other.group:
            raise TypeError("only FreeGroup elements of same FreeGroup can "
                             "be compared")
        a = self.letter_form
        b = other.letter_form
        l = len(self)
        m = len(other)
        # implement lenlex order
        if l < m:
            return True
        elif l > m:
            return False
        for i in range(l):
            p = abs(a[i])
            q = abs(b[i])
            if p < q:
                return True
            elif p > q:
                return False
            elif a[i] < b[i]:
                return True
            elif a[i] > b[i]:
                return False
        return False

    def __le__(self, other):
        return (self == other or self < other)

    def __gt__(self, other):
        """

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> b = FreeGroup( 3 )
        >>> b[1]**2 > b[0]**2
        True
        >>> b[1]*b[2] > b[2]*b[1]
        False
        >>> b[0] > b[0].inverse()
        True

        """
        if not isinstance(other, FreeGroupElm) or self.group != other.group:
            raise TypeError("only FreeGroup elements of same FreeGroup can "
                             "be compared")
        return not self <= other

    def __ge__(self, other):
        return not self < other

    def exponent_sum_word(self, gen):
        """
        For an associative word `w` and a generator `gen`, ``exponent_sum_word``
        returns the number of times `gen` appears in `w` minus the number of
        times its inverse appears in `w`. If both `gen` and its inverse do not
        occur in `w` then 0 is returned. `gen` may also be the inverse of a
        generator.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> a, b = f[0], f[1]
        >>> w = a**5*b*a**2*b**-4*a
        >>> w.exponent_sum_word(a)
        8
        >>> w.exponent_sum_word(b)
        3
        >>> ((a*b*a**-1)**3).exponent_sum_word(a)
        0
        >>> w.exponent_sum_word(b**-1)
        3

        """
        w = self.letter_form
        gen = gen.letter_form
        if len(gen) != 1:
            raise ValueError("<gen> must be a generator")
        n = 0
        g = abs(gen[0])
        for i in w:
            if i == g:
                n = n + 1
            elif i == -g:
                n = n - 1

        if gen[0] < 0:
            n = -n
        return n

    def subword(self, i, j):
        """
        For  an associative word w and two positive integers from and to, Subword
        returns the subword of w that begins at position from and ends at position to.
        Indexing is done with origin 1.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> w = f[0]**2*f[1]**4
        >>> w.subword(1, 4)
        f0*f1**2

        """
        if i < 0 or j >= len(self) or j < i:
            raise ValueError("`i`, `j` must be within bounds, and `i` should "
                            "be less than or equal to `j`")
        if i == j:
            return IdentityElm()
        else:
            letter_form = self.letter_form[i: j]
            array_form = letter_form_to_array_form(letter_form)
            return FreeGroupElm(self.group, array_form, self.str_expr)

    def AssocWordByLetterRep(self, lrep):
        """
        """
        pass

    def number_syllables(self):
        """Returns the number of syllables of the associative word w.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4, "swapnil" )
        >>> (f[1]**3*f[0]*f[1]**-1).number_syllables()
        3

        """
        return len(self.array_form)

    def exponent_syllable(self, i):
        """
        Returns the exponent of the i-th syllable of the associative word w.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> (f[0]**2*f[1]**-3*f[1]).exponent_syllable(1)
        -3

        """
        return self.array_form[i][1]

    def generator_syllable(self, i):
        """
        Returns the number of the generator that is involved in the
        i-th syllable of the associative word w.
        """
        pass


class IdentityElm(Basic):
    """

    """
    is_Identity = True

    def __new__(cls, *args):
        obj = Basic.__new__(cls, *args)
        return obj

    def __mul__(self, other):
        return other

    __rmul__ = __mul__

    def __pow__(self, other):
        return self

    def __div__(self, other):
        return other.inverse()

    def __rdiv__(self, other):
        return other

    def __str__(self):
        return '<identity ...>'

    __repr__ = __str__


def letter_form_to_array_form(array_form):
    a = array_form[:]
    new_array = []
    n = 1
    for i in range(len(a)):
        if i == len(a) - 1:
            if a[i] == a[i - 1]:
                new_array.append((a[i] -1, n))
            else:
                new_array.append((a[i] - 1, 1))
            return new_array
        elif a[i] == a[i + 1]:
            n += 1
        else:
            new_array.append((a[i] - 1, n))
            n = 1

#class Group(Basic):
#    pass
#
#
#class GroupElm(Basic):
#    """
#    Group element definition
#
#    """
#
#    def __init__(self, elem, group=None):
#        if group is None:
#            raise ValueError("")
#        if not elem in group.Set:
#            raise TypeError("elem is not an element of group")
#        self.elem = elem
#        self.group = group
#
#    def __str__(self):
#        return str(self.elem)
#
#    def __eq__(self, other):
#        """
#        Two GroupElems are equal if they represent the same element,
#        regardless of the Groups they belong to
#        """
#
#        if not isinstance(other, GroupElem):
#            raise TypeError("other is not a GroupElem")
#        return self.elem == other.elem
#
#    def __ne__(self, other):
#        return not self == other
#
#    def __hash__(self):
#        return hash(self.elem)
#
#    def __mul__(self, other):
#        """
#        If other is a group element, returns self * other.
#        If other = n is an int, and self is in an abelian group, returns self**n
#        """
#        if self.group.is_abelian() and isinstance(other, (int, long)):
#            return self ** other
#
#        if not isinstance(other, GroupElem):
#            raise TypeError("other must be a GroupElem, or an int " \
#                            "(if self's group is abelian)")
#        try:
#            return GroupElem(self.group.bin_op((self.elem, other.elem)), \
#                             self.group)
#        # This can return a TypeError in Funcion.__call__ if self and other
#        # belong to different Groups. So we see if we can make sense of this
#        # operation the other way around.
#        except TypeError:
#            return other.__rmul__(self)
#
#    def __rmul__(self, other):
#        """
#        If other is a group element, returns other * self.
#        If other = n is an int, and self is in an abelian group, returns self**n
#        """
#        if self.group.is_abelian() and isinstance(other, (int, long)):
#            return self ** other
#
#        if not isinstance(other, GroupElem):
#            raise TypeError("other must be a GroupElem, or an int " \
#                            "if self's group is abelian")
#
#        return GroupElem(self.group.bin_op((other.elem, self.elem)), self.group)
#
#    def __add__(self, other):
#        """Returns self + other for Abelian groups"""
#        if self.group.is_abelian():
#            return self * other
#        raise TypeError("not an element of an abelian group")
#
#    def __pow__(self, n, modulo=None):
#        """
#        Returns self**n
#
#        modulo is included as an argument to comply with the API, and ignored
#        """
#        if not isinstance(n, (int, long)):
#            raise TypeError("n must be an int or a long")
#
#        if n == 0:
#            return self.group.e
#        elif n < 0:
#            return self.group.inverse(self) ** -n
#        elif n % 2 == 1:
#            return self * (self ** (n - 1))
#        else:
#            return (self * self) ** (n / 2)
#
#    def __neg__(self):
#        """Returns self ** -1 if self is in an abelian group"""
#        if not self.group.is_abelian():
#            raise TypeError("self must be in an abelian group")
#        return self ** (-1)
#
#    def __sub__(self, other):
#        """Returns self * (other ** -1) if self is in an abelian group"""
#        if not self.group.is_abelian():
#            raise TypeError("self must be in an abelian group")
#        return self * (other ** -1)
#
#    def order(self):
#        """Returns the order of self in the Group"""
#        return len(self.group.generate([self]))
#
#class GroupIdentity(Basic):
#    """
#    Represents the idenity for any Group.
#    """
#    is_Identity = True
#    is_Abelian = True
#
#    def __repr__(self):
#        return "<identity ...>"
#
#    def order(self):
#        return 1
#
#    def __mul__(self, other):
#        return other
#
#    __rmul__ = __mul__
#
#    def __pow__(self, other):
#        return self
#
#
#def One(group):
#    """
#    Returns the Identity Element of Group `group`.
#    """
#    pass
#
#
#class InfiniteGenerators(Basic):
#    """
#    Returns an InfiniteList of Generators.
#    It prints the Infinite list of generators in the form
#    of a limited with two elements [ f0, f1 ... ]
#    """
#    pass
