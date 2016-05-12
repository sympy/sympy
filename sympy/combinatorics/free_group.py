from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.core.compatibility import as_int
from sympy.core.sympify import CantSympify
from collections.abc import Sequence
from sympy.core import S
from sympy.utilities import public
from sympy.utilities.iterables import flatten
from sympy.core.power import Pow

@public
def free_group(rank, str_expr="f"):
    _free_group = FreeGroup(rank, str_expr)
    return (_free_group,) + tuple(_free_group.generators)


##############################################################################
#                          FREE GROUP                                        #
##############################################################################


class FreeGroup(Basic):
    """
    Called with a positive integer `rank`, ``FreeGroup`` returns a free group
    on `rank` generators. If the optional argument `name` is given then the
    generators are printed as `name0`, `name1` etc., that is, each name is the
    concatenation of the string `name` and an integer from `0` to `range-1`.
    The default for `name` is the string "f".

    ``FreeGroup( rank )``     ..........................................  (1)

    Called in the second form, ``FreeGroup`` returns a free group on as many
    generators as arguments, printed as `name0`, `name1` etc.

    ``FreeGroup( rank, "name" )`` ......................................  (2)

    Called in the third form, ``FreeGroup`` returns a free group on as many
    generators as the length of the list `names`, the i-th generator being
    printed as `names[i]`.

    ``FreeGroup( "string0", "string1", "string2", ... )`` ..............  (3)

    Called in the fourth form, ``FreeGroup`` returns a free group on
    infinitely many generators, where the first generators are printed by the
    names in the list `init`, and the other generators by `name` and an appended
    number. Like ``FreeGroup( S.Infinity )``.

    ```FreeGroup( S.Infinity, name, init )```   ........................  (4)

    References
    ==========

    [1] https://en.wikipedia.org/wiki/Free_group

    [2] https://www.gap-system.org

    """
    is_associative = True
    is_group = True
    is_FreeGroup = True
    is_PermutationGroup = False

    def __new__(cls, *args, **kwargs):
        """
        Since any number of arguments can be passed in the form of string,
        hence `rank` is "not" used as argument.
        """

        obj = Basic.__new__(cls, *args, **kwargs)

        if as_int(args[0]) and args[0] >= 0:

            # (1) form of the API used here
            if len(args) == 1:
                obj._str = "f"

            # (2) form of the API used here
            elif len(args) == 2:
                if not isinstance(args[1], str):
                    raise ValueError("second argument should be of type `str`, "
                            "not of type %s" % type(args[1]))
                obj._str = args[1]

            else:
                raise ValueError("Invalid arguments")

            obj._rank = args[0]

        # (3) form of API, where all the generators of `FreeGroup`
        # given as strings.
        elif all([isinstance(i, str) for i in args]):
            pass

        # (4) form of API, not sure how this should be implemented
        # right now.
        elif args[0] is S.Infinity:
            pass

        obj.dtype = type("FreeGroupElm", (FreeGroupElm,), {"group": obj,
                        "str_expr": obj._str})
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
        _gens = []
        for i in range(self.rank):
            elm = self.identity
            elm.append((i, 1))
            _gens.append(elm)
        return _gens

    def __getitem__(self, i):
        return self.generators[i]

    def __contains__(self, i):
        """
        Return True if `i` is contained in FreeGroup.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4, "swap" )
        >>> g = FreeGroup( 4, "swapnil" )

        >>> f[0]**2*f[3] in f
        True
        >>> f[0]**2*f[3] in g
        False

        """
        if not isinstance(i, FreeGroupElm):
            raise TypeError("FreeGroup contains elements of type "
                    "`FreeGroupElm`")
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
        """No ``FreeGroup`` is equal to any "other" ``FreeGroup``.

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
        return self is other

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
            # A set containing Identity element of `FreeGroup` self is returned
            return FreeGroupElm(self, [])
        else:
            raise ValueError("Group contains infinitely many elements"
                            ", hence can't be represented")

    @property
    def rank(self):
        r"""
        In group theory, the `rank` of a group `G`, denoted `G.rank`,
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
        if self.rank == 0 or self.rank == 1:
            return True
        else:
            return False

    @property
    def identity(self):
        return self.dtype()

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

    def is_subgroup(self, F):
        """Return True if all elements of `self` belong to `F`.

        Examples
        ========
        """
        return F.is_group and all([self.contains(gen) for gen in F.generators])

    def assign_variables(self):
        """
        If  self is a group, whose generators are represented by symbols
        (for example a free group, a finitely presented group or a pc group)
        this function assigns these generators to global variables with the
        same names.
        """
        pass


############################################################################
#                          FreeGroupElm                                    #
############################################################################


class FreeGroupElm(CantSympify, list):
    """
    ``FreeGroupElm`` is actually not usable as public import, since it is
    always associated with a ``FreeGroup``. It is always called only by the
    ``FreeGroup`` class. Called with a ``FreeGroup`` as the first argument and
    the second argument called the `array_form` for the ``FreeGroupElm``, and
    the third argument is a string representation for the element.

    """
    is_identity = None
    is_AssocWord = True

    @property
    def expt(self):
        return self._expt

    @property
    def is_identity(self):
        if self.array_form == list():
            return True
        else:
            return False

    @property
    def array_form(self):
        """
        SymPy provides two different internal kinds of representation
        of associative words. The first one is called the `array_form`
        which is a list containing `tuples` as its elements, where the
        size of each tuple is two. At the first position the tuple
        contains the `generator-index`, while at the second position
        of tuple contains the exponent of that generator at the position.
        Since elements (i.e. words) don't commute, the indexing of list
        makes that property to stay.

        The structure in `array_form` of `FreeGroupElm` is shown below,

        [  ( index_of_gen , exponent ), ( , ), ... ( , )  ]

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
        return list([i for i in self])

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
        # ** Warning **
        # Note here that the representation adds 1 or -1 to make it represent
        # since 0 removes the `-` sign from it, hence it making it non-usable,
        # a non-pythonic way. But no other option is there.
        return flatten([[i + 1]*j if j > 0 else [-i - 1]*(-j)
                        for i, j in self.array_form])

    def ext_rep(self):
        """This is called the External Representation of `FreeGroupElm`
        """
        return flatten(self.array_form)

    def __str__(self):
        if self.is_identity:
            return "<identity>"

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

    def __pow__(self, n):
        n = as_int(n)
        group = self.group
        if n == 0:
            return group.identity

        if n < 0:
            n = -n
            return (self.inverse())**n

        result = self
        for i in range(n - 1):
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
        >>> (f[0]**2*f[1]*f[1]**-1*f[0]**-2)
        <identity>

        """
        group = self.group
        r = group.identity
        if not isinstance(other, group.dtype):
            raise TypeError("only FreeGroup elements of same FreeGroup can "
                    "be multiplied")
        if self.is_identity:
            return other

        if self.array_form[-1][0] == other.array_form[0][0]:
            r.extend(self.array_form[:-1])
            r.append((self.array_form[-1][0], self.array_form[-1][1] +
                other.array_form[0][1]))
            r.extend(other.array_form[1:])
        else:
            r.extend(self.array_form + other.array_form)

        return r

    def __div__(self, other):
        return self*(other.inverse())

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
        group = self.group
        r = group.identity
        r.extend([(i, -j) for i, j in self.array_form])
        return r

    def order(self):
        """Find the order of a `FreeGroupElm`.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> (f[0]**2*f[1]*f[1]**-1*f[0]**-2).order()
        1

        """
        if self.is_identity:
            return 1
        else:
            return S.Infinity

    def commutator(self, other):
        """Returns the commutator of self and x: ``~x*~self*x*self``
        """
        group = self.group
        if not isinstance(other, group.dtype):
            raise ValueError("commutator of only `FreeGroupElm` of the same "
                    "`FreeGroup` exists")
        else:
            return self.inverse()*other.inverse()*self*other

    def eliminate_word(self, gen, by):
        """
        For an associative word `self`, a generator `gen`, and an associative word
        by, `eliminate_word` returns the associative word obtained by replacing
        each occurrence of `gen` in `self` by `by`.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> w = f[0]**5*f[1]*f[0]**2*f[1]**-4*f[0]
        >>> w.eliminate_word( f[0], f[0]**2 )
        f0**10*f1*f0**4*f1**-4*f0**2

        >>> w.eliminate_word( f[0], f[1]**-1 )
        f1**-11

        """
        e = self.ext_rep()
        gen = gen.generator_syllable(0)
        l = []
        for i in range(0, len(e) - 1, 2):
            if e[i] == gen:
                app = (by**e[i + 1]).ext_rep()
            else:
                app = e[i: i + 2]
            j = len(l) - 1
            while j > 0 and len(app) > 0 and l[j - 1] == app[0]:
                s = l[j] + app[1]
                if s == 0:
                    j = j - 2
                else:
                    l[j] = s
                app = app[2: len(app)]

            if j + 1 < len(l):
                l = l[0: j + 1]

            if len(app) > 0:
                l.append(tuple(app))
        zero_simp(l)
        mult_simp(l)
        zero_simp(l)
        return FreeGroupElm(self.group, l, self.str_expr)

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
        >>> f = FreeGroup( 2, "swapnil" )
        >>> f
        <free group on the generators [swapnil0, swapnil1]>
        >>> g = FreeGroup( 2 , "swapnil" )
        >>> g
        <free group on the generators [swapnil0, swapnil1]>

        >>> f[0] == f[1]
        False
        >>> f[0]*f[1] == f[1]/f[1]*f[0]*f[1]
        True
        >>> f[0]*f[1] == f[1]*f[0]
        False

        >>> f[1]**0 == g[0]**0
        False

        """
        group = self.group
        if not isinstance(other, group.dtype):
            return False
        return list.__eq__(self, other)

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
        >>> a = FreeGroup( 4, "swapnil" )
        >>> a[1] < a[0]
        False
        >>> a[0] < a[0].inverse()
        False
        """
        group = self.group
        if not isinstance(other, group.dtype):
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
        >>> b = FreeGroup( 3, "swapnil" )
        >>> b[1]**2 > b[0]**2
        True
        >>> b[1]*b[2] > b[2]*b[1]
        False
        >>> b[0] > b[0].inverse()
        True

        """
        group = self.group
        if not isinstance(other, group.dtype):
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

    def subword(self, from_i, to_j):
        """
        For an associative word `self` and two positive integers `from_i` and
        `to_j`, subword returns the subword of `self` that begins at position
        `from_to` and ends at `to_j`, indexing is done with origin 0.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> w = f[0]**5*f[1]*f[0]**2*f[1]**-4*f[0]
        >>> w.subword(2, 6)
        f0**3*f1

        """
        if from_i < 0 or to_j >= len(self):
            raise ValueError("`from_i`, `to_j` must be positive and less than "
                    "the length of associative word")
        if to_j <= from_i:
            return FreeGroupElm(self.group, [])
        else:
            letter_form = self.letter_form[from_i: to_j]
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
        Returns the exponent of the `i`-th syllable of the associative word
        `self`.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4, "swap" )
        >>> w = f[0]**5*f[1]*f[0]**2*f[1]**-4*f[0]
        >>> w.exponent_syllable( 2 )
        2

        """
        return self.array_form[i][1]

    def generator_syllable(self, i):
        """
        Returns the number of the generator that is involved in the
        i-th syllable of the associative word w.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 2 )
        >>> w = f[0]**5*f[1]*f[0]**2*f[1]**-4*f[0]
        >>> w.generator_syllable( 3 )
        1

        """
        return self.array_form[i][0]

    def sub_syllables(self, from_i, to_j):
        """
        `sub_syllables` returns the subword of the associative word `self` that
        consists of syllables from positions `from_to` to `to_j`, where
        `from_to` and `to_j` must be positive integers and indexing is done
        with origin 0.

        Examples
        ========

        >>> from sympy import FreeGroup
        >>> f = FreeGroup( 4 )
        >>> w = f[0]**5*f[1]*f[0]**2*f[1]**-4*f[0]
        >>> w.sub_syllables(1, 2)
        f1
        >>> w.sub_syllables(3, 3)
        <identity ...>

        """
        if not isinstance(from_i, int) or not isinstance(to_j, int):
            raise ValueError("both arguments should be integers")
        if to_j <= from_i:
            return FreeGroupElm(self.group, [])
        else:
            return FreeGroupElm(self.group, self.array_form[from_i: to_j],
                    self.str_expr)

    def substitute_word(self, from_i, to_j, by):
        """
        """
        lw = len(self)
        if from_i > to_j or from_i > lw or to_j > lw:
            raise ValueError("values should be within bounds")

        # otherwise there are four possibilities

        # first if from=1 and to=Length(w) then
        if from_i == 0 and to_j == lw - 1:
            return by
        elif from_i == 0:  # second if from_i=1 (and to_j < lw) then
            return by*self.subword(to_j, lw - 1)
        elif to_j == lw:   # third if to_j=1 (and fromi_i > 1) then
            return self.subword(0, from_i - 1)*by;
        else:              # finally
            return self.subword(0, from_i - 1)*by*self.subword(to_j + 1, lw)


def letter_form_to_array_form(array_form):
    """
    This method converts a list given with possible repetitions of elements in
    it. It returns a new list such that repetitions of consecutive elements is
    removed and replace with a tuple element of size two such that the first
    index contains `value` and the second index contains the number of
    consecutive repetitions of `value`.

    Examples
    ========

    >>> l = [1, 1, 1, -2, -2, 5, 1, 1, 1, 1, 4, 4, 4, 4]
    >>> letter_form_to_array_form(l)
    [(0, 3), (-1, 2), (4, 1), (0, 3), (3, 4)]

    """
    a = array_form[:]
    new_array = []
    n = 1
    for i in range(len(a)):
        if i == len(a) - 1:
            if a[i] == a[i - 1]:
                # 1 has been subtracted in accordance with the meaning that in
                # letter-form `-ve` sign indicates presence of inverse which is
                # not possi
                new_array.append((a[i] -1, n))
            else:
                new_array.append((a[i] - 1, 1))
            return new_array
        elif a[i] == a[i + 1]:
            n += 1
        else:
            new_array.append((a[i] - 1, n))
            n = 1


def zero_simp(array_form):
    for i in range(len(array_form) - 1, -1, -1):
        if array_form[i][1] == 0:
            del array_form[i]


def mult_simp(array_form):
    for i in range(len(array_form) - 1, 0, -1):
        if array_form[i][0] == array_form[i - 1][0]:
            array_form[i] = (array_form[i][0],
                    array_form[i][1] + array_form[i - 1][1])
            del array_form[i - 1]
