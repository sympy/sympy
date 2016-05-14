from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.core.compatibility import as_int, string_types
from sympy.core.symbol import Symbol, symbols as _symbols
from sympy.core.sympify import CantSympify
from mpmath import isint
from collections.abc import Sequence
from sympy.core import S
from sympy.utilities import public
from sympy.utilities.iterables import flatten
from sympy.core.power import Pow
from sympy.utilities.magic import pollute


@public
def free_group(string_gens):
    symbols = tuple(_parse_symbols(string_gens))
    _free_group = FreeGroup(symbols)
    return (_free_group,) + tuple(_free_group.generators)

@public
def xfree_group(string_gens):
    symbols = tuple(_parse_symbols(string_gens))
    _free_group = FreeGroup(symbols)
    return (_free_group, _free_group.generators)

@public
def vfree_group(string_gens):
    symbols = tuple(_parse_symbols(string_gens))
    _free_group = FreeGroup(symbols)
    pollute([sym.name for sym in _free_group.symbols], _free_group.generators)
    return _free_group


def _parse_symbols(symbols):
    if isinstance(symbols, string_types):
        return _symbols(symbols, seq=True)
    elif isinstance(symbols, Expr):
        return (symbols,)
    elif is_sequence(symbols):
        if all(isinstance(s, string_types) for s in symbols):
            return _symbols(symbols)
        elif all(isinstance(s, Expr) for s in symbols):
            return symbols
    raise ValueError("")


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

    ``FreeGroup(rank)``                                      ............ (1)

    Called in the second form, ``FreeGroup`` returns a free group on as many
    generators as arguments, printed as `name0`, `name1` etc.

    ``FreeGroup(rank, "name")``                              ............ (2)

    Called in the third form, ``FreeGroup`` returns a free group on as many
    generators as the length of the list `names`, the i-th generator being
    printed as `names[i]`.

    ``FreeGroup("string0", "string1", "string2", .. , "string_n-1" )`` ... (3)

    Called in the fourth form, ``FreeGroup`` returns a free group on
    infinitely many generators, where the first generators are printed by the
    names in the list `init`, and the other generators by `name` and an appended
    number. Like ``FreeGroup( S.Infinity )``.

    ``FreeGroup( S.Infinity, name, init )``                   ............ (4)

    References
    ==========

    [1] https://en.wikipedia.org/wiki/Free_group

    [2] https://www.gap-system.org

    """
    is_associative = True
    is_group = True
    is_FreeGroup = True
    is_PermutationGroup = False

    def __new__(cls, symbols):
        """
        Since any number of arguments can be passed in the form of string,
        hence `rank` is "not" used as argument.
        """

        obj = Basic.__new__(cls, symbols)
        rank = len(symbols)
        obj._rank = rank
        obj.dtype = type("FreeGroupElm", (FreeGroupElm,), {"group": obj})
        obj.symbols = symbols
        obj.generators = obj._generators()
        #obj._gens_set = set(obj.generators)
        for symbol, generator in zip(obj.symbols, obj.generators):
            if isinstance(symbol, Symbol):
                name = symbol.name
                if hasattr(obj, name):
                    setattr(obj, name, generator)
        return obj

    def _generators(group):
        """Returns the generators of the FreeGroup

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y, z = free_group("x, y, z")
        >>> f.generators
        [x, y, z]

        """
        gens = []
        for i in range(group.rank):
            elm = ((group.symbols[i], 1),)
            gens.append(group.dtype(elm))
        return gens

    def __getitem__(self, i):
        return self.generators[i]

    def __contains__(self, i):
        """
        Return True if `i` is contained in FreeGroup.

        Examples
        ========

        """
        if not isinstance(i, FreeGroupElm):
            raise TypeError("FreeGroup contains only FreeGroupElm as elements "
                        ", not elements of type %s" % type(i))
        group = i.group
        return self == group

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

        """
        return self is other

    def order(self):
        if self.rank == 0:
            return 1
        else:
            return S.Infinity

    @property
    def elements(self):
        if self.rank == 0:
            # A set containing Identity element of `FreeGroup` self is returned
            return set([self.identity])
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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y, z = free_group("x y z")
        >>> f.is_abelian
        False

        """
        # this needs to be tested
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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y, z = free_group("x y z")
        >>> f.contains(x**3*y**2)
        True

        """
        if not isinstance(g, FreeGroupElm):
            return False
        elif self != g.group:
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


class FreeGroupElm(CantSympify, tuple):
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
    def is_identity(self):
        if self.array_form == tuple():
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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y, z = free_group("x y z")
        >>> (x*z).array_form
        [(0, 1), (2, 1)]
        >>> (x**2*z*y*x**2).array_form
        [(0, 2), (2, 1), (1, 1), (0, 2)]

        See Also
        ========

        letter_repr

        """
        return tuple(self)

    @property
    def letter_form(self):
        """
        The  letter  representation  of an `FreeGroupElm` is as a
        list of integers, each entry corresponding to a group
        generator. Inverses of the generators are represented by
        negative numbers.

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, a, b, c, d = free_group("a b c d")
        >>> (a**3).letter_form
        [1, 1, 1]
        >>> (a**2*d**-2*a*b**-4).letter_form
        [1, 1, -4, -4, 1, -2, -2, -2, -2]
        >>> (a**-2*b**3*d).letter_form
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

        symbols = self.group.symbols
        str_form = ""
        array_form = self.array_form
        for i in range(len(array_form)):
            if i == len(array_form) - 1:
                if array_form[i][1] == 1:
                    str_form += str(array_form[i][0])
                else:
                    str_form += str(array_form[i][0]) + \
                                    "**" + str(array_form[i][1])
            else:
                if array_form[i][1] == 1:
                    str_form += str(array_form[i][0]) + "*"
                else:
                    str_form += str(array_form[i][0]) + \
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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y, z = free_group("x y z")
        >>> x*y**2*y**-4
        x*y**-2
        >>> z*y**-2
        z*y**-2
        >>> x**2*y*y**-1*x**-2
        <identity>

        """
        group = self.group
        if not isinstance(other, group.dtype):
            raise TypeError("only FreeGroup elements of same FreeGroup can "
                    "be multiplied")
        if self.is_identity:
            return other
        if other.is_identity:
            return self
        r = tuple(zero_mul_simp(list(self.array_form + other.array_form),
                        len(self.array_form) - 1))
        return group.dtype(r)

    def __div__(self, other):
        return self*(other.inverse())

    def __rdiv__(self, other):
        return other*(self.inverse())

    def inverse(self):
        """
        Returns the inverse of a `FreeGroupElm` element

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y, z = free_group("x y z")
        >>> x.inverse()
        x**-1
        >>> (x*y).inverse()
        y**-1*x**-1

        """
        group = self.group
        r = tuple([(i, -j) for i, j in self.array_form])
        return group.dtype(r)

    def order(self):
        """Find the order of a `FreeGroupElm`.

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y = free_group("x y")
        >>> (x**2*y*y**-1*x**-2).order()
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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y = free_group("x y")
        >>> w = x**5*y*x**2*y**-4*x
        >>> w.eliminate_word( x, x**2 )
        x**10*y*x**4*y**-4*x**2
        >>> w.eliminate_word( x, y**-1 )
        y**-5*y*y**-2*y**-4*y**-1

        """
        group = self.group
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
        # zero_mul_simp to be used
        return group.dtype(l)

    def __len__(self):
        """
        For an associative word `self`, this returns the number of letters in it.

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, a, b = free_group("a b")
        >>> w = a**5*b*a**2*b**-4*a
        >>> len(w)
        13
        >>> len(a**17)
        17
        >>> len(w**0)
        0

        """
        return sum([abs(j) for (i, j) in self])

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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, swapnil0, swapnil1 = free_group("swapnil0 swapnil1")
        >>> f
        <free group on the generators [swapnil0, swapnil1]>
        >>> g, swap0, swap1 = free_group("swap0 swap1")
        >>> g
        <free group on the generators [swapl0, swap1]>

        >>> swapnil0 == swapnil1
        False
        >>> swapnil0*swapnil1 == swapnil1/swapnil1*swapnil0*swapnil1
        True
        >>> swapnil0*swapnil1 == swapnil1*swapnil0
        False

        >>> swapnil1**0 == swap0**0
        False

        """
        group = self.group
        if not isinstance(other, group.dtype):
            return False
        return tuple.__eq__(self, other)

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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, a, b = free_group("a b")
        >>> b < a
        False
        >>> a < a.inverse()
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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, x, y, z = free_group("x y z")
        >>> y**2 > x**2
        True
        >>> y*z > z*y
        False
        >>> x > x.inverse()
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
        For an associative word `self` and a generator `gen`, ``exponent_sum_word``
        returns the number of times `gen` appears in `self` minus the number of
        times its inverse appears in `self`. If both `gen` and its inverse do
        not occur in `self` then 0 is returned. `gen` may also be the inverse of
        a generator.

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, a, b = free_group("a b")
        >>> w = a**5*b*a**2*b**-4*a
        >>> w.exponent_sum_word(a)
        8
        >>> w.exponent_sum_word(b)
        -3
        >>> ((a*b*a**-1)**3).exponent_sum_word(a)
        0
        >>> w.exponent_sum_word(b**-1)
        3

        """
        w = self.letter_form
        gen = gen.letter_form
        if len(gen) != 1:
            raise ValueError("gen must be a generator or inverse of a generator")
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
        
        >>> from sympy.combinatorics.free_group import free_group
        >>> f, a, b = free_group("a b")
        >>> w = a**5*b*a**2*b**-4*a
        >>> w.subword(2, 6)
        x**3*y

        """
        group = self.group
        if from_i < 0 or to_j >= len(self):
            raise ValueError("`from_i`, `to_j` must be positive and less than "
                    "the length of associative word")
        if to_j <= from_i:
            return group.identity
        else:
            letter_form = self.letter_form[from_i: to_j]
            array_form = letter_form_to_array_form(letter_form)
            return group.dtype(array_form)

    def assoc_word_by_letter_rep(self, lrep):
        """
        """
        pass


    def number_syllables(self):
        """Returns the number of syllables of the associative word `self`.

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f = free_group("swapnil0 swapnil1")
        >>> swapnil0, swapnil1 = f[0][0], f[0][1]
        >>> (swapnil1**3*swapnil0*swapnil1**-1).number_syllables()
        3

        """
        return len(self.array_form)

    def exponent_syllable(self, i):
        """
        Returns the exponent of the `i`-th syllable of the associative word
        `self`.

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, a, b = free_group("a b")
        >>> w = a**5*b*a**2*b**-4*a
        >>> w.exponent_syllable( 2 )
        2

        """
        return self.array_form[i][1]

    def generator_syllable(self, i):
        """
        Returns the number of the generator that is involved in the
        i-th syllable of the associative word `self`.

        Examples
        ========

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, a, b = free_group("a b")
        >>> w = a**5*b*a**2*b**-4*a
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

        >>> from sympy.combinatorics.free_group import free_group
        >>> f, a, b = free_group("a b")
        >>> w = a**5*b*a**2*b**-4*a
        >>> w.sub_syllables(1, 2)
        f1
        >>> w.sub_syllables(3, 3)
        <identity>

        """
        if not isinstance(from_i, int) or not isinstance(to_j, int):
            raise ValueError("both arguments should be integers")
        group = self.group
        if to_j <= from_i:
            return group.identity
        else:
            r = tuple(self.array_form[from_i: to_j])
            return group.dtype(r)

    def substitute_word(self, from_i, to_j, by):
        """
        """
        lw = len(self)
        if from_i > to_j or from_i > lw or to_j > lw:
            raise ValueError("values should be within bounds")

        # otherwise there are four possibilities

        # first if from=1 and to=lw then
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
    a = list(array_form[:])
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


def zero_mul_simp(array_form, index):
    while index >= 0 and index < len(array_form) - 1:
        if array_form[index][0] == array_form[index + 1][0]:
            updated_exp = array_form[index][1] + array_form[index + 1][1]
            updated_base = array_form[index][0]
            if updated_exp == 0:
                del array_form[index], array_form[index]
                return zero_mul_simp(array_form, index - 1)
            else:
                array_form[index] = (updated_base, updated_exp)
                del array_form[index + 1]
        else:
            return array_form
    return array_form
