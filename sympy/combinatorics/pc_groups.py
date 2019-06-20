from sympy.core import Basic
from sympy import sieve
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.printing.defaults import DefaultPrinting

class PolycyclicGroup(Basic):

    is_group = True
    is_solvable = True

    def __init__(self, pcgs_):
        self.perm_group = PermutationGroup(pcgs_)
        self.pc_series = self._pc_series()
        self.pcgs = self._compute_pcgs()

    def _pc_series(self):
        return self.perm_group.composition_series()

    def _compute_pcgs(self):
        # computes the generating sequence for polycyclic groups.
        series = self.pc_series
        pcgs = []
        for i in range(len(series)-1):
            for g in series[i].generators:
                if not g in series[i+1]:
                    pcgs.append(g)
        return pcgs

    def relative_orders(self):
        rel_orders = []
        for i in range(len(self.pc_series)-1):
            G = self.pc_series[i]
            H = self.pc_series[i+1]
            rel_orders.append(G.order()//H.order())
        return rel_orders

    def is_prime_order(self):
        for order in self.relative_orders():
            if order not in sieve:
                return False
        return True

    def length(self):
        return len(self.pcgs)

    def pc_element_exponent(self, element):
        series = self.pc_series
        pcgs = self.pcgs
        exponent = [0]*len(series)
        for i in range(len(series)):
            exp = 0
            if not element in series[i]:
                for j in range(len(pcgs)):
                    element = (pcgs[j]**-1)*element
                    exp = exp + 1
                    if element in series[i]:
                        exponent[i] = exp
                        break
        return exponent


class Collector(DefaultPrinting):

    def __init__(self, pc_relators, word):
        self.pc_relators = pc_relators
        self.word = word
        self.group = word.group

    def minimal_uncollected_subwords(self):
        """
        Returns the minimal uncollected subwords.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group

        Example 8.7 Pg. 281 from Handbook
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : 1, x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> word = x2**2*x1**7
        >>> collector = Collector(pc_relators, word)
        >>> collector.minimal_uncollected_subwords()
        {((x1, 7),): 1, ((x2, 2), (x1, 1)): 0}

        """
        group = self.word.group
        gens = group.symbols
        index = {}
        i = 1
        for g in gens:
            index[g] = i
            i += 1
        array = self.word.array_form
        uncollected_subwords = {}
        for i in range(len(array)-1):
            if array[i+1][1] > 0 and index[array[i][0]] > index[array[i+1][0]]:
                # case-1:  v = x[i]**a*x[i+1]
                uncollected_subwords[((array[i][0], array[i][1]), (array[i+1][0], 1))] = 0

            elif array[i+1][1] < 0 and index[array[i][0]] > index[array[i+1][0]]:
                # case-2: v = x[i]**a*x[i+1]*-1
                uncollected_subwords[((array[i][0], array[i][1]), (array[i+1][0], -1))] = 0

            if array[i][1] < 0 or array[i][1] > index[array[i][0]]:
                # case-3: v = x[i]**a
                uncollected_subwords[((array[i][0], array[1][1]), )] = 1

        i = len(array)-1
        if array[i][1] < 0 or array[i][1] > index[array[i][0]]:
            # case-3: v = x[i]**a
            uncollected_subwords[((array[i][0], array[1][1]), )] = 1

        return uncollected_subwords

    def relations(self):
        """
        Separates the given relators of pc presentation in power and
        conjugate relations.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : 1, x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> word = x2**2*x1**7
        >>> collector = Collector(pc_relators, word)
        >>> power_rel, conj_rel = collector.relations()
        >>> power_rel
        {x1**2: 1}
        >>> conj_rel
        {x1**-1*x2*x1: x2**-1, x1*x2*x1**-1: x2**-1}

        """
        power_relators = {}
        conjugate_relators = {}
        for key, value in self.pc_relators.items():
            if len(key.array_form) == 1:
                power_relators[key] = value
            else:
                conjugate_relators[key] = value
        return power_relators, conjugate_relators

    def reduce_word(self):
        """
        Reduce the given word with the help of power
        relators only.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : 1, x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> word = x1**7
        >>> collector = Collector(pc_relators, word)
        >>> collector.reduce_word()
        x1

        """
        group = self.word.group
        power_relators = self.relations()[0]
        dividend = self.word.array_form[0][1]
        rem = dividend
        for w, v in power_relators.items():
            divisor = w.array_form[0][1]
            if w.letter_form_elm[0] == self.word.letter_form_elm[0] and dividend % divisor < rem:
                rem = dividend % divisor
                quo = dividend // divisor
                value = power_relators[w]
                reduced_word = (value**quo)*(self.word.letter_form_elm[0]**rem)
                return group.dtype(reduced_word)
        return None

    def conjugate_word(self):
        """
        Returns the reduced word with the help of
        conjugate relators.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group

        Example 8.7 Pg. 282 from Handbook
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : 1, x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> word = x2**2*x1
        >>> collector = Collector(pc_relators, word)
        >>> collector.conjugate_word()
        x1*x2**-2

        """
        conjugate_relators = self.relations()[1]
        exp = self.word.array_form[0][1]
        w = (self.word[len(self.word)-1].inverse()*self.word.letter_form_elm[0]*self.word.subword(exp, len(self.word)))
        if w in conjugate_relators:
            conj_word = self.word[len(self.word)-1]*(conjugate_relators[w]**exp)
        else:
            loop = self.word.array_form[len(self.word.array_form)-1][1]
            if loop%2 == 0:
                conj_word = self.word[len(self.word)-1]**loop*self.word.subword(0, len(self.word)-loop)
            else:
                low, high = abs(self.word.array_form[0][1]), abs(self.word.array_form[1][1])
                high = high + low
                conj_word = self.word.eliminate_word(self.word.subword(low, high), self.word.subword(low, high).inverse())
                word = conj_word[len(conj_word)-1]**loop*conj_word.subword(0, len(conj_word)-loop)
        return conj_word

    def _index(self, w):
        """
        Returns the start and ending index of a given
        subword in a word.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : 1, x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> word = x2**2*x1**7
        >>> collector = Collector(pc_relators, word)
        >>> w = x2**2*x1
        >>> collector._index(w)
        (0, 3)
        >>> w = x1**7
        >>> collector._index(w)
        (2, 9)

        """
        low = -1
        high = -1
        for i in range(len(self.word)-len(w)+1):
            if self.word.subword(i, i+len(w)) == w:
                low = i
                high = i+len(w)
                break
        if low == high == -1:
            raise ValueError("Given word is not a subword")
        return low, high

def collected_word(pc_relators, word):
    """
    Examples
    ========
    >>> from sympy.combinatorics.pc_groups import collected_word
    >>> from sympy.combinatorics.pc_groups import Collector
    >>> from sympy.combinatorics.free_groups import free_group
    >>> F, x1, x2 = free_group("x1, x2")
    >>> pc_relators = {x1**2 : 1, x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
    >>> word = x2**2*x1**7
    >>> collector = Collector(pc_relators, word)
    >>> collected_word(pc_relators, word)
    x1*x2**-2

    """
    collector = Collector(pc_relators, word)
    group = collector.word.group
    uncollected_subwords = collector.minimal_uncollected_subwords()
    power_relators, conjugate_relators = collector.relations()

    for w, case in uncollected_subwords.items():
        w = group.dtype(w)
        if case == 1:
            _word = collector.reduce_word()
            collector.word = collector.word.eliminate_word(w, _word)
            if collector.word.array_form[len(collector.word.array_form)-1][1] > 0:
                collector.word = collector.conjugate_word()
        else:
            low, high = collector._index(w)
            collector.word = w
            w = collector.conjugate_word()
            collector.word = word
            collector.word = collector.word.substituted_word(low, high, w)
            collector.word = collector.conjugate_word()
            _word = collector.word.subword(0, collector.word.array_form[0][1])
            collector.word = collector.word.eliminate_word( _word, collector.reduce_word())
    return collector.word
