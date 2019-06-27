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

    """
    References
    ==========

    .. [1] Holt, D., Eick, B., O'Brien, E.
           "Handbook of Computational Group Theory"
           Section 8.1.3
    """

    def __init__(self, pc_relators, word, relative_order):
        self.pc_relators = pc_relators
        self.word = word
        self.group = word.group
        self.relative_order = relative_order
        self.index = {s: i+1 for i, s in enumerate(self.group.symbols)}

    def minimal_uncollected_subwords(self):
        """
        Returns the minimal uncollected subwords.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group

        Example 8.7 Pg. 281 from [1]
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : (), x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> word = x2**2*x1**7
        >>> relative_order = {x1: 2, x2: 3}
        >>> collector = Collector(pc_relators, word, relative_order)
        >>> collector.minimal_uncollected_subwords()
        {((x1, 7),): 2, ((x2, 2), (x1, 1)): 0}

        """
        group = self.group
        index = self.index
        array = self.word.array_form
        re = self.relative_order
        uncollected_subwords = {}

        for i in range(len(array)-1):
            s1, e1 = array[i]
            s2, e2 = array[i+1]
            if e2 > 0 and index[s1] > index[s2]:
                # case-0:  v = x[i]**a*x[i+1], where index[x[i]] > index[x[i+1]]
                uncollected_subwords[((s1, e1), (s2, 1))] = 0

            elif e2 < 0 and index[s1] > index[s2]:
                # case-1: v = x[i]**a*x[i+1]*-1, where index[x[i]] > index[x[i+1]]
                uncollected_subwords[((s1, e1), (s2, -1))] = 1
            s = ((s1, 1), )
            s = group.dtype(s)
            if e1 > re[s]-1:
                # case-2: v = x[i]**a
                uncollected_subwords[((s1, e1), )] = 2

        i = len(array)-1
        s1, e1 = array[i]
        s = ((s1, 1), )
        s = group.dtype(s)
        if e1 > re[s]-1:
            # case-2: v = x[i]**a
            uncollected_subwords[((s1, e1), )] = 2

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
        >>> relative_order = {x1: 2, x2: 3}
        >>> collector = Collector(pc_relators, word, relative_order)
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

    def subword_index(self, w):
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
        >>> relative_order = {x1: 2, x2: 3}
        >>> collector = Collector(pc_relators, word, relative_order)
        >>> w = x2**2*x1
        >>> collector.subword_index(w)
        (0, 3)
        >>> w = x1**7
        >>> collector.subword_index(w)
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
            return -1, -1
        return low, high

    def map_relation(self, w):
        """
        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x0, x1, x2 = free_group("x0, x1, x2")
        >>> pc_relators = {x0**-1*x1*x0: x1**2, x1**-1*x2*x1:x2, x0**-1*x2*x0:x2*x1}
        >>> relative_order = {x0: 2, x1: 2, x2: 3}
        >>> word = x2*x1*x0
        >>> collector = Collector(pc_relators, word, relative_order)
        >>> w = x2*x1
        >>> collector.map_relation(w)
        x2
        >>> w = x1*x0
        >>> collector.map_relation(w)
        x1**2

        """
        group = w.group
        gens = list(sorted(w.contains_generators()))
        key = gens[0]**-1*gens[1]*gens[0]
        return self.pc_relators[key]


def collected_word(pc_relators, word, relative_order):
    """
    Examples
    ========
    >>> from sympy.combinatorics.pc_groups import collected_word
    >>> from sympy.combinatorics.pc_groups import Collector
    >>> from sympy.combinatorics.free_groups import free_group
    >>> F, x1, x2 = free_group("x1, x2")
    >>> pc_relators = {x1**2 : (), x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
    >>> word = x2**2*x1**7
    >>> relative_order = {x1: 2, x2: 3}
    >>> collected_word(pc_relators, word, relative_order)
    x1*x2**-2

    """
    collector = Collector(pc_relators, word, relative_order)
    power_relators, conjugate_relators = collector.relations()
    group = collector.group
    flag = 0
    while True:
        uncollected_subwords = collector.minimal_uncollected_subwords()
        if not uncollected_subwords:
            break
        for w, case in uncollected_subwords.items():
            w = group.dtype(w)
            low, high = collector.subword_index(w)
            if low == -1:
                continue
            if case == 0:
                gens = list(sorted(w.contains_generators()))
                array = w.array_form
                s, e = array[0]
                word_ = collector.map_relation(w)
                word_ = gens[0]*word_**e
                word_ = group.dtype(word_)
                collector.word = collector.word.substituted_word(low, high, word_)

            elif case == 1:
                gens = list(sorted(w.contains_generators()))
                array = w.array_form
                s, e = array[0]
                word_ = collector.map_relation(w)
                word_ = (gens[0]**-1)*word_**e
                word_ = group.dtype(word_)
                collector.word = collector.word.substituted_word(low, high, word_)

            else:
                array = w.array_form
                s, e = array[0]
                s1 = ((s, 1), )
                s = group.dtype(s1)
                re = relative_order[s]
                q = e//re
                r = e-q*re
                key = w[0]**re
                if pc_relators[key]:
                    word_ = ((w[0], r), (pc_relators[key], re))
                    word_ = group.dtype(word_)
                else:
                    if r != 0:
                        word_ = w[0]**r
                    else:
                        word_ = None
                collector.word = collector.word.eliminate_word(w, word_)
    return collector.word
