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

    def __init__(self, pc_relators, relative_order, group):
        self.pc_relators = pc_relators
        self.relative_order = relative_order
        self.group = group

    def minimal_uncollected_subword(self, word):
        """
        Returns the minimal uncollected subwords.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group

        Example 8.7 Pg. 281 from [1]
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : (), x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> relative_order = {x1: 2}
        >>> word = x2**2*x1**7
        >>> group = word.group
        >>> collector = Collector(pc_relators, relative_order, group)
        >>> collector.minimal_uncollected_subword(word)
        ((x2, 2), (x1, 1))

        """
        # To handle the case word = <identity>
        if not word:
            return None

        group = self.group
        index = {s: i+1 for i, s in enumerate(group.symbols)}
        array = word.array_form
        re = self.relative_order

        for i in range(len(array)-1):
            s1, e1 = array[i]
            s2, e2 = array[i+1]
            s = ((s1, 1), )
            s = group.dtype(s)
            if s in re and (e1 < 0 or e1 > re[s]-1):
                return ((s1, e1), )

            if index[s1] > index[s2]:
                e = 1 if e2 > 0 else -1
                return ((s1, e1), (s2, e))

        i = len(array)-1
        s1, e1 = array[i]
        s = ((s1, 1), )
        s = group.dtype(s)
        if s in re and (e1 < 0 or e1 > re[s]-1):
            return ((s1, e1), )
        return None

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
        >>> relative_order = {x1: 2}
        >>> word = x2**2*x1**7
        >>> group = word.group
        >>> collector = Collector(pc_relators, relative_order, group)
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

    def subword_index(self, word, w):
        """
        Returns the start and ending index of a given
        subword in a word.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : 1, x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> relative_order = {x1: 2}
        >>> word = x2**2*x1**7
        >>> group = word.group
        >>> collector = Collector(pc_relators, relative_order, group)
        >>> w = x2**2*x1
        >>> collector.subword_index(word, w)
        (0, 3)
        >>> w = x1**7
        >>> collector.subword_index(word, w)
        (2, 9)

        """
        low = -1
        high = -1
        for i in range(len(word)-len(w)+1):
            if word.subword(i, i+len(w)) == w:
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
        >>> word = x2**2*x1**7
        >>> group = word.group
        >>> collector = Collector(pc_relators, relative_order, group)
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


    def collected_word(self, word):
        """
        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : (), x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> relative_order = {x1: 2}
        >>> word = x2**2*x1**7
        >>> group = word.group
        >>> collector = Collector(pc_relators, relative_order, group)
        >>> collector.collected_word(word)
        x1*x2**-2

        """
        group = self.group
        while True:
            w = self.minimal_uncollected_subword(word)
            if not w:
                break
            low, high = self.subword_index(word, group.dtype(w))
            if low == -1:
                continue
            s, e = w[0]
            if len(w) == 2 and w[1][1] > 0:
                gens = list(sorted(group.dtype(w).contains_generators()))
                word_ = self.map_relation(group.dtype(w))
                word_ = gens[0]*word_**e
                word_ = group.dtype(word_)
                word = word.substituted_word(low, high, word_)

            elif len(w) == 2 and w[1][1] < 0:
                gens = list(sorted(group.dtype(w).contains_generators()))
                word_ = self.map_relation(group.dtype(w))
                word_ = gens[0]**-1*word_**e
                word_ = group.dtype(word_)
                word = word.substituted_word(low, high, word_)

            s1 = ((s, 1), )
            s = group.dtype(s1)
            if not s in self.relative_order:
                continue
            re = self.relative_order[s]
            q = e//re
            r = e-q*re
            if r < 0 or r > re-1:
                continue

            if len(w) == 1 and (e < 0 or e > re-1):
                key = ((w[0][0], re), )
                key = group.dtype(key)
                if self.pc_relators[key]:
                    word_ = ((w[0][0], r), (self.pc_relators[key], q))
                    word_ = group.dtype(word_)
                else:
                    if r != 0:
                        word_ = ((w[0][0], r), )
                        word_ = group.dtype(word_)
                    else:
                        word_ = None
                word = word.eliminate_word(group.dtype(w), word_)
        return word
