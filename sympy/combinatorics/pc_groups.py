from sympy.core import Basic
from sympy import isprime
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.printing.defaults import DefaultPrinting

class PolycyclicGroup(DefaultPrinting):

    is_group = True
    is_solvable = True

    def __init__(self, pc_sequence, pc_series):
        self.pcgs = pc_sequence
        self.pc_series = pc_series

    def relative_order(self):
        rel_orders = []
        for i in range(len(self.pc_series)-1):
            G = self.pc_series[i]
            H = self.pc_series[i+1]
            rel_orders.append(G.order() // H.order())
        return rel_orders

    def is_prime_order(self):
        return all(isprime(order) for order in self.relative_order())

    def length(self):
        return len(self.pcgs)

    def exponent_vector(self, element, free_group):
        """
        Return the exponent vector of length equal to the
        length of polycyclic generating sequence.

        For a given generator/element `g` of the polycyclic group,
        it can be represented as `g = x{1}**e{1}....x{n}**e{n}`,
        where `x{i}` represents polycyclic generators and `n` is
        the number of generators in the free_group equal to the length
        of pcgs.

        Examples
        ========
        >>> from sympy.combinatorics.named_groups import SymmetricGroup
        >>> from sympy.combinatorics.free_groups import free_group
        >>> from sympy.combinatorics.permutations import Permutation
        >>> G = SymmetricGroup(4)
        >>> PcGroup = G.polycyclic_group()
        >>> pcgs = PcGroup.pcgs
        >>> free_group, x0, x1, x2, x3 = free_group("x0, x1, x2, x3")
        >>> PcGroup.exponent_vector(G[0], free_group)
        [1, 0, 0, 0]
        >>> exp = PcGroup.exponent_vector(G[1], free_group)
        >>> g = Permutation()
        >>> for i in range(len(exp)):
        ...     g = g*pcgs[i] if exp[i] else g
        >>> assert g == G[1]

        References
        ==========

        .. [1] Holt, D., Eick, B., O'Brien, E.
               "Handbook of Computational Group Theory"
                Section 8.1.1, Definition 8.4

        """
        G = PermutationGroup()
        for g in self.pcgs:
            G = PermutationGroup([g] + G.generators)
        gens = G.generator_product(element, original = True)
        gens.reverse()

        perm_to_free = {}
        for sym, g in zip(free_group.generators, self.pcgs):
            perm_to_free[g**-1] = sym**-1
            perm_to_free[g] = sym
        w = free_group.identity
        for g in gens:
            w = w*perm_to_free[g]

        pc_presentation = self.pc_presentation(free_group)
        collector = Collector(pc_presentation, self.relative_order(), free_group)
        word = collector.collected_word(w)

        index = {s: i for i, s in enumerate(free_group.symbols)}
        exp_vector = [0]*len(free_group)
        word = word.array_form
        for t in word:
            exp_vector[index[t[0]]] = t[1]
        return exp_vector

    def depth(self, element, free_group):
        """
        Return the depth of a given element.

        The depth of a given element `g` is defined by
        `dep{g} = i if e{1} = e{2} = ... = e{i-1} = 0`
        and `e{i} != 0`, where `e` represents the exponent-vector.

        Examples
        ========
        >>> from sympy.combinatorics.named_groups import SymmetricGroup
        >>> from sympy.combinatorics.free_groups import free_group
        >>> free_group, x, y = free_group("x, y")
        >>> G = SymmetricGroup(3)
        >>> PcGroup = G.polycyclic_group()
        >>> PcGroup.depth(G[0], free_group)
        2
        >>> PcGroup.depth(G[1], free_group)
        1

        References
        ==========

        .. [1] Holt, D., Eick, B., O'Brien, E.
               "Handbook of Computational Group Theory"
                Section 8.1.1, Definition 8.5

        """
        exp_vector = self.exponent_vector(element, free_group)
        return next((i+1 for i, x in enumerate(exp_vector) if x), len(self.pcgs)+1)

    def leading_exponent(self, element, free_group):
        """
        Return the leading non-zero exponent.

        The leading exponent for a given element `g` is defined
        by `leading_exponent{g} = e{i}`, if `depth{g} = i`.

        Examples
        ========
        >>> from sympy.combinatorics.named_groups import SymmetricGroup
        >>> from sympy.combinatorics.free_groups import free_group
        >>> free_group, x, y = free_group("x, y")
        >>> G = SymmetricGroup(3)
        >>> PcGroup = G.polycyclic_group()
        >>> PcGroup.leading_exponent(G[1], free_group)
        1

        """
        exp_vector = self.exponent_vector(element, free_group)
        depth = self.depth(element, free_group)
        if depth != len(self.pcgs)+1:
            return exp_vector[depth-1]
        return None

    def pc_presentation(self, free_group):
        """
        Return the polycyclic presentation.

        There are two types of relations used in polycyclic
        presentation.
        i) Power relations of the form `x{i}^re{i} = R{i}{i}`,
        `for 0 <= i < length(pcgs)` where `x` represents polycyclic
        generator and `re` is the corresponding relative order.

        ii) Conjugate relations of the form `x{j}^-1*x{i}*x{j}`,
        `for 0 <= j < i <= length(pcgs)`.

        Examples
        ========
        >>> from sympy.combinatorics.named_groups import SymmetricGroup
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.free_groups import free_group
        >>> S = SymmetricGroup(49).sylow_subgroup(7)
        >>> der = S.derived_series()
        >>> G = der[len(der)-2]
        >>> PcGroup = G.polycyclic_group()
        >>> pcgs = PcGroup.pcgs
        >>> len(pcgs)
        6
        >>> free_group, x0, x1, x2, x3, x4, x5 = free_group("x0, x1, x2, x3, x4, x5")
        >>> pc_resentation = PcGroup.pc_presentation(free_group)
        >>> free_to_perm = {}
        >>> for s, g in zip(free_group.symbols, pcgs):
        ...     free_to_perm[s] = g

        >>> for k, v in pc_resentation.items():
        ...     k_array = k.array_form
        ...     if v != ():
        ...        v_array = v.array_form
        ...     lhs = Permutation()
        ...     for gen in k_array:
        ...         s = gen[0]
        ...         e = gen[1]
        ...         lhs = lhs*free_to_perm[s]**e
        ...     if v == ():
        ...         assert lhs.is_identity
        ...         continue
        ...     rhs = Permutation()
        ...     for gen in v_array:
        ...         s = gen[0]
        ...         e = gen[1]
        ...         rhs = rhs*free_to_perm[s]**e
        ...     assert lhs == rhs

        """
        rel_order = self.relative_order()
        pc_relators = {}
        perm_to_free = {}
        pcgs = self.pcgs

        for gen, s in zip(pcgs, free_group.generators):
            perm_to_free[gen**-1] = s**-1
            perm_to_free[gen] = s

        collector = Collector(pc_relators, rel_order, free_group)
        pcgs.reverse()
        series = self.pc_series
        series.reverse()
        collected_gens = []

        for i, gen in enumerate(pcgs):
            re = rel_order[len(rel_order)-i-1]
            relation = perm_to_free[gen]**re
            G = series[i] if i != 0 else series[i+1]

            l = G.generator_product(gen**re, original = True)
            l.reverse()

            word = free_group.identity
            for g in l:
                word = word*perm_to_free[g]

            word = collector.collected_word(word)
            pc_relators[relation] = word if word else ()
            collector.pc_relators = pc_relators

            collected_gens.append(gen)
            if len(collected_gens) > 1:
                conj = collected_gens[len(collected_gens)-1]
                conjugator = perm_to_free[conj]

                for j in range(len(collected_gens)-1):
                    conjugated = perm_to_free[collected_gens[j]]

                    relation = conjugator**-1*conjugated*conjugator
                    gens = conj**-1*collected_gens[j]*conj

                    l = G.generator_product(gens, original = True)
                    l.reverse()
                    word = free_group.identity
                    for g in l:
                        word = word*perm_to_free[g]

                    word = collector.collected_word(word)
                    pc_relators[relation] = word if word else ()
                    collector.pc_relators = pc_relators

        series.reverse()
        pcgs.reverse()
        return pc_relators


class Collector(DefaultPrinting):

    """
    References
    ==========

    .. [1] Holt, D., Eick, B., O'Brien, E.
           "Handbook of Computational Group Theory"
           Section 8.1.3
    """

    def __init__(self, pc_relators, relative_order, free_group):
        self.pc_relators = pc_relators
        self.relative_order = relative_order
        self.free_group = free_group
        self.index = {s: i for i, s in enumerate(free_group.symbols)}

    def minimal_uncollected_subword(self, word):
        """
        Returns the minimal uncollected subwords.

        A word `v` defined on generators in `X` is a minimal
        uncollected subword of the word `w` if `v` is a subword
        of `w` and it has one of the following form

        i) `v = x[i+1]**a_j*x[i]`

        ii) `v = x[i+1]**a_j*x[i]**-1`

        iii) `v = x[i]**a_j` for relative_order of `x[i] != infinity`
        and `a_j` is not in `{1, ..., s-1}`. Where, s is the power
        exponent of the corresponding generator.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group

        Example 8.7 Pg. 281 from [1]
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : (), x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> relative_order = [2, None]
        >>> word = x2**2*x1**7
        >>> free_group = word.group
        >>> collector = Collector(pc_relators, relative_order, free_group)
        >>> collector.minimal_uncollected_subword(word)
        ((x1, 7),)

        """
        # To handle the case word = <identity>
        if not word:
            return None

        array = word.array_form
        re = self.relative_order
        index = self.index

        for i in range(len(array)):
            s1, e1 = array[i]

            if re[index[s1]] and (e1 < 0 or e1 > re[index[s1]]-1):
                return ((s1, e1), )

        for i in range(len(array)-1):
            s1, e1 = array[i]
            s2, e2 = array[i+1]

            if index[s1] > index[s2]:
                e = 1 if e2 > 0 else -1
                return ((s1, e1), (s2, e))

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
        >>> relative_order = [2, None]
        >>> word = x2**2*x1**7
        >>> free_group = word.group
        >>> collector = Collector(pc_relators, relative_order, free_group)
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
        >>> relative_order = [2, None]
        >>> word = x2**2*x1**7
        >>> free_group = word.group
        >>> collector = Collector(pc_relators, relative_order, free_group)
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
        Return a conjugate relation.
        Given a word formed by two free group elements, the
        corresponding conjugate relation with those free
        group elements is formed and mapped with the collected
        word in the polycyclic presentation.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x0, x1, x2 = free_group("x0, x1, x2")
        >>> pc_relators = {x0**-1*x1*x0: x1**2, x1**-1*x2*x1: x2, x0**-1*x2*x0: x2*x1}
        >>> relative_order = [2, 2, 3]
        >>> word = x2**2*x1**7
        >>> free_group = word.group
        >>> collector = Collector(pc_relators, relative_order, free_group)
        >>> w = x2*x1
        >>> collector.map_relation(w)
        x2
        >>> w = x1*x0
        >>> collector.map_relation(w)
        x1**2

        """
        array = w.array_form
        s1 = array[0][0]
        s2 = array[1][0]
        key = ((s2, -1), (s1, 1), (s2, 1))
        key = self.free_group.dtype(key)
        return self.pc_relators[key]


    def collected_word(self, word):
        """
        Return the collected form of a word.

        A word `w` is called collected, if `w = x{i_1}**a_1*...*x{i_r}**a_r`
        with `i_1 < i_2< ... < i_r` and `a_j` is in `{1, ..., s_j-1}`
        if `s_j != infinity`.
        Otherwise w is uncollected.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import Collector
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x1, x2 = free_group("x1, x2")
        >>> pc_relators = {x1**2 : (), x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
        >>> relative_order = [2, None]
        >>> word = x2**2*x1**7
        >>> free_group = word.group
        >>> collector = Collector(pc_relators, relative_order, free_group)
        >>> collector.collected_word(word)
        x1*x2**-2

        """
        free_group = self.free_group
        while True:
            w = self.minimal_uncollected_subword(word)
            if not w:
                break

            low, high = self.subword_index(word, free_group.dtype(w))
            if low == -1:
                continue
            s1, e1 = w[0]
            if len(w) == 2 and w[1][1] > 0:
                s2, e2 = w[1]
                s2 = ((s2, 1), )
                s2 = free_group.dtype(s2)
                word_ = self.map_relation(free_group.dtype(w))
                word_ = s2*word_**e1
                word_ = free_group.dtype(word_)
                word = word.substituted_word(low, high, word_)

            if not self.relative_order[self.index[s1]]:
                continue
            re = self.relative_order[self.index[s1]]
            q = e1 // re
            r = e1-q*re
            if r < 0 or r > re-1:
                continue

            if len(w) == 1 and (e1 < 0 or e1 > re-1):
                key = ((w[0][0], re), )
                key = free_group.dtype(key)
                if self.pc_relators[key]:
                    word_ = ((w[0][0], r), (self.pc_relators[key], q))
                    word_ = free_group.dtype(word_)
                else:
                    if r != 0:
                        word_ = ((w[0][0], r), )
                        word_ = free_group.dtype(word_)
                    else:
                        word_ = None
                word = word.eliminate_word(free_group.dtype(w), word_)

            elif len(w) == 2 and w[1][1] < 0:
                s2, e2 = w[1]
                s2 = ((s2, 1), )
                s2 = free_group.dtype(s2)
                word_ = self.map_relation(free_group.dtype(w))
                word_ = s2**-1*word_**e1
                word_ = free_group.dtype(word_)
                word = word.substituted_word(low, high, word_)

        return word
