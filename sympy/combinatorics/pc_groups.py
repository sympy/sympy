from sympy.core import Basic
from sympy import sieve
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
        for order in self.relative_order():
            if order not in sieve:
                return False
        return True

    def length(self):
        return len(self.pcgs)

    def exponent_vector(self, element, group):
        G = PermutationGroup()
        for g in self.pcgs:
            G = PermutationGroup(G.generators + [g])
        gens = G.generator_product(element, original = True)

        perm_to_free = {}
        for sym, g in zip(group.generators, self.pcgs):
            perm_to_free[g] = sym
        w = group.identity
        for g in gens:
            w = w*perm_to_free[g]

        pc_presentation = self.pc_presentation(group)
        collector = Collector(pc_presentation, self.relative_order(), group)
        word = collector.collected_word(w)

        index = {s: i for i, s in enumerate(group.symbols)}
        exp_vector = [0]*len(group)
        word = word.array_form
        for t in word:
            exp_vector[index[t[0]]] = t[1]
        return exp_vector

    def power_relations(self, group):
        """
        Return a list of power relations.
        If `x` is a free group element and `re` represents
        the relative order of `x` then power relation of a
        single generator `x` is `x**re`.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import PolycyclicGroup
        >>> from sympy.combinatorics.named_groups import SymmetricGroup
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x0, x1 = free_group("x0, x1")
        >>> G = SymmetricGroup(4)
        >>> pc_group = G.polycyclic_group()
        >>> group = F
        >>> pc_group.power_relations(group)
        [x1**3, x0**2]

        """
        power_relators = []
        for i, s in enumerate(group.symbols):
            gen = ((s, self.relative_order()[i]), )
            gen = group.dtype(gen)
            power_relators.append(gen)
        power_relators.reverse()
        return power_relators

    def conjugate_relations(self, group):
        """
        Return a list of conjugate relations.
        Conjugate relations of a polycyclic group are of the
        form `x[i]**-1*x[i+1]*x[i]`.

        Examples
        ========
        >>> from sympy.combinatorics.pc_groups import PolycyclicGroup
        >>> from sympy.combinatorics.named_groups import SymmetricGroup
        >>> from sympy.combinatorics.free_groups import free_group
        >>> F, x0, x1, x2, x3 = free_group("x0, x1, x2, x3")
        >>> G = SymmetricGroup(4)
        >>> pc_group = G.polycyclic_group()
        >>> group = F
        >>> pc_group.conjugate_relations(group)
        [x2**-1*x3*x2, x1**-1*x3*x1, x1**-1*x2*x1, x0**-1*x3*x0,
         x0**-1*x2*x0, x0**-1*x1*x0]

        """
        conjugate_relations = []
        g = group.symbols
        for i in range(len(g)):
            for j in range(i+1, len(g)):
                x1 = g[i]
                x2 = g[j]
                relation = ((x1, -1), (x2, 1), (x1, 1))
                conjugate_relations.append(group.dtype(relation))
        conjugate_relations.reverse()
        return conjugate_relations

    def pc_presentation(self, group):
        pc_relators = {}
        sym = list(group.generators)
        # start from the bottom that is, `x2->x1->x0`
        sym.reverse()

        pc_sequence = self.pcgs
        pc_sequence.reverse()

        # store the mapping of polycyclic sequence with the
        # free group elements and vice-versa
        perm_to_free = {}
        free_to_perm = {}

        for gen, s in zip(pc_sequence, sym):
            # since `s**-1` can also be produced by the generator_product
            perm_to_free[gen**-1] = s**-1
            perm_to_free[gen] = s
            free_to_perm[s**-1] = gen**-1
            free_to_perm[s] = gen

        # the LHS of the power and conjugate relations
        power_relators = self.power_relations(group)
        conjugate_relators = self.conjugate_relations(group)

        # form the group with only the last element in polycyclic
        # sequence which will be incremented as moving up
        G = PermutationGroup([pc_sequence[0]])

        # compute the RHS of power relations
        for i, rel in enumerate(power_relators):
            array = rel.array_form
            s, e = array[0]

            l = G.generator_product(free_to_perm[rel[0]]**e, original = True)
            word = group.identity

            # concatenate the generators in `l` to form the word
            # to be collected.
            for gens in l:
                word = word*perm_to_free[gens]
            collector = Collector(pc_relators, self.relative_order(), group)

            word = collector.collected_word(word)
            pc_relators[rel] = word if word else ()

            if i < len(pc_sequence) and pc_sequence[i] not in G:
                G = PermutationGroup(G.generators + [pc_sequence[i]])

        G = PermutationGroup([pc_sequence[0]])
        index = {s: i for i, s in enumerate(group.generators)}
        index1 = {i: s for i, s in enumerate(group.generators)}
        pc_sequence.reverse()

        # compute the RHS of conjugate relations
        for rel in conjugate_relators:
            gens = sorted(rel.contains_generators())
            if free_to_perm[index1[index[gens[0]]+1]] not in G:
                G = PermutationGroup(G.generators + [free_to_perm[index1[index[gens[0]]+1]]])

            # map free group generators to the pc_sequence elements
            # and compute the conjugate relation `x[i]**-1*x[i+1]*x[i]`
            gens = [free_to_perm[gens[0]], free_to_perm[gens[1]]]
            relation = gens[0]**-1*gens[1]*gens[0]

            l = G.generator_product(relation, original = True)
            l.reverse()
            word = group.identity

            # concatenate the generators in `l` to form the word
            # to be collected.
            for gen in l:
                word = word*perm_to_free[gen]
            collector = Collector(pc_relators, self.relative_order(), group)
            word = collector.collected_word(word)
            pc_relators[rel] = word if word else ()

        return pc_relators


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
        self.index = {s: i for i, s in enumerate(group.symbols)}

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
        >>> group = word.group
        >>> collector = Collector(pc_relators, relative_order, group)
        >>> collector.minimal_uncollected_subword(word)
        [((x2, 2), (x1, 1)), ((x1, 7),)]

        """
        # To handle the case word = <identity>
        if not word:
            return None
        l = []
        group = self.group
        array = word.array_form
        re = self.relative_order
        index = self.index

        for i in range(len(array)-1):
            s1, e1 = array[i]
            s2, e2 = array[i+1]

            if re[index[s1]] and (e1 < 0 or e1 > re[index[s1]]-1):
                l.append(((s1, e1), ))

            if index[s1] > index[s2]:
                e = 1 if e2 > 0 else -1
                l.append(((s1, e1), (s2, e)))

        i = len(array)-1
        s1, e1 = array[i]

        if re[index[s1]] and (e1 < 0 or e1 > re[index[s1]]-1):
            l.append(((s1, e1), ))

        return l

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
        >>> relative_order = [2, None]
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
        gens = sorted(w.contains_generators())
        key = gens[0]**-1*gens[1]*gens[0]
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
        >>> group = word.group
        >>> collector = Collector(pc_relators, relative_order, group)
        >>> collector.collected_word(word)
        x1*x2**-2

        """
        group = self.group
        while True:
            l = self.minimal_uncollected_subword(word)
            if not l:
                break
            l.sort(key = len)
            for w in l:
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

                if not self.relative_order[self.index[s]]:
                    continue
                re = self.relative_order[self.index[s]]
                q = e // re
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
