from sympy.combinatorics.free_groups import free_group, FreeGroupElement
from sympy.combinatorics.fp_groups import FpGroup, FpSubgroup, subgroup_quotient
from sympy import sieve

class PcGroup(FpGroup):
    '''
    Class implementing Polycyclic groups.
    To-do : Docstring completion.
            Add error and exception handling
    The generators will be converted to equivalent forms on
    another free group.
    The implementation for this methods will be similar to
    the coin change problem in DP.
    '''
    def __init__(self, pc_series, power_exponents=None):
        self.pc_series = pc_series
        self.fp_group = self.pc_series[0]
        self.pc_gens = self.compute_pcgs()
        self.rel_orders = self.relative_orders()
        if not power_exponents:
            self.power_exponents = self.relative_orders
        self.power_relations = {}
        self.conjugate_relations = {}
        self.relations = self.power_relations + self.conjugate_relations
        # All the elements of the `pc_gens` will be replaces
        # by the elements of another `FreeGroup`.
        self.F, self.x, self.pc_gens_dict = self._pc_gens_dict()

    def _pc_gens_dict(self):
        '''
        Map the `pc_gens` to a new set of generators of
        another FreeGroup.
        This is useful in the implementation of the collection
        algorithm where a word is converted to a form
        where each element in the word is a distince generator.

        '''
        len_pc_gens = len(self.pc_gens)
        fr_grp_str = ""
        for i in range(0, 3):
            fr_grp_str = fr_grp_str + "x_" + str(i) + ", "
        fr_grp_str = fr_grp_str[:-2]

        # Define a few group on a new set of generators.
        free_group = free_group(fr_grp_str)
        F = free_group[0]
        # New set of generators defined as an array.
        x = [free_group[i] for i in range(1, len(free_group))]
        for i in range(0, len_pc_gens):
            dict[x[i]] = self.pc_gens[i]

        return F, x, dict

    def compute_pcgs(self):
        pcgs = []
        for i in range(0, len(self.pc_series) - 1):
            quotient = subgroup_quotient(self.fp_grp.free_group, self.fp_grp,self.pc_series[i], self.pc_series[i+1])
            pcgs.append(quotient.generators[0])
        self.pc_gens = pcgs

    def relative_orders(self):
        rel_orders = []
        for i in range(0, len(self.pc_series) - 1):
            first_grp = self.pc_series[i]
            sec_grp = self.pc_series[i+1]
            if isinstance(self.pc_series[i], list):
                first_grp = fp_group.subgroup(self.pc_series[i])
            if isinstance(self.pc_series[i+1], list):
                sec_grp = fp_group.subgroup(self.pc_series[i+1])
            rel_orders.append(first_grp.order()/sec_grp.order())
        self.rel_orders = rel_orders

    def is_prime_order(self):
        for elem in self.relative_orders:
            if elem not in sieve:
                return False
        return True

    def is_consistent(self):
        if self.power_exponents.sort() == self.relative_orders.sort():
            return True
        return False

    def _reorder_word(self, w):
        '''
        To find the minimal uncollected subword,
        the word has to be reordered in the same order
        as that of the polycyclic generating sequence.
        '''
        if not isinstance(w, FreeGroupElement):
            raise ValueError("The word must be an instance of FreeGroupElement")
        k = w.letter_form_elm
        k.sort()
        reordered_word = self.fp_group.identity
        for x in k:
            reordered_word *= x
        return reordered_word

    def minimal_uncollected_subwords(self, w):
        # This store the subword along with the
        # type and the exponent of the subword.
        subword_dict = {}
        arr_form = list(w.array_form)
        word_index = 0

        for i in range(0, len(arr_form)):
            if i != (len(arr_form) - 1):
                # Stores the first, last index and the
                # type of the subword
                index_arr = []
                current_index = word_index
                if arr_form[i+1][1] > 0:
                    subword = arr_form[i][0]**arr_form[i][1]*arr_form[i+1]
                    # subword of the type `x**k*y`.
                    index_arr.append(1)
                else:
                    subword = arr_form[i][0]**arr_form[i][1]*arr_form[i+1]**-1
                    # subword of the type `x**k*y**-1`.
                    index_arr.append(2)
                # Append the exponent `k`.
                index_arr.append(arr_form[i][1])
                word_index += abs(arr_form[i][1])
                subword_dict[subword] = index_arr

            # Check for subwords of the type `x**k`.
            if not arr_form[i][1] in self.power_exponents:
                index_arr = []
                subword = arr_form[i][0]**arr_form[i][1]
                # Subword of the type `x**k`.
                index_arr.append(0)
                current_index += abs(arr_form[i][1])
                subword_dict[subword] = index_arr

        return subword_dict

    def _find_index_range(self, elem, w):
        # Returns the index range of a subword in a word
        len_elem = len(elem)
        for i in range(0, len(w)-len):
            if w.subword(i, i+len_elem) == elem:
                return i, i+len_elem
        raise ValueError("Subword not found")

    def _find_relation(w, type):
        '''Find the approproate relation'''
        if type == 0:
            # It's a power relation
            elem = w.letter_form_elm[0]
            elem_index = self.pc_gens.index(elem)
            exp = self.relative_orders(elem_index)
            return elem**exp
        else:
            # It's a conjugate relation
            elem_1 = w.letter_form_elm[0]
            elem_2 = w.letter_form_elm[len(w)-1]
            rel = elem_2**(elem_2.arr_form[0][1])*elem_1*elem_2
            if rel in self.conjugate_relations:
                return rel
        return None

    def collected_word(self, w):
        uncollected_list = minimal_uncollected_subwords(w)
        # Choose a minimal uncollected subword.
        for elem in uncollected_list:
            exp = uncollected_list[elem][1]
            l_index, h_index = self._find_index_range(elem, w)
            # if the subword is of the form `x**k`.
            if uncollected_list[elem][0] == 0:
                relation = self._find_relation(elem, 0)
                rel_val = self.power_relations[relation]
                # the exponent of the relation
                rel_base = rel_val.array_form[0][0]
                rel_exp = rel_val.array_form[0][1]
                # the subword has to be replaced by x**r*(R[r]**q)
                # where R[r] is the appropriate value of the realtion
                q = exp/rel_exp
                r = exp % rel_exp
                new_subword = rel_base**r*rel_val**q
            # the subword is of the form x**k*y or x**k*y**-1
            else:
                relation = self._find_relation(elem, 1)
                rel_val = self.conjugate_relations[relation]
                word_arr = elem.letter_form_elm
                len_word_arr = len(word_arr)
                new_subword = word_arr[len_word_arr-1]*rel_val**exp
            w.substitute_word(l_index, h_index, new_subword)
        return w

    ##########
    #  To-Do #
    ##########

    # 1. Documentation
    # 2. Following methods.
    # 3. Convert the words and relations using DP.

    def compute_presentation(self):
        '''
        Compute the presentation of the PcGroup
        '''
        return
