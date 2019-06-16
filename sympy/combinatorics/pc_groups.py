from sympy.core import Basic
from sympy import sieve
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup

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


def minimal_uncollected_subwords(word, relative_order):
    """
    Returns the minimal uncollected subwords.

    Examples
    ========
    >>> from sympy.combinatorics.pc_groups import minimal_uncollected_subwords
    >>> from sympy.combinatorics.free_groups import free_group

    Example 8.7 Pg. 281 from Handbook
    >>> F, x1, x2 = free_group("x1, x2")
    >>> word = x2**2*x1**7
    >>> relative_order = [3, 2]
    >>> minimal_uncollected_subwords(word, relative_order)
    {x1**7: 1, x1*x2**2: 0}

    """
    array = word.array_form
    uncollected_subwords = {}
    for i in range(len(array)-1):
        if array[i+1][1] > 0:
            # case-1:  v = x[i]**a*x[i+1]
            uncollected_subwords[array[i][0]**array[i][1]*array[i+1][0]] = 0
        else:
            # case-2: v = x[i]**a*x[i+1]*-1
            uncollected_subwords[array[i][0]**array[i][1]*array[i+1][0]**-1] = 0

        if all(array[i][1]!=exp for exp in range(relative_order[i])):
            # case-3: v = x[i]**a
            uncollected_subwords[array[i][0]**array[i][1]] = 1

    i = len(array)-1
    if all(array[i][1]!=exp for exp in range(relative_order[i])):
        # case-3: v = x[i]**a
        uncollected_subwords[array[i][0]**array[i][1]] = 1

    return uncollected_subwords

def _relations(relators):
    """
    Separates the given relators of pc presentation in power and
    conjugate relations.

    Examples
    ========
    >>> from sympy.combinatorics.pc_groups import _relations
    >>> from sympy.combinatorics.free_groups import free_group
    >>> F, x1, x2 = free_group("x1, x2")
    >>> pc_relators = {x1**2 : 1, x1*x2*x1**-1 : x2**-1, x1**-1*x2*x1 : x2**-1}
    >>> power_rel, conj_rel = _relations(pc_relators)
    >>> power_rel
    {x1**2: 1}
    >>> conj_rel
    {x1**-1*x2*x1: x2**-1, x1*x2*x1**-1: x2**-1}

    """
    power_relators = {}
    conjugate_relators = {}
    for key, value in relators.items():
        if len(key.array_form) == 1:
            power_relators[key] = value
        else:
            conjugate_relators[key] = value
    return power_relators, conjugate_relators

def _reduce_word(word, power_relators):
    """
    Reduce the given word with the help of power
    relators only.

    Examples
    ========
    >>> from sympy.combinatorics.pc_groups import _reduce_word
    >>> from sympy.combinatorics.free_groups import free_group
    >>> F, x1, x2 = free_group("x1, x2")
    >>> power_relators = {x1**2: 1}
    >>> word = x1**7
    >>> _reduce_word(word, power_relators)
    x1

    """
    group = word.group
    dividend = word.array_form[0][1]
    rem = dividend
    for w, v in power_relators.items():
        divisor = w.array_form[0][1]
        if w.letter_form_elm[0] == word.letter_form_elm[0] and dividend % divisor < rem:
            rem = dividend % divisor
            quo = dividend // divisor
            value = power_relators[w]
    word = (value**quo)*(word.letter_form_elm[0]**rem)
    return group.dtype(word)

def _conjugate_word(word, conjugate_relators):
    """
    Returns the reduced word with the help of
    conjugate relators.

    Examples
    ========
    >>> from sympy.combinatorics.pc_groups import _conjugate_word
    >>> from sympy.combinatorics.free_groups import free_group

    Example 8.7 Pg. 282 from Handbook
    >>> F, x1, x2 = free_group("x1, x2")
    >>> word = x2**2*x1
    >>> conjugate_relators = {x1*x2*x1**-1: x2**-1, x1**-1*x2*x1: x2**-1}
    >>> _conjugate_word(word, conjugate_relators)
    x1*x2**-2

    """
    exp = word.array_form[0][1]
    w = (word[len(word)-1].inverse()*word.letter_form_elm[0]*word.subword(exp, len(word)))
    if w in conjugate_relators:
        word = word[len(word)-1]*(conjugate_relators[w]**exp)
    else:
        loop = word.array_form[len(word.array_form)-1][1]
        if loop%2 == 0:
            word = word[len(word)-1]**loop*word.subword(0, len(word)-loop)
        else:
            low, high = abs(word.array_form[0][1]), abs(word.array_form[1][1])
            high = high + low
            print(word.subword(low, high))
            word = word.eliminate_word(word.subword(low, high), word.subword(low, high).inverse())
            word = word[len(word)-1]**loop*word.subword(0, len(word)-loop)
    return word

def _index(word, w):
    """
    Returns the start and ending index of a given
    subword in a word.

    Examples
    ========
    >>> from sympy.combinatorics.pc_groups import _index
    >>> from sympy.combinatorics.free_groups import free_group
    >>> F, x1, x2 = free_group("x1, x2")
    >>> word = x2**2*x1**7
    >>> w = x2**2*x1
    >>> _index(word, w)
    (0, 3)
    >>> w = x1**7
    >>> _index(word, w)
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
        raise ValueError("Given word is not a subword")
    return low, high

def collected_word(pc_relators, word, relative_order):
    """

    """
    group = word.group
    uncollected_subwords = minimal_uncollected_subwords(word, relative_order)
    power_relators, conjugate_relators = _relations(pc_relators)
    for w, case in uncollected_subwords.items():
        if case == 1:
            _word = _reduce_word(w, power_relators)
            word = word.eliminate_word(w, _word)
            word = _conjugate_word(word, conjugate_relators)

        else:
            low, high = _index(word, w)
            w = _conjugate_word(w, conjugate_relators)
            word = word.substituted_word(low, high, w)
            word = _conjugate_word(word, conjugate_relators)
            _word = word.subword(0, word.array_form[0][1])
            word = word.eliminate_word( _word, _reduce_word(_word, power_relators))
    return word
