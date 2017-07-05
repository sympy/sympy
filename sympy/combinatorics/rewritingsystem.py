from __future__ import print_function, division

from sympy import S
from sympy.combinatorics.free_groups import FreeGroupElement

class RewritingSystem(object):
    '''
    A class implementing rewriting systems for `FpGroup`s.

    References
    ==========
    [1] Epstein, D., Holt, D. and Rees, S. (1991).
        The use of Knuth-Bendix methods to solve the word problem in automatic groups.
        Journal of Symbolic Computation, 12(4-5), pp.397-414.

    [2] GAP's Manual on its KBMAG package
        https://www.gap-system.org/Manuals/pkg/kbmag-1.5.3/doc/manual.pdf

    '''
    def __init__(self, group):
        self.group = group
        self.alphabet = group.generators
        self._is_confluent = None

        # these values are taken from [2]
        self.maxeqns = 32767 # max rules
        self.tidyint = 100 # rules before tidying

        # dictionary of reductions
        self.rules = {}
        self._init_rules()

    @property
    def is_confluent(self):
        '''
        Return `True` if the system is confluent

        '''
        if self._is_confluent is None:
            self._is_confluent = self._check_confluence()
        return self._is_confluent

    def _init_rules(self):
        rels = self.group.relators[:]
        identity = self.group.free_group.identity
        for r in rels:
            self.add_rule(r, identity)
        self._remove_redundancies()
        return

    def add_rule(self, w1, w2):
        new_keys = set()
        if len(self.rules) + 1 > self.maxeqns:
            return ValueError("Too many rules were defined.")

        if w1 < w2:
            s1 = w1
            w1 = w2
            w2 = s1

        s1 = w1
        s2 = w2

        # The following is the equivalent of checking
        # s1 for overlaps with the implicit reductions
        # {g*g**-1 -> <identity>} and {g**-1*g -> <identity>}
        # for any generator g without installing the
        # redundant rules that would result from processing
        # the overlaps. See [1], Section 3 for details.

        # overlaps on the right
        while len(s1) - len(s2) > 0:
            g = s1[len(s1)-1]
            s1 = s1.subword(0, len(s1)-1)
            s2 = s2*g**-1
            if len(s1) - len(s2) in [0, 1, 2]:
                if s2 < s1:
                    self.rules[s1] = s2
                    new_keys.add(s1)
                elif s1 < s2:
                    self.rules[s2] = s1
                    new_keys.add(s2)

        # overlaps on the left
        while len(w1) - len(w2) > 1:
            g = w1[0]
            w1 = w1.subword(1, len(w1))
            w2 = g**-1*w2
            if len(w1) - len(w2) in [0, 1, 2]:
                if w2 < w1:
                    self.rules[w1] = w2
                    new_keys.add(w1)
                elif w1 < w2:
                    self.rules[w2] = w1
                    new_keys.add(w2)

        return new_keys

    def _remove_redundancies(self):
        '''
        Reduce left- and right-hand sides of reduction rules
        and remove redundant equations (i.e. those for which
        lhs == rhs)

        '''
        rules = self.rules.copy()
        for r in rules:
            v = self.reduce(r, exclude=r)
            w = self.reduce(rules[r])
            if v != r:
                del self.rules[r]
                if v != w:
                    self.add_rule(v, w)
            else:
                self.rules[v] = w
        return

    def make_confluent(self, check=False):
        '''
        Try to make the system confluent using the Knuth-Bendix
        completion algorithm

        '''
        lhs = list(self.rules.keys())

        def _overlaps(r1, r2):
            len1 = len(r1)
            len2 = len(r2)
            result = []
            j = 0
            while j < len1 + len2 - 1:
                j += 1
                if (r1.subword(len1 - j, len1 + len2 - j, strict=False)
                       == r2.subword(j - len1, j, strict=False)):
                    a = r1.subword(0, len1-j, strict=False)
                    a = a*r2.subword(0, j-len1, strict=False)
                    b = r2.subword(j-len1, j, strict=False)
                    c = r2.subword(j, len2, strict=False)
                    c = c*r1.subword(len1 + len2 - j, len1, strict=False)
                    result.append(a*b*c)
            return result

        def _process_overlap(w, r1, r2):
                s = w.eliminate_word(r1, self.rules[r1])
                t = w.eliminate_word(r2, self.rules[r2])
                if self.reduce(s) != self.reduce(t):
                    new_keys = self.add_rule(t, s)
                    return new_keys
                return

        added = 0
        for i, r1 in enumerate(lhs):
            # this could be lhs[i+1:] to not
            # check each pair twice but lhs
            # is extended in the loop and the new
            # elements have to be checked with the
            # preceding ones. there is probably a better way
            # to handle this
            for r2 in lhs:
                overlaps = _overlaps(r1, r2)
                if not overlaps:
                    continue
                for w in overlaps:
                    new_keys = _process_overlap(w, r1, r2)
                    if new_keys:
                        if check:
                            return False
                        added += 1
                        if added > self.tidyint:
                            # tidy up
                            self._remove_redundancies()
                        lhs.extend([n for n in new_keys if n in self.rules])

        self._is_confluent = True
        if check:
            return True
        else:
            self._remove_redundancies()
        return

    def _check_confluence(self):
        return self.make_confluent(check=True)

    def reduce(self, word, exclude=None):
        '''
        Apply reduction rules to `word` excluding the reduction rule
        for the lhs equal to `exclude`

        '''
        rules = {r: self.rules[r] for r in self.rules if r != exclude}
        # the following is essentially `eliminate_words()` code from the
        # `FreeGroupElement` class, the only difference being the first
        # "if" statement
        again = True
        new = word
        while again:
            again = False
            for r in rules:
                prev = new
                if len(r) == len(rules[r]):
                    new = new.eliminate_word(r, rules[r], _all=True, inverse=False)
                else:
                    new = new.eliminate_word(r, rules[r], _all=True)
                if new != prev:
                    again = True
        return new
