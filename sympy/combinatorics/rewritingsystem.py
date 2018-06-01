from __future__ import print_function, division

from sympy import S
from sympy.combinatorics.free_groups import FreeGroupElement
from sympy.combinatorics.rewritingsystem_fsm import State, StateMachine

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
        from collections import deque
        self.group = group
        self.alphabet = group.generators
        self._is_confluent = None

        # these values are taken from [2]
        self.maxeqns = 32767 # max rules
        self.tidyint = 100 # rules before tidying

        # _max_exceeded is True if maxeqns is exceeded
        # at any point
        self._max_exceeded = False

        # dictionary of reductions
        self.rules = {}
        self.rules_cache = deque([], 50)
        self._init_rules()

        # inverse rules - used in automaton
        self.inverse_rules = {}

        # Automaton variables
        self.reduction_automaton = self.construct_automaton()

    def set_max(self, n):
        '''
        Set the maximum number of rules that can be defined

        '''
        if self._max_exceeded and n > self.maxeqns:
            self._max_exceeded = False
        self.maxeqns = n
        return

    @property
    def is_confluent(self):
        '''
        Return `True` if the system is confluent

        '''
        if self._is_confluent is None:
            self._is_confluent = self._check_confluence()
        return self._is_confluent

    def _init_rules(self):
        identity = self.group.free_group.identity
        for r in self.group.relators:
            self.add_rule(r, identity)
        self._remove_redundancies()
        return

    def _add_rule(self, r1, r2):
        '''
        Add the rule r1 -> r2 with no checking or further
        deductions

        '''
        if len(self.rules) + 1 > self.maxeqns:
            self._is_confluent = self._check_confluence()
            self._max_exceeded = True
            raise RuntimeError("Too many rules were defined.")
        self.rules[r1] = r2

    def add_rule(self, w1, w2, check=False):
        new_keys = set()

        if w1 == w2:
            return new_keys

        if w1 < w2:
            w1, w2 = w2, w1

        if (w1, w2) in self.rules_cache:
            return new_keys
        self.rules_cache.append((w1, w2))

        s1, s2 = w1, w2

        # The following is the equivalent of checking
        # s1 for overlaps with the implicit reductions
        # {g*g**-1 -> <identity>} and {g**-1*g -> <identity>}
        # for any generator g without installing the
        # redundant rules that would result from processing
        # the overlaps. See [1], Section 3 for details.

        if len(s1) - len(s2) < 3:
            if s1 not in self.rules:
                new_keys.add(s1)
                if not check:
                    self._add_rule(s1, s2)
            if s2**-1 > s1**-1 and s2**-1 not in self.rules:
                new_keys.add(s2**-1)
                if not check:
                    self._add_rule(s2**-1, s1**-1)

        # overlaps on the right
        while len(s1) - len(s2) > -1:
            g = s1[len(s1)-1]
            s1 = s1.subword(0, len(s1)-1)
            s2 = s2*g**-1
            if len(s1) - len(s2) < 0:
                if s2 not in self.rules:
                    if not check:
                        self._add_rule(s2, s1)
                    new_keys.add(s2)
            elif len(s1) - len(s2) < 3:
                new = self.add_rule(s1, s2, check)
                new_keys.update(new)

        # overlaps on the left
        while len(w1) - len(w2) > -1:
            g = w1[0]
            w1 = w1.subword(1, len(w1))
            w2 = g**-1*w2
            if len(w1) - len(w2) < 0:
                if w2 not in self.rules:
                    if not check:
                        self._add_rule(w2, w1)
                    new_keys.add(w2)
            elif len(w1) - len(w2) < 3:
                new = self.add_rule(w1, w2, check)
                new_keys.update(new)

        return new_keys

    def _remove_redundancies(self, changes=False):
        '''
        Reduce left- and right-hand sides of reduction rules
        and remove redundant equations (i.e. those for which
        lhs == rhs). If `changes` is `True`, return a set
        containing the removed keys and a set containing the
        added keys

        '''
        removed = set()
        added = set()
        rules = self.rules.copy()
        for r in rules:
            v = self.reduce(r, exclude=r)
            w = self.reduce(rules[r])
            if v != r:
                del self.rules[r]
                removed.add(r)
                if v > w:
                    added.add(v)
                    self.rules[v] = w
                elif v < w:
                    added.add(w)
                    self.rules[w] = v
            else:
                self.rules[v] = w
        if changes:
            return removed, added
        return

    def make_confluent(self, check=False):
        '''
        Try to make the system confluent using the Knuth-Bendix
        completion algorithm

        '''
        if self._max_exceeded:
            return self._is_confluent
        lhs = list(self.rules.keys())

        def _overlaps(r1, r2):
            len1 = len(r1)
            len2 = len(r2)
            result = []
            for j in range(1, len1 + len2):
                if (r1.subword(len1 - j, len1 + len2 - j, strict=False)
                       == r2.subword(j - len1, j, strict=False)):
                    a = r1.subword(0, len1-j, strict=False)
                    a = a*r2.subword(0, j-len1, strict=False)
                    b = r2.subword(j-len1, j, strict=False)
                    c = r2.subword(j, len2, strict=False)
                    c = c*r1.subword(len1 + len2 - j, len1, strict=False)
                    result.append(a*b*c)
            return result

        def _process_overlap(w, r1, r2, check):
                s = w.eliminate_word(r1, self.rules[r1])
                s = self.reduce(s)
                t = w.eliminate_word(r2, self.rules[r2])
                t = self.reduce(t)
                if s != t:
                    if check:
                        # system not confluent
                        return [0]
                    try:
                        new_keys = self.add_rule(t, s, check)
                        return new_keys
                    except RuntimeError:
                        return False
                return

        added = 0
        i = 0
        while i < len(lhs):
            r1 = lhs[i]
            i += 1
            # j could be i+1 to not
            # check each pair twice but lhs
            # is extended in the loop and the new
            # elements have to be checked with the
            # preceding ones. there is probably a better way
            # to handle this
            j = 0
            while j < len(lhs):
                r2 = lhs[j]
                j += 1
                if r1 == r2:
                    continue
                overlaps = _overlaps(r1, r2)
                overlaps.extend(_overlaps(r1**-1, r2))
                if not overlaps:
                    continue
                for w in overlaps:
                    new_keys = _process_overlap(w, r1, r2, check)
                    if new_keys:
                        if check:
                            return False
                        lhs.extend(new_keys)
                        added += len(new_keys)
                    elif new_keys == False:
                        # too many rules were added so the process
                        # couldn't complete
                        return self._is_confluent

                if added > self.tidyint and not check:
                    # tidy up
                    r, a = self._remove_redundancies(changes=True)
                    added = 0
                    if r:
                        # reset i since some elements were removed
                        i = min([lhs.index(s) for s in r])
                    lhs = [l for l in lhs if l not in r]
                    lhs.extend(a)
                    if r1 in r:
                        # r1 was removed as redundant
                        break

        self._is_confluent = True
        if not check:
            self._remove_redundancies()
        return True

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
                if rules[r]**-1 > r**-1:
                    new = new.eliminate_word(r, rules[r], _all=True, inverse=False)
                else:
                    new = new.eliminate_word(r, rules[r], _all=True)
                if new != prev:
                    again = True
        return new

    def compute_inverse_rules(self):
        '''
        This computes the inverse rules and stores them in the inverse_rules dictionary.
        The inverse rules are used in the automaton for word reduction.

        '''
        for r in self.rules:
            rule_key_inverse = r**-1
            rule_value_inverse = (self.rules[r])**-1
            if (rule_value_inverse < rule_key_inverse):
                self.inverse_rules[rule_key_inverse] = rule_value_inverse
            else:
                self.inverse_rules[rule_value_inverse] = rule_key_inverse

    def construct_automaton(self):
        '''
        Construct the automaton based on the set of reduction rules of the system.

        Automata Design:
        1. The accept states of the automaton are the proper prefixes of the left hand side of the rules.
        2. The complete left hand side of the rules are the dead states of the automaton

        '''
        automaton_alphabet = []
        left_hand_rules = []
        proper_prefixes = {}

        # compute all_rules when the automaton is constucted.
        all_rules = self.rules
        self.compute_inverse_rules()
        all_rules.update(self.inverse_rules)

        generators = list(self.alphabet)
        generators += [gen**-1 for gen in generators]

        # Contains the alphabets that will be used for state transitions.
        automaton_alphabet = generators

        # Store the complete left hand side of the rules - dead states.
        left_hand_rules = list(all_rules)

        # Compute the proper prefixes for every rule.
        for r in all_rules:
            proper_prefixes[r] = []
            letter_word_array = [s for s in r.letter_form_elm]
            for i in range (1, len(letter_word_array)):
                letter_word_array[i] = letter_word_array[i-1]*letter_word_array[i]
            proper_prefixes[r] = letter_word_array

        # Create the states in the automaton.
        # The left-hand side of the rules are the dead states.
        # The proper left-hand side of the rules are the accept states.

        # Create a finite state machine as an instance of the StateMachine object
        fsm = StateMachine('fsm')
        fsm.add_state('start', is_start=True)

        # Add dead states.
        for rule in left_hand_rules:
            if not rule in fsm.states:
                fsm.add_state(rule, is_dead=True)

        # Add accept states.
        for r in all_rules:
            prop_prefix = proper_prefixes[r]
            for elem in prop_prefix:
                if not elem in fsm.states:
                    fsm.add_state(elem, is_accept=True)

        # Add transitions for every state
        for state in fsm.states:
            current_state_name = state
            current_state_type = fsm.states[state].state_type
            if current_state_type == "start":
                for letter in automaton_alphabet:
                    if letter in fsm.states:
                        fsm.states[state].add_transition(letter, letter)
                    else:
                        fsm.states[state].add_transition(letter, current_state_name)
            elif current_state_type == "accept":
                for letter in automaton_alphabet:
                    next = current_state_name*letter
                    len_next_word = len(next)
                    while True:
                        if len(next) <= 1:
                            if next in fsm.states:
                                fsm.states[state].add_transition(letter, next)
                            else:
                                fsm.states[state].add_transition(letter, 'start')
                            break
                        else:
                            if next in fsm.states:
                                fsm.states[state].add_transition(letter, next)
                                break
                            next = next.subword(1, len(next))
        return fsm

    def reduce_using_automaton(self, word):
        '''
        The method for word reduction using automaton is mentioned in the section 13.1.3 of the Handbook.
        All the elements of the automaton are stored in an array and are given as the input to the automaton.
        If the automaton reaches a dead state that subword is replaced and the automaton is run from the beginning.
        This is repeated until the word reaches the end and the automaton stays in the accept state.

        '''
        # compute all_rules.
        all_rules = self.rules
        self.compute_inverse_rules()
        all_rules.update(self.inverse_rules)

        flag = 1
        while flag:
            flag = 0
            current_state = self.reduction_automaton.states['start']
            word_array = [s for s in word.letter_form_elm]
            for i in range (0, len(word_array)):
                next_state_name = current_state.transitions[word_array[i]]
                next_state = None
                for state in self.reduction_automaton.states:
                    if state == next_state_name:
                        next_state = self.reduction_automaton.states[state]
                if next_state.state_type == "dead":
                    subst = all_rules[next_state_name]
                    word = word.substituted_word(i - len(next_state_name) + 1, i+1, subst)
                    flag = 1
                    break
                current_state = next_state
        return word
