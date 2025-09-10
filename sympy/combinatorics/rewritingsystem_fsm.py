from __future__ import annotations
from typing import TypeVar, Mapping, Collection, DefaultDict


class State:
    '''
    A representation of a state managed by a ``StateMachine``.

    Attributes:
        name (instance of FreeGroupElement or string) -- State name which is also assigned to the Machine.
        transisitons (OrderedDict) -- Represents all the transitions of the state object.
        state_type (string) -- Denotes the type (accept/start/dead) of the state.
        rh_rule (instance of FreeGroupElement) -- right hand rule for dead state.
        state_machine (instance of StateMachine object) -- The finite state machine that the state belongs to.
    '''

    def __init__(self, name, state_machine, state_type=None, rh_rule=None):
        self.name = name
        self.transitions = {}
        self.state_machine = state_machine
        self.state_type = state_type[0]
        self.rh_rule = rh_rule

    def add_transition(self, letter, state):
        '''
        Add a transition from the current state to a new state.

        Keyword Arguments:
            letter -- The alphabet element the current state reads to make the state transition.
            state -- This will be an instance of the State object which represents a new state after in the transition after the alphabet is read.

        '''
        self.transitions[letter] = state

class StateMachine:
    '''
    Representation of a finite state machine the manages the states and the transitions of the automaton.

    Attributes:
        states (dictionary) -- Collection of all registered `State` objects.
        name (str) -- Name of the state machine.
    '''

    def __init__(self, name, automaton_alphabet):
        self.name = name
        self.automaton_alphabet = automaton_alphabet
        self.states = {} # Contains all the states in the machine.
        self.add_state('start', state_type='s')

    def add_state(self, state_name, state_type=None, rh_rule=None):
        '''
        Instantiate a state object and stores it in the 'states' dictionary.

        Arguments:
            state_name (instance of FreeGroupElement or string) -- name of the new states.
            state_type (string) -- Denotes the type (accept/start/dead) of the state added.
            rh_rule (instance of FreeGroupElement) -- right hand rule for dead state.

        '''
        new_state = State(state_name, self, state_type, rh_rule)
        self.states[state_name] = new_state

    def __repr__(self):
        return "%s" % (self.name)


_N = TypeVar('_N')
_T = TypeVar('_T')


def epsilon_closure(
    transition: Mapping[_N, Mapping[_T, Collection[_N]]],
    epsilon: _T,
    states: Collection[_N]
) -> set[_N]:
    r"""Compute the $\epsilon$ closure for a nondeterministic finite state
    machine.

    Explanation
    ===========

    Given a transition table for nondeterministic finite state machine with
    $\epsilon$ is defined as a mapping:

    .. math::
        N \to ((T | \epsilon) \to 2^N)

    And given the initial set $S$ of states in $2^N$,
    the $\epsilon$ closure is defined as a set of all states reachable by the
    elements of $S$ by $\epsilon$ transitions in zero or finite number of
    steps.

    Parameters
    ==========

    transition: Mapping[N, Mapping[T, Collection[N]]]
        A transition table of nondeterministic finite machine in the type
        $N \to ((T | \epsilon) \to 2^N)$ represented as a dictionary of
        dictionaries of sets.

    epsilon: T
        An empty transition symbol $\epsilon$.

    states: Collection[N]
        A set of initial states of type $2^N$.

    Returns
    =======

    set[N]
        A set of states of type $2^N$ which is reflexively-transitively closed
        under $\epsilon$ transition.

    Examples
    ========

    Computing the $\epsilon$ closure of a nondeterministic finite state machine
    for a single state:

    >>> from sympy.combinatorics.rewritingsystem_fsm import epsilon_closure

    >>> transition = {
    ...     0: {'e': {1, 7}},
    ...     1: {'e': {2, 4}},
    ...     2: {'a': {3}},
    ...     3: {'e': {6}},
    ...     4: {'b': {5}},
    ...     5: {'e': {6}},
    ...     6: {'e': {1, 7}},
    ...     7: {'a': {8}},
    ...     8: {'b': {9}},
    ...     9: {'b': {10}},
    ... }
    >>> epsilon_closure(transition, 'e', {0})
    {0, 1, 2, 4, 7}

    Computing the $\epsilon$ closure of a nondeterministic finite state machine
    for multiple states:

    >>> epsilon_closure(transition, 'e', {3, 8})
    {1, 2, 3, 4, 6, 7, 8}

    References
    ==========

    .. [*] Aho, Alfred V., Monica S. Lam, Ravi Sethi and Jeffrey D. Ullman.
       "Compilers: Principles, Techniques, and Tools (2nd Edition)." (2006).
    """
    stack: list[_N] = list(states)
    closure = set(stack)

    while stack:
        state = stack.pop()
        if state not in transition:
            continue

        if epsilon not in transition[state]:
            continue

        for new_state in transition[state][epsilon]:
            if new_state not in closure:
                closure.add(new_state)
                stack.append(new_state)

    return closure


def instructions(
    transition: Mapping[_N, Mapping[_T, Collection[_N]]],
    states: Collection[_N],
) -> DefaultDict[_T, set[_N]]:
    r"""Compute the union of instructions from the given set of input states.

    Explanation
    ===========

    Given a transition table for nondeterministic finite state machine defined
    as a mapping:

    .. math::
        \delta: N \to (T \to 2^N)

    And given the initial set $S$ of states in $2^N$,
    the instructions from $S$ is defined as a mapping:

    .. math::
        \delta': T \to 2^N

    where

    .. math::
        \delta'(a) = \bigcup_{A \in S} \delta(A)(a)

    Parameters
    ==========

    transition: Mapping[N, Mapping[T, Collection[N]]]
        A transition table of nondeterministic finite machine in the type
        $N \to (T \to 2^N)$ represented as a dictionary of
        dictionaries of sets.

    states: Collection[N]
        A set of initial states of type $2^N$.

    Returns
    =======

    defaultdict[T, set[N]]
        A instruction table of type $T \to 2^N$ represented as a dictionary
        of sets.

    Examples
    ========

    >>> from sympy.combinatorics.rewritingsystem_fsm import instructions

    >>> transition = {
    ...     0: {'e': {1, 7}},
    ...     1: {'e': {2, 4}},
    ...     2: {'a': {3}},
    ...     3: {'e': {6}},
    ...     4: {'b': {5}},
    ...     5: {'e': {6}},
    ...     6: {'e': {1, 7}},
    ...     7: {'a': {8}},
    ...     8: {'b': {9}},
    ...     9: {'b': {10}},
    ... }

    Computing all possible transitions that is reachable by the state $2$ for
    each input symbols:

    >>> instructions(transition, {2})
    {'a': {3}}

    Computing all possible transitions that is reachable by the states
    $\{0, 1, 2, 4, 7\}$ for each input symbols:

    >>> instructions(transition, {0, 1, 2, 4, 7})
    {'a': {3, 8}, 'b': {5}, 'e': {1, 2, 4, 7}}

    References
    ==========

    .. [*] Johnson, J.H., Wood, D. (1997). Instruction computation in subset
       construction. In: Raymond, D., Wood, D., Yu, S. (eds) Automata
       Implementation. WIA 1996. Lecture Notes in Computer Science, vol 1260.
       Springer, Berlin, Heidelberg. https://doi.org/10.1007/3-540-63174-7_6
    """
    out: DefaultDict[_T, set[_N]] = DefaultDict(set)
    for state in states:
        if state not in transition:
            continue
        for enter in transition[state]:
            out[enter] |= set(transition[state][enter])
    return out


def subset_construction(
    transition: Mapping[_N, Mapping[_T, Collection[_N]]],
    start: Collection[_N]
) -> dict[frozenset[_N], dict[_T, frozenset[_N]]]:
    r"""Convert the nondeterministic finite state machine without $\epsilon$
    to a deterministic finite state machine.

    Explanation
    ===========

    Given a transition table for nondeterministic finite state machine defined
    as a mapping:

    .. math::
        N \to (T \to 2^N)

    An equivalent deterministic finite state machine can be computed as a
    mapping from the set of instructions to the set of instructions:

    .. math::
        2^N \to (T \to 2^N)

    Parameters
    ==========

    transition: Mapping[N, Mapping[T, Collection[N]]]
        A transition table of nondeterministic finite state machine of the type
        $N \to (T \to 2^N)$ represented as a dictionary of
        dictionaries of sets.

    start: Collection[N]
        The initial $\epsilon$ closure of the set of start symbols.

    Returns
    =======

    dict[frozenset[N], dict[T, frozenset[N]]]
        The transition table of the deterministic finite state machine of the
        type $2^N \to (T \to 2^N)$ represented as a dictionary of
        dictionaries of frozen sets.

    Examples
    ========

    Compute the table for deterministic finite state machine without $\epsilon$
    from the initial state $1$:

    >>> from sympy.combinatorics.rewritingsystem_fsm import subset_construction

    >>> transition = {
    ...     0: {'e': {1, 7}},
    ...     1: {'e': {2, 4}},
    ...     2: {'a': {3}},
    ...     3: {'e': {6}},
    ...     4: {'b': {5}},
    ...     5: {'e': {6}},
    ...     6: {'e': {1, 7}},
    ...     7: {'a': {8}},
    ...     8: {'b': {9}},
    ...     9: {'b': {10}},
    ... }
    >>> dfa = subset_construction(transition, {0})
    >>> for A in dfa:
    ...     for x in dfa[A]:
    ...         B = dfa[A][x]
    ...         print(f"{set(A)}, {x} -> {set(B)}")
    {0}, e -> {1, 7}
    {1, 7}, e -> {2, 4}
    {1, 7}, a -> {8}
    {8}, b -> {9}
    {9}, b -> {10}
    {2, 4}, a -> {3}
    {2, 4}, b -> {5}
    {5}, e -> {6}
    {6}, e -> {1, 7}
    {3}, e -> {6}

    Notes
    =====

    Although accepting the set of start symbols is has computationally
    well-defined behavior, it may be ambiguous generalization for users,
    and in most applications, supplying a single element set like ``{'A'}``
    is sufficient.

    The interface is designed such that it aligns with the inteface of
    :func:`subset_construction_epsilon`.

    References
    ==========

    .. [*] Gertjan van Noord; Treatment of Epsilon Moves in Subset
       Construction. Computational Linguistics 2000; 26 (1): 61-76.
       doi: https://doi.org/10.1162/089120100561638

    .. [*] Aho, Alfred V., Monica S. Lam, Ravi Sethi and Jeffrey D. Ullman.
       "Compilers: Principles, Techniques, and Tools (2nd Edition)." (2006).
    """
    _start = frozenset(start)
    stack: list[frozenset[_N]] = [_start]
    out: dict[frozenset[_N], dict[_T, frozenset[_N]]] = {}

    while stack:
        state = stack.pop()
        for enter, _new_state in instructions(transition, state).items():
            new_state = frozenset(_new_state)
            if state not in out:
                out[state] = {}
            out[state][enter] = new_state
            if new_state not in out:
                stack.append(new_state)

    return out


def subset_construction_epsilon(
    transition: Mapping[_N, Mapping[_T, Collection[_N]]],
    epsilon: _T,
    start: Collection[_N]
) -> dict[frozenset[_N], dict[_T, frozenset[_N]]]:
    r"""Convert the nondeterministic finite state machine without $\epsilon$
    to a deterministic finite state machine.

    Explanation
    ===========

    Given a transition table for nondeterministic finite state machine defined
    as a mapping:

    .. math::
        N \to ((T | \epsilon) \to 2^N)

    An equivalent deterministic finite state machine can be computed as a
    mapping from the set of instructions to the set of instructions:

    .. math::
        2^N \to (T \to 2^N)

    Parameters
    ==========

    transition: Mapping[N, Mapping[T, Collection[N]]]
        A transition table of nondeterministic finite state machine of the type
        $N \to ((T | \epsilon) \to 2^N)$ represented as a dictionary of
        dictionaries of sets.

    start: Collection[N]
        The initial $\epsilon$ closure of the set of start symbols.
        We warn that asserting the $\epsilon$ closure is assumed for users
        composing the function :func:`epsilon_closure` in priori, and not done
        here automatically to avoid some redundancy of computation and bring
        some symmetry of implementation.

    epsilon: T
        The symbol for $\epsilon$ transition.

    Returns
    =======

    dict[frozenset[N], dict[T, frozenset[N]]]
        The transition table of the deterministic finite state machine

    Examples
    ========

    Compute the table for deterministic finite state machine with $\epsilon$
    from the initial state $1$:

    >>> from sympy.combinatorics.rewritingsystem_fsm import (
    ...     epsilon_closure, subset_construction_epsilon)

    >>> transition = {
    ...     0: {'e': {1, 7}},
    ...     1: {'e': {2, 4}},
    ...     2: {'a': {3}},
    ...     3: {'e': {6}},
    ...     4: {'b': {5}},
    ...     5: {'e': {6}},
    ...     6: {'e': {1, 7}},
    ...     7: {'a': {8}},
    ...     8: {'b': {9}},
    ...     9: {'b': {10}},
    ... }
    >>> start = epsilon_closure(transition, 'e', {0})
    >>> dfa = subset_construction_epsilon(transition, 'e', start)
    >>> for A in dfa:
    ...     for x in dfa[A]:
    ...         B = dfa[A][x]
    ...         print(f"{set(A)}, {x} -> {set(B)}")
    {0, 1, 2, 4, 7}, a -> {1, 2, 3, 4, 6, 7, 8}
    {0, 1, 2, 4, 7}, b -> {1, 2, 4, 5, 6, 7}
    {1, 2, 4, 5, 6, 7}, a -> {1, 2, 3, 4, 6, 7, 8}
    {1, 2, 4, 5, 6, 7}, b -> {1, 2, 4, 5, 6, 7}
    {1, 2, 3, 4, 6, 7, 8}, a -> {1, 2, 3, 4, 6, 7, 8}
    {1, 2, 3, 4, 6, 7, 8}, b -> {1, 2, 4, 5, 6, 7, 9}
    {1, 2, 4, 5, 6, 7, 9}, a -> {1, 2, 3, 4, 6, 7, 8}
    {1, 2, 4, 5, 6, 7, 9}, b -> {1, 2, 4, 5, 6, 7, 10}
    {1, 2, 4, 5, 6, 7, 10}, a -> {1, 2, 3, 4, 6, 7, 8}
    {1, 2, 4, 5, 6, 7, 10}, b -> {1, 2, 4, 5, 6, 7}

    Notes
    =====

    Deterministic finite state machine shouldn't use $\epsilon$ transition,
    so the type for $\epsilon$ can be removed from the output, however, because
    of the limitations of type checking libraries of python, it may not be
    able to automatically infer such type.

    References
    ==========

    .. [*] Gertjan van Noord; Treatment of Epsilon Moves in Subset
       Construction. Computational Linguistics 2000; 26 (1): 61-76.
       doi: https://doi.org/10.1162/089120100561638

    .. [*] Aho, Alfred V., Monica S. Lam, Ravi Sethi and Jeffrey D. Ullman.
       "Compilers: Principles, Techniques, and Tools (2nd Edition)." (2006).
    """
    _start = frozenset(start)
    stack: list[frozenset[_N]] = [_start]
    out: dict[frozenset[_N], dict[_T, frozenset[_N]]] = {}

    while stack:
        state = stack.pop()
        for enter, _new_state in instructions(transition, state).items():
            if enter == epsilon:
                continue

            new_state = frozenset(epsilon_closure(transition, epsilon, _new_state))
            if state not in out:
                out[state] = {}
            out[state][enter] = new_state
            if new_state not in out:
                stack.append(new_state)

    return out
