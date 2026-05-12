from __future__ import annotations
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

    def accepts(self, words):
        '''
        Check whether the state machine accepts each word in a list.

        Arguments:
            words -- a list of strings, each string is checked against
                the machine character by character.

        Returns
        =======
        dict
            A dictionary mapping each word to ``True`` if accepted,
            ``False`` otherwise.
        '''
        results = {}
        for word in words:
            state = self.states['start']
            accepted = True
            for symbol in word:
                if symbol not in self.automaton_alphabet:
                    accepted = False
                    break
                if symbol not in state.transitions:
                    accepted = False
                    break
                state = state.transitions[symbol]
            if accepted:
                results[word] = state.state_type == 'a'
            else:
                results[word] = False
        return results

    def validate(self):
        '''
        Validate that the machine is a well-formed DFA by checking all five
        components of the formal 5-tuple (Q, Sigma, delta, q0, F).

        - Q     -- at least one state exists; no empty state names.
        - Sigma -- the alphabet is non-empty; no empty symbols.
        - delta -- every state has a transition for every symbol in the
                   alphabet; no transitions outside the alphabet; all
                   transition targets exist in Q.
        - q0    -- exactly one start state exists (state_type ``'s'``);
                   a state named ``'start'`` must be present in Q.
        - F     -- at least one accept state exists (state_type ``'a'``).

        Returns
        =======
        True
            If all checks pass.

        Raises
        ======
        ValueError
            With a descriptive message for the first violated check.
        '''
        # Q: at least one state
        if not self.states:
            raise ValueError(
                "Q is empty: the machine must have at least one state.")

        # Q: no empty state names
        if any(name == '' for name in self.states):
            raise ValueError(
                "Q violation: state name cannot be an empty string.")

        # Sigma: non-empty alphabet
        if not self.automaton_alphabet:
            raise ValueError(
                "Sigma is empty: the alphabet must contain at least one symbol.")

        # Sigma: no empty symbols
        if any(s == '' for s in self.automaton_alphabet):
            raise ValueError(
                "Sigma violation: alphabet symbol cannot be an empty string.")

        # q0: 'start' state must exist in Q
        if 'start' not in self.states:
            raise ValueError(
                "q0 violation: no state named 'start' found in Q.")

        alphabet = set(self.automaton_alphabet)
        valid_types = {'s', 'a', 'd'}

        # state_type validity
        for name, state in self.states.items():
            if state.state_type not in valid_types:
                raise ValueError(
                    f"State '{name}' has unrecognized state_type "
                    f"'{state.state_type}'; must be one of {valid_types}.")

        # q0: exactly one start state
        start_states = [s for s in self.states.values() if s.state_type == 's']
        if len(start_states) != 1:
            raise ValueError(
                f"q0 violation: expected exactly 1 start state, "
                f"found {len(start_states)}.")

        # F: at least one accept state
        accept_states = [s for s in self.states.values() if s.state_type == 'a']
        if not accept_states:
            raise ValueError(
                "F is empty: the machine must have at least one accept state.")

        # delta: missing, extra, and invalid targets
        for name, state in self.states.items():
            defined = set(state.transitions.keys())

            missing = alphabet - defined
            if missing:
                raise ValueError(
                    f"delta violation: state '{name}' is missing transitions "
                    f"for: {missing}.")

            extra = defined - alphabet
            if extra:
                raise ValueError(
                    f"delta violation: state '{name}' has transitions outside "
                    f"the alphabet: {extra}.")

            for symbol, dest in state.transitions.items():
                if dest.name not in self.states:
                    raise ValueError(
                        f"delta violation: state '{name}' on '{symbol}' "
                        f"transitions to unknown state '{dest.name}'.")

        return True

    def __repr__(self):
        return "%s" % (self.name)
