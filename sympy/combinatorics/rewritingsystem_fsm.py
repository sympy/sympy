class State(object):
    """A representation of a state managed by a ``StateMachine``.

    Attributes:
        name (instance of FreeGroupElement or string): State name which is also assigned to the Machine.
        transisitons (OrderedDict): Represents all the transitions of the state object.
        is_start (boolean): To mark if the state is a start state.
        is_dead (boolean): To mark if the state is a dead state.
        is_accept (boolean): To mark if the state is an accept state.
    """


    def __init__(self, name, is_start=False, is_dead=False, is_accept=False):
        self.name = name
        self.is_start = is_start
        self.is_dead = is_dead
        self.is_accept = is_accept
        self.transitions = {}
        self.state_type = None
        if is_start:
            self.state_type = "start"
        elif is_accept:
            self.state_type = "accept"
        elif is_dead:
            self.state_type = "dead"

    def add_transition(self, letter, state):
        """
        This method is to add a transition from the current state to a new state.

        Arguments:
            letter: The alphabet element the current state reads to make the state transition.
            state: This will be an instance of the State object which represents a new state after in the transition after the alphabet is read.
        """
        self.transitions[letter] = state

class StateMachine(object):
    """
    The state machine class which manages the states and their transisitons of the automata.

    Attributes:
        states (dictionary): Collection of all registered.
        name (str): Name of the state machine.
    """

    def __init__(self, name):
        self.name = name
        self.states = {} # Contains all the states in the machine.

    def add_state(self, state_name, is_start=False, is_dead=False, is_accept=False):
        """
        This instantiates a state object and stores this in the 'states' list.

        Arguments:
            state_name (instance of FreeGroupElement or string): name of the new states.
            is_start (boolean): To mark if the state is a start state.
            is_dead (boolean): To mark if the state is a dead state.
            is_accept (boolean): To mark if the state is an accept state.
        """
        new_state = State(state_name, is_start, is_dead, is_accept)
        self.states[state_name] = new_state

    def __repr__(self):
        return "%s" % (self.name)
