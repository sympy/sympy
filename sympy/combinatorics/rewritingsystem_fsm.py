class State(object):
    '''
    A representation of a state managed by a ``StateMachine``.

    Attributes:
        name (instance of FreeGroupElement or string): State name which is also assigned to the Machine.
        transisitons (OrderedDict): Represents all the transitions of the state object.
        is_start (boolean): To mark if the state is a start state.
        is_dead (boolean): To mark if the state is a dead state.
        is_accept (boolean): To mark if the state is an accept state.
        rh_rule (instance of FreeGroupElement): right hand rule for dead state.
    '''


    def __init__(self, name, is_start=False, is_dead=False, is_accept=False, subst=None):
        self.name = name
        self.is_start = is_start
        self.is_dead = is_dead
        self.is_accept = is_accept
        self.transitions = {}
        self.rh_rule = subst
        self.state_type = None
        if is_start:
            self.state_type = "start"
        elif is_accept:
            self.state_type = "accept"
        elif is_dead:
            self.state_type = "dead"

    def add_transition(self, letter, state):
        '''
        Add a transition from the current state to a new state.

        Keyword Arguments:
            letter -- The alphabet element the current state reads to make the state transition.
            state -- This will be an instance of the State object which represents a new state after in the transition after the alphabet is read.

        '''
        self.transitions[letter] = state

class StateMachine(object):
    '''
    Represents the attributes and methods of a finite state machine.

    Attributes:
        states (dictionary) -- Collection of all registered `State` objects.
        name (str) -- Name of the state machine.
    '''

    def __init__(self, name):
        self.name = name
        self.states = {} # Contains all the states in the machine.

    def add_state(self, state_name, is_start=False, is_dead=False, is_accept=False, subst=None):
        '''
        Instantiate a state object and stores it in the 'states' dictionary.

        Arguments:
            state_name (instance of FreeGroupElement or string) -- name of the new states.
            is_start (boolean) -- To mark if the state is a start state.
            is_dead (boolean) -- To mark if the state is a dead state.
            is_accept (boolean) -- To mark if the state is an accept state.
            subt (instance of FreeGroupElement) -- right hand rule for dead state.

        '''
        new_state = State(state_name, is_start, is_dead, is_accept, subst)
        self.states[state_name] = new_state

    def __repr__(self):
        return "%s" % (self.name)
