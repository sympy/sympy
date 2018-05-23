class State(object):
    """A representation of a state managed by a ``StateMachine``.

    Attributes:
        name (str): State name which is also assigned to the Machine.
        transisitons (OrderedDict): Represents all the transitions of the state object.
        is_start (boolean): To mark if the state is a start state.
        is_dead (boolean): To mark if the state is a dead state.
        is_accept (boolean): To mark if te state is an accept state.
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

    def add_transition(self, alphabet, state):
        """
        This method is to add a transition from the current state to a new state. 

        Arguments:
            alphabet: The alphabet the current state reads to make the state transition. 
            state: This will be an instance of the State object which represents a new state after in the transition after the alphabet is read. 
        """
        self.transitions[alphabet] = state

    def set_start(self):
        self.is_start = True
        self.state_type = "start"
    
    def set_dead(self):
        self.is_accept = False
        self.is_dead = True
        self.transisitons = {} # empty the transitions if it is a dead state.
        self.state_type = "dead"

    def set_accept(self):
        self.is_dead = False
        self.is_accept = True
        self.state_type = "accept"

class StateMachine(object):
    """ The state machine class which manages the states and their transisitons of the automata.

    Attribute:
        states (list of States): Collection of all registered.
        name (str): Name of the state machine.
    """

    def __init__(self, name):
        self.name = name
        self.states = [] # Contains all the states in the machine.
        self.state_names = []

    def add_state(self, state_name, is_start=False, is_dead=False, is_accept=False):
        """
        This instantiates a state object and stores this in the 'states' list.

        Arguments: 
            Same as the __init__ function of the State class.
        """
        new_state = State(state_name, is_start, is_dead, is_accept)
        self.states.append(new_state)
        self.state_names.append(new_state.name)

    def __repr__(self):
        return "%s" % (self.name)