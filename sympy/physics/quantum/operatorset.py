from sympy.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet
from sympy.physics.quantum.operator import Operator
from sympy.physics.quantum.state import StateBase, KetBase, BraBase

__all__ = [
    'operators_to_state',
    'state_to_operators'
]

#state_mapping stores the mappings between states and their associated operators or tuples
#of operators. This should be updated when new classes are written! Entries are of the form
# PxKet : PxOp or something like 3DKet : (ROp, ThetaOp, PhiOp)

# TODO: Update dict with full list of state : operator pairs

state_mapping = { PxKet : PxOp,
                  XKet : XOp }

rev_mapping = dict((v,k) for k,v in state_mapping.iteritems())

def operators_to_state(operators, **options):
    """ operators_to_state

    A global function for mapping operator classes to their associated states.
    It takes either an Operator or a set of operators and returns the state associated
    with these.

    This function can handle both instances of a given operator or just the class itself
    (i.e. both XOp() and XOp)

    There are multiple use cases to consider:

    1) A class or set of classes is passed: First, we try to instantiate default instances for these operators. If this fails,
    then the class is simply returned. If we succeed in instantiating default instances, then we try to call
    state._operators_to_state on the operator instances. If this fails, the class is returned. Otherwise, the instance returned
    by _operators_to_state is returned.

    2) An instance or set of instances is passed: In this case, state._operators_to_state is called on the instances passed. If
    this fails, a state class is returned. If the method returns an instance, that instance is returned.

    In both cases, if the operator class or set does not exist in the state_mapping dictionary, None is returned.

    Parameters
    ===========

    arg: Operator or set
         The class or instance of the operator or set of operators to be mapped to a state

    Examples
    ===========
    >>> from sympy.physics.quantum.cartesian import XOp, PxOp
    >>> from sympy.physics.quantum.operatorset import operators_to_state
    >>> operators_to_state(XOp)
    |x>
    >>> operators_to_state(XOp())
    |x>
    >>> operators_to_state(PxOp)
    |px>
    >>> operators_to_state(PxOp())
    |px>
    """

    if not (isinstance(operators, Operator) or issubclass(operators, Operator) or isinstance(operators, set)):
        raise NotImplementedError("Argument is not an Operator or a set!")

    if isinstance(operators, set):
        for s in operators:
            if not isinstance(operators, Operator) or issubclass(operators, Operator):
                raise NotImplementedError("Members of given set are not all Operators!")

        ops = tuple(operators)

        if ops in rev_mapping: #ops is a list of classes in this case
            #Try to get an object from default instances of the operators...if this fails, return the class
            try:
                op_instances = [op() for op in ops]
                ret = _get_state(rev_mapping[ops], op_instances, **options)
            except NotImplementedError:
                ret = rev_mapping[ops]

            return ret
        else:
            tmp = [type(o) for o in ops]
            classes = tuple(tmp)

            if classes in rev_mapping:
                ret = _get_state(rev_mapping[classes], ops, **options)
            else:
                ret = None

            return ret
    else:
        if operators in rev_mapping:
            try:
                op_instance = operators()
                ret = _get_state(rev_mapping[operators], op_instance, **options)
            except NotImplementedError:
                ret = rev_mapping[operators]

            return ret
        elif type(operators) in rev_mapping:
            return _get_state(rev_mapping[type(operators)], operators, **options)
        else:
            return None

def state_to_operators(state, **options):
    """ state_to_operators

    A global function for mapping state classes to their associated operators or sets of operators.
    It takes either a state class or instance.

    This function can handle both instances of a given state or just the class itself
    (i.e. both XKet() and XKet)

    There are multiple use cases to consider:

    1) A state class is passed: In this case, we first try instantiating a default instance of the class. If this succeeds,
    then we try to call state._state_to_operators on that instance. If the creation of the default instance or if the calling of
    _state_to_operators fails, then either an operator class or set of operator classes is returned. Otherwise, the appropriate
    operator instances are returned.

    2) A state instance is returned: Here, state._state_to_operators is called for the instance. If this fails, then a class
    or set of operator classes is returned. Otherwise, the instances are returned.

    In either case, if the state's class does not exist in state_mapping, None is returned.

    Parameters
    ===========

    arg: StateBase class or instance (or subclasses)
         The class or instance of the state to be mapped to an operator or set of operators

    Examples
    ===========
    >>> from sympy.physics.quantum.cartesian import XKet, PxKet
    >>> from sympy.physics.quantum.operatorset import state_to_operators
    >>> state_to_operators(XKet)
    X
    >>> state_to_operators(XKet())
    X
    >>> state_to_operators(PxKet)
    Px
    >>> state_to_operators(PxKet())
    Px
    """

    if not (isinstance(state, StateBase) or issubclass(state, StateBase)):
        raise NotImplementedError("Argument is not a state!")

    if state in state_mapping: #state is a class
        try:
            state_inst = state()
            ret = _get_ops(state_inst, _tuple_to_set(state_mapping[state]), **options)
        except NotImplementedError:
            ret = state_mapping[state]
    elif type(state) in state_mapping:
        ret = _get_ops(state, _tuple_to_set(state_mapping[type(state)]), **options)
    elif (isinstance(state, BraBase) or issubclass(state, BraBase)) and state.dual_class() in state_mapping:
        ret = _get_ops(state, _tuple_to_set(state_mapping[state.dual_class()]))
    else:
        ret = None

    return _tuple_to_set(ret)

def _get_state(state_class, ops, **options):
    #Try to get a state instance from the operator INSTANCES. If this fails, get the class
    try:
        ret = state_class._operators_to_state(ops, **options)
    except NotImplementedError:
        ret = state_class

    return ret

def _get_ops(state_inst, op_classes, **options):
    #Try to get operator instances from the state INSTANCE. If this fails, just return the classes
    try:
        ret = state_inst._state_to_operators(op_classes, **options)
    except NotImplementedError:
        ret = op_classes

    return ret

def _tuple_to_set(ops):
    if isinstance(ops, tuple):
        return set(ops)
    else:
        return ops
