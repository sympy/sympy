from sympy.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet
from sympy.physics.quantum.operator import Operator
from sympy.physics.quantum.state import StateBase, KetBase, BraBase

__all__ = [
    'operator_to_state',
    'state_to_operator'
]

#state_mapping stores the mappings between states and their associated operators or tuples
#of operators. This should be updated when new classes are written! Entries are of the form
# PxKet : PxOp or something like 3DKet : (ROp, ThetaOp, PhiOp)

# TODO: Update dict with full list of state : operator pairs

state_mapping = { PxKet : PxOp,
                  XKet : XOp }

rev_mapping = dict((v,k) for k,v in state_mapping.iteritems())

def operator_to_state(arg, **options):
    """ operator_to_state

    A global function for mapping operator classes to their associated states.
    It takes either an Operator or a set of operators and returns the state associated
    with these.

    This function can handle both instances of a given operator or just the class itself
    (i.e. both XOp() and XOp)

    Parameters
    ===========

    arg: Operator or set
         The class or instance of the operator or set of operators to be mapped to a state

    Examples
    ===========
    >>> from sympy.physics.quantum.cartesian import XOp, PxOp
    >>> from sympy.physics.quantum.operatorset import operator_to_state
    >>> operator_to_state(XOp)
    <class 'sympy.physics.quantum.cartesian.XKet'>
    >>> operator_to_state(XOp())
    <class 'sympy.physics.quantum.cartesian.XKet'>
    >>> operator_to_state(PxOp)
    <class 'sympy.physics.quantum.cartesian.PxKet'>
    >>> operator_to_state(PxOp())
    <class 'sympy.physics.quantum.cartesian.PxKet'>
    """

    if not (isinstance(arg, Operator) or issubclass(arg, Operator) or isinstance(arg, set)):
        raise NotImplementedError("Argument is not an Operator or a set!")

    if isinstance(arg, set):
        for s in arg:
            if not isinstance(arg, Operator) or issubclass(arg, Operator):
                raise NotImplementedError("Members of given set are not all Operators!")

        ops = tuple(arg)

        if ops in rev_mapping:
            return rev_mapping[ops]
        else:
            tmp = [type(o) for o in ops]
            classes = tuple(tmp)

            if classes in rev_mapping:
                return rev_mapping[classes]
            else:
                return None
    else:
        if arg in rev_mapping:
            return rev_mapping[arg]
        elif type(arg) in rev_mapping:
            return rev_mapping[type(arg)]
        else:
            return None

def state_to_operator(arg, **options):
    """ state_to_operator

    A global function for mapping state classes to their associated operators or sets of operators.
    It takes either a state class or instance.

    This function can handle both instances of a given state or just the class itself
    (i.e. both XKet() and XKet)

    Parameters
    ===========

    arg: StateBase class or instance (or subclasses)
         The class or instance of the state to be mapped to an operator or set of operators

    Examples
    ===========
    >>> from sympy.physics.quantum.cartesian import XKet, PxKet
    >>> from sympy.physics.quantum.operatorset import state_to_operator
    >>> state_to_operator(XKet)
    <class 'sympy.physics.quantum.cartesian.XOp'>
    >>> state_to_operator(XKet())
    <class 'sympy.physics.quantum.cartesian.XOp'>
    >>> state_to_operator(PxKet)
    <class 'sympy.physics.quantum.cartesian.PxOp'>
    >>> state_to_operator(PxKet())
    <class 'sympy.physics.quantum.cartesian.PxOp'>
    """

    if not (isinstance(arg, StateBase) or issubclass(arg, StateBase)):
        raise NotImplementedError("Argument is not a state!")

    if arg in state_mapping:
        ret = state_mapping[arg]
    elif type(arg) in state_mapping:
        ret = state_mapping[type(arg)]
    elif isinstance(arg, BraBase) and arg.dual_class in state_mapping:
        ret = state_mapping[arg.dual_class]
    else:
        ret = None

    if isinstance(ret, tuple):
        return set(ret)
    else:
        return ret
