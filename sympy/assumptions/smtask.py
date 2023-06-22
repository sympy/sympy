"""
Module to evaluate the proposition with assumptions using SMT algorithm.
Can handle relational propositions and assumptions.
"""

from sympy.core.relational import Eq, Ne, Gt, Lt, Ge, Le
from sympy.assumptions.ask import ask, Q
from sympy.logic.boolalg import And, Or, Not
from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.printing.smtlib import smtlib_code
from sympy.core.power import Pow
from sympy.external.importtools import import_module


def smtask(proposition, assumptions=True, context=global_assumptions):
    """
    An alternative to satask that handles relational assumptions. Will pick an
    external SMT solver depending on the problem and what is installed.
    For exmaple, if a problem involves non-linear real functions that Z3 can't
    handle, dreal can be used. Could also call an internal SMT solver
    if one is ever built.


     Parameters
     ==========

     proposition : Any boolean expression.
         Proposition which will be evaluated to boolean value.

     assumptions : Any boolean expression, optional.
         Local assumptions to evaluate the *proposition*.

     Returns
     -------

     ``True``, ``False``, or ``None``

     Examples
     ========

     >>> from sympy.assumptions.smtask import z3ask
     >>> from sympy.abc import z, w
     >>> from sympy import Q
     >>> z3ask(z**2 + w**2 > 0, Q.positive(z) & Q.positive(w))
     True

    """
    assumptions = And(assumptions, And(*list(context)))
    try:
        return z3ask(proposition, assumptions)
    except:
        return None


def z3ask(proposition, assumptions=True):
    """
    Checks the satisfiability of a proposition using the Z3 SMT solver which
    can handle inequalities. Note that it currently lacks many features satask
    has; for example, it does not yet have a function equivalent to satask's
    get_all_relevant_facts().

    Parameters
    ==========

    proposition : Any boolean expression.
        Proposition which will be evaluated to boolean value.

    assumptions : Any boolean expression, optional.
        Local assumptions to evaluate the *proposition*.


    Returns
    -------

    ``True``, ``False``, or ``None``

    Examples
    ========

    >>> from sympy.assumptions.smtask import z3ask
    >>> from sympy.abc import x, y
    >>> from sympy import Q
    >>> z3ask(x**2 + y**2 > 0, Q.positive(x) & Q.positive(y))
    True
    >>> z3ask(y > 0, (y < -x**2 +.5) & Q.integer(y))
    False
    >>> print(z3ask(y > 0, (y < -x**2 +.5))) # y may be positive or not
    None

    """
    z3 = import_module('z3')
    if z3 is None:
        return None


    props = proposition
    _props = Not(proposition)

    s_true = z3.Solver()
    s_true.from_string(smtlib_code(And(props,assumptions)).replace("Pow","^"))
    can_be_true = str(s_true.check())

    s_false = z3.Solver()
    s_false.from_string(smtlib_code(And(_props,assumptions)).replace("Pow","^"))
    can_be_false = str(s_false.check())

    if can_be_true == "unknown" or can_be_false == "unknown":
        return None
    else:
        can_be_true = can_be_true == "sat"
        can_be_false = can_be_false == "sat"

    if can_be_true and can_be_false:
        return None

    if can_be_true and not can_be_false:
        return True

    if not can_be_true and can_be_false:
        return False

    if not can_be_true and not can_be_false:
        raise ValueError("Inconsistent assumptions")
