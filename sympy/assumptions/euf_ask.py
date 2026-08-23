"""
Ask interface backed by the theory of Equality with Uninterpreted Functions.

This mirrors ``lra_satask``: the proposition and its negation are each added
to the assumptions, and the EUF theory solver decides both sides.  What EUF
contributes over plain ``satask`` is congruence -- once ``x = y`` is asserted,
any predicate or function application over ``x`` and ``y`` has to agree.
"""
from __future__ import annotations

from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.logic.inference import satisfiable


def euf_ask(proposition, assumptions=True, context=global_assumptions):
    """
    Evaluate the proposition under the assumptions with the EUF theory solver.

    Returns ``True`` if the assumptions force the proposition, ``False`` if
    they force its negation, and ``None`` if neither.

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.assumptions.ask import Q
    >>> from sympy.assumptions.euf_ask import euf_ask
    >>> x, y = symbols('x y')
    >>> euf_ask(Q.positive(x), Q.positive(y) & Q.eq(x, y))
    True
    >>> euf_ask(Q.positive(x), Q.positive(y) & Q.ne(x, y)) is None
    True
    """
    props = CNF.from_prop(proposition)
    _props = CNF.from_prop(~proposition)

    factbase = EncodedCNF()
    factbase.from_cnf(CNF.from_prop(assumptions))
    if context:
        factbase.add_from_cnf(CNF().extend(context))

    return check_satisfiability(props, _props, factbase)


def check_satisfiability(prop, _prop, factbase):
    """
    Decide whether the assumptions in ``factbase`` leave room for ``prop``,
    for ``_prop``, or only for one of them.
    """
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(prop)
    sat_false.add_from_cnf(_prop)

    can_be_true = satisfiable(sat_true, use_euf_theory=True) is not False
    can_be_false = satisfiable(sat_false, use_euf_theory=True) is not False

    if can_be_true and can_be_false:
        return None
    if can_be_true and not can_be_false:
        return True
    if not can_be_true and can_be_false:
        return False
    raise ValueError("Inconsistent assumptions")
