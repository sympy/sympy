from sympy.assumptions.ask import _ask_single_fact, Q
from sympy.assumptions.cnf import CNF

def test_ask_single_fact():
    assert _ask_single_fact(Q.real, CNF()) is None
    assert _ask_single_fact(Q.even, CNF.from_prop(Q.zero)) is True
    assert _ask_single_fact(Q.even, CNF.from_prop(Q.odd)) is False
    assert _ask_single_fact(Q.even, CNF.from_prop(Q.real)) is None
    assert _ask_single_fact(Q.integer, CNF.from_prop(Q.even | Q.odd)) is None
    assert _ask_single_fact(Q.integer, CNF.from_prop(Q.prime)) is True
    assert _ask_single_fact(Q.prime,   CNF.from_prop(Q.composite)) is False
    assert _ask_single_fact(Q.zero, CNF.from_prop(~Q.even)) is False

    # previously _ask_single_fact gave None for this type of case
    assert _ask_single_fact(Q.zero, CNF.from_prop(~Q.even & Q.real)) is False
