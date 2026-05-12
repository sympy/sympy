from __future__ import annotations
from sympy.combinatorics.fp_groups import FpGroup
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.rewritingsystem_fsm import StateMachine
from sympy.testing.pytest import raises

def test_rewriting():
    F, a, b = free_group("a, b")
    G = FpGroup(F, [a*b*a**-1*b**-1])
    a, b = G.generators
    R = G._rewriting_system
    assert R.is_confluent

    assert G.reduce(b**-1*a) == a*b**-1
    assert G.reduce(b**3*a**4*b**-2*a) == a**5*b
    assert G.equals(b**2*a**-1*b, b**4*a**-1*b**-1)

    assert R.reduce_using_automaton(b*a*a**2*b**-1) == a**3
    assert R.reduce_using_automaton(b**3*a**4*b**-2*a) == a**5*b
    assert R.reduce_using_automaton(b**-1*a) == a*b**-1

    G = FpGroup(F, [a**3, b**3, (a*b)**2])
    R = G._rewriting_system
    R.make_confluent()
    # R._is_confluent should be set to True after
    # a successful run of make_confluent
    assert R.is_confluent
    # but also the system should actually be confluent
    assert R._check_confluence()
    assert G.reduce(b*a**-1*b**-1*a**3*b**4*a**-1*b**-15) == a**-1*b**-1
    # check for automaton reduction
    assert R.reduce_using_automaton(b*a**-1*b**-1*a**3*b**4*a**-1*b**-15) == a**-1*b**-1

    G = FpGroup(F, [a**2, b**3, (a*b)**4])
    R = G._rewriting_system
    assert G.reduce(a**2*b**-2*a**2*b) == b**-1
    assert R.reduce_using_automaton(a**2*b**-2*a**2*b) == b**-1
    assert G.reduce(a**3*b**-2*a**2*b) == a**-1*b**-1
    assert R.reduce_using_automaton(a**3*b**-2*a**2*b) == a**-1*b**-1
    # Check after adding a rule
    R.add_rule(a**2, b)
    assert R.reduce_using_automaton(a**2*b**-2*a**2*b) == b**-1
    assert R.reduce_using_automaton(a**4*b**-2*a**2*b**3) == b

    R.set_max(15)
    raises(RuntimeError, lambda:  R.add_rule(a**-3, b))
    R.set_max(20)
    R.add_rule(a**-3, b)

    assert R.add_rule(a, a) == set()

def test_state_machine_accepts():

    # Build odd-length DFA over {"0", "1"}
    M = StateMachine("odd_length", ["0", "1"])
    M.add_state("q1", state_type="a")
    start = M.states["start"]
    q1 = M.states["q1"]
    start.add_transition("0", q1)
    start.add_transition("1", q1)
    q1.add_transition("0", start)
    q1.add_transition("1", start)

    results = M.accepts(["0", "010", "01", "", "2"])
    assert results["0"] is True       # length 1
    assert results["010"] is True     # length 3
    assert results["01"] is False     # length 2
    assert results[""] is False       # empty
    assert results["2"] is False      # not in alphabet


def test_state_machine_validate():

    # Valid machine passes
    M = StateMachine("odd_length", ["0", "1"])
    M.add_state("q1", state_type="a")
    start = M.states["start"]
    q1 = M.states["q1"]
    start.add_transition("0", q1)
    start.add_transition("1", q1)
    q1.add_transition("0", start)
    q1.add_transition("1", start)
    assert M.validate() is True

    # No accept states
    M2 = StateMachine("no_accept", ["0", "1"])
    raises(ValueError, lambda: M2.validate())

    # Missing transition
    M3 = StateMachine("missing", ["0", "1"])
    M3.add_state("q1", state_type="a")
    M3.states["start"].add_transition("0", M3.states["q1"])
    raises(ValueError, lambda: M3.validate())

    # Empty alphabet
    M4 = StateMachine("empty_alpha", [])
    raises(ValueError, lambda: M4.validate())

    # Extra transition outside alphabet
    M5 = StateMachine("extra_trans", ["0", "1"])
    M5.add_state("q1", state_type="a")
    M5.states["start"].add_transition("0", M5.states["q1"])
    M5.states["start"].add_transition("1", M5.states["q1"])
    M5.states["start"].add_transition("2", M5.states["q1"])
    M5.states["q1"].add_transition("0", M5.states["start"])
    M5.states["q1"].add_transition("1", M5.states["start"])
    raises(ValueError, lambda: M5.validate())
