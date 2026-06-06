from __future__ import annotations
from sympy.core.symbol import Symbol
from sympy.logic.boolalg import And, Or, Not, Implies
from sympy.core.relational import (
    Eq, LessThan, GreaterThan, StrictLessThan, StrictGreaterThan
)
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.functions.elementary.piecewise import Piecewise
from sympy.core.function import Function
from sympy.testing.pytest import raises
from sympy.parsing.smtlib.smtlib_parser import parse_smtlib, SMTLibSyntaxError

def test_parse_qf_uf_equalities_or_chain():
    source = """
    (set-option :print-success false)
    (set-logic QF_UF)
    (declare-sort x5 0)
    (declare-fun x4 () x5)
    (declare-fun x2 () Bool)
    (declare-fun x1 () x5)
    (declare-fun x3 () x5)
    (assert x2)
    (assert (not (= x4 x3)))
    (assert (not (or (= x4 x1) (= x1 x3))))
    (check-sat)
    (exit)
    """
    symbols, assertions = parse_smtlib(source)

    x4 = symbols['x4']
    x2 = symbols['x2']
    x1 = symbols['x1']
    x3 = symbols['x3']

    assert len(assertions) == 3
    assert assertions[0] == x2
    assert assertions[1] == Not(Eq(x4, x3))
    assert assertions[2] == Not(Or(Eq(x4, x1), Eq(x1, x3)))


def test_parse_qf_lra_implies_and_lessthan():
    source = """
    (set-option :print-success false)
    (set-logic QF_LRA)
    (declare-fun x2 () Real)
    (declare-fun x1 () Real)
    (assert (not (=> (and (= 0 x1) (= 16 x2)) (<= x2 16))))
    (check-sat)
    (exit)
    """
    symbols, assertions = parse_smtlib(source)

    x1 = symbols['x1']
    x2 = symbols['x2']

    assert len(assertions) == 1
    assert assertions[0] == Not(
        Implies(And(Eq(0, x1), Eq(16, x2)), LessThan(x2, 16))
    )


def test_parse_qf_lra_implies_and_greaterthan():
    source = """
    (set-option :print-success false)
    (set-logic QF_LRA)
    (declare-fun x2 () Real)
    (declare-fun x1 () Real)
    (assert (not (=> (and (= 0 x1) (= 16 x2)) (>= x1 0))))
    (check-sat)
    (exit)
    """
    symbols, assertions = parse_smtlib(source)

    x1 = symbols['x1']
    x2 = symbols['x2']

    assert len(assertions) == 1
    assert assertions[0] == Not(
        Implies(And(Eq(0, x1), Eq(16, x2)), GreaterThan(x1, 0))
    )


def test_parse_qf_lra_nested_conditions():
    source = """
    (set-option :print-success false)
    (set-logic QF_LRA)
    (declare-fun x2 () Real)
    (declare-fun x1 () Real)
    (assert (not (=>
        (and (and (>= x2 1) (and (= 10 x2) (= 0 x1))) (<= x2 12))
        (or (or (= 3 x1) (= x1 1)) (>= 10 x2)))))
    (check-sat)
    (exit)
    """
    symbols, assertions = parse_smtlib(source)

    x1 = symbols['x1']
    x2 = symbols['x2']

    assert len(assertions) == 1

    expected_condition = Not(
        Implies(
            And(
                And(GreaterThan(x2, 1), And(Eq(10, x2), Eq(0, x1))),
                LessThan(x2, 12)
            ),
            Or(
                Or(Eq(3, x1), Eq(x1, 1)),
                GreaterThan(10, x2)
            )
        )
    )
    assert assertions[0] == expected_condition


def test_parse_qf_uf_double_negations():
    source = """
    (set-option :print-success false)
    (set-logic QF_UF)
    (declare-fun x1 () Bool)
    (declare-fun x2 () Bool)
    (declare-fun x4 () Bool)
    (declare-fun x3 () Bool)
    (assert (not x1))
    (assert x3)
    (assert x2)
    (assert (not (not x4)))
    (check-sat)
    (exit)
    """
    symbols, assertions = parse_smtlib(source)

    x1 = symbols['x1']
    x2 = symbols['x2']
    x3 = symbols['x3']
    x4 = symbols['x4']

    assert len(assertions) == 4
    assert assertions[0] == Not(x1)
    assert assertions[1] == x3
    assert assertions[2] == x2
    assert assertions[3] == Not(Not(x4))


def test_parse_qf_uf_implies_double_negations_equalities():
    source = """
    (set-option :print-success false)
    (set-logic QF_UF)
    (declare-sort x3 0)
    (declare-fun x2 () x3)
    (declare-fun x4 () x3)
    (declare-fun x6 () Bool)
    (declare-fun x1 () Bool)
    (declare-fun x5 () Bool)
    (assert (not (not x1)))
    (assert (=> x1 (not x5)))
    (assert (not (= x2 x4)))
    (assert x6)
    (check-sat)
    (exit)
    """
    symbols, assertions = parse_smtlib(source)

    x1 = symbols['x1']
    x2 = symbols['x2']
    x4 = symbols['x4']
    x5 = symbols['x5']
    x6 = symbols['x6']

    assert len(assertions) == 4
    assert assertions[0] == Not(Not(x1))
    assert assertions[1] == Implies(x1, Not(x5))
    assert assertions[2] == Not(Eq(x2, x4))
    assert assertions[3] == x6


def test_parse_smtlib_syntax_errors():
    # Missing closing parenthesis
    with raises(SMTLibSyntaxError):
        parse_smtlib("(assert (= x 1)")

    # Unmatched closing parenthesis
    with raises(SMTLibSyntaxError):
        parse_smtlib("(assert (= x 1)))")


def test_parse_let_bindings():
    source = """
    (set-logic LRA)
    (declare-fun x () Real)
    (declare-fun y () Real)
    (assert (let ((a (+ x 1)) (b (- y 1))) (and (> a 0) (< b 0))))
    """
    symbols, assertions = parse_smtlib(source)
    x = symbols['x']
    y = symbols['y']
    assert len(assertions) == 1
    assert assertions[0] == And(
        StrictGreaterThan(Add(x, 1), 0), StrictLessThan(y - 1, 0)
    )


def test_parse_quantifiers():
    source = """
    (set-logic LRA)
    (declare-fun x () Real)
    (assert (forall ((a Real) (b Real)) (=> (> a b) (> (+ a x) (+ b x)))))
    (assert (exists ((c Real)) (= c (* x 2))))
    """
    symbols, assertions = parse_smtlib(source)
    x = symbols['x']
    a = Symbol('a')
    b = Symbol('b')
    c = Symbol('c')
    assert len(assertions) == 2
    assert assertions[0] == Function('forall')(
        (a, b),
        Implies(
            StrictGreaterThan(a, b), StrictGreaterThan(Add(a, x), Add(b, x))
        )
    )
    assert assertions[1] == Function('exists')((c,), Eq(c, Mul(x, 2)))


def test_parse_macro_definition():
    source = """
    (set-logic LRA)
    (define-fun max ((a Real) (b Real)) Real (ite (> a b) a b))
    (declare-fun x () Real)
    (declare-fun y () Real)
    (assert (= (max x y) 10))
    """
    symbols, assertions = parse_smtlib(source)
    x = symbols['x']
    y = symbols['y']
    assert len(assertions) == 1
    assert assertions[0] == Eq(
        Piecewise((x, StrictGreaterThan(x, y)), (y, True)), 10
    )


def test_parse_bitvector_fallback():
    source = """
    (set-logic QF_BV)
    (declare-fun a () (_ BitVec 32))
    (declare-fun b () (_ BitVec 32))
    (assert (= (bvadd a b) #x00000000))
    """
    symbols, assertions = parse_smtlib(source)
    a = symbols['a']
    b = symbols['b']
    assert len(assertions) == 1
    assert assertions[0] == Eq(Function('bvadd')(a, b), 0)
