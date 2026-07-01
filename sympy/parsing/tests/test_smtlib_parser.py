from __future__ import annotations

from sympy.assumptions import Q
from sympy.core.add import Add
from sympy.core.function import Function
from sympy.core.mul import Mul
from sympy.core.relational import (
    Eq, GreaterThan, LessThan, Ne, StrictGreaterThan, StrictLessThan
)
from sympy.core.symbol import Symbol
from sympy.external import import_module
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.piecewise import Piecewise
from sympy.logic.boolalg import And, Implies, Not, Or, Xor
from sympy.parsing.smtlib.lark.smtlib_parser import (
    SMTLibSyntaxError, UnknownSMTLibCommandError, parse_smtlib
)
from sympy.testing.pytest import raises

lark = import_module('lark')

# disable tests if lark is not present
disabled = lark is None


def test_parse_equalities_or_chain():
    source = """
    (set-option :print-success false)
    (set-logic QF_UF)
    (declare-sort x5 0)
    (declare-fun x4 () x5)
    (declare-fun x6 () x5)
    (declare-fun x7 () x5)
    (declare-fun x1 () x5)
    (declare-fun x3 () x5)
    (assert (= x6 x7))
    (assert (not (= x4 x3)))
    (assert (not (or (= x4 x1) (= x1 x3))))
    (check-sat)
    (exit)
    """
    parsed_symbols, assertions = parse_smtlib(source)

    x1 = parsed_symbols['x1']
    x6 = parsed_symbols['x6']
    x7 = parsed_symbols['x7']
    x3 = parsed_symbols['x3']
    x4 = parsed_symbols['x4']

    assert len(assertions) == 3
    assert assertions[0] == Eq(x6, x7)
    assert assertions[1] == Not(Eq(x4, x3))
    assert assertions[2] == Not(Or(Eq(x4, x1), Eq(x1, x3)))


def test_parse_implies_double_negations_equalities():
    source = """
    (set-option :print-success false)
    (set-logic QF_UF)
    (declare-sort x0 0)
    (declare-fun x2 () x0)
    (declare-fun x4 () x0)
    (declare-fun x7 () x0)
    (declare-fun x8 () x0)
    (declare-fun x1 () x0)
    (declare-fun x3 () x0)
    (declare-fun x5 () x0)
    (declare-fun x6 () x0)
    (assert (not (not (= x1 x3))))
    (assert (=> (= x1 x3) (not (= x5 x6))))
    (assert (not (= x2 x4)))
    (assert (= x7 x8))
    (check-sat)
    (exit)
    """
    parsed_symbols, assertions = parse_smtlib(source)

    x1 = parsed_symbols['x1']
    x3 = parsed_symbols['x3']
    x2 = parsed_symbols['x2']
    x4 = parsed_symbols['x4']
    x5 = parsed_symbols['x5']
    x6 = parsed_symbols['x6']
    x7 = parsed_symbols['x7']
    x8 = parsed_symbols['x8']

    assert len(assertions) == 4
    assert assertions[0] == Not(Not(Eq(x1, x3)))
    assert assertions[1] == Implies(Eq(x1, x3), Not(Eq(x5, x6)))
    assert assertions[2] == Not(Eq(x2, x4))
    assert assertions[3] == Eq(x7, x8)


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
    parsed_symbols, assertions = parse_smtlib(source)
    x = parsed_symbols['x']
    y = parsed_symbols['y']
    assert len(assertions) == 3
    assert assertions[0] == Q.real(x)
    assert assertions[1] == Q.real(y)
    assert assertions[2] == And(
        StrictGreaterThan(Add(x, 1), 0), StrictLessThan(y - 1, 0)
    )


def test_parse_quantifiers():
    source = """
    (set-logic LRA)
    (declare-fun x () Real)
    (assert (forall ((a Real) (b Real)) (=> (> a b) (> (+ a x) (+ b x)))))
    (assert (exists ((c Real)) (= c (* x 2))))
    """
    parsed_symbols, assertions = parse_smtlib(source)
    x = parsed_symbols['x']
    a = Symbol('a')
    b = Symbol('b')
    c = Symbol('c')
    assert len(assertions) == 3
    assert assertions[0] == Q.real(x)
    assert assertions[1] == Function('forall')(
        (a, b),
        Implies(
            StrictGreaterThan(a, b), StrictGreaterThan(Add(a, x), Add(b, x))
        )
    )
    assert assertions[2] == Function('exists')((c,), Eq(c, Mul(x, 2)))


def test_parse_macro_definition():
    source = """
    (set-logic LRA)
    (define-fun max ((a Real) (b Real)) Real (ite (> a b) a b))
    (declare-fun x () Real)
    (declare-fun y () Real)
    (assert (= (max x y) 10))
    """
    parsed_symbols, assertions = parse_smtlib(source)
    x = parsed_symbols['x']
    y = parsed_symbols['y']
    assert len(assertions) == 3
    assert assertions[0] == Q.real(x)
    assert assertions[1] == Q.real(y)
    assert assertions[2] == Eq(
        Piecewise((x, StrictGreaterThan(x, y)), (y, True)), 10
    )


def test_parse_bitvector_fallback():
    source = """
    (set-logic QF_BV)
    (declare-fun a () (_ BitVec 32))
    (declare-fun b () (_ BitVec 32))
    (assert (= (bvadd a b) #x00000000))
    """
    parsed_symbols, assertions = parse_smtlib(source)
    a = parsed_symbols['a']
    b = parsed_symbols['b']
    assert len(assertions) == 1
    assert assertions[0] == Eq(Function('bvadd')(a, b), 0)


def test_parse_lexer_edge_cases():
    # Test strings
    source = '''
    (declare-fun |a b| () String)
    (assert (= |a b| "hello world"))
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    assert assertions[0] == Eq(parsed_symbols['|a b|'], Symbol('"hello world"'))

    # Unclosed string error
    with raises(SMTLibSyntaxError):
        parse_smtlib('(assert (= x "unclosed))')

    # Unclosed quoted symbol error
    with raises(SMTLibSyntaxError):
        parse_smtlib('(declare-fun |unclosed () S)')


def test_parse_hex_and_binary():
    source = '''
    (assert (= #x0A 10))
    (assert (= #b1010 10))
    '''
    _, assertions = parse_smtlib(source)
    assert assertions[0] == Eq(10, 10)
    assert assertions[1] == Eq(10, 10)


def test_parse_type_and_errors():
    source = '''
    (declare-sort A 0)
    (declare-fun y () A)
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    assert 'y' in parsed_symbols

    # Just coverage for invalid arity
    with raises(SMTLibSyntaxError):
        parse_smtlib('(declare-sort B 1 2)')


def test_parse_smtlib_file():
    import os

    source = '''
    (declare-const x Int)
    (assert (= x 10))
    '''

    temp_path = "test.smt2"
    with open(temp_path, "w", encoding="utf-8") as f:
        f.write(source)

    parsed_symbols, assertions = parse_smtlib(temp_path)
    assert len(assertions) == 2
    assert assertions[-1] == Eq(parsed_symbols['x'], 10)

    os.remove(temp_path)


def test_parse_recursive_functions():
    source = '''
    (define-fun-rec f ((x Int)) Int (ite (> x 0) 1 0))
    '''
    with raises(NotImplementedError):
        parse_smtlib(source)


def test_parse_unknown_and_syntax_errors():
    with raises(UnknownSMTLibCommandError):
        parse_smtlib('(assert (unknown-op x y))')

    # Extra or missing parentheses
    with raises(SMTLibSyntaxError):
        parse_smtlib('(assert (= x 1)))')


def test_parse_extra_coverage():
    # Test division with multiple params and 'false' bool
    source1 = '''
    (declare-fun x () Real)
    (declare-fun y () Real)
    (declare-fun z () Real)
    (assert (= (/ x y z) 1.0))
    (assert false)
    '''
    parsed_symbols, assertions = parse_smtlib(source1)
    x = parsed_symbols['x']
    y = parsed_symbols['y']
    z = parsed_symbols['z']
    assert len(assertions) == 5
    assert assertions[0] == Q.real(x)
    assert assertions[1] == Q.real(y)
    assert assertions[2] == Q.real(z)
    assert assertions[3] == Eq(x / (y * z), 1.0)
    assert assertions[4] is False

    # Test and with 1 param and => with nested
    source2 = '''
    (assert (and true))
    '''
    _, assertions = parse_smtlib(source2)
    assert assertions[0] == And(True)


def test_parse_multi_arg_and_unary():
    source = '''
    (declare-fun x () Int)
    (declare-fun y () Int)
    (declare-fun z () Int)
    (assert (= x y z))
    (assert (distinct x y z))
    (assert (xor x y z))
    (assert (= (- x) (abs y)))
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    x, y, z = parsed_symbols['x'], parsed_symbols['y'], parsed_symbols['z']

    assert len(assertions) == 7
    assert assertions[0] == Q.integer(x)
    assert assertions[1] == Q.integer(y)
    assert assertions[2] == Q.integer(z)
    assert assertions[3] == And(Eq(x, y), Eq(y, z))
    assert assertions[4] == And(Ne(x, y), Ne(x, z), Ne(y, z))
    assert assertions[5] == Xor(x, y, z)
    assert assertions[6] == Eq(-x, Abs(y))


def test_parse_state_and_config():
    source = '''
    (set-info :smt-lib-version 2.6)
    (set-logic QF_LIA)
    (declare-fun x () Int)
    (push 1)
    (assert (> x 0))
    (pop 1)
    (get-info :name)
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    x = parsed_symbols['x']
    assert len(assertions) == 1
    assert assertions[0] == Q.integer(x)


def test_parse_reset():
    source = '''
    (declare-fun x () Int)
    (assert (= x 10))
    (reset-assertions)
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    assert len(assertions) == 0  # reset-assertions clears assertions
    assert 'x' in parsed_symbols # but keeps parsed_symbols

    source_reset = '''
    (declare-fun y () Int)
    (reset)
    '''
    parsed_symbols_r, assertions_r = parse_smtlib(source_reset)
    assert len(assertions_r) == 0
    assert 'y' not in parsed_symbols_r  # reset clears everything


def test_parse_match_fallback():
    with raises(NotImplementedError):
        parse_smtlib('(assert (match 10 ((x 20))))')


def test_parse_semantic_variadic_implies():
    source = '''
    (declare-sort S 0)
    (declare-fun a () S)
    (declare-fun b () S)
    (declare-fun c () S)
    (assert (=> a b c))
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    a, b, c = parsed_symbols['a'], parsed_symbols['b'], parsed_symbols['c']
    assert assertions[0] == Implies(a, Implies(b, c))


def test_parse_semantic_chained_comparisons():
    source = '(assert (< 1 2 0))'
    parsed_symbols, assertions = parse_smtlib(source)
    assert assertions[0] == And(StrictLessThan(1, 2), StrictLessThan(2, 0))


def test_parse_semantic_simultaneous_macro():
    source = '''
    (declare-fun a () Int)
    (declare-fun b () Int)
    (define-fun f ((a Int)(b Int)) Int (- a b))
    (assert (= (f b a) 0))
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    a, b = parsed_symbols['a'], parsed_symbols['b']
    assert len(assertions) == 3
    assert assertions[0] == Q.integer(a)
    assert assertions[1] == Q.integer(b)
    assert assertions[2] == Eq(b - a, 0)


def test_parse_semantic_simultaneous_let():
    source = '''
    (declare-fun a () Int)
    (declare-fun b () Int)
    (assert (let ((b a) (a 1)) (= b 1)))
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    a, b = parsed_symbols['a'], parsed_symbols['b']

    assert len(assertions) == 3
    assert assertions[0] == Q.integer(a)
    assert assertions[1] == Q.integer(b)
    assert assertions[2] == Eq(a, 1)


def test_parse_indexed_and_as_identifiers():
    source = '''
    (declare-sort S 0)
    (declare-fun x () S)
    (assert (= ((_ extract 7 0) x) 0))
    '''
    parsed_symbols, assertions = parse_smtlib(source)
    x = parsed_symbols['x']
    extract_func = Function("_extract_7_0")
    assert assertions[0] == Eq(extract_func(x), 0)

    source_as = '''
    (assert (= ((as const Array) 1) 1))
    '''
    parsed_symbols, assertions = parse_smtlib(source_as)
    const_func = Function("const")
    assert assertions[0] == Eq(const_func(1), 1)


def test_parse_complex_nested_logic():
    source = """
    (set-option :print-success false)
    (set-logic QF_LRA)
    (declare-fun x1 () Real)
    (declare-fun x3 () Real)
    (declare-fun x2 () Real)
    (assert (not (=> (and (<= x3 (+ 10 x2)) (and (and (>= x3 1) (and (= 10 x3) (= 0 x1))) (<= x3 12))) (or (>= 10 x3) (= x1 3)))))
    (check-sat)
    (exit)
    """
    parsed_symbols, assertions = parse_smtlib(source)
    x1 = parsed_symbols['x1']
    x2 = parsed_symbols['x2']
    x3 = parsed_symbols['x3']

    assert len(assertions) == 4
    assert assertions[0] == Q.real(x1)
    assert assertions[1] == Q.real(x3)
    assert assertions[2] == Q.real(x2)
    assert assertions[3] == Not(Implies(And(Eq(0, x1), Eq(10, x3), GreaterThan(x3, 1), LessThan(x3, 12), LessThan(x3, x2 + 10)), Or(Eq(x1, 3), GreaterThan(10, x3))))


def test_parse_complex_nested_arithmetic():
    source = """
    (set-option :print-success false)
    (set-logic QF_LRA)
    (declare-fun x1 () Real)
    (declare-fun x3 () Real)
    (declare-fun x2 () Real)
    (assert (not (=> (and (and (and (>= x3 1) (and (= 5 x3) (= 2 x2))) (<= x3 12)) (>= x3 (- 5 (* 2 x1)))) (or (>= x3 5) (= x2 1)))))
    (check-sat)
    (exit)
    """
    parsed_symbols, assertions = parse_smtlib(source)
    x1 = parsed_symbols['x1']
    x2 = parsed_symbols['x2']
    x3 = parsed_symbols['x3']

    assert len(assertions) == 4
    assert assertions[0] == Q.real(x1)
    assert assertions[1] == Q.real(x3)
    assert assertions[2] == Q.real(x2)
    assert assertions[3] == Not(Implies(And(Eq(2, x2), Eq(5, x3), GreaterThan(x3, 1), LessThan(x3, 12), GreaterThan(x3, 5 - 2*x1)), Or(Eq(x2, 1), GreaterThan(x3, 5))))


def test_parse_2_args():
    source = """
    (set-logic QF_LRA)
    (declare-fun x () Real)
    (declare-fun y () Real)
    (assert (distinct x y))
    (assert (= (/ x y) 10.0))
    """
    parsed_symbols, assertions = parse_smtlib(source)
    x = parsed_symbols['x']
    y = parsed_symbols['y']

    assert len(assertions) == 4
    assert assertions[0] == Q.real(x)
    assert assertions[1] == Q.real(y)
    assert assertions[2] == Ne(x, y)
    assert assertions[3] == Eq(x / y, 10.0)


def test_parse_variadic_inequalities():
    source = """
    (set-logic QF_LRA)
    (declare-fun x () Real)
    (declare-fun y () Real)
    (declare-fun z () Real)
    (assert (<= x y z))
    (assert (> x y z))
    (assert (>= x y z))
    """
    parsed_symbols, assertions = parse_smtlib(source)
    x = parsed_symbols['x']
    y = parsed_symbols['y']
    z = parsed_symbols['z']

    assert len(assertions) == 6
    assert assertions[0] == Q.real(x)
    assert assertions[1] == Q.real(y)
    assert assertions[2] == Q.real(z)
    assert assertions[3] == And(LessThan(x, y), LessThan(y, z))
    assert assertions[4] == And(StrictGreaterThan(x, y), StrictGreaterThan(y, z))
    assert assertions[5] == And(GreaterThan(x, y), GreaterThan(y, z))
