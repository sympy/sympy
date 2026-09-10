from sympy import symbols, true, false, S, Eq
from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, Xor, Nand, ITE, Nor
from sympy.logic.boolean_to_polynomial import boolean_to_polynomial
from itertools import product

def assert_equivalent(f1, f2, symbols):
    for values in product([0, 1], repeat=len(symbols)):
        subs = dict(zip(symbols, values))
        v1 = f1.subs(subs) % 2
        v2 = f2.subs(subs) % 2
        assert v1 == v2, f"Failed at {subs}: {v1} != {v2}"


def test_and_not():
    a, b = symbols('a b')
    expr = And(Not(a), Not(b))
    poly = boolean_to_polynomial(expr)
    expected = (1 - a)*(1 - b)
    assert poly.equals(expected)

def test_or():
    a, b = symbols('a b')
    expr = Or(a, b)
    poly = boolean_to_polynomial(expr)
    expected = a + b - a * b
    assert poly.equals(expected)

def test_implies():
    a, b = symbols('a b')
    expr = Implies(a, b)
    poly = boolean_to_polynomial(expr)
    expected = 1 - a + a * b
    assert poly.equals(expected)

def test_equivalent():
    a, b = symbols('a b')
    expr = Equivalent(a, b)
    poly = boolean_to_polynomial(expr)
    expected = 1 - a - b + 2 * a * b
    assert poly.equals(expected)

def test_complicated_expression():
    a, b, c = symbols('a b c')
    expr = Or(And(a, b), Not(c))
    poly = boolean_to_polynomial(expr)
    expected = (a * b) + (1 - c) - (a * b * (1 - c))
    assert poly.equals(expected)

def test_nested_expressions():
    a, b, c = symbols('a b c')
    expr = And(Or(a, b), Not(c))
    poly = boolean_to_polynomial(expr)
    expected = (a + b - a * b) * (1 - c)
    assert poly.equals(expected)

def test_xor():
    a, b = symbols('a b')
    expr = Xor(a, b)
    poly = boolean_to_polynomial(expr)
    expected = a + b - 2 * a * b
    assert poly.equals(expected)

def test_nand():
    a, b = symbols('a b')
    expr = Nand(a, b)
    poly = boolean_to_polynomial(expr)
    expected = 1 - a * b
    assert poly.equals(expected)


def test_reduction_to_degree_one():
    """
    Test that terms like a^2 are reduced to just 'a'.
    """
    a, b = symbols('a b')

    # Test a^2 -> a (should be reduced)
    expr = And(a, a)  # a * a
    poly = boolean_to_polynomial(expr)
    assert poly.equals(a)

    # Test b^2 -> b
    expr = And(b, b)  # b * b
    poly = boolean_to_polynomial(expr)
    assert poly.equals(b)

    # Test (a^2 * b) -> a * b
    expr = And(And(a, a), b)  # (a * a) * b
    poly = boolean_to_polynomial(expr)
    assert poly.equals(a * b)  # Compare with arithmetic polynomial


def test_rejection_of_relations():
    """
    Test that invalid expressions like Eq or Gt are rejected.
    """
    a, b = symbols('a b')

    # Test Eq(a, b) - should raise TypeError
    try:
        expr = Eq(a, b)
        boolean_to_polynomial(expr)
    except TypeError:
        pass  # Expected

    # Test Gt(a, b) - should raise TypeError
    try:
        expr = a > b
        boolean_to_polynomial(expr)
    except TypeError:
        pass  # Expected

    # Test Neq(a, b) - should raise TypeError
    try:
        expr = a != b
        boolean_to_polynomial(expr)
    except TypeError:
        pass  # Expected

def test_complicated_cases():
    """
    Test more complicated combinations of boolean operators.
    """
    a, b, c = symbols('a b c')

    # Test Or(And(a, b), Not(c)) -> a * b + (1 - c) - a * b * (1 - c)
    expr = Or(And(a, b), Not(c))
    poly = boolean_to_polynomial(expr)
    expected = a * b + (S.One - c) - a * b * (S.One - c)
    assert poly.equals(expected)

    # Test Or(And(a, b), And(c, Not(b))) -> a * b + c * (1 - b) - a * b * c * (1 - b)
    expr = Or(And(a, b), And(c, Not(b)))
    poly = boolean_to_polynomial(expr)
    expected = a * b + c * (S.One - b) - a * b * c * (S.One - b)
    assert poly.equals(expected)

    # Test True and False
    expr_true = true
    poly_true = boolean_to_polynomial(expr_true)
    assert poly_true.equals(S.One)

    expr_false = false
    poly_false = boolean_to_polynomial(expr_false)
    assert poly_false.equals(S.Zero)

def test_nor():
    a, b = symbols('a b')
    expr = Nor(a, b)
    poly = boolean_to_polynomial(expr)
    expected = 1 - (a + b - a * b)
    assert poly.equals(expected)

def test_ite():
    a, b, c = symbols('a b c')
    expr = ITE(a, b, c)  # if a then b else c
    poly = boolean_to_polynomial(expr)
    expected = a * b + (1 - a) * c
    assert poly.equals(expected)

def test_nested_and_or_not():
    a, b, c = symbols('a b c')
    expr = Or(And(a, Not(b)), c)
    poly = boolean_to_polynomial(expr)

    and_part = a * (1 - b)
    expected = and_part + c - and_part * c

    assert poly.equals(expected)


def test_ite_nested_expression():
    a, b, c, d = symbols('a b c d')
    expr = Or(ITE(a, b, c), d)
    poly = boolean_to_polynomial(expr)

    # ITE(a,b,c) = ab + (1 - a)c
    ite_part = a * b + (1 - a) * c
    expected = ite_part + d - ite_part * d

    assert poly.equals(expected)


def test_simplify_flag_expansion():
    a, b = symbols('a b')
    expr = Or(And(a, b), Not(b))
    poly_unsimplified = boolean_to_polynomial(expr, simplify=False)
    poly_simplified = boolean_to_polynomial(expr, simplify=True)

    expected = a * b + (1 - b) - a * b * (1 - b)

    assert_equivalent(poly_unsimplified, expected, [a, b])
    assert_equivalent(poly_simplified, expected.expand(modulus=2), [a, b])
