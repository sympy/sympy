from sympy import symbols, Function, Eq, Wild
from sympy_invariant_density import Presentation

def test_unit_redundancy():
    op = Function("op", commutative=False)
    e, x = symbols("e x")
    a, b = Wild("a"), Wild("b")

    relations = [
        Eq(op(e, x), x),
        Eq(op(x, e), x),
        Eq(op(a, b), op(b, a)),
    ]

    P = Presentation(generators=[op, e], relations=relations)

    assert P.is_redundant(1) is True
