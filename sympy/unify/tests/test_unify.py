from sympy.unify.core import Compound, Variable, CondVariable, allcombinations, subst
from sympy.unify import core
from typing import Any


a, b, c = 'a', 'b', 'c'
w, x, y, z = map(Variable, 'wxyz')


class Function:
    """A helper for building compound term.
    This should be added in the core for convenience
    """
    def __init__(self, name: str):
        self.name = name

    def __call__(self, *args: Any) -> Compound:
        return Compound(self.name, args)


Add = Function('Add')
Mul = Function('Mul')
A = Function('A')
C = Function('C')
AC = Function('AC')


def is_associative(x) -> bool:
    return isinstance(x, Compound) and (x.op in ('A', 'AC'))


def is_commutative(x) -> bool:
    return isinstance(x, Compound) and (x.op in ('C', 'AC'))


def unify(a: Any, b: Any, s={}):
    return core.unify(
        a, b, s=s,
        is_associative=is_associative,
        is_commutative=is_commutative
    )


def test_basic():
    assert list(unify(a, x, {})) == [{x: a}]
    assert list(unify(a, x, {x: 10})) == []
    assert list(unify(1, x, {})) == [{x: 1}]
    assert list(unify(a, a, {})) == [{}]
    assert list(unify(x, (a, b), {})) == [{x: (a, b)}]
    assert list(unify((a, (b, c)), (a, (x, y)), {})) == [{x: b, y: c}]

    assert list(unify((a, b), (x, x), {})) == []
    assert list(unify((y, z), (x, x), {})) == [{y: x, z: x}]

    assert list(unify((w, x), (y, z), {})) == [{w: y, x: z}]
    assert list(unify((w, x), (w, x), {})) == [{}]


def test_operators():
    assert list(unify(Add(a, b, c), Add(a, x, y), {})) == [{x: b, y: c}]
    assert list(unify(Add(Mul(1, 2), b, c), Add(x, y, c), {})) == [{x: Mul(1, 2), y: b}]


def test_associative():
    assert list(unify(A(1, 2, 3), A(x, y), {})) == [
        {x: 1, y: A(2, 3)},
        {x: A(1, 2), y: 3}
    ]


def test_commutative():
    assert list(unify(C(1, 2, 3), C(x, y, z), {})) == [
        {x: 1, y: 2, z: 3},
        {x: 1, y: 3, z: 2},
        {x: 2, y: 1, z: 3},
        {x: 2, y: 3, z: 1},
        {x: 3, y: 1, z: 2},
        {x: 3, y: 2, z: 1},
    ]


def test_associative_commutative():
    assert list(unify(AC(1, 2, 3), AC(x, y), {})) == [
        {x: 1, y: AC(2, 3)},
        {x: AC(1, 2), y: 3},
        {x: 1, y: AC(3, 2)},
        {x: AC(1, 3), y: 2},
        {x: 2, y: AC(1, 3)},
        {x: AC(2, 1), y: 3},
        {x: 2, y: AC(3, 1)},
        {x: AC(2, 3), y: 1},
        {x: 3, y: AC(1, 2)},
        {x: AC(3, 1), y: 2},
        {x: 3, y: AC(2, 1)},
        {x: AC(3, 2), y: 1}
    ]


def test_allcombinations():
    assert set(allcombinations((1, 2), (1, 2), 'commutative')) ==\
        {(((1,), (2,)), ((1,), (2,))), (((1,), (2,)), ((2,), (1,)))}


def test_conditional():
    expr = Add(1, 2)
    x = Variable('x')
    y = CondVariable('y', lambda a: a % 2 == 0)
    z = CondVariable('z', lambda a: a > 3)

    pattern = Add(x, y)
    assert list(unify(expr, pattern, {})) == [{x: 1, y: 2}]
    pattern = Add(z, y)
    assert list(unify(expr, pattern, {})) == []


def test_defaultdict():
    assert next(unify(Variable('x'), 'foo')) == {Variable('x'): 'foo'}


def test_issue_21917():
    x = Variable('x')
    y = Variable('y')
    a = Compound('a', ())

    A = Compound('f', (x, y))
    B = Compound('f', (Compound('h', (a,)), x))
    for s in unify(A, B):
        assert subst(s, A) == subst(s, B)

    A = Compound('f', (x, y))
    B = Compound('f', (y, x))
    for s in unify(A, B):
        assert subst(s, A) == subst(s, B)

    A = Compound('f', (x, y))
    B = Compound('f', (Compound('f', (y,)), Compound('f', (x,))))
    for s in unify(A, B):
        assert subst(s, A) == subst(s, B)

    z = Variable('z')
    g_a = Compound('g', (a,))
    g_x = Compound('g', (x,))
    g_y = Compound('g', (y,))
    g_z = Compound('g', (z,))
    g_g_x = Compound('g', (g_x,))

    A = Compound('f', (x, g_a, g_z))
    B = Compound('f', (g_y, g_y, g_g_x))
    for s in unify(A, B):
        assert subst(s, A) == subst(s, B)
