from sympy.logpy.core import refine_one, asko
from logpy import Relation, facts, fact
from logpy.unify import unify, reify
from logpy.variables import variables
from sympy import Basic, Q, Dummy, Symbol, Abs, Mul, S
from functools import partial
from sympy.assumptions import assuming

from sympy.strategies.traverse import top_down

x, k = map(Dummy, 'xk')

reduces = Relation('reduces')

facts(reduces, *[
        (Basic(5), Basic(3), True),
        (Abs(x), x, Q.positive(x)),
        (Abs(x)**k, x**k, Q.real(x) & Q.even(k)),
    ])

vars = [x, k]

refine = partial(refine_one, vars=vars, reduces=reduces)

refine_deep = top_down(refine)


y = Symbol('y')
z = Symbol('z')

def test_op_args():
    from logpy.assoccomm import op_args
    from sympy import Add, S
    op, args = op_args(y+z)
    assert op == Add
    assert set(args) == set((y, z))

    assert op_args(S.One) == op_args(y) == (None, None)


def test_matching():
    with variables(*vars):
        assert unify(4, k, {})
    assert asko(Q.even(4), True)

def test_basic():
    assert refine(Basic(5)) == Basic(3)

def test_exprs_simple():
    assert refine(Abs(y)) == Abs(y)
    assert refine(Abs(y), Q.positive(y)) == y

def rebuild(obj):
    return obj.func(*obj.args)

def test_exprs_pow_of_abs():
    assert refine(Abs(y)**4, Q.real(y)) == y**4
    assert rebuild(refine(Abs(y*z)**4, Q.real(y) & Q.real(z))) == (y*z)**4

def test_pow_of_abs_deep():
    with assuming(Q.real(y)):
        assert refine_deep(2*Abs(y)**4) == 2*y**4
        assert refine_deep(Abs(2*y)**4) == (2*y)**4

def test_commutativity():
    from logpy.assoccomm import eq_assoccomm as eqac
    from logpy.core import goaleval
    from logpy.variables import variables

    y = Symbol('_tmp')
    x = Symbol('x')

    reduces = Relation('reduces')
    fact(reduces, x * 2, x, True)
    refine = partial(refine_one, reduces=reduces, vars=[x])
    mul = partial(Mul, evaluate=False)

    assert refine(3) == 3
    assert refine(2*y) == y
    assert refine(mul(2, y)) == y
    assert refine(mul(y, 2)) == y
    assert refine(z*2*y) == z*y
    assert refine(z*y*2) == z*y
    assert refine(y*z*2) == z*y

    with variables(x):
        assert reify(x, next(goaleval(eqac(y+x, y+1))({}))) == 1
