from sympy import Basic, Q, Dummy, Symbol, Abs, Add, Mul, S, assuming
from sympy.strategies.traverse import top_down
from functools import partial

from sympy.external import import_module
logpy = import_module('logpy')
if logpy:
    from sympy.logpy.core import refine_one, asko
    from logpy import Relation, facts, fact, run, unify, reify
    from logpy.variables import variables
    from logpy.core import goaleval, eq
    from logpy.assoccomm import op_args, build, buildo, commutative
else:
    #bin/test will not execute any tests now
    disabled = True

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


x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

def test_as_from_tuple():
    examples = [x, S(2), S(3.0), S.One, x+5]
    assert all(e._from_tuple(e._as_tuple()) == e for e in examples)

def test_op_args():
    op, args = op_args(y+z)
    assert op == Add
    assert set(args) == set((y, z))

    op, args = op_args(y)
    assert op == Symbol
    assert 'y' in args

def test_build():
    assert build(Add, (y, z)) == Add(y, z, evaluate=False)

def test_buildo():
    with variables(x):
        result = next(goaleval(buildo(Add, (y, z), x))({}))
        expected = {x: Add(y, z, evaluate=False)}
        assert result == expected

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

add = partial(Add, evaluate=False)
def test_commutativity():
    from logpy.assoccomm import eq_comm
    assert (Add,) in commutative.facts

    with variables(x):
        assert set(run(0, x, eq_comm(x, add(y, 1)))) == {add(y, 1), add(1, y)}

def test_commutativity_refine():
    from logpy.assoccomm import eq_assoccomm as eqac
    from logpy.variables import variables

    y = Symbol('y')
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
        assert reify(x, next(goaleval(eq(x, y))({}))) == y
        assert reify(x, next(goaleval(eq(x, 1))({}))) == 1
        assert reify(x, next(goaleval(eqac(add(x, y), add(1, y)))({}))) == 1
        assert reify(x, next(goaleval(eqac(add(x, y), add(y, 1)))({}))) == 1
