from sympy.logpy import refine_one
from logpy import Relation, facts, fact
from sympy import Basic, Q, Dummy, Symbol, Abs
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

def test_basic():
    assert refine(Basic(5)) == Basic(3)

def test_exprs_simple():
    assert refine(Abs(y)) == Abs(y)
    assert refine(Abs(y), Q.positive(y)) == y

def test_exprs_pow_of_abs():
    assert refine(Abs(y)**4, Q.real(y)) == y**4
    assert refine(Abs(y*z)**4, Q.real(y) & Q.real(z)) == (y*z)**4

def test_pow_of_abs_deep():
    with assuming(Q.real(y)):
        assert refine_deep(2*Abs(y)**4) == 2*y**4
        assert refine_deep(Abs(2*y)**4) == (2*y)**4

def test_commutativity():
    from logpy.assoccomm import op_registry
    from logpy.assoccomm import eq_assoccomm as eqac
    from logpy.core import goaleval
    from logpy.variables import variables
    op_registry.append((lambda x: isinstance(x, type) and issubclass(x, Basic),
                        lambda x: isinstance(x, Basic),
                        type,
                        lambda x: tuple(x.args),
                        lambda op, args: op.func(*args)))

    y = Symbol('_tmp')

    with variables(x):
        print next(goaleval(eqac(y+x, y+1))({}))
        assert next(goaleval(eqac(y+x, y+1))({})) == {x: 1}

    reduces = Relation('reduces')
    fact(reduces, x * y, x, True)
    refine = partial(refine_one, reduces=reduces, vars=[x])

    assert refine(3) == 3
    print refine(3*y)
    assert refine(3*y) == 3
    print refine(z*3*y)
    assert refine(z*3*y) == z*3
    assert refine(z*y*3) == z*3
    assert refine(y*z*3) == z*3

    del op_registry[-1]
