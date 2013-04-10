from sympy.logpy import refine_one
from logpy import Relation, facts
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
