from sympy.logpy import refine_one
from logpy import Relation, facts
from sympy import Basic, Q, Dummy, Symbol, Abs
from functools import partial

x = Dummy('x')

reduces = Relation('reduces')

facts(reduces, *[
        (Basic(5), Basic(3), True),
        (Abs(x), x, Q.positive(x))

    ])

vars = [x]

refine = partial(refine_one, vars=vars, reduces=reduces)



def test_basic():
    assert refine(Basic(5), [True]) == Basic(3)

def test_exprs_simple():
    y = Symbol('y')
    assert refine(Abs(y), []) == Abs(y)
    assert refine(Abs(y), [Q.positive(y)]) == y
