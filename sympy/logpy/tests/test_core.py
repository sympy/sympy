from sympy.logpy import refine_one
from logpy import Relation, facts
from sympy import Basic
from functools import partial

reduces = Relation('reduces')

facts(reduces, *[(Basic(5), Basic(3), True)])

vars = []

refine = partial(refine_one, vars=vars, reduces=reduces)

def test_basic():
    assert refine(Basic(5), [True]) == Basic(3)
