from sympy import SymTuple, symbols
from sympy.core.symtuple import tuple_wrapper

def test_SymTuple():
    t = (1, 2, 3, 4)
    st =  SymTuple(*t)
    assert set(t) == set(st)
    assert len(t) == len(st)
    assert set(t[:2]) == set(st[:2])
    assert isinstance(st[:], SymTuple)
    assert st == SymTuple(1, 2, 3, 4)
    assert st.func(*st.args) == st
    p, q, r, s = symbols('p q r s')
    t2 = (p, q, r, s)
    st2 = SymTuple(*t2)
    assert st2.atoms() == set(t2)
    assert st == st2.subs({p:1, q:2, r:3, s:4})

def test_tuple_wrapper():

    @tuple_wrapper
    def wrap_tuples_and_return(*t):
        return t

    p = symbols('p')
    assert wrap_tuples_and_return(p, 1) == (p, 1)
    assert wrap_tuples_and_return((p, 1)) == (SymTuple(p, 1),)
    assert wrap_tuples_and_return(1, (p, 2), 3) == (1, SymTuple(p, 2), 3)
