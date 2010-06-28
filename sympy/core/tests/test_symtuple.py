from sympy import SymTuple, symbols

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

