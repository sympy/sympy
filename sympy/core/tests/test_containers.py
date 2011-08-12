from sympy import Matrix, Tuple, symbols, sympify, Basic
from sympy.core.containers import tuple_wrapper
from sympy.utilities.pytest import raises
from sympy.core.compatibility import is_sequence, iterable

def test_Tuple():
    t = (1, 2, 3, 4)
    st =  Tuple(*t)
    assert set(sympify(t)) == set(st)
    assert len(t) == len(st)
    assert set(sympify(t[:2])) == set(st[:2])
    assert isinstance(st[:], Tuple)
    assert st == Tuple(1, 2, 3, 4)
    assert st.func(*st.args) == st
    p, q, r, s = symbols('p q r s')
    t2 = (p, q, r, s)
    st2 = Tuple(*t2)
    assert st2.atoms() == set(t2)
    assert st == st2.subs({p:1, q:2, r:3, s:4})
    # issue 2406
    assert all([ isinstance(arg, Basic) for arg in st.args ])
    assert Tuple(p, 1).subs(p, 0) == Tuple(0, 1)
    assert Tuple(p, Tuple(p, 1)).subs(p, 0) == Tuple(0, Tuple(0, 1))

    assert Tuple(t2) == Tuple(Tuple(*t2))
    assert Tuple.fromiter(t2) == Tuple(*t2)
    assert Tuple.fromiter(x for x in xrange(4)) == Tuple(0, 1, 2, 3)
    assert st2.fromiter(st2.args) == st2

def test_Tuple_contains():
    t1, t2 = Tuple(1), Tuple(2)
    assert t1 in Tuple(1, 2, 3, t1, Tuple(t2))
    assert t2 not in Tuple(1, 2, 3, t1, Tuple(t2))

def test_Tuple_concatenation():
    assert Tuple(1, 2) + Tuple(3, 4) == Tuple(1, 2, 3, 4)
    assert (1, 2) + Tuple(3, 4) == Tuple(1, 2, 3, 4)
    assert Tuple(1, 2) + (3, 4) == Tuple(1, 2, 3, 4)
    raises(TypeError, 'Tuple(1, 2) + 3')
    raises(TypeError, '1 + Tuple(2, 3)')

    #the Tuple case in __radd__ is only reached when a subclass is involved
    class Tuple2(Tuple):
        def __radd__(self, other):
            return Tuple.__radd__(self, other + other)
    assert Tuple(1, 2) + Tuple2(3, 4) == Tuple(1, 2, 1, 2, 3, 4)
    assert Tuple2(1, 2) + Tuple(3, 4) == Tuple(1, 2, 3, 4)

def test_Tuple_equality():
    assert Tuple(1, 2) is not (1, 2)
    assert (Tuple(1, 2) == (1, 2)) is True
    assert (Tuple(1, 2) != (1, 2)) is False
    assert (Tuple(1, 2) == (1, 3)) is False
    assert (Tuple(1, 2) != (1, 3)) is True
    assert (Tuple(1, 2) == Tuple(1, 2)) is True
    assert (Tuple(1, 2) != Tuple(1, 2)) is False
    assert (Tuple(1, 2) == Tuple(1, 3)) is False
    assert (Tuple(1, 2) != Tuple(1, 3)) is True

def test_tuple_wrapper():

    @tuple_wrapper
    def wrap_tuples_and_return(*t):
        return t

    p = symbols('p')
    assert wrap_tuples_and_return(p, 1) == (p, 1)
    assert wrap_tuples_and_return((p, 1)) == (Tuple(p, 1),)
    assert wrap_tuples_and_return(1, (p, 2), 3) == (1, Tuple(p, 2), 3)

def test_iterable_is_sequence():
    ordered = [list(), tuple(), Tuple(), Matrix([[]])]
    unordered = [set()]
    not_sympy_iterable = [{}, '', u'']
    assert all(is_sequence(i) for i in ordered)
    assert all(not is_sequence(i) for i in unordered)
    assert all(iterable(i) for i in ordered + unordered)
    assert all(not iterable(i) for i in not_sympy_iterable)
    assert all(iterable(i, exclude=None) for i in not_sympy_iterable)
