from sympy.core.logic import fuzzy_not, name_not, Logic, And, Or, Not
from sympy.utilities.pytest import raises

T = True
F = False
U = None

def test_fuzzy_not():
    assert fuzzy_not(T) == F
    assert fuzzy_not(F) == T
    assert fuzzy_not(U) == U


def test_name_not():
    assert name_not('zero')  == '!zero'
    assert name_not('!zero') == 'zero'


def test_logic_cmp():
    l1 = And('a', Not('b'))
    l2 = And('a', Not('b'))

    assert hash(l1) == hash(l2)
    assert (l1==l2) == True
    assert (l1!=l2) == False
    assert cmp(l1, l2) == 0

    assert And('a','b','c') == And('b','a','c')
    assert And('a','b','c') == And('c','b','a')
    assert And('a','b','c') == And('c','a','b')


def test_logic_onearg():
    raises(TypeError, 'And()')
    raises(TypeError, 'Or ()')

    assert And(T)   == T
    assert And(F)   == F
    assert Or (T)   == T
    assert Or (F)   == F

    assert And('a') == 'a'
    assert Or ('a') == 'a'



def test_logic_xnotx():
    assert And('a', '!a') == F
    assert Or ('a', '!a') == T


def test_logic_eval_TF():
    assert And(F, F)   == F
    assert And(F, T)   == F
    assert And(T, F)   == F
    assert And(T, T)   == T

    assert Or (F, F)   == F
    assert Or (F, T)   == T
    assert Or (T, F)   == T
    assert Or (T, T)   == T

    assert And('a', T) == 'a'
    assert And('a', F) == F
    assert Or ('a', T) == T
    assert Or ('a', F) == 'a'


def test_logic_combine_args():
    assert And('a', 'b', 'a')   == And('a', 'b')
    assert Or ('a', 'b', 'a')   == Or ('a', 'b')

    assert And( And('a','b'), And('c','d') )    == And('a','b','c','d')
    assert Or ( Or ('a','b'), Or ('c','d') )    == Or ('a','b','c','d')

    assert Or( 't', And('n','p','r'), And('n','r'), And('n','p','r'), 't', And('n','r') ) == \
                    Or('t', And('n','p','r'), And('n','r'))


def test_logic_expand():
    t = And(Or('a','b'), 'c')
    assert t.expand()  == Or(And('a','c'), And('b','c'))

    t = And(Or('a','!b'),'b')
    assert t.expand()  == And('a','b')

    t = And(Or('a','b'), Or('c','d'))
    assert t.expand()  == Or(And('a','c'), And('a','d'), And('b','c'), And('b','d'))



def test_logic_fromstring():
    S = Logic.fromstring

    assert S('a')           == 'a'
    assert S('!a')          == '!a'
    assert S('a & b')       == And('a','b')
    assert S('a | b')       == Or ('a','b')
    assert S('a | b & c')   == And(Or ('a','b'), 'c')
    assert S('a & b | c')   == Or (And('a','b'), 'c')
    assert S('a & b & c')   == And('a','b','c')
    assert S('a | b | c')   == Or ('a','b','c')

    raises(ValueError, "S('| a')")
    raises(ValueError, "S('& a')")
    raises(ValueError, "S('a | | b')")
    raises(ValueError, "S('a | & b')")
    raises(ValueError, "S('a & & b')")
    raises(ValueError, "S('a |')")



def test_logic_not():
    assert Not('a')     == '!a'
    assert Not('!a')    == 'a'

    # NOTE: we may want to change default Not behaviour and put this
    # functionality into some method.
    assert Not(And('a','b'))  == Or ('!a','!b')
    assert Not(Or ('a','b'))  == And('!a','!b')
