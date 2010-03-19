"""Tests for useful utilities for higher level polynomial classes. """

from sympy import S, I, Integer, sqrt, symbols, pi
from sympy.utilities.pytest import raises

from sympy.polys.polyutils import (
    _sort_gens,
    _unify_gens,
    _analyze_gens,
    _sort_factors,
    _analyze_power,
    _dict_from_basic_if_gens,
    _dict_from_basic_no_gens,
    dict_from_basic,
    basic_from_dict,
)

from sympy.polys.polyerrors import (
    GeneratorsNeeded,
    PolynomialError,
)

from sympy.polys.algebratools import ZZ, QQ, EX

x,y,z,p,q,r,s,t,u,v,w = symbols('x,y,z,p,q,r,s,t,u,v,w')

def test__sort_gens():
    assert _sort_gens([]) == ()

    assert _sort_gens([x]) == (x,)
    assert _sort_gens([p]) == (p,)
    assert _sort_gens([q]) == (q,)

    assert _sort_gens([x,p]) == (x,p)
    assert _sort_gens([p,x]) == (x,p)
    assert _sort_gens([q,p]) == (p,q)

    assert _sort_gens([q,p,x]) == (x,p,q)

    assert _sort_gens([x,p,q], wrt='x') == (x,p,q)
    assert _sort_gens([x,p,q], wrt='p') == (p,x,q)
    assert _sort_gens([x,p,q], wrt='q') == (q,x,p)

    assert _sort_gens([x,p,q], sort='x > p > q') == (x, p, q)
    assert _sort_gens([x,p,q], sort='p > x > q') == (p, x, q)
    assert _sort_gens([x,p,q], sort='p > q > x') == (p, q, x)

    assert _sort_gens([x,p,q], wrt='x', sort='q > p') == (x, q, p)
    assert _sort_gens([x,p,q], wrt='p', sort='q > x') == (p, q, x)
    assert _sort_gens([x,p,q], wrt='q', sort='p > x') == (q, p, x)

def test__unify_gens():
    assert _unify_gens([], []) == ()

    assert _unify_gens([x], [x]) == (x,)
    assert _unify_gens([y], [y]) == (y,)

    assert _unify_gens([x,y], [x]) == (x,y)
    assert _unify_gens([x], [x,y]) == (x,y)

    assert _unify_gens([x,y], [x,y]) == (x,y)
    assert _unify_gens([y,x], [y,x]) == (y,x)

    assert _unify_gens([x], [y]) == (x,y)
    assert _unify_gens([y], [x]) == (y,x)

    assert _unify_gens([x], [y,x]) == (y,x)
    assert _unify_gens([y,x], [x]) == (y,x)

    assert _unify_gens([x,y,z], [x,y,z]) == (x,y,z)
    assert _unify_gens([z,y,x], [x,y,z]) == (z,y,x)
    assert _unify_gens([x,y,z], [z,y,x]) == (x,y,z)
    assert _unify_gens([z,y,x], [z,y,x]) == (z,y,x)

    assert _unify_gens([x,y,z], [t,x,p,q,z]) == (t,x,y,p,q,z)

def test__analyze_gens():
    assert _analyze_gens((x,y,z)) == (x,y,z)
    assert _analyze_gens([x,y,z]) == (x,y,z)

    assert _analyze_gens(([x,y,z],)) == (x,y,z)
    assert _analyze_gens(((x,y,z),)) == (x,y,z)

def test__sort_factors():
    assert _sort_factors([], multiple=True) == []
    assert _sort_factors([], multiple=False) == []

    F = [[1,2,3], [1,2], [1]]
    G = [[1], [1,2], [1,2,3]]

    assert _sort_factors(F, multiple=False) == G

    F = [[1,2], [1,2,3], [1,2], [1]]
    G = [[1], [1,2], [1,2], [1,2,3]]

    assert _sort_factors(F, multiple=False) == G

    F = [[2,2], [1,2,3], [1,2], [1]]
    G = [[1], [1,2], [2,2], [1,2,3]]

    assert _sort_factors(F, multiple=False) == G

    F = [([1,2,3], 1), ([1,2], 1), ([1], 1)]
    G = [([1], 1), ([1,2], 1), ([1,2,3], 1)]

    assert _sort_factors(F, multiple=True) == G

    F = [([1,2], 1), ([1,2,3], 1), ([1,2], 1), ([1], 1)]
    G = [([1], 1), ([1,2], 1), ([1,2], 1), ([1,2,3], 1)]

    assert _sort_factors(F, multiple=True) == G

    F = [([2,2], 1), ([1,2,3], 1), ([1,2], 1), ([1], 1)]
    G = [([1], 1), ([1,2], 1), ([2,2], 1), ([1,2,3], 1)]

    assert _sort_factors(F, multiple=True) == G

    F = [([2,2], 1), ([1,2,3], 1), ([1,2], 2), ([1], 1)]
    G = [([1], 1), ([2,2], 1), ([1,2], 2), ([1,2,3], 1)]

    assert _sort_factors(F, multiple=True) == G

def test__analyze_power():
    assert _analyze_power(x, S(1)) == (x, S(1))
    assert _analyze_power(x, S(2)) == (x, S(2))
    assert _analyze_power(x, -S(1)) == (x**(-1), S(1))
    assert _analyze_power(x, -S(2)) == (x**(-1), S(2))

    assert _analyze_power(x, S(1)/3) == (x**(S(1)/3), S(1))
    assert _analyze_power(x, S(2)/3) == (x**(S(1)/3), S(2))
    assert _analyze_power(x, -S(1)/3) == (x**(-S(1)/3), S(1))
    assert _analyze_power(x, -S(2)/3) == (x**(-S(1)/3), S(2))

    assert _analyze_power(x, y) == (x**y, S(1))
    assert _analyze_power(x, -y) == (x**(-y), S(1))
    assert _analyze_power(x, 2*y) == (x**y, S(2))
    assert _analyze_power(x, -2*y) == (x**(-y), S(2))

    assert _analyze_power(x, y/3) == (x**(y/3), S(1))
    assert _analyze_power(x, -y/3) == (x**(-y/3), S(1))
    assert _analyze_power(x, 2*y/3) == (x**(y/3), S(2))
    assert _analyze_power(x, -2*y/3) == (x**(-y/3), S(2))

    assert _analyze_power(x, S(1.0)) == (x**S(1.0), S(1))
    assert _analyze_power(x, S(2.0)) == (x**S(2.0), S(1))
    assert _analyze_power(x, -S(1.0)) == (x**(-S(1.0)), S(1))
    assert _analyze_power(x, -S(2.0)) == (x**(-S(2.0)), S(1))

    assert _analyze_power(x, S(1.0)*y) == (x**(S(1.0)*y), S(1))
    assert _analyze_power(x, S(2.0)*y) == (x**(S(2.0)*y), S(1))
    assert _analyze_power(x, -S(1.0)*y) == (x**(-S(1.0)*y), S(1))
    assert _analyze_power(x, -S(2.0)*y) == (x**(-S(2.0)*y), S(1))

def test__dict_from_basic_if_gens():
    assert _dict_from_basic_if_gens(Integer(17), (x,)) == {(0,): Integer(17)}
    assert _dict_from_basic_if_gens(Integer(17), (x,y)) == {(0,0): Integer(17)}
    assert _dict_from_basic_if_gens(Integer(17), (x,y,z)) == {(0,0,0): Integer(17)}

    assert _dict_from_basic_if_gens(Integer(-17), (x,)) == {(0,): Integer(-17)}
    assert _dict_from_basic_if_gens(Integer(-17), (x,y)) == {(0,0): Integer(-17)}
    assert _dict_from_basic_if_gens(Integer(-17), (x,y,z)) == {(0,0,0): Integer(-17)}

    assert _dict_from_basic_if_gens(Integer(17)*x, (x,)) == {(1,): Integer(17)}
    assert _dict_from_basic_if_gens(Integer(17)*x, (x,y)) == {(1,0): Integer(17)}
    assert _dict_from_basic_if_gens(Integer(17)*x, (x,y,z)) == {(1,0,0): Integer(17)}

    assert _dict_from_basic_if_gens(Integer(17)*x**7, (x,)) == {(7,): Integer(17)}
    assert _dict_from_basic_if_gens(Integer(17)*x**7*y, (x,y)) == {(7,1): Integer(17)}
    assert _dict_from_basic_if_gens(Integer(17)*x**7*y*z**12, (x,y,z)) == {(7,1,12): Integer(17)}

    assert _dict_from_basic_if_gens(x+2*y+3*z, (x,)) == \
        {(1,): Integer(1), (0,): 2*y+3*z}
    assert _dict_from_basic_if_gens(x+2*y+3*z, (x,y)) == \
        {(1,0): Integer(1), (0,1): Integer(2), (0,0): 3*z}
    assert _dict_from_basic_if_gens(x+2*y+3*z, (x,y,z)) == \
        {(1,0,0): Integer(1), (0,1,0): Integer(2), (0,0,1): Integer(3)}

    assert _dict_from_basic_if_gens(x*y+2*x*z+3*y*z, (x,)) == \
        {(1,): y+2*z, (0,): 3*y*z}
    assert _dict_from_basic_if_gens(x*y+2*x*z+3*y*z, (x,y)) == \
        {(1,1): Integer(1), (1,0): 2*z, (0,1): 3*z}
    assert _dict_from_basic_if_gens(x*y+2*x*z+3*y*z, (x,y,z)) == \
        {(1,1,0): Integer(1), (1,0,1): Integer(2), (0,1,1): Integer(3)}

    assert _dict_from_basic_if_gens(2**y*x, (x,)) == {(1,): 2**y}
    raises(PolynomialError, "_dict_from_basic_if_gens(2**y*x, (x,y))")

def test__dict_from_basic_no_gens():
    raises(GeneratorsNeeded, "_dict_from_basic_no_gens(Integer(17))")

    assert _dict_from_basic_no_gens(x) == ({(1,): Integer(1)}, (x,))
    assert _dict_from_basic_no_gens(y) == ({(1,): Integer(1)}, (y,))

    assert _dict_from_basic_no_gens(x*y) == ({(1,1): Integer(1)}, (x,y))
    assert _dict_from_basic_no_gens(x+y) == ({(1,0): Integer(1), (0,1): Integer(1)}, (x,y))

    assert _dict_from_basic_no_gens(sqrt(2)) == ({(1,): Integer(1)}, (sqrt(2),))
    raises(GeneratorsNeeded, "_dict_from_basic_no_gens(sqrt(2), greedy=False)")

    assert _dict_from_basic_no_gens(x*y, domain=ZZ[x]) == ({(1,): x}, (y,))
    assert _dict_from_basic_no_gens(x*y, domain=ZZ[y]) == ({(1,): y}, (x,))

    assert _dict_from_basic_no_gens(3*sqrt(2)*pi*x*y, extension=False) == ({(1,1,1,1): 3}, (x,y,sqrt(2),pi))
    assert _dict_from_basic_no_gens(3*sqrt(2)*pi*x*y, extension=True) == ({(1,1,1): 3*sqrt(2)}, (x,y,pi))

