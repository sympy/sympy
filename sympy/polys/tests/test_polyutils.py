"""Tests for useful utilities for higher level polynomial classes. """

from sympy import S, I, sqrt, symbols

from sympy.polys.polyutils import (
    _sort_gens,
    _unify_gens,
    _analyze_power,
    _dict_from_basic_if_gens,
    _dict_from_basic_no_gens,
    dict_from_basic,
    basic_from_dict,
)

from sympy.polys.polyerrors import (
    GeneratorsNeeded,
)

x,y,z,p,q,r,s,t,u,v,w = symbols('x,y,z,p,q,r,s,t,u,v,w')

def test__sort_gens():
    assert _sort_gens([]) == []

    assert _sort_gens([x]) == [x]
    assert _sort_gens([p]) == [p]
    assert _sort_gens([q]) == [q]

    assert _sort_gens([x,p]) == [x,p]
    assert _sort_gens([p,x]) == [x,p]
    assert _sort_gens([q,p]) == [p,q]

    assert _sort_gens([q,p,x]) == [x,p,q]

    assert _sort_gens([x,p,q], wrt='x') == [x,p,q]
    assert _sort_gens([x,p,q], wrt='p') == [p,x,q]
    assert _sort_gens([x,p,q], wrt='q') == [q,x,p]

    assert _sort_gens([x,p,q], sort='x < p < q') == [x, p, q]
    assert _sort_gens([x,p,q], sort='p < x < q') == [p, x, q]
    assert _sort_gens([x,p,q], sort='p < q < x') == [p, q, x]

    assert _sort_gens([x,p,q], wrt='x', sort='q < p') == [x, q, p]
    assert _sort_gens([x,p,q], wrt='p', sort='q < x') == [p, q, x]
    assert _sort_gens([x,p,q], wrt='q', sort='p < x') == [q, p, x]

def test__unify_gens():
    assert _unify_gens([], []) == []

    assert _unify_gens([x], [x]) == [x]
    assert _unify_gens([y], [y]) == [y]

    assert _unify_gens([x,y], [x]) == [x,y]
    assert _unify_gens([x], [x,y]) == [x,y]

    assert _unify_gens([x,y], [x,y]) == [x,y]
    assert _unify_gens([y,x], [y,x]) == [y,x]

    assert _unify_gens([x], [y]) == [x,y]
    assert _unify_gens([y], [x]) == [y,x]

    assert _unify_gens([x], [y,x]) == [y,x]
    assert _unify_gens([y,x], [x]) == [y,x]

    assert _unify_gens([x,y,z], [x,y,z]) == [x,y,z]
    assert _unify_gens([z,y,x], [x,y,z]) == [z,y,x]
    assert _unify_gens([x,y,z], [z,y,x]) == [x,y,z]
    assert _unify_gens([z,y,x], [z,y,x]) == [z,y,x]

    assert _unify_gens([x,y,z], [t,x,p,q,z]) == [t,x,y,p,q,z]

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

