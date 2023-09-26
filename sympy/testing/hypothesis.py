from hypothesis import strategies as st
from hypothesis import assume
from sympy.abc import x
from sympy.polys.domains import ZZ, QQ
from sympy.polys.polytools import Poly, lcm, gcd
from typing import Any


@st.composite
def coefficients(draw: st.DrawFn, empty=True, degree=None, rational=False):
    min_size = 0 if empty else 1
    elements_strat = Any
    if rational:
        group = QQ
        elements_strat = st.fractions()
    else:
        group = ZZ
        elements_strat = st.integers()
    # converting to the underlying type
    raw_l = draw(st.lists(elements_strat, min_size=min_size, max_size=degree))
    l = [group(i) for i in raw_l]
    # ensuring there is no leading zero or fully zero list.
    if len(l) > 0:
        assume(any(l))
        assume(l[0] != 0)
    return l


@st.composite
def polys(
    draw: st.DrawFn, empty=True, sparse=False, domain="ZZ", gens=(x,), degree=None
):
    if sparse:
        raise NotImplementedError("Sparse polynomials are not yet supported.")
    else:
        return Poly(draw(coefficients(empty, degree=degree)), *gens, domain=domain)


def lattice_axioms_singular(polysList, func=lcm):
    if not polysList:
        return ValueError("Empty list of polynomials.")
    f, g, h = polysList[0], polysList[1], polysList[2]
    # Associativity
    assert func(func(f, g), h) == func(f, func(g, h))
    # Commutativity
    assert func(f, g) == func(g, f)
    # Idempotence
    assert func(f, f) == f  # when f is a negative integer
    return True


def lattice_axioms_dual(polysList, funcs=[lcm, gcd]):
    if not polysList:
        return ValueError("Empty list of polynomials.")
    f, g, h = polysList[0], polysList[1], polysList[2]
    # Absorption
    absorb1 = lambda func1, func2: func1(f, func2(f, g)) == f
    absorb2 = lambda func1, func2: func2(f, func1(f, g)) == f
    assert absorb1(funcs[0], funcs[1])
    assert absorb2(funcs[0], funcs[1])
    # Identity
    assert funcs[0](f, 0) == 0
    assert funcs[0](f, 1) == f  # lcm(0,1) = 0
    assert funcs[1](f, 1) == 1
    assert funcs[1](f, 0) == f
    # Distributivity
    dist1 = lambda func1, func2: func2(f, func1(g, h)) == func1(
        func2(f, g), func2(f, h)
    )
    dist2 = lambda func1, func2: func1(f, func2(g, h)) == func2(
        func1(f, g), func1(f, h)
    )
    assert dist1(funcs[0], funcs[1])
    assert dist2(funcs[0], funcs[1])
    return True
