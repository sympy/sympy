from __future__ import annotations

from sympy.polys import sparsemonomials as smm
from sympy.polys.domains import GF, QQ, ZZ
from sympy.testing.pytest import raises


def test_sparse_monomial_conversion_and_degree():
    mon = ((0, 3), (2, 1), (5, 4))

    assert smm.to_dense(mon, 6) == (3, 0, 1, 0, 0, 4)
    assert smm.from_dense((3, 0, 1, 0, 0, 4)) == mon

    assert smm.degree(mon, 0) == 3
    assert smm.degree(mon, 1) == 0
    assert smm.degree(mon, 5) == 4
    assert smm.degree(mon, 7) == 0


def test_sparse_monomial_set_exp():
    mon = ((1, 2), (3, 4))

    assert smm.set_exp(mon, 0, 5) == ((0, 5), (1, 2), (3, 4))
    assert smm.set_exp(mon, 2, 5) == ((1, 2), (2, 5), (3, 4))
    assert smm.set_exp(mon, 4, 5) == ((1, 2), (3, 4), (4, 5))
    assert smm.set_exp(mon, 1, 5) == ((1, 5), (3, 4))
    assert smm.set_exp(mon, 1, 0) == ((3, 4),)


def test_sparse_monomial_multiplication_and_lex():
    a = ((0, 2), (3, 1))
    b = ((1, 4), (3, 2), (5, 1))

    assert smm.mul_monom(a, b) == ((0, 2), (1, 4), (3, 3), (5, 1))

    assert smm._lex_gt(((0, 2),), ((0, 1),))
    assert smm._lex_gt(((0, 1),), ((1, 9),))
    assert not smm._lex_gt(((1, 9),), ((0, 1),))
    assert smm._lex_gt(((0, 1), (2, 1)), ((0, 1),))


def test_sparse_polynomial_arithmetic():
    x = ((0, 1),)
    y = ((1, 1),)
    f = {x: ZZ(2), y: ZZ(3), (): ZZ(4)}
    g = {x: ZZ(-2), y: ZZ(1), (): ZZ(-4)}

    assert smm.add(f, g, ZZ) == {y: ZZ(4)}
    assert smm.sub(f, f, ZZ) == {}
    assert smm.neg(f) == {x: ZZ(-2), y: ZZ(-3), (): ZZ(-4)}

    assert smm.add_ground(f, ZZ(-4), ZZ) == {x: ZZ(2), y: ZZ(3)}
    assert smm.sub_ground(f, ZZ(4), ZZ) == {x: ZZ(2), y: ZZ(3)}
    assert smm.mul_ground(f, ZZ(0)) == {}
    assert smm.mul_ground(f, ZZ(2)) == {
        x: ZZ(4), y: ZZ(6), (): ZZ(8)}

    assert smm.mul({x: ZZ(1), (): ZZ(1)},
                   {x: ZZ(1), (): ZZ(-1)}, ZZ) == {
        ((0, 2),): ZZ(1), (): ZZ(-1)}

    h = {x: ZZ(1), (): ZZ(1)}
    assert smm.square(h, ZZ) == {
        ((0, 2),): ZZ(1), x: ZZ(2), (): ZZ(1)}
    assert smm.pow_generic(h, 3, ZZ) == {
        ((0, 3),): ZZ(1),
        ((0, 2),): ZZ(3),
        x: ZZ(3),
        (): ZZ(1),
    }

    assert smm.pow_generic({}, 2, ZZ) == {}
    assert smm.pow_generic({}, 0, ZZ) == {(): ZZ.one}
    raises(ValueError, lambda: smm.pow_generic(h, -1, ZZ))

    # In characteristic 2 the cross term in a square vanishes.
    assert smm.square({x: GF(2).one, (): GF(2).one}, GF(2)) == {
        ((0, 2),): GF(2).one, (): GF(2).one}


def test_sparse_polynomial_degrees_and_calculus():
    d = {
        ((0, 3), (2, 1)): ZZ(2),
        ((1, 4),): ZZ(5),
        (): ZZ(7),
    }

    assert smm.poly_degree(d, 0) == 3
    assert smm.poly_degree(d, 1) == 4
    assert smm.poly_degree({}, 0) == -1
    assert smm.poly_degrees(d, 3) == (3, 4, 1)
    assert smm.poly_degrees({}, 3) == (-1, -1, -1)
    assert smm.total_degree(d) == 4
    assert smm.total_degree({}) == -1

    f = {((0, 3),): ZZ(2), ((0, 1),): ZZ(5), (): ZZ(7)}
    assert smm.diff(f, 0, ZZ) == {
        ((0, 2),): ZZ(6), (): ZZ(5)}

    g = {((0, 2),): QQ(6), (): QQ(5)}
    assert smm.integrate(g, 0, QQ) == {
        ((0, 3),): QQ(2), ((0, 1),): QQ(5)}


def test_sparse_polynomial_coefficients():
    x2 = ((0, 2),)
    y = ((1, 1),)
    d = {y: ZZ(6), x2: ZZ(4), (): ZZ(8)}

    assert smm.LC(d, ZZ) == ZZ(4)
    assert smm.LC({}, ZZ) == ZZ.zero
    assert smm.content(d, ZZ) == ZZ(2)
    assert smm.content({}, ZZ) == ZZ.zero
    assert smm.primitive(d, ZZ) == (
        ZZ(2), {y: ZZ(3), x2: ZZ(2), (): ZZ(4)})
    assert smm.primitive({}, ZZ) == (ZZ.zero, {})

    q = {x2: QQ(1, 2), (): QQ(2, 3)}
    assert smm.clear_denoms(q, QQ) == (
        ZZ(6), {x2: QQ(3), (): QQ(4)})
    assert smm.clear_denoms(d, ZZ) == (ZZ.one, d)


def test_sparse_polynomial_truncation():
    x = ((0, 1),)

    assert smm.trunc_ground({x: ZZ(7), (): ZZ(-6)}, ZZ(5), ZZ) == {
        x: ZZ(2), (): ZZ(-1)}
    assert smm.trunc_ground({x: ZZ(5), (): ZZ(10)}, ZZ(5), ZZ) == {}


def test_sparse_polynomial_substitution_and_composition():
    x = ((0, 1),)
    y = ((1, 1),)
    z = ((2, 1),)

    f = {
        ((0, 1), (2, 1)): ZZ(2),
        ((1, 1), (2, 1)): ZZ(3),
        z: ZZ(-7),
    }

    # Substituting x = 2 and y = 1 drops both coordinates.  All terms
    # become multiples of the remaining z coordinate and cancel.
    assert smm.subs_drop(f, {0: ZZ(2), 1: ZZ(1)}, 3, ZZ) == {}

    f = {((0, 2),): ZZ(1), ((1, 1),): ZZ(1)}
    reps = {0: {x: ZZ(1), (): ZZ(1)}}
    assert smm.compose(f, reps, ZZ) == {
        ((0, 2),): ZZ(1),
        x: ZZ(2),
        y: ZZ(1),
        (): ZZ(1),
    }

    # A zero replacement removes every term containing that generator.
    assert smm.compose({x: ZZ(1), y: ZZ(2)}, {0: {}}, ZZ) == {
        y: ZZ(2)}


def test_sparse_polynomial_predicates():
    x = ((0, 1),)
    x2 = ((0, 2),)
    xy = ((0, 1), (1, 1))

    assert smm.is_zero({})
    assert not smm.is_zero({x: ZZ.one})

    assert smm.is_one({(): ZZ.one}, ZZ)
    assert not smm.is_one({(): ZZ(2)}, ZZ)

    assert smm.is_ground({})
    assert smm.is_ground({(): ZZ(2)})
    assert not smm.is_ground({x: ZZ.one})

    assert smm.is_monic({x: ZZ.one, (): ZZ(2)}, ZZ)
    assert not smm.is_monic({x: ZZ(2)}, ZZ)

    assert smm.is_primitive({x: ZZ(2), (): ZZ(3)}, ZZ)
    assert not smm.is_primitive({x: ZZ(2), (): ZZ(4)}, ZZ)

    assert smm.is_linear({x: ZZ.one, (): ZZ.one})
    assert not smm.is_linear({xy: ZZ.one})
    assert smm.is_quadratic({x2: ZZ.one, xy: ZZ.one})
    assert not smm.is_quadratic({((0, 3),): ZZ.one})
