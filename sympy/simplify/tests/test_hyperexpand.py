from sympy.simplify.hyperexpand import (ShiftA, ShiftB, UnShiftA, UnShiftB,
                       ReduceOrder, reduce_order, apply_operators,
                       devise_plan, make_derivative_operator, Formula,
                       hyperexpand, IndexPair)
from sympy import hyper, I, S
from sympy.utilities.pytest import raises
from sympy.abc import z
from sympy.utilities.randtest import test_numerically as tn
from sympy.utilities.pytest import XFAIL
from random import randrange

from sympy import cos, sin, log, exp, asin

def test_hyperexpand():
    # Luke, Y. L. (1969), The Special Functions and Their Approximations,
    # Volume 1, section 6.2

    assert hyperexpand(hyper([], [], z)) == exp(z)
    assert hyperexpand(hyper([1, 1], [2], -z)*z) == log(1 + z)
    assert hyperexpand(hyper([], [S.Half], -z**2/4)) == cos(z)
    assert hyperexpand(z*hyper([], [S('3/2')], -z**2/4)) == sin(z)
    assert hyperexpand(hyper([S('1/2'), S('1/2')], [S('3/2')], z**2)*z) \
           == asin(z)

    # TODO
    # More from above
    # Kelly' Roach's Papers
    # "Integrals and Series"

def randrat():
    """ Steer clear of integers. """
    return S(randrange(25) + 10)/50

def randcplx():
    """ Polys is not good with real coefficients. """
    return randrat() + I*randrat()

def test_formulae():
    from sympy.simplify.hyperexpand import FormulaCollection
    formulae = FormulaCollection().formulae
    for formula in formulae:
        h = hyper(formula.indices.ap, formula.indices.bq, formula.z)
        rep = {}
        for sym in formula.symbols:
            rep[sym] = randcplx()

        #print h, closed_form

        # first test if the closed-form is actually correct
        h = h.subs(rep)
        closed_form = formula.closed_form.subs(rep)
        z = formula.z
        assert tn(h, closed_form, z)

        # now test the computed matrix
        cl = (formula.C * formula.B)[0].subs(rep)
        assert tn(closed_form, cl, z)
        deriv1 = z*formula.B.diff(z)
        deriv2 = formula.M * formula.B
        for d1, d2 in zip(deriv1, deriv2):
            assert tn(d1.subs(rep), d2.subs(rep), z)

def op(f): return z*f.diff(z)

def test_plan():
    raises(ValueError, 'devise_plan(IndexPair([0], ()), IndexPair([0], ()), z)')
    raises(ValueError, 'devise_plan(IndexPair([1], ()), IndexPair((), ()), z)')
    raises(ValueError, 'devise_plan(IndexPair([2], [1]), IndexPair([2], [2]), z)')
    raises(KeyError,
           'devise_plan(IndexPair([2], []), IndexPair([S("1/2")], []), z)')

    a1, a2, b1 = map(lambda _: randcplx(), range(3))
    b1 += 2*I
    h = hyper([a1], [b1], z)

    h2 = hyper((a1 + 1, a2), [b1], z)
    tn(apply_operators(h, devise_plan(IndexPair((a1 + 1, a2), [b1]),
                                      IndexPair((a1, a2), [b1]), z), op),
       h2, z)

    h2 = hyper((a1 + 1, a2 - 1), [b1], z)
    tn(apply_operators(h, devise_plan(IndexPair((a1 + 1, a2 - 1), [b1]),
                                      IndexPair((a1, a2), [b1]), z), op),
       h2, z)

def test_plan_derivatives():
    a1, a2, a3 = 1, 2, S('1/2')
    b1, b2 = 3, S('5/2')
    h = hyper((a1, a2, a3), (b1, b2), z)
    h2 = hyper((a1 + 1, a2 + 1, a3 + 2), (b1 + 1, b2 + 1), z)
    ops = devise_plan(IndexPair((a1 + 1, a2 + 1, a3 + 2), (b1 + 1, b2 + 1)),
                      IndexPair((a1, a2, a3), (b1, b2)), z)
    f = Formula((a1, a2, a3), (b1, b2), z, h, [])
    deriv = make_derivative_operator(f.M, z)
    tn(apply_operators(f.C, ops, deriv)[0], h2, z)

    h2 = hyper((a1, a2 - 1, a3 - 2), (b1 - 1, b2 - 1), z)
    ops = devise_plan(IndexPair((a1, a2 - 1, a3 - 2), (b1 - 1, b2 - 1)),
                      IndexPair((a1, a2, a3), (b1, b2)), z)
    tn(apply_operators(f.C, ops, deriv)[0], h2, z)

def test_reduction_operators():
    a1, a2, b1 = map(lambda _: randcplx(), range(3))
    h = hyper([a1], [b1], z)

    assert ReduceOrder(2, 0) is None
    assert ReduceOrder(2, -1) is None
    assert ReduceOrder(0, -1) is None
    assert ReduceOrder(1, S('1/2')) is None

    h2 = hyper((a1, a2), (b1, a2), z)
    assert tn(ReduceOrder(a2, a2).apply(h, op), h2, z)

    h2 = hyper((a1, a2 + 1), (b1, a2), z)
    assert tn(ReduceOrder(a2 + 1, a2).apply(h, op), h2, z)

    h2 = hyper((a2 + 4, a1), (b1, a2), z)
    assert tn(ReduceOrder(a2 + 4, a2).apply(h, op), h2, z)

    # test several step order reduction
    ap = (a2 + 4, a1, b1 + 1)
    bq = (a2, b1, b1)
    nip, ops = reduce_order(IndexPair(ap, bq))
    assert nip.ap == (a1,)
    assert nip.bq == (b1,)
    assert tn(apply_operators(h, ops, op), hyper(ap, bq, z), z)

def test_shift_operators():
    a1, a2, b1, b2, b3 = map(lambda _: randcplx(), range(5))
    h = hyper((a1, a2), (b1, b2, b3), z)

    raises(ValueError, 'ShiftA(0)')
    raises(ValueError, 'ShiftB(1)')

    assert tn(ShiftA(a1).apply(h, op), hyper((a1 + 1, a2), (b1, b2, b3), z), z)
    assert tn(ShiftA(a2).apply(h, op), hyper((a1, a2 + 1), (b1, b2, b3), z), z)
    assert tn(ShiftB(b1).apply(h, op), hyper((a1, a2), (b1 - 1, b2, b3), z), z)
    assert tn(ShiftB(b2).apply(h, op), hyper((a1, a2), (b1, b2 - 1, b3), z), z)
    assert tn(ShiftB(b3).apply(h, op), hyper((a1, a2), (b1, b2, b3 - 1), z), z)

def test_ushift_operators():
    a1, a2, b1, b2, b3 = map(lambda _: randcplx(), range(5))
    h = hyper((a1, a2), (b1, b2, b3), z)

    raises(ValueError, 'UnShiftA((1,), (), 0, z)')
    raises(ValueError, 'UnShiftB((), (-1,), 0, z)')
    raises(ValueError, 'UnShiftA((1,), (0, -1, 1), 0, z)')
    raises(ValueError, 'UnShiftB((0, 1), (1,), 0, z)')

    s = UnShiftA((a1, a2), (b1, b2, b3), 0, z)
    assert tn(s.apply(h, op), hyper((a1 - 1, a2), (b1, b2, b3), z), z)
    s = UnShiftA((a1, a2), (b1, b2, b3), 1, z)
    assert tn(s.apply(h, op), hyper((a1, a2 - 1), (b1, b2, b3), z), z)

    s = UnShiftB((a1, a2), (b1, b2, b3), 0, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1 + 1, b2, b3), z), z)
    s = UnShiftB((a1, a2), (b1, b2, b3), 1, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1, b2 + 1, b3), z), z)
    s = UnShiftB((a1, a2), (b1, b2, b3), 2, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1, b2, b3 + 1), z), z)
