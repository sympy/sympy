from sympy.simplify.hyperexpand import (ShiftA, ShiftB, UnShiftA, UnShiftB,
                       ReduceOrder, reduce_order, apply_operators,
                       devise_plan, make_derivative_operator, Formula,
                       hyperexpand, IndexPair)
from sympy import hyper, I, S
from sympy.utilities.pytest import raises
from sympy.abc import z, a, b, c
from sympy.utilities.randtest import test_numerically as tn
from sympy.utilities.pytest import XFAIL, skip
from random import randrange

from sympy import cos, sin, log, exp, asin, lowergamma, atanh, besseli, gamma

# whether to veryify that we can indeed do everything in the tables
# beware: this takes a *long* time
do_tables = False

def test_hyperexpand():
    # Luke, Y. L. (1969), The Special Functions and Their Approximations,
    # Volume 1, section 6.2

    assert hyperexpand(hyper([], [], z)) == exp(z)
    assert hyperexpand(hyper([1, 1], [2], -z)*z) == log(1 + z)
    assert hyperexpand(hyper([], [S.Half], -z**2/4)) == cos(z)
    assert hyperexpand(z*hyper([], [S('3/2')], -z**2/4)) == sin(z)
    assert hyperexpand(hyper([S('1/2'), S('1/2')], [S('3/2')], z**2)*z) \
           == asin(z)

def tables(fn):
    def wrapper():
        skip("This is too slow.")
    wrapper.__name__ = fn.__name__
    if do_tables:
        return fn
    else:
        return wrapper

def can_do(ap, bq, numerical=True):
    r = hyperexpand(hyper(ap, bq, z))
    if r.has(hyper):
        return False

    if not numerical:
        return True

    repl = {}
    for a in r.free_symbols - set([z]):
        repl[a] = randcplx()
    return tn(hyper(ap, bq, z).subs(repl), r.subs(repl), z)

def test_roach():
    # Kelly B. Roach.  Meijer G Function Representations.
    # Section "Gallery"
    assert can_do([S(1)/2], [S(9)/2])
    assert can_do([], [1, S(5)/2, 4])
    assert can_do([-S.Half, 1, 2], [3, 4])
    assert can_do([S(1)/3], [-S(2)/3, -S(1)/2, S(1)/2, 1])
    assert can_do([-S(3)/2, -S(1)/2], [-S(5)/2, 1])

@XFAIL
def test_roach_fail():
    assert can_do([-S(3)/2,], [-S(1)/2, S(1)/2]) # shine-integral
    assert can_do([-S(3)/2, -S(1)/2], [2]) # elliptic integrals
    assert can_do([-S(1)/2, 1], [S(1)/4, S(1)/2, S(3)/4]) # PFDD
    assert can_do([S(3)/2], [S(5)/2, 5]) # polylog
    assert can_do([-S(1)/2, S(1)/2, 1], [S(3)/2, S(5)/2]) # polylog, pfdd
    assert can_do([1, 2, 3], [S(1)/2, 4]) # XXX ?
    assert can_do([S(1)/2], [-S(1)/3, -S(1)/2, -S(2)/3]) # PFDD ?

# For the long table tests, see end of file

def test_polynomial():
    from sympy import oo
    assert hyperexpand(hyper([], [-1], z)) == oo
    assert hyperexpand(hyper([-2], [-1], z)) == oo
    assert hyperexpand(hyper([0, 0], [-1], z)) == 1
    assert can_do([-5, -2, randcplx(), randcplx()], [-10, randcplx()])

def test_hyperexpand_bases():
    assert hyperexpand(hyper([2], [a], z)) == \
  -1 + a + z**(1 - a)*(-2 - z + 3*a + a*z - a**2)*exp(z)*lowergamma(-1 + a, z)
    # TODO [a+1, a-S.Half], [2*a]
    assert hyperexpand(hyper([1, 2], [3], z)) == -2/z - 2*log(-z + 1)/z**2
    assert hyperexpand(hyper([S.Half, 2], [S(3)/2], z)) == \
      1/(-2*z + 2) + log((z**(S(1)/2) + 1)/(-z**(S(1)/2) + 1))/(4*z**(S(1)/2))
    assert hyperexpand(hyper([S(1)/2, S(1)/2], [S(5)/2], z)) == \
     3*(-z + 1)**(S(1)/2)/(4*z) + (2*z - 1)*asin(z**(S(1)/2))*3/(4*z**(S(3)/2))
    assert hyperexpand(hyper([1, 2], [S(3)/2], z)) == 1/(-2*z + 2) \
            + asin(z**(S(1)/2))/(z**(S(1)/2)*(-2*z + 2)*(-z + 1)**(S(1)/2))
    assert hyperexpand(hyper([-S.Half - 1, 1, 2], [S.Half, 3], z)) == \
             z**(S(1)/2)*(6*z/7 - S(6)/5)*atanh(z**(S(1)/2)) \
           + (-15*z**2 + 16*z - 3)/(35*z)*2 - 6*log(-z + 1)/(35*z**2)
    assert hyperexpand(hyper([1+S.Half, 1, 1], [2, 2], z)) == \
           -4*log((-z + 1)**(S(1)/2)/2 + S(1)/2)/z
    # TODO hyperexpand(hyper([a], [2*a + 1], z))
    # TODO [S.Half, a], [S(3)/2, a+1]
    assert hyperexpand(hyper([2], [b, 1], z)) == \
             z**(-b/2 + S(1)/2)*besseli(b - 1, 2*z**(S(1)/2))*gamma(b) \
           + z**(-b/2 + 1)*besseli(b, 2*z**(S(1)/2))*gamma(b)
    # TODO [a], [a - S.Half, 2*a]

def test_hyperexpand_parametric():
    assert hyperexpand(hyper([a, S(1)/2 + a], [S(1)/2], z)) \
        == (1 + z**(S(1)/2))**(-2*a)/2 + (1 - z**(S(1)/2))**(-2*a)/2
    assert hyperexpand(hyper([a, -S(1)/2 + a], [2*a], z)) \
        == 2**(2*a - 1)*((-z + 1)**(S(1)/2) + 1)**(-2*a + 1)

def test_shifted_sum():
    from sympy import simplify
    assert simplify(hyperexpand(z**4*hyper([2], [3, S('3/2')], -z**2))) \
           == -S(1)/2 + cos(2*z)/2 + z*sin(2*z) - z**2*cos(2*z)

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
    assert devise_plan(IndexPair([0], ()), IndexPair([0], ()), z) == []
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

@tables
def test_prudnikov_misc():
    assert can_do([1, (3 + I)/2, (3 - I)/2], [S(3)/2, 2])
    assert can_do([S.Half, a - 1], [S(3)/2, a + 1])
    assert can_do([], [b + 1])
    assert can_do([a], [a - 1, b + 1])

    assert can_do([a], [a - S.Half, 2*a])
    assert can_do([a], [a - S.Half, 2*a + 1])
    assert can_do([a], [a - S.Half, 2*a - 1])
    assert can_do([a], [a + S.Half, 2*a])
    assert can_do([a], [a + S.Half, 2*a + 1])
    assert can_do([a], [a + S.Half, 2*a - 1])
    assert can_do([S.Half], [b, 2-b])
    assert can_do([S.Half], [b, 3-b])
    assert can_do([1], [2, b])

    assert can_do([a, a+S.Half], [2*a, b, 2*a - b + 1])
    assert can_do([a, a+S.Half], [S.Half, 2*a, 2*a + S.Half])

@tables
def test_prudnikov_1():
    # A. P. Prudnikov, Yu. A. Brychkov and O. I. Marichev (1990).
    # Integrals and Series: More Special Functions, Vol. 3,.
    # Gordon and Breach Science Publisher

    # 7.3.1
    assert can_do([a, -a], [S.Half])
    assert can_do([a, 1 - a], [S.Half])
    assert can_do([a, 1 - a], [S(3)/2])
    assert can_do([a, 2 - a], [S.Half])
    assert can_do([a, 2 - a], [S(3)/2])
    assert can_do([a, 2 - a], [S(3)/2])
    assert can_do([a, a + S(1)/2], [2*a - 1])
    assert can_do([a, a + S(1)/2], [2*a])
    assert can_do([a, a + S(1)/2], [2*a + 1])
    assert can_do([a, a + S(1)/2], [S(1)/2])
    assert can_do([a, a + S(1)/2], [S(3)/2])
    assert can_do([a, a/2 + 1], [a/2])
    assert can_do([1, b], [2])

    assert can_do([a], [2*a])
    assert can_do([a], [2*a + 1])
    assert can_do([a], [2*a - 1])

@tables
def test_prudnikov_2():
    h = S.Half
    assert can_do([-h, -h], [h])
    assert can_do([-h, h], [3*h])
    assert can_do([-h, h], [5*h])
    assert can_do([-h, h], [7*h])
    assert can_do([-h, 1], [h])

    for p in [-h, h]:
      for n in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]:
          for m in [-h, h, 3*h, 5*h, 7*h]:
              assert can_do([p, n], [m])
      for n in [1, 2, 3, 4]:
          for m in [1, 2, 3, 4]:
              assert can_do([p, n], [m])

@tables
def test_prudnikov_3():
    h = S.Half
    assert can_do([S(1)/4, S(3)/4], [h])
    assert can_do([S(1)/4, S(3)/4], [3*h])
    assert can_do([S(1)/3, S(2)/3], [3*h])
    assert can_do([S(3)/4, S(5)/4], [h])
    assert can_do([S(3)/4, S(5)/4], [3*h])

    for p in [1, 2, 3, 4]:
      for n in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4, 9*h]:
          for m in [1, 3*h, 2, 5*h, 3, 7*h, 4]:
              assert can_do([p, m], [n])


@tables
def test_prudnikov_4():
    h = S.Half
    for p in [3*h, 5*h, 7*h]:
      for n in [-h, h, 3*h, 5*h, 7*h]:
          for m in [3*h, 2, 5*h, 3, 7*h, 4]:
              assert can_do([p, m], [n])
      for n in [1, 2, 3, 4]:
          for m in [2, 3, 4]:
              assert can_do([p, m], [n])

@tables
def test_prudnikov_5():
    h = S.Half

    for p in [1, 2, 3]:
        for q in range(p, 4):
            for r in [1, 2, 3]:
                for s in range(r, 4):
                    assert can_do([-h, p, q], [r, s])

    for p in [h, 1, 3*h, 2, 5*h, 3]:
        for q in [h, 3*h, 5*h]:
            for r in [h, 3*h, 5*h]:
                for s in [h, 3*h, 5*h]:
                    if s <= q and s <= r:
                        assert can_do([-h, p, q], [r, s])

    for p in [h, 1, 3*h, 2, 5*h, 3]:
        for q in [1, 2, 3]:
            for r in [h, 3*h, 5*h]:
                for s in [1, 2, 3]:
                    assert can_do([-h, p, q], [r, s])

@tables
def test_prudnikov_6():
    h = S.Half

    for m in [3*h, 5*h]:
        for n in [1, 2, 3]:
            for q in [h, 1, 2]:
                for p in [1, 2, 3]:
                     assert can_do([h, q, p], [m, n])
            for q in [1, 2, 3]:
                for p in [3*h, 5*h]:
                     assert can_do([h, q, p], [m, n])

    for q in [1, 2]:
      for p in [1, 2, 3]:
         for m in [1, 2, 3]:
             for n in [1, 2, 3]:
                 assert can_do([h, q, p], [m, n])

    assert can_do([h, h, 5*h], [3*h, 3*h])
    assert can_do([h, 1, 5*h], [3*h, 3*h])
    assert can_do([h, 2, 2], [1, 3])

    # pages 435 to 457 contain more PFDD and stuff like this

@tables
def test_prudnikov_7():
    assert can_do([3], [6])

    h = S.Half
    for n in [h, 3*h, 5*h, 7*h]:
        assert can_do([-h], [n])
    for m in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]: # HERE
        for n in [-h, h, 3*h, 5*h, 7*h, 1, 2, 3, 4]:
            assert can_do([m], [n])

@tables
def test_prudnikov_8():
    h = S.Half

    # 7.12.2
    for a in [1, 2, 3]:
        for b in [1, 2, 3]:
            for c in range(1, a+1):
                for d in [h, 1, 3*h, 2, 5*h, 3]:
                    assert can_do([a, b], [c, d])
        for b in [3*h, 5*h]:
            for c in [h, 1, 3*h, 2, 5*h, 3]:
                for d in [1, 2, 3]:
                    assert can_do([a, b], [c, d])

    for a in [-h, h, 3*h, 5*h]:
        for b in [1, 2, 3]:
            for c in [h, 1, 3*h, 2, 5*h, 3]:
                for d in [1, 2, 3]:
                    assert can_do([a, b], [c, d])
        for b in [h, 3*h, 5*h]:
            for c in [h, 3*h, 5*h, 3]:
                for d in [h, 1, 3*h, 2, 5*h, 3]:
                    if c <= b:
                        assert can_do([a, b], [c, d])

@tables
def test_prudnikov_9():
    # 7.13.1 [we have a general formula ... so this is a bit pointless]
    for i in range(9):
        assert can_do([], [(S(i) + 1)/2])
    for i in range(5):
        assert can_do([], [-(2*S(i) + 1)/2])

@tables
def test_prudnikov_10():
    # 7.14.2
    h = S.Half
    for p in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]:
      for m in [1, 2, 3, 4]:
          for n in range(m, 5):
              assert can_do([p], [m, n])

    for p in [1, 2, 3, 4]:
      for n in [h, 3*h, 5*h, 7*h]:
          for m in [1, 2, 3, 4]:
              assert can_do([p], [n, m])

    for p in [3*h, 5*h, 7*h]:
      for m in [h, 1, 2, 5*h, 3, 7*h, 4]:
         assert can_do([p], [h, m])
         assert can_do([p], [3*h, m])

    for m in [h, 1, 2, 5*h, 3, 7*h, 4]:
       assert can_do([7*h], [5*h, m])

@tables
def test_prudnikov_11():
    # 7.15
    assert can_do([a, a+S.Half], [2*a, b, 2*a - b])
    assert can_do([a, a+S.Half], [S(3)/2, 2*a, 2*a - S(1)/2])

    assert can_do([S(1)/4, S(3)/4], [S(1)/2, S(1)/2, 1])
    assert can_do([S(5)/4, S(3)/4], [S(3)/2, S(1)/2, 2])
    assert can_do([S(5)/4, S(3)/4], [S(3)/2, S(3)/2, 1])
    assert can_do([S(5)/4, S(7)/4], [S(3)/2, S(5)/2, 2])

@tables
def test_prudnikov_12():
    # 7.16
    assert can_do([], [a, a + S.Half, 2*a], False) # branches only agree for some z!
    assert can_do([], [a, a + S.Half, 2*a+1], False) # dito
    assert can_do([], [S.Half, a, a+S.Half])
    assert can_do([], [S(3)/2, a, a+S.Half])

    assert can_do([], [S(1)/4, S(1)/2, S(3)/4])
    assert can_do([], [S(1)/2, S(1)/2, 1])
    assert can_do([], [S(1)/2, S(3)/2, 1])
    assert can_do([], [S(3)/4, S(3)/2, S(5)/4])
    assert can_do([], [1, 1, S(3)/2])
    assert can_do([], [1, 2, S(3)/2])
    assert can_do([], [1, S(3)/2, S(3)/2])
    assert can_do([], [S(5)/4, S(3)/2, S(7)/4])
    assert can_do([], [2, S(3)/2, S(3)/2])

@XFAIL
def test_prudnikov_fail_2F1():
    assert can_do([a, b], [b + 1]) # incomplete beta function
    assert can_do([1, b], [b + 1]) # Lerch Phi
    assert can_do([-1, b], [c])    # Poly. also -2, -3 etc

    # TODO polys

    # Legendre functions:
    assert can_do([a, b], [a + b + S.Half])
    assert can_do([a, b], [a + b - S.Half])
    assert can_do([a, b], [a + b + S(3)/2])
    assert can_do([a, b], [(a + b + 1)/2])
    assert can_do([a, b], [(a + b)/2 + 1])
    assert can_do([a, b], [a - b + 1])
    assert can_do([a, b], [a - b + 2])
    assert can_do([a, b], [2*b])
    assert can_do([a, b], [S.Half])
    assert can_do([a, b], [S(3)/2])
    assert can_do([a, 1 - a], [c])
    assert can_do([a, 2 - a], [c])
    assert can_do([a, 3 - a], [c])
    assert can_do([a, a + S(1)/2], [c])
    assert can_do([1, b], [c])
    assert can_do([1, b], [S(3)/2])

    h = S.Half
    # Elliptic integrals
    for p in [-h, h]:
        for m in [h, 3*h, 5*h, 7*h]:
            for n in [1, 2, 3, 4]:
                assert can_do([p, m], [n])

    assert can_do([S(1)/4, S(3)/4], [1])

    # PFDD
    o = S(1)
    assert can_do([o/8, 1], [o/8*9])
    assert can_do([o/6, 1], [o/6*7])
    assert can_do([o/6, 1], [o/6*13])
    assert can_do([o/5, 1], [o/5*6])
    assert can_do([o/5, 1], [o/5*11])
    assert can_do([o/4, 1], [o/4*5])
    assert can_do([o/4, 1], [o/4*9])
    assert can_do([o/3, 1], [o/3*4])
    assert can_do([o/3, 1], [o/3*7])
    assert can_do([o/8*3, 1], [o/8*11])
    assert can_do([o/5*2, 1], [o/5*7])
    assert can_do([o/5*2, 1], [o/5*12])
    assert can_do([o/5*3, 1], [o/5*8])
    assert can_do([o/5*3, 1], [o/5*13])
    assert can_do([o/8*5, 1], [o/8*13])
    assert can_do([o/4*3, 1], [o/4*7])
    assert can_do([o/4*3, 1], [o/4*11])
    assert can_do([o/3*2, 1], [o/3*5])
    assert can_do([o/3*2, 1], [o/3*8])
    assert can_do([o/5*4, 1], [o/5*9])
    assert can_do([o/5*4, 1], [o/5*14])
    assert can_do([o/6*5, 1], [o/6*11])
    assert can_do([o/6*5, 1], [o/6*17])
    assert can_do([o/8*7, 1], [o/8*15])

@XFAIL
def test_prudnikov_fail_3F2():
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(1)/3, S(2)/3])
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(2)/3, S(4)/3])
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(4)/3, S(5)/3])

    # page 421
    assert can_do([a, a + S(1)/3, a + S(2)/3], [3*a/2, (3*a+1)/2])

    # pages 422 ...
    assert can_do([-S.Half, S.Half, S.Half], [1, 1]) # elliptic integrals
    assert can_do([-S.Half, S.Half, 1], [S(3)/2, S(3)/2])
    # TODO LOTS more

    # PFDD
    assert can_do([S(1)/8, S(3)/8, 1], [S(9)/8, S(11)/8])
    assert can_do([S(1)/8, S(5)/8, 1], [S(9)/8, S(13)/8])
    assert can_do([S(1)/8, S(7)/8, 1], [S(9)/8, S(15)/8])
    assert can_do([S(1)/6, S(1)/3, 1], [S(7)/6, S(4)/3])
    assert can_do([S(1)/6, S(2)/3, 1], [S(7)/6, S(5)/3])
    assert can_do([S(1)/6, S(2)/3, 1], [S(5)/3, S(13)/6])
    assert can_do([S.Half, 1, 1], [S(1)/4, S(3)/4])
    # LOTS more


@XFAIL
def test_prudnikov_fail_other():
    # 7.11.2
    assert can_do([a], [a+1]) # lowergamma ... why??

    # 7.12.1
    assert can_do([1, a], [b, 1 - 2*a + b]) # ???

    # 7.14.2
    assert can_do([-S(1)/2], [S(1)/2, S(1)/2]) # shine-integral shi
    assert can_do([-S(1)/2], [S(1)/2, 1]) # poly-log
    assert can_do([1], [S(1)/2, S(1)/2])  # poly-log
    assert can_do([S(1)/4], [S(1)/2, S(5)/4]) # PFDD
    assert can_do([S(3)/4], [S(3)/2, S(7)/4]) # PFDD
    assert can_do([1], [S(1)/4, S(3)/4]) # PFDD
    assert can_do([1], [S(3)/4, S(5)/4]) # PFDD
    assert can_do([1], [S(5)/4, S(7)/4]) # PFDD
    # TODO LOTS more

    # 7.15.2
    assert can_do([S(1)/2, 1], [S(3)/4, S(5)/4, S(3)/2]) # PFDD
    assert can_do([S(1)/2, 1], [S(7)/4, S(5)/4, S(3)/2]) # PFDD
    assert can_do([1, 1], [S(3)/2, 2, 2]) # cosh-integral chi

    # 7.16.1
    assert can_do([], [S(1)/3, S(2/3)]) # PFDD
    assert can_do([], [S(2)/3, S(4/3)]) # PFDD
    assert can_do([], [S(5)/3, S(4/3)]) # PFDD

    # XXX this does not *evaluate* right??
    assert can_do([], [a, a + S.Half, 2*a-1])
