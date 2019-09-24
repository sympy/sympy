from sympy import *
del test
del C

def solveset_multi(eqs, syms, domains):
    '''Basic implementation of a multivariate solveset'''

    rep = {}
    for sym, dom in zip(syms, domains):
        if dom is S.Reals:
            rep[sym] = Symbol(sym.name, real=True)
    eqs = [eq.subs(rep) for eq in eqs]
    syms = [sym.subs(rep) for sym in syms]

    syms = tuple(syms)

    if len(eqs) == 0:
        return ProductSet(*domains)

    if len(syms) == 1:
        sym = syms[0]
        domain = domains[0]
        solsets = [solveset(eq, sym, domain) for eq in eqs]
        solset = Intersection(*solsets)
        return ImageSet(Lambda((sym,), (sym,)), solset).doit()

    eqs = sorted(eqs, key=lambda eq: len(eq.free_symbols & set(syms)))

    for n, eq in enumerate(eqs):
        sols = []
        all_handled = True
        for sym in syms:
            if sym not in eq.free_symbols:
                continue
            sol = solveset(eq, sym, domains[syms.index(sym)])

            if isinstance(sol, FiniteSet):
                i = syms.index(sym)
                symsp = syms[:i] + syms[i+1:]
                domainsp = domains[:i] + domains[i+1:]
                eqsp = eqs[:n] + eqs[n+1:]
                for s in sol:
                    eqsp_sub = [eq.subs(sym, s) for eq in eqsp]
                    sol_others = solveset_multi(eqsp_sub, symsp, domainsp)
                    fun = Lambda((symsp,), symsp[:i] + (s,) + symsp[i:])
                    sols.append(ImageSet(fun, sol_others).doit())
            else:
                all_handled = False
        if all_handled:
            return Union(*sols)


def test_basic():
    from sympy.abc import x
    assert solveset_multi([x**2-1], [x], [S.Reals]) == FiniteSet((1,), (-1,))

def test_two():
    from sympy.abc import x, y
    assert solveset_multi([x+y, x+1], [x, y], [Reals, Reals]) == FiniteSet((-1, 1))
    assert solveset_multi([x+y, x+1], [y, x], [Reals, Reals]) == FiniteSet((1, -1))
    assert solveset_multi([x+y, x-y-1], [x, y], [Reals, Reals]) == FiniteSet((S(1)/2, -S(1)/2))
    assert solveset_multi([x-1, y-2], [x, y], [Reals, Reals]) == FiniteSet((1, 2))
    #assert solveset_multi([x+y], [x, y], [Reals, Reals]) == ImageSet(Lambda(x, (x, -x)), Reals)
    assert solveset_multi([x+y], [x, y], [Reals, Reals]) == Union(ImageSet(Lambda(((x,),), (x, -x)), ProductSet(Reals)), ImageSet(Lambda(((y,),), (-y, y)), ProductSet(Reals)))
    assert solveset_multi([x+y, x+y+1], [x, y], [Reals, Reals]) == S.EmptySet
    assert solveset_multi([x+y, x-y, x-1], [x, y], [Reals, Reals]) == S.EmptySet
    assert solveset_multi([x+y, x-y, x-1], [y, x], [Reals, Reals]) == S.EmptySet

def test_three():
    from sympy.abc import x, y, z
    assert solveset_multi([x+y+z-1, x+y-z-2, x-y-z-3], [x, y, z], [Reals,
        Reals, Reals]) == FiniteSet((2, -S.Half, -S.Half))

def test_nonlin():
    from sympy.abc import r, theta, z, x, y
    assert solveset_multi([x**2+y**2-2, x+y], [x, y], [Reals, Reals]) == FiniteSet((-1, 1), (1, -1))
    assert solveset_multi([x**2-1, y], [x, y], [Reals, Reals]) == FiniteSet((1, 0), (-1, 0))
    assert solveset_multi([r*cos(theta)-1, r*sin(theta)], [theta, r],
            [Interval(0, pi), Interval(-1, 1)]) == FiniteSet((0, 1), (pi, -1))
    assert solveset_multi([r*cos(theta)-1, r*sin(theta)], [r, theta],
            [Interval(0, 1), Interval(0, pi)]) == FiniteSet((1, 0))
    #assert solveset_multi([r*cos(theta)-r, r*sin(theta)], [r, theta],
    #        [Interval(0, 1), Interval(0, pi)]) == ?
    assert solveset_multi([r*cos(theta)-r, r*sin(theta)], [r, theta],
            [Interval(0, 1), Interval(0, pi)]) == Union(ImageSet(Lambda(((r,),), (r, 0)), ImageSet(Lambda(r, (r,)), Interval(0, 1))), ImageSet(Lambda(((theta,),), (0, theta)), ImageSet(Lambda(theta, (theta,)), Interval(0, pi))))
