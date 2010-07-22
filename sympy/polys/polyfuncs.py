"""High--level polynomials manipulation functions. """

from sympy.polys.polytools import poly_from_expr, Poly
from sympy.polys.specialpolys import symmetric_poly
from sympy.polys.polyoptions import allowed_flags

from sympy.polys.polyerrors import (
    PolificationFailed, ComputationFailed,
)

from sympy.utilities import (
    all, any, numbered_symbols,
)

from sympy.core import S, Add, Mul

def symmetrize(f, *gens, **args):
    """Rewrite a polynomial in terms of elementary symmetric polynomials. """
    allowed_flags(args, ['formal'])

    try:
        f, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        return ComputationFailed('symmetrize', 1, exc)

    polys, symbols = [], numbered_symbols('s', start=1)

    gens, dom = f.gens, f.get_domain()

    for i in range(0, len(f.gens)):
        poly = symmetric_poly(i+1, gens, polys=True)
        polys.append((symbols.next(), poly.set_domain(dom)))

    indices = range(0, len(gens) - 1)
    weights = range(len(gens), 0, -1)

    symmetric = []

    if not f.is_homogeneous:
        symmetric.append(f.TC())
        f -= f.TC()

    while f:
        _height, _monom, _coeff = -1, None, None

        for i, (monom, coeff) in enumerate(f.terms()):
            if all(monom[i] >= monom[i+1] for i in indices):
                height = max([ n*m for n, m in zip(weights, monom) ])

                if height > _height:
                    _height, _monom, _coeff = height, monom, coeff

        if _height != -1:
            monom, coeff = _monom, _coeff
        else:
            break

        exponents = []

        for m1, m2 in zip(monom, monom[1:] + (0,)):
            exponents.append(m1 - m2)

        term = [ s**n for (s, _), n in zip(polys, exponents) ]
        poly = [ p**n for (_, p), n in zip(polys, exponents) ]

        symmetric.append(Mul(coeff, *term))

        product = poly[0].mul(coeff)

        for p in poly[1:]:
            product = product.mul(p)

        f -= product

    polys = [ (s, p.as_basic()) for s, p in polys ]

    if opt.formal:
        return (Add(*symmetric), f.as_basic(), dict(polys))
    else:
        return (Add(*symmetric).subs(polys), f.as_basic())

def horner(f, *gens, **args):
    """Rewrite a polynomial in Horner form. """
    allowed_flags(args, [])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed, exc:
        return exc.expr

    form, gen = S.Zero, F.gen

    if F.is_univariate:
        for coeff in F.all_coeffs():
            form = form*gen + coeff
    else:
        F, gens = Poly(F, gen), gens[1:]

        for coeff in F.all_coeffs():
            form = form*gen + horner(coeff, *gens, **args)

    return form
