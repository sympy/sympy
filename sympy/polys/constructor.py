"""Tools for constructing domains for expressions. """

from sympy.polys.polyutils import parallel_dict_from_basic
from sympy.polys.polyoptions import build_options
from sympy.polys.domains import ZZ, QQ, RR, EX
from sympy.assumptions import ask
from sympy.core import S, sympify

def _construct_simple(coeffs, opt):
    """Handle simple domains, e.g.: ZZ, QQ, RR and algebraic domains. """
    result, rationals, reals, algebraics = {}, False, False, False

    if opt.extension is True:
        is_algebraic = lambda coeff: ask(coeff, 'algebraic')
    else:
        is_algebraic = lambda coeff: False

    # XXX: add support for a + b*I coefficients
    for coeff in coeffs:
        if coeff.is_Rational:
            if not coeff.is_Integer:
                rationals = True
        elif coeff.is_Real:
            if not algebraics:
                reals = True
            else:
                # there are both reals and algebraics -> EX
                return False
        elif is_algebraic(coeff):
            if not reals:
                algebraics = True
            else:
                # there are both algebraics and reals -> EX
                return False
        else:
            # this is a composite domain, e.g. ZZ[X], EX
            return None

    if algebraics:
        domain, result = _construct_algebraic(coeffs, opt)
    else:
        if reals:
            domain = RR
        else:
            if opt.field or rationals:
                domain = QQ
            else:
                domain = ZZ

        result = []

        for coeff in coeffs:
            result.append(domain.from_sympy(coeff))

    return domain, result

def _construct_algebraic(coeffs, opt):
    """We know that coefficients are algebraic so construct the extension. """
    from sympy.polys.numberfields import primitive_element

    result, exts = [], set([])

    for coeff in coeffs:
        if coeff.is_Rational:
            coeff = (None, 0, QQ.from_sympy(coeff))
        else:
            a = coeff.as_coeff_add()[0]
            coeff -= a

            b = coeff.as_coeff_mul()[0]
            coeff /= b

            exts.add(coeff)

            a = QQ.from_sympy(a)
            b = QQ.from_sympy(b)

            coeff = (coeff, b, a)

        result.append(coeff)

    exts = list(exts)

    g, span, H = primitive_element(exts, ex=True, polys=True)
    root = sum([ s*ext for s, ext in zip(span, exts) ])

    domain, g = QQ.algebraic_field((g, root)), g.rep.rep

    for i, (coeff, a, b) in enumerate(result):
        if coeff is not None:
            coeff = a*domain.dtype.from_list(H[exts.index(coeff)], g, QQ) + b
        else:
            coeff = domain.dtype.from_list([b], g, QQ)

        result[i] = coeff

    return domain, result

def _construct_composite(coeffs, opt):
    """Handle composite domains, e.g.: ZZ[X], QQ[X], ZZ(X), QQ(X). """
    numers, denoms = [], []

    for coeff in coeffs:
        numer, denom = coeff.as_numer_denom()

        numers.append(numer)
        denoms.append(denom)

    polys, gens = parallel_dict_from_basic(numers + denoms) # XXX: sorting

    if any(gen.is_number for gen in gens):
        return None # generators are number-like so lets better use EX

    n = len(gens)
    k = len(polys)//2

    numers = polys[:k]
    denoms = polys[k:]

    if opt.field:
        fractions = True
    else:
        fractions, monom = False, (0,)*n

        for denom in denoms:
            if len(denom) > 1 or monom not in denom:
                fractions = True
                break

    result = []

    if not fractions:
        rationals = False

        for i, denom in enumerate(denoms):
            coeff = denom[monom]
            denoms[i] = coeff

            if coeff is not S.One:
                rationals = True

        if not rationals:
            domain = ZZ.poly_ring(*gens)

            for numer in numers:
                for monom, coeff in numer.iteritems():
                    numer[monom] = ZZ.from_sympy(coeff)

                result.append(domain(numer))
        else:
            domain = QQ.poly_ring(*gens)

            for numer, denom in zip(numers, denoms):
                for monom, coeff in numer.iteritems():
                    numer[monom] = QQ.from_sympy(coeff/denom)

                result.append(domain(numer))
    else:
        domain = ZZ.frac_field(*gens)

        for numer, denom in zip(numers, denoms):
            for monom, coeff in numer.iteritems():
                numer[monom] = ZZ.from_sympy(coeff)

            for monom, coeff in denom.iteritems():
                denom[monom] = ZZ.from_sympy(coeff)

            result.append(domain((numer, denom)))

    return domain, result

def _construct_expression(coeffs, opt):
    """The last resort case, i.e. use the expression domain. """
    domain, result = EX, []

    for coeff in coeffs:
        result.append(domain.from_sympy(coeff))

    return domain, result

def construct_domain(obj, **args):
    """Construct a minimal domain for the list of coefficients. """
    opt = build_options(**args)

    if isinstance(obj, dict):
        monoms, coeffs = zip(*obj.items())
    else:
        coeffs = obj

    coeffs = map(sympify, coeffs)
    result = _construct_simple(coeffs, opt)

    if result is not None:
        if result is not False:
            domain, coeffs = result
        else:
            domain, coeffs = _construct_expression(coeffs, opt)
    else:
        result = _construct_composite(coeffs, opt)

        if result is not None:
            domain, coeffs = result
        else:
            domain, coeffs = _construct_expression(coeffs, opt)

    if isinstance(obj, dict):
        return domain, dict(zip(monoms, coeffs))
    else:
        return domain, coeffs
