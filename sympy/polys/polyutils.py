"""Useful utilities for higher level polynomial classes. """

from sympy import S, Basic, sympify, Integer, Rational, Symbol, Add, Mul, Pow, ask
from sympy.polys.polyerrors import PolynomialError, GeneratorsNeeded
from sympy.polys.polyoptions import build_options

import re

_gens_order = {
    'a': 301, 'b': 302, 'c': 303, 'd': 304,
    'e': 305, 'f': 306, 'g': 307, 'h': 308,
    'i': 309, 'j': 310, 'k': 311, 'l': 312,
    'm': 313, 'n': 314, 'o': 315, 'p': 216,
    'q': 217, 'r': 218, 's': 219, 't': 220,
    'u': 221, 'v': 222, 'w': 223, 'x': 124,
    'y': 125, 'z': 126,
}

_max_order = 1000
_re_gen = re.compile(r"^(.+?)(\d*)$")

def _sort_gens(gens, **args):
    """Sort generators in a reasonably intelligent way. """
    opt = build_options(args)

    gens_order, wrt = {}, None

    if opt is not None:
        gens_order, wrt = {}, opt.wrt

        for i, gen in enumerate(opt.sort):
            gens_order[gen] = i+1

    def order_key(gen):
        gen = str(gen)

        if wrt is not None:
            try:
                return (-len(wrt) + wrt.index(gen), gen, 0)
            except ValueError:
                pass

        name, index = _re_gen.match(gen).groups()

        if index:
            index = int(index)
        else:
            index = 0

        try:
            return ( gens_order[name], name, index)
        except KeyError:
            pass

        try:
            return (_gens_order[name], name, index)
        except KeyError:
            pass

        return (_max_order, name, index)

    try:
        gens = sorted(gens, key=order_key)
    except TypeError: # pragma: no cover
        pass

    return tuple(gens)

def _unify_gens(f_gens, g_gens):
    """Unify generators in a reasonably intelligent way. """
    f_gens = list(f_gens)
    g_gens = list(g_gens)

    if f_gens == g_gens:
        return tuple(f_gens)

    gens, common, k = [], [], 0

    for gen in f_gens:
        if gen in g_gens:
            common.append(gen)

    for i, gen in enumerate(g_gens):
        if gen in common:
            g_gens[i], k = common[k], k+1

    for gen in common:
        i = f_gens.index(gen)

        gens.extend(f_gens[:i])
        f_gens = f_gens[i+1:]

        i = g_gens.index(gen)

        gens.extend(g_gens[:i])
        g_gens = g_gens[i+1:]

        gens.append(gen)

    gens.extend(f_gens)
    gens.extend(g_gens)

    return tuple(gens)

def _analyze_gens(gens):
    """Support for passing generators as `*gens` and `[gens]`. """
    if len(gens) == 1 and hasattr(gens[0], '__iter__'):
        return tuple(gens[0])
    else:
        return tuple(gens)

def _sort_factors(factors, **args):
    """Sort low-level factors in increasing 'complexity' order. """
    def order_if_multiple_key((f, n)):
        return (len(f), n, f)

    def order_no_multiple_key(f):
        return (len(f), f)

    if args.get('multiple', True):
        return sorted(factors, key=order_if_multiple_key)
    else:
        return sorted(factors, key=order_no_multiple_key)

def _analyze_power(base, exp):
    """Extract non-integer part of `exp` to the `base`. """
    if exp.is_Number:
        if exp.is_Rational:
            if not exp.is_Integer:
                base = Pow(base, Rational(1, exp.q))

            exp = exp.p
        else:
            base, exp = Pow(base, exp), 1
    else:
        exp, tail = exp.as_coeff_terms()

        if exp.is_Number:
            if exp.is_Rational:
                if not exp.is_Integer:
                    tail += (Rational(1, exp.q),)

                exp = exp.p
            else:
                tail, exp = (exp,) + tail, 1
        else: # pragma: no cover
            raise PolynomialError("got invalid polynomial term")

        base = Pow(base, Mul(*tail))

    if exp < 0:
        exp, base = -exp, Pow(base, -S.One)

    return base, exp

def _parallel_dict_from_basic_if_gens(exprs, opt):
    """Transform expressions into a multinomial form given generators. """
    k, indices = len(opt.gens), {}

    for i, g in enumerate(opt.gens):
        indices[g] = i

    polys = []

    for expr in exprs:
        poly = {}

        for term in expr.as_Add():
            coeff, monom = [], [0]*k

            for factor in term.as_Mul():
                if factor.is_Number:
                    coeff.append(factor)
                else:
                    try:
                        base, exp = _analyze_power(*factor.as_Pow())
                        monom[indices[base]] = exp
                    except KeyError:
                        if not factor.has(*opt.gens):
                            coeff.append(factor)
                        else:
                            raise PolynomialError("%s contains an element of the generators set" % factor)

            monom = tuple(monom)

            if monom in poly:
                poly[monom] += Mul(*coeff)
            else:
                poly[monom] = Mul(*coeff)

        polys.append(poly)

    return polys, opt.gens

def _parallel_dict_from_basic_no_gens(exprs, opt):
    """Transform expressions into a multinomial form and figure out generators. """
    if opt.domain is not None:
        def _is_coeff(factor):
            return factor in opt.domain
    elif opt.extension is True:
        def _is_coeff(factor):
            return ask(factor, 'algebraic')
    elif opt.greedy is not False:
        def _is_coeff(factor):
            return False
    else:
        def _is_coeff(factor):
            return factor.is_number

    gens, reprs = set([]), []

    for expr in exprs:
        terms = []

        for term in expr.as_Add():
            coeff, elements = [], {}

            for factor in term.as_Mul():
                if factor.is_Number or _is_coeff(factor):
                    coeff.append(factor)
                else:
                    base, exp = _analyze_power(*factor.as_Pow())

                    elements[base] = exp
                    gens.add(base)

            terms.append((coeff, elements))

        reprs.append(terms)

    if not gens:
        raise GeneratorsNeeded("specify generators to give %s a meaning" % (exprs,))

    gens = _sort_gens(gens, opt=opt)
    k, indices = len(gens), {}

    for i, g in enumerate(gens):
        indices[g] = i

    polys = []

    for terms in reprs:
        poly = {}

        for coeff, term in terms:
            monom = [0]*k

            for base, exp in term.iteritems():
                monom[indices[base]] = exp

            monom = tuple(monom)

            if monom in poly:
                poly[monom] += Mul(*coeff)
            else:
                poly[monom] = Mul(*coeff)

        polys.append(poly)

    return polys, tuple(gens)

def _dict_from_basic_if_gens(expr, opt):
    """Transform an expression into a multinomial form given generators. """
    (poly,), gens = _parallel_dict_from_basic_if_gens((expr,), opt)
    return poly, gens

def _dict_from_basic_no_gens(expr, opt):
    """Transform an expression into a multinomial form and figure out generators. """
    (poly,), gens = _parallel_dict_from_basic_no_gens((expr,), opt)
    return poly, gens

def parallel_dict_from_basic(exprs, **args):
    """Transform expressions into a multinomial form. """
    opt = build_options(args)

    if opt.expand is not False:
        exprs = [ expr.expand() for expr in exprs ]

    if opt.gens:
        return _parallel_dict_from_basic_if_gens(exprs, opt)
    else:
        return _parallel_dict_from_basic_no_gens(exprs, opt)

def dict_from_basic(expr, **args):
    """Transform an expression into a multinomial form. """
    opt = build_options(args)

    if opt.expand is not False:
        expr = expr.expand()

    if opt.gens:
        return _dict_from_basic_if_gens(expr, opt)
    else:
        return _dict_from_basic_no_gens(expr, opt)

def basic_from_dict(rep, *gens):
    """Convert a multinomial form into an expression. """
    result = []

    for monom, coeff in rep.iteritems():
        term = [coeff]

        for g, m in zip(gens, monom):
            term.append(Pow(g, m))

        result.append(Mul(*term))

    return Add(*result)

def _dict_reorder(rep, gens, new_gens):
    """Reorder levels using dict representation. """
    gens = list(gens)

    monoms = rep.keys()
    coeffs = rep.values()

    new_monoms = [ [] for _ in xrange(len(rep)) ]

    for gen in new_gens:
        try:
            j = gens.index(gen)

            for M, new_M in zip(monoms, new_monoms):
                new_M.append(M[j])
        except ValueError:
            for new_M in new_monoms:
                new_M.append(0)

    return map(tuple, new_monoms), coeffs

def _parallel_dict_from_expr(exprs, opt):
    """Transform expressions into a multinomial form. """
    if opt.expand is not False:
        exprs = [ expr.expand() for expr in exprs ]

    if opt.gens:
        reps, gens = _parallel_dict_from_basic_if_gens(exprs, opt)
    else:
        reps, gens = _parallel_dict_from_basic_no_gens(exprs, opt)

    return reps, opt.clone({'gens': gens})

def _dict_from_expr(expr, opt):
    """Transform an expression into a multinomial form. """
    if opt.expand is not False:
        expr = expr.expand()

    if opt.gens:
        rep, gens = _dict_from_basic_if_gens(expr, opt)
    else:
        rep, gens = _dict_from_basic_no_gens(expr, opt)

    return rep, opt.clone({'gens': gens})

