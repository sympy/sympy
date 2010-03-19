"""Useful utilities for higher level polynomial classes. """

from sympy import S, sympify, Integer, Rational, Symbol, Add, Mul, Pow, ask

from sympy.polys.polyerrors import PolynomialError, GeneratorsNeeded

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

def _sort_gens(gens, **args):
    """Sort generators in a reasonably intelligent way. """
    sort = args.get('sort')
    wrt = args.get('wrt')

    gens_order = {}

    if sort is not None:
        for i, elt in enumerate(sort.split('>')):
            gens_order[elt.strip()] = i+1

    if wrt is not None:
        wrt = str(wrt)

    def order_key(x):
        x = str(x)

        if x == wrt:
            return (0, x)

        try:
            return ( gens_order[x], x)
        except KeyError:
            pass

        try:
            return (_gens_order[x], x)
        except KeyError:
            pass

        return (_max_order, x)

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
        exp, tail = exp.as_coeff_mul()

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

def _dict_from_basic_if_gens(ex, gens, **args):
    """Convert `ex` to a multinomial given a generators list. """
    k, indices = len(gens), {}

    for i, g in enumerate(gens):
        indices[g] = i

    result = {}

    for term in Add.make_args(ex):
        coeff, monom = [], [0]*k

        for factor in Mul.make_args(term):
            if factor.is_Number:
                coeff.append(factor)
            else:
                try:
                    base, exp = _analyze_power(*factor.as_base_exp())
                    monom[indices[base]] = exp
                except KeyError:
                    if not factor.has(*gens):
                        coeff.append(factor)
                    else:
                        raise PolynomialError("%s contains an element of the generators set" % factor)

        monom = tuple(monom)

        if monom in result:
            result[monom] += Mul(*coeff)
        else:
            result[monom] = Mul(*coeff)

    return result

def _dict_from_basic_no_gens(ex, **args):
    """Figure out generators and convert `ex` to a multinomial. """
    domain = args.get('domain')

    if domain is not None:
        def _is_coeff(factor):
            return factor in domain
    else:
        extension = args.get('extension')

        if extension is True:
            def _is_coeff(factor):
                return ask(factor, 'algebraic')
        else:
            greedy = args.get('greedy', True)

            if greedy is True:
                def _is_coeff(factor):
                    return False
            else:
                def _is_coeff(factor):
                    return factor.is_number

    gens, terms = set([]), []

    for term in Add.make_args(ex):
        coeff, elements = [], {}

        for factor in Mul.make_args(term):
            if factor.is_Number or _is_coeff(factor):
                coeff.append(factor)
            else:
                base, exp = _analyze_power(*factor.as_base_exp())

                elements[base] = exp
                gens.add(base)

        terms.append((coeff, elements))

    if not gens:
        raise GeneratorsNeeded("specify generators to give %s a meaning" % ex)

    gens = _sort_gens(gens, **args)

    k, indices = len(gens), {}

    for i, g in enumerate(gens):
        indices[g] = i

    result = {}

    for coeff, term in terms:
        monom = [0]*k

        for base, exp in term.iteritems():
            monom[indices[base]] = exp

        monom = tuple(monom)

        if monom in result:
            result[monom] += Mul(*coeff)
        else:
            result[monom] = Mul(*coeff)

    return result, tuple(gens)

def dict_from_basic(ex, gens=None, **args):
    """Converts a SymPy expression to a multinomial. """
    if args.get('expand', True):
        ex = ex.expand()

    if gens is not None:
        return _dict_from_basic_if_gens(ex, gens, **args)
    else:
        return _dict_from_basic_no_gens(ex, **args)

def basic_from_dict(rep, *gens):
    """Converts a multinomial to a SymPy expression. """
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

