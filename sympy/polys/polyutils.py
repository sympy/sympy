"""Useful utilities for higher level polynomial classes. """

from sympy import S, sympify, Integer, Rational, Symbol, Add, Mul, Pow

from sympy.polys.polyerrors import PolynomialError, GeneratorsNeeded

def _update_args(args, key, value):
    """Add a new `(key, value)` pair to arguments dict. """
    args = dict(args)

    if not args.has_key(key):
        args[key] = value

    return args

def _analyze_gens(gens):
    """Support for passing generators as `*gens` and `[gens]`. """
    if len(gens) == 1 and hasattr(gens[0], '__iter__'):
        return tuple(gens[0])
    else:
        return gens

def _analyze_modulus(args):
    """Convert `modulus` to an internal representation. """
    modulus = args.get('modulus')

    if modulus is not None:
        modulus = sympify(modulus)

        if modulus.is_Integer and modulus >= 2:
            modulus = int(modulus)
        else:
            raise PolynomialError("modulus must be an integer >= 2, got %s" % modulus)

    return modulus

def _analyze_extension(args):
    """Convert `extension` to an internal representation. """
    extension = args.get('extension')

    if extension is not None:
        if not hasattr(extension, '__iter__'):
            return set([extension])
        else:
            return set(extension)
    else:
        return set([])

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
        else:
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

    for term in ex.as_Add():
        coeff, monom = [], [0]*k

        for factor in term.as_Mul():
            if factor.is_Number:
                coeff.append(factor)
            else:
                try:
                    base, exp = _analyze_power(*factor.as_Pow())
                    monom[indices[base]] = exp
                except KeyError:
                    if not factor.has(*gens):
                        coeff.append(factor)
                    else:
                        raise PolynomialError("%s contains an element of the generators set" % factor)

        monom = tuple(monom)

        if result.has_key(monom):
            result[monom] += Mul(*coeff)
        else:
            result[monom] = Mul(*coeff)

    return result

def _dict_from_basic_no_gens(ex, **args):
    """Figure out generators and convert `ex` to a multinomial. """
    greedy = args.get('greedy', True)
    gens, terms = set([]), []

    for term in ex.as_Add():
        coeff, elements = [], {}

        for factor in term.as_Mul():
            if factor.is_Number:
                coeff.append(factor)
            else:
                if not greedy and factor.is_number:
                    coeff.append(factor)
                else:
                    base, exp = _analyze_power(*factor.as_Pow())

                    elements[base] = exp
                    gens.add(base)

        terms.append((coeff, elements))

    if not gens:
        raise GeneratorsNeeded("specify generators to give %s meaning" % ex)

    try:
        gens = sorted(gens)
    except TypeError:
        gens = list(gens)

    k, indices = len(gens), {}

    for i, g in enumerate(gens):
        indices[g] = i

    result = {}

    for coeff, term in terms:
        monom = [0]*k

        for base, exp in term.iteritems():
            monom[indices[base]] = exp

        monom = tuple(monom)

        if result.has_key(monom):
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

