"""Computational algebraic number field theory. """

from sympy.core import (
    S, Basic, I, Integer, Symbol, Add, Mul, sympify,
)

from sympy.polys.polytools import (
    Poly, factor_list, groebner,
)

from sympy.polys.polyerrors import (
    NotAlgebraic,
)

from sympy.utilities import (
    any, all, numbered_symbols,
)

def minpoly(ex, x=None, **args):
    """Computes the minimal polynomial of an algebraic number. """
    generator = numbered_symbols('a', dummy=True)
    mapping, symbols = {}, {}

    ex = sympify(ex)

    if x is not None:
        x = sympify(x)
    else:
        x = Symbol('x', dummy=True)

    def update_mapping(ex, exp, base):
        a = generator.next()

        symbols[ex] = a
        mapping[ex] = a**exp + base

        return a

    def bottom_up_scan(ex):
        if ex.is_Atom:
            if ex is S.ImaginaryUnit:
                if ex not in mapping:
                    return update_mapping(ex, 2, 1)
                else:
                    return symbols[ex]
            elif ex.is_Rational:
                return ex
        elif ex.is_Add:
            return Add(*[ bottom_up_scan(g) for g in ex.args ])
        elif ex.is_Mul:
            return Mul(*[ bottom_up_scan(g) for g in ex.args ])
        elif ex.is_Pow:
            if ex.exp.is_Rational:
                base = bottom_up_scan(ex.base)

                if base != ex.base:
                    power = base**ex.exp
                else:
                    power = ex

                if power not in mapping:
                    return update_mapping(power, 1/ex.exp, -base)
                else:
                    return symbols[power]

                return a

        raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)

    if ex.is_AlgebraicNumber:
        if args.get('polys', False):
            return ex.minpoly
        else:
            return ex.minpoly.as_basic()

    F = [x - bottom_up_scan(ex)] + mapping.values()
    G = groebner(F, *(symbols.values() + [x]), order='lex')

    _, factors = factor_list(G[-1])

    if len(factors) == 1:
        ((result, _),) = factors
    else:
        for result, _ in factors:
            if result.subs(x, ex).expand().is_zero:
                break
        else: # pragma: no cover
            raise NotImplementedError("multiple candidates for the minimal polynomial of %s" % ex)

    if args.get('polys', False):
        return Poly(result, x)
    else:
        return result

class AlgebraicNumber(Basic):
    """Class for representing algebraic numbers in SymPy. """

    __slots__ = ['expr', 'alias', 'minpoly']

    is_AlgebraicNumber = True

    def __new__(cls, ex, alias=None):
        obj = Basic.__new__(cls)

        obj.expr = sympify(ex)
        obj.minpoly = minpoly(obj.expr)

        if alias is not None:
            if isinstance(alias, basestring):
                obj.alias = Symbol(alias)
            else:
                obj.alias = alias
        else:
            obj.alias = None

        return obj

    @property
    def is_aliased(self):
        return self.alias is not None

