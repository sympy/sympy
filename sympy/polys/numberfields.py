"""Computational algebraic number field theory. """

from sympy.core import (
    S, Basic, I, Integer, Symbol, Add, Mul, sympify,
)

from sympy.polys.polytools import (
    Poly, sqf_norm, factor_list, groebner,
)

from sympy.polys.polyutils import (
    basic_from_dict,
)

from sympy.polys.polyclasses import (
    ANP, DMP,
)

from sympy.polys.polyerrors import (
    IsomorphismFailed,
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

    def update_mapping(ex, exp, base=None):
        a = generator.next()
        symbols[ex] = a

        if base is not None:
            mapping[ex] = a**exp + base
        else:
            mapping[ex] = exp.as_basic(a)

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
        elif ex.is_AlgebraicNumber:
            if ex.root not in mapping:
                return update_mapping(ex.root, ex.minpoly)
            else:
                return symbols[ex.root]

        raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)

    polys = args.get('polys', False)

    if ex.is_AlgebraicNumber:
        if not polys:
            return ex.minpoly.as_basic()
        else:
            return ex.minpoly
    elif ex.is_Rational:
        result = ex.q*x - ex.p
    else:
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

    if polys:
        return Poly(result, x, field=True)
    else:
        return result

def primitive_element(*extension, **args):
    """Construct a common number field for all extensions.  """
    extension = map(AlgebraicNumber, extension)

    minpoly = extension[0].minpoly
    root    = extension[0].root

    for ext in extension[1:]:
        s, _, minpoly = sqf_norm(minpoly, extension=ext)
        root = s*root + ext.root

    return minpoly, root

def field_isomorphism(a, b, **args):
    """Construct an isomorphism between two number fields. """
    a, b = sympify(a), sympify(b)

    if not a.is_AlgebraicNumber:
        a = AlgebraicNumber(a)

    if not b.is_AlgebraicNumber:
        b = AlgebraicNumber(b)

    if a == b:
        return a.coeffs()

    n = a.minpoly.degree()
    m = b.minpoly.degree()

    if n == 1:
        return [a.root]

    if m % n != 0:
        return None

    A = a.minpoly.LC()
    B = b.minpoly.LC()

    a = a.to_algebraic_integer()
    b = b.to_algebraic_integer()

    _, factors = factor_list(a.minpoly, extension=b)

    for f, _ in factors:
        if f.degree() == 1:
            isomorphism = f.rep.TC()

            if isomorphism.LC() < 0:
                continue

            if a.root.is_negative != b.root.is_negative:
                isomorphism = -isomorphism

            coeffs = isomorphism.to_sympy_list()

            for i in xrange(0, len(coeffs)):
                coeffs[i] = B**i*coeffs[i]/A

            return coeffs
    else:
        return None

minimal_polynomial = minpoly

class AlgebraicNumber(Basic):
    """Class for representing algebraic numbers in SymPy. """

    __slots__ = ['rep', 'root', 'alias', 'minpoly']

    is_AlgebraicNumber = True

    def __new__(cls, expr, theta=None, **args):
        """Construct a new algebraic number. """
        if not isinstance(expr, (list, ANP)):
            expr = sympify(expr)

        if theta is not None:
            theta = sympify(theta)

            if not theta.is_AlgebraicNumber:
                theta = AlgebraicNumber(theta)

            dom = theta.minpoly.get_domain()

            root = theta.root
            alias = theta.alias
            minpoly = theta.minpoly

            if isinstance(expr, ANP):
                expr = expr.rep

            if isinstance(expr, list):
                rep = DMP.from_list(expr, 0, dom)
            else:
                coeffs = field_isomorphism(expr, theta)

                if coeffs is not None:
                    rep = DMP.from_sympy_list(coeffs, 0, dom)
                else:
                    raise IsomorphismFailed("%s is not in a subfield of %s" % (expr, theta.root))
        else:
            if not hasattr(expr, '__iter__'):
                minpoly, root = minimal_polynomial(expr, polys=True), expr
            else:
                minpoly, root = primitive_element(*expr)

            rep = DMP.from_list([1, 0], 0, minpoly.get_domain())

        alias = args.get('alias')

        if alias is not None:
            alias = sympify(alias)

        obj = Basic.__new__(cls)

        obj.rep = rep
        obj.root = root
        obj.alias = alias
        obj.minpoly = minpoly

        return obj

    def __eq__(a, b):
        if not b.is_AlgebraicNumber:
            try:
                b = AlgebraicNumber(b, a)
            except (NotAlgebraic, IsomorphismFailed):
                return False

        return a.rep == b.rep and \
            a.minpoly.all_coeffs() == b.minpoly.all_coeffs()

    @property
    def is_aliased(self):
        """Returns `True` if `alias` was set. """
        return self.alias is not None

    def as_poly(self, x=None):
        """Create a Poly instance from `self`. """
        if x is not None:
            return Poly(self.rep, x)
        else:
            if self.alias is not None:
                return Poly(self.rep, self.alias)
            else:
                return Poly(self.rep, Symbol('x', dummy=True))

    def as_basic(self, x=None):
        """Create a Basic expression from `self`. """
        return self.as_poly(x or self.root).as_basic().expand()

    def coeffs(self):
        """Returns all coefficients of an algebraic number. """
        return [ self.rep.dom.to_sympy(c) for c in self.rep.all_coeffs() ]

    def to_algebraic_integer(self):
        """Convert `self` to an algebraic integer. """
        if self.minpoly.LC() == 1:
            return self
        else:
            raise NotImplementedError

