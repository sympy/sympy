"""Computational algebraic number field theory. """

from sympy import (
    S, Basic, I, Integer, Rational, Real,
    Symbol, Add, Mul, sympify, Q, ask,
)

from sympy.polys.polytools import (
    Poly, sqf_norm, invert, factor_list, groebner,
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

def minimal_polynomial(ex, x=None, **args):
    """Computes the minimal polynomial of an algebraic number. """
    generator = numbered_symbols('a', dummy=True)
    mapping, symbols, replace = {}, {}, []

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
                if ex.exp < 0 and ex.base.is_Add:
                    coeff, terms = ex.base.as_coeff_factors()
                    elt, _ = primitive_element(terms, polys=True)

                    alg = ex.base - coeff

                    inverse = invert(elt.gen + coeff, elt)
                    base = inverse.subs(elt.gen, alg).expand()

                    if ex.exp == -1:
                        return bottom_up_scan(base)
                    else:
                        ex = base**(-ex.exp)

                if not ex.exp.is_Integer:
                    base, exp = (ex.base**ex.exp.p).expand(), Rational(1, ex.exp.q)
                else:
                    base, exp = ex.base, ex.exp

                base = bottom_up_scan(base)
                expr = base**exp

                if expr not in mapping:
                    return update_mapping(expr, 1/exp, -base)
                else:
                    return symbols[expr]
        elif ex.is_AlgebraicNumber:
            if ex.root not in mapping:
                return update_mapping(ex.root, ex.minpoly)
            else:
                return symbols[ex.root]

        raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)

    polys = args.get('polys', False)

    if ex.is_AlgebraicNumber:
        if not polys:
            return ex.minpoly.as_basic(x)
        else:
            return ex.minpoly.replace(x)
    elif ex.is_Rational:
        result = ex.q*x - ex.p
    else:
        F = [x - bottom_up_scan(ex)] + mapping.values()
        G = groebner(F, symbols.values() + [x], order='lex')

        _, factors = factor_list(G[-1])

        if len(factors) == 1:
            ((result, _),) = factors
        else:
            for result, _ in factors:
                if result.subs(x, ex).evalf(chop=True) == 0:
                    break
            else: # pragma: no cover
                raise NotImplementedError("multiple candidates for the minimal polynomial of %s" % ex)

    if polys:
        return Poly(result, x, field=True)
    else:
        return result

minpoly = minimal_polynomial

def primitive_element(extension, x=None, **args):
    """Construct a common number field for all extensions. """
    if not extension:
        raise ValueError("can't compute primitive element for nothing")

    extension = [ AlgebraicNumber(ext, gen=x) for ext in extension ]

    minpoly = extension[0].minpoly
    root    = extension[0].root

    for ext in extension[1:]:
        s, _, minpoly = sqf_norm(minpoly, extension=ext)
        root = s*root + ext.root

    if not args.get('polys', False):
        return minpoly.as_basic(), root
    else:
        return minpoly, root

primelt = primitive_element

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

    _, factors = factor_list(a.minpoly, extension=b)

    for f, _ in factors:
        if f.degree() == 1:
            coeffs = f.rep.TC().to_sympy_list()
            d, terms = len(coeffs)-1, []

            for i, coeff in enumerate(coeffs):
                terms.append(coeff*b.root**(d-i))

            root = Add(*terms)

            if (a.root - root).evalf(chop=True) == 0:
                return coeffs

            if (a.root + root).evalf(chop=True) == 0:
                return [ -c for c in coeffs ]
    else:
        return None

def to_number_field(extension, theta=None, **args):
    """Express `extension` in the field generated by `theta`. """
    gen = args.get('gen')

    if hasattr(extension, '__iter__'):
        extension = list(extension)
    else:
        extension = [extension]

    minpoly, root = primitive_element(extension, gen, polys=True)

    if theta is None:
        return AlgebraicNumber((minpoly, root))
    else:
        theta = sympify(theta)

        if not theta.is_AlgebraicNumber:
            theta = AlgebraicNumber(theta, gen=gen)

        coeffs = field_isomorphism(root, theta)

        if coeffs is not None:
            return AlgebraicNumber(theta, coeffs)
        else:
            raise IsomorphismFailed("%s is not in a subfield of %s" % (root, theta.root))

class AlgebraicNumber(Basic):
    """Class for representing algebraic numbers in SymPy. """

    __slots__ = ['rep', 'root', 'alias', 'minpoly']

    is_AlgebraicNumber = True

    def __new__(cls, expr, coeffs=None, **args):
        """Construct a new algebraic number. """
        expr = sympify(expr)

        if type(expr) is tuple:
            minpoly, root = expr

            if not minpoly.is_Poly:
                minpoly = Poly(minpoly)
        elif expr.is_AlgebraicNumber:
            minpoly, root = expr.minpoly, expr.root
        else:
            minpoly, root = minimal_polynomial(expr, args.get('gen'), polys=True), expr

        dom = minpoly.get_domain()

        if coeffs is not None:
            if not isinstance(coeffs, ANP):
                rep = DMP.from_sympy_list(sympify(coeffs), 0, dom)
            else:
                rep = DMP.from_list(coeffs.to_list(), 0, dom)

            if rep.degree() >= minpoly.degree():
                rep = rep.rem(minpoly.rep)
        else:
            rep = DMP.from_list([1, 0], 0, dom)

            if ask(root, Q.negative):
                rep = -rep

        alias = args.get('alias')

        if alias is not None:
            if not isinstance(alias, Symbol):
                alias = Symbol(alias)

        obj = Basic.__new__(cls)

        obj.rep = rep
        obj.root = root
        obj.alias = alias
        obj.minpoly = minpoly

        return obj

    def __eq__(a, b):
        if not b.is_AlgebraicNumber:
            try:
                b = to_number_field(b, a)
            except (NotAlgebraic, IsomorphismFailed):
                return False

        return a.rep == b.rep and \
            a.minpoly.all_coeffs() == b.minpoly.all_coeffs()

    def __ne__(a, b):
        if not b.is_AlgebraicNumber:
            try:
                b = to_number_field(b, a)
            except (NotAlgebraic, IsomorphismFailed):
                return True

        return a.rep != b.rep or \
            a.minpoly.all_coeffs() != b.minpoly.all_coeffs()

    def _eval_evalf(self, prec):
        return self.as_basic()._evalf(prec)

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
        """Returns all SymPy coefficients of an algebraic number. """
        return [ self.rep.dom.to_sympy(c) for c in self.rep.all_coeffs() ]

    def native_coeffs(self):
        """Returns all native coefficients of an algebraic number. """
        return self.rep.all_coeffs()

    def to_algebraic_integer(self):
        """Convert `self` to an algebraic integer. """
        f = self.minpoly

        if f.LC() == 1:
            return self

        coeff = f.LC()**(f.degree()-1)
        poly  = f.compose(Poly(f.gen/f.LC()))

        minpoly = poly*coeff
        root    = f.LC()*self.root

        return AlgebraicNumber((minpoly, root), self.coeffs())

