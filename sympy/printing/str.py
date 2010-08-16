"""
A Printer for generating readable representation of most sympy classes.
"""

from printer import Printer
from sympy.printing.precedence import precedence, PRECEDENCE
from sympy.core.basic import S
from sympy.core.numbers import Rational
from sympy.core.power import Pow
from sympy.core.basic import Basic

import sympy.mpmath.libmp as mlib
from sympy.mpmath.libmp import prec_to_dps

from sympy.polys.polyerrors import PolynomialError

class StrPrinter(Printer):
    printmethod = "_sympystr"
    _default_settings = {
        "order": None,
        "full_prec": "auto",
    }

    def parenthesize(self, item, level):
        if precedence(item) <= level:
            return "(%s)"%self._print(item)
        else:
            return self._print(item)

    def stringify(self, args, sep, level=0):
        return sep.join([self.parenthesize(item, level) for item in args])

    def emptyPrinter(self, expr):
        if isinstance(expr, str):
            return expr
        elif isinstance(expr, Basic):
            if hasattr(expr, "args"):
                return repr(expr)
            else:
                raise
        else:
            return str(expr)

    def _print_Add(self, expr, order=None):
        # Now we need to sort the factors in Add, which are in "rest". Any
        # ordering is fine, but some ordering looks better and some looks bad.
        # This particular solution is slow, but it ensures a sane ordering. It
        # can of course be improved:

        if order is None and self.order is None:
            terms = sorted(expr.args, Basic._compare_pretty)
        else:
            terms = [ elt[-1] for elt in self.analyze(expr, order) ]

        PREC = precedence(expr)
        l = []
        for term in terms:
            t = self._print(term)
            if t.startswith('-'):
                sign = "-"
                t = t[1:]
            else:
                sign = "+"
            if precedence(term) < PREC:
                l.extend([sign, "(%s)"%t])
            else:
                l.extend([sign, t])
        sign = l.pop(0)
        if sign=='+':
            sign = ""
        return sign + ' '.join(l)

    def _print_Assume(self, expr):
        return 'Assume(%s, %r)' % (expr.expr, expr.key)

    def _print_Basic(self, expr):
        l = [self._print(o) for o in expr.args]
        return expr.__class__.__name__ + "(%s)"%", ".join(l)

    def _print_Catalan(self, expr):
        return 'Catalan'

    def _print_ComplexInfinity(self, expr):
        return 'zoo'

    def _print_Derivative(self, expr):
        return 'D(%s)'%", ".join(map(self._print, expr.args))

    def _print_dict(self, expr):
        keys = expr.keys()
        keys.sort( Basic.compare_pretty )

        items = []
        for key in keys:
            item = "%s: %s" % (self._print(key), self._print(expr[key]))
            items.append(item)

        return "{%s}"%", ".join(items)

    def _print_Dummy(self, expr):
        return '_' + expr.name

    def _print_EulerGamma(self, expr):
        return 'EulerGamma'

    def _print_Exp1(self, expr):
        return 'E'

    def _print_ExprCondPair(self, expr):
        return '(%s, %s)' % (expr.expr, expr.cond)

    def _print_Factorial(self, expr):
        return "%s!" % self.parenthesize(expr.args[0], PRECEDENCE["Pow"])

    def _print_Function(self, expr):
        return expr.func.__name__ + "(%s)"%self.stringify(expr.args, ", ")

    def _print_GeometryEntity(self, expr):
        # GeometryEntity is special -- it's base is tuple
        return str(expr)

    def _print_GoldenRatio(self, expr):
        return 'GoldenRatio'

    def _print_ImaginaryUnit(self, expr):
        return 'I'

    def _print_Infinity(self, expr):
        return 'oo'

    def _print_Integral(self, expr):
        def _xab_tostr(xab):
            x, ab = xab
            if ab is None:
                return self._print(x)
            else:
                return self._print((x,) + ab)
        L = ', '.join([_xab_tostr(l) for l in expr.limits])
        return 'Integral(%s, %s)' % (self._print(expr.function), L)

    def _print_Interval(self, i):
        if i.left_open:
            left = '('
        else:
            left = '['

        if i.right_open:
            right = ')'
        else:
            right = ']'

        return "%s%s, %s%s" % \
               (left, self._print(i.start), self._print(i.end), right)

    def _print_LatticeOp(self, expr):
        args = sorted(expr.args, cmp=expr._compare_pretty)
        return expr.func.__name__ + "(%s)"%", ".join(self._print(arg) for arg in args)

    def _print_Limit(self, expr):
        e, z, z0, dir = expr.args
        if dir == "+":
            return "Limit(%s, %s, %s)" % (e, z, z0)
        else:
            return "Limit(%s, %s, %s, dir='%s')" % (e, z, z0, dir)

    def _print_list(self, expr):
        return "[%s]"%self.stringify(expr, ", ")

    def _print_Matrix(self, expr):
        return expr._format_str(lambda elem: self._print(elem))

    def _print_DeferredVector(self, expr):
        return expr.name

    def _print_Mul(self, expr):
        coeff, terms = expr.as_coeff_terms()
        if coeff.is_negative:
            coeff = -coeff
            if coeff is not S.One:
                terms = (coeff,) + terms
            sign = "-"
        else:
            terms = (coeff,) + terms
            sign = ""

        a = [] # items in the numerator
        b = [] # items that are in the denominator (if any)

        # Gather terms for numerator/denominator
        for item in terms:
            if item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                b.append(Pow(item.base, -item.exp))
            elif item.is_Rational and item is not S.Infinity:
                if item.p != 1:
                    a.append(Rational(item.p))
                if item.q != 1:
                    b.append(Rational(item.q))
            else:
                a.append(item)

        if len(a)==0:
            a = [S.One]

        a_str = map(lambda x:self.parenthesize(x, precedence(expr)), a)
        b_str = map(lambda x:self.parenthesize(x, precedence(expr)), b)

        if len(b)==0:
            return sign + '*'.join(a_str)
        elif len(b)==1:
            if len(a)==1 and not (a[0].is_Atom or a[0].is_Add):
                return sign + "%s/"%a_str[0] + '*'.join(b_str)
            else:
                return sign + '*'.join(a_str) + "/%s"%b_str[0]
        else:
            return sign + '*'.join(a_str) + "/(%s)"%'*'.join(b_str)

    def _print_NaN(self, expr):
        return 'nan'

    def _print_NegativeInfinity(self, expr):
        return '-oo'

    def _print_Normal(self, expr):
        return "Normal(%s, %s)"%(expr.mu, expr.sigma)

    def _print_Order(self, expr):
        if len(expr.symbols) <= 1:
            return 'O(%s)'%self._print(expr.expr)
        else:
            return 'O(%s)'%self.stringify(expr.args, ', ', 0)

    def _print_PDF(self, expr):
        return 'PDF(%s, (%s, %s, %s))' % \
            (self._print(expr.pdf.args[1]), self._print(expr.pdf.args[0]), \
            self._print(expr.domain[0]), self._print(expr.domain[1]))

    def _print_Pi(self, expr):
        return 'pi'

    def _print_Poly(self, expr):
        terms, gens = [], [ self._print(s) for s in expr.gens ]

        for monom, coeff in expr.terms():
            s_monom = []

            for i, exp in enumerate(monom):
                if exp > 0:
                    if exp == 1:
                        s_monom.append(gens[i])
                    else:
                        s_monom.append(gens[i] + "**%d" % exp)

            s_monom = "*".join(s_monom)

            if coeff.is_Add:
                if s_monom:
                    s_coeff = "(" + self._print(coeff) + ")"
                else:
                    s_coeff = self._print(coeff)
            else:
                if s_monom:
                    if coeff is S.One:
                        terms.extend(['+', s_monom])
                        continue

                    if coeff is S.NegativeOne:
                        terms.extend(['-', s_monom])
                        continue

                s_coeff = self._print(coeff)

            if not s_monom:
                s_term = s_coeff
            else:
                s_term = s_coeff + "*" + s_monom

            if s_term.startswith('-'):
                terms.extend(['-', s_term[1:]])
            else:
                terms.extend(['+', s_term])

        if terms[0] in ['-', '+']:
            modifier = terms.pop(0)

            if modifier == '-':
                terms[0] = '-' + terms[0]

        format = expr.__class__.__name__ + "(%s, %s"

        try:
            format += ", modulus=%s" % expr.get_modulus()
        except PolynomialError:
            format += ", domain='%s'" % expr.get_domain()

        format += ")"

        return format % (' '.join(terms), ', '.join(gens))

    def _print_AlgebraicNumber(self, expr):
        if expr.is_aliased:
            return self._print(expr.as_poly().as_basic())
        else:
            return self._print(expr.as_basic())

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp is S.NegativeOne:
            return '1/%s'%(self.parenthesize(expr.base, PREC))
        else:
            return '%s**%s'%(self.parenthesize(expr.base, PREC),
                             self.parenthesize(expr.exp, PREC))

    def _print_Integer(self, expr):
        return str(expr.p)

    def _print_int(self, expr):
        return str(expr)

    def _print_mpz(self, expr):
        return str(expr)

    def _print_Rational(self, expr):
        return '%s/%s' % (expr.p, expr.q)

    def _print_Fraction(self, expr):
        return '%s/%s' % (expr.numerator, expr.denominator)

    def _print_mpq(self, expr):
        return '%s/%s' % (expr.numer(), expr.denom())

    def _print_Real(self, expr):
        prec = expr._prec
        if prec < 5:
            dps = 0
        else:
            dps = prec_to_dps(expr._prec)
        if self._settings["full_prec"] == True:
            strip = False
        elif self._settings["full_prec"] == False:
            strip = True
        elif self._settings["full_prec"] == "auto":
            strip = self._print_level > 1
        return mlib.to_str(expr._mpf_, dps, strip_zeros=strip)

    def _print_Relational(self, expr):
        return '%s %s %s'%(self.parenthesize(expr.lhs, precedence(expr)),
                           expr.rel_op,
                           self.parenthesize(expr.rhs, precedence(expr)))

    def _print_DMP(self, expr):
        return "%s(%s, %s)" % (expr.__class__.__name__, self._print(expr.rep), self._print(expr.dom))

    def _print_RootOf(self, expr):
        return "RootOf(%s, %d)" % (self._print_Add(expr.expr, order='lex'), expr.index)

    def _print_RootSum(self, expr):
        args = [self._print_Add(expr.expr, order='lex')]

        if not (isinstance(expr.func, Basic) and expr.func.is_identity):
            args.append(self._print(expr.func))

        return "RootSum(%s)" % ", ".join(args)

    def _print_Sample(self, expr):
        return "Sample([%s])"%self.stringify(expr, ", ", 0)

    def __print_set(self, expr):
        items = list(expr)
        items.sort( Basic.compare_pretty )

        args = ', '.join(self._print(item) for item in items)
        if args:
            args = '[%s]' % args
        return '%s(%s)' % (type(expr).__name__, args)

    _print_set       = __print_set
    _print_frozenset = __print_set

    def _print_SMatrix(self, expr):
        return self._print(expr.toMatrix())

    def _print_Sum(self, expr):
        limits = str.join(', ', \
            tuple([ '(%s)' % self.stringify(l, ', ') for l in expr.limits ]))

        return 'Sum(%s, %s)' % (self._print(expr.function), limits)

    def _print_Sum2(self, expr):
        return "Sum2(%r, (%r, %r, %r))" % (expr.f, expr.i, expr.a, expr.b)

    def _print_Symbol(self, expr):
        return expr.name

    def _print_Predicate(self, expr):
        return "Q.%s" % expr.name

    def _print_str(self, expr):
        return expr

    def _print_tuple(self, expr):
        if len(expr)==1:
            return "(%s,)"%self._print(expr[0])
        else:
            return "(%s)"%self.stringify(expr, ", ")

    def _print_Uniform(self, expr):
        return "Uniform(%s, %s)"%(expr.a, expr.b)

    def _print_Unit(self, expr):
        return expr.abbrev

    def _print_Wild(self, expr):
        return expr.name + '_'

    def _print_WildFunction(self, expr):
        return expr.name + '_'

    def _print_Zero(self, expr):
        return "0"


def sstr(expr, **settings):
    """return expr in str form"""

    p = StrPrinter(settings)
    s = p.doprint(expr)

    return s


class StrReprPrinter(StrPrinter):
    """(internal) -- see sstrrepr"""

    def _print_str(self, s):
        return repr(s)

def sstrrepr(expr, **settings):
    """return expr in mixed str/repr form

       i.e. strings are returned in repr form with quotes, and everything else
       is returned in str form.

       This function could be useful for hooking into sys.displayhook
    """

    p = StrReprPrinter(settings)
    s = p.doprint(expr)

    return s

