"""
A Printer for generating executable code.

The most important function here is srepr that returns a string so that the
relation eval(srepr(expr))=expr holds in an apropriate environment.
"""

from printer import Printer
from sympy.printing.precedence import precedence
import sympy.mpmath.lib as mlib

class ReprPrinter(Printer):
    def reprify(self, args, sep):
        return sep.join([self.doprint(item) for item in args])

    def emptyPrinter(self, expr):
        if isinstance(expr, str):
            return expr
        elif hasattr(expr, "__srepr__"):
            return expr.__srepr__()
        elif hasattr(expr, "args") and hasattr(expr.args, "__iter__"):
            l = []
            for o in expr.args:
                l.append(self._print(o))
            return expr.__class__.__name__ + '(%s)'%', '.join(l)
        elif hasattr(expr, "__module__") and hasattr(expr, "__name__"):
            return "<'%s.%s'>"%(expr.__module__, expr.__name__)
        else:
            return str(expr)

    def _print_Function(self, expr):
        r = '%s(%r)' % (expr.func.__base__.__name__, expr.func.__name__)
        r+= '(%s)' % ', '.join([self._print(a) for a in expr.args])
        return r

    def _print_Integer(self, expr):
        return '%s(%s)' % (expr.__class__.__name__, self._print(expr.p))

    def _print_list(self, expr):
        return "[%s]"%self.reprify(expr, ", ")

    def _print_Matrix(self, expr):
        s = expr._format_str(lambda elem: self._print(elem), ',\n  ')
        return '%s([\n  %s,\n])' % (expr.__class__.__name__, s)

    def _print_Poly(self, expr):
        terms = []

        for coeff, monom in expr.iter_terms():
            terms.append("(%s, %s)" % (self._print(coeff), monom))

        format = expr.__class__.__name__ + "([%s], %s, order='%s')"

        symbols = [ self._print(s) for s in expr.symbols ]

        return format % (', '.join(terms),
            ', '.join(symbols), expr.order)

    def _print_Polynomial(self, expr):
        return "Polynomial(%s, %s, %s, '%s')" % (self._print(expr.sympy_expr),
                  self._print(expr.coeffs), self._print(expr.var), self._print(expr.order))

    def _print_Rational(self, expr):
        return '%s(%s, %s)' % (expr.__class__.__name__, self._print(expr.p), self._print(expr.q))

    def _print_Real(self, expr):
        dps = mlib.prec_to_dps(expr._prec)
        r = mlib.to_str(expr._mpf_, dps+3)
        return '%s(%s, prec=%i)' % (expr.__class__.__name__, r, dps)

    def _print_Sum2(self, expr):
        return "Sum2(%s, (%s, %s, %s))" % (self._print(expr.f), self._print(expr.i),
                                           self._print(expr.a), self._print(expr.b))

    def _print_Symbol(self, expr):
        return "%s('%s')" % (expr.__class__.__name__, self._print(expr.name))

    def _print_tuple(self, expr):
        if len(expr)==1:
            return "(%s,)"%self._print(expr[0])
        else:
            return "(%s)"%self.reprify(expr, ", ")

    def _print_WildFunction(self, expr):
        return "%s('%s')" % (expr.__class__.__name__, expr.name)

RPrinter = ReprPrinter()

def srepr(s):
    return RPrinter.doprint(s)
