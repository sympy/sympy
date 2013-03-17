from sympy.utilities import default_sort_key

from sympy.printing.printer import Printer

import theano
Scalar = theano.scalar.Scalar


class TheanoPrinter(Printer):
    """ Code printer for Theano computations """
    printmethod = "_theano"

    def _print_Symbol(self, s, dtypes={}):
        dtype = dtypes.get(s, 'floatX')
        return Scalar(dtype)(s.name)

    def doprint(self, expr):
        """Returns printer's representation for expr (as a string)"""
        return self._print(expr)

def theano_code(expr, **settings):
    return TheanoPrinter(settings).doprint(expr)
