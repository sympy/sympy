# -*- coding: utf-8 -*-

from sympy.core import C
from sympy import oo
from printer import Printer

class PythonPrinter(Printer):
    """A printer which converts an expression into its Python interpretation."""

    def __init__(self):
        Printer.__init__(self)
        self.symbols = []
        self.functions = []

    # Functions for printing exact type of expression
    def _print_Add(self, expr):
        for value in expr.args:
            self._print(value)

    def _print_Derivative(self, expr):
        for arg in expr.args:
            self._print(arg)

    # procedure (!) for defining functions which have to be defined in print_python()
    def _print_Function(self, expr, exp=None):
        func = expr.func.__name__
        args = ''.join([self._print(arg) for arg in expr.args])
        import sympy
        if not hasattr(sympy, func) and func not in self.functions:
            self.functions.append(func)

    def _print_Integral(self, expr):
        for arg in expr._args:
            self._print(arg)

    def _print_Mul(self, expr):
        terms = expr.as_coeff_terms()[1]
        for term in terms:
            self._print(term)

    def _print_Pow(self, expr):
        if expr.exp.is_Rational and expr.exp.q == 2:
            self._print(expr.base)
        else:
            if expr.base.is_Function or expr.base.is_Symbol:
                self._print(expr.base)
                self._print(expr.exp)

    def _print_Rational(self, expr):
        if expr == oo:
            return oo
        if expr == 0:
            return 0
        return 'Rational(' + self._str(expr.p) + ', ' + self._str(expr.q) + ')'

    def _print_Relational(self, expr):
        self._print(expr.lhs) + self._print(expr.rhs)

    # procedure (!) for defining symbols which have be defined in print_python()
    def _print_Symbol(self, expr):
        symbol = self._str(expr)
        if symbol not in self.symbols:
            self.symbols.append(symbol)

    # Main function for printing
    def doprint(self, expr):
        return Printer.doprint(self, expr)

    # Function definitions for clearer output
    # (it's possible they don't have to be defined)
    def _print_int(self, expr):
        pass

    def _print_module(self, expr):
        raise Exception('Modules in the expression are unacceptable')

    def _print_str(self, expr):
        raise Exception('Strings in the expression are unacceptable')

    def emptyPrinter(self, expr):
        pass

def python(expr):
    """Return Python interpretation of passed expression
    (can be passed to the exec() function without any modifications)"""

    printer = PythonPrinter()
    expr = printer.doprint(expr)

    result = ''
    # Returning found symbols and functions
    for symbol in printer.symbols:
        result += symbol + ' = Symbol(\'' + symbol + '\')\n'
    for function in printer.functions:
        result += function + ' = Function(\'' + function + '\')\n'

    result += 'e = ' + printer._str(expr)
    return result

def print_python(expr):
    """Print output of python() function"""
    print python(expr)
