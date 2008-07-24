# -*- coding: utf-8 -*-

import sympy
from repr import ReprPrinter
from str import StrPrinter

class PythonPrinter(ReprPrinter, StrPrinter):
    """A printer which converts an expression into its Python interpretation."""

    def __init__(self):
        ReprPrinter.__init__(self)
        StrPrinter.__init__(self)
        self.symbols = []
        self.functions = []

    def _print_Add(self, expr):
        return StrPrinter._print_Add(self, expr)

    def _print_Function(self, expr):
        func = expr.func.__name__
        if not hasattr(sympy, func) and not func in self.functions:
            self.functions.append(func)
        return StrPrinter._print_Function(self, expr)

    def _print_Infinity(self, expr):
        return StrPrinter._print_Infinity(self, expr)

    def _print_Integer(self, expr):
        return StrPrinter._print_Integer(self, expr)

    def _print_Mul(self, expr):
        return StrPrinter._print_Mul(self, expr)

    def _print_Pow(self, expr):
        return StrPrinter._print_Pow(self, expr)

    # procedure (!) for defining symbols which have be defined in print_python()
    def _print_Symbol(self, expr):
        symbol = self._str(expr)
        if symbol not in self.symbols:
            self.symbols.append(symbol)
        return StrPrinter._print_Symbol(self, expr)

    def _print_module(self, expr):
        raise Exception('Modules in the expression are unacceptable')

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
