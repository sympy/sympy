from __future__ import (absolute_import, division, print_function)
"""
C++ code printer
"""
from functools import wraps
from .ccode import C89CodePrinter, C99CodePrinter


# from http://en.cppreference.com/w/cpp/keyword
reserved = {
    'C++98': [
        'and', 'and_eq', 'asm', 'auto', 'bitand', 'bitor', 'bool', 'break',
        'case', 'catch,', 'char', 'class', 'compl', 'const', 'const_cast',
        'continue', 'default', 'delete', 'do', 'double', 'dynamic_cast',
        'else', 'enum', 'explicit', 'export', 'extern', 'false', 'float',
        'for', 'friend', 'goto', 'if', 'inline', 'int', 'long', 'mutable',
        'namespace', 'new', 'not', 'not_eq', 'operator', 'or', 'or_eq',
        'private', 'protected', 'public', 'register', 'reinterpret_cast',
        'return', 'short', 'signed', 'sizeof', 'static', 'static_cast',
        'struct', 'switch', 'template', 'this', 'throw', 'true', 'try',
        'typedef', 'typeid', 'typename', 'union', 'unsigned', 'using',
        'virtual', 'void', 'volatile', 'wchar_t', 'while', 'xor', 'xor_eq'
    ]
}

reserved['C++11'] = reserved['C++98'][:] + [
    'alignas', 'alignof', 'char16_t', 'char32_t', 'constexpr', 'decltype',
    'noexcept', 'nullptr', 'static_assert', 'thread_local'
]
reserved['C++17'] = []
# TM TS: atomic_cancel, atomic_commit, atomic_noexcept, synchronized
# concepts TS: concept, requires
# module TS: import, module


_math_functions = {
    'C++98': {
        'Mod': 'fmod',
        'ceiling': 'ceil',
    },
    'C++11': {
        'gamma': 'tgamma',
    },
    'C++17': {
        'beta': 'beta',
        'Ei': 'expint',
        'zeta': 'riemann_zeta',
    }
}

# from http://en.cppreference.com/w/cpp/header/cmath
for k in ('Abs', 'exp', 'log', 'log10', 'sqrt', 'sin', 'cos', 'tan',  # 'Pow'
          'asin', 'acos', 'atan', 'atan2', 'sinh', 'cosh', 'tanh', 'floor'):
    _math_functions['C++98'][k] = k.lower()


for k in ('asinh', 'acosh', 'atanh', 'erf', 'erfc'):
    _math_functions['C++11'][k] = k.lower()


def _attach_print_method(cls, sympy_name, func_name):
    meth_name = '_print_%s' % sympy_name
    if hasattr(cls, meth_name):
        raise ValueError("Edit method (or subclass) instead of overwriting.")
    def _print_method(self, expr):
        return '{0}{1}({2})'.format(self._ns, func_name, ', '.join(map(self._print, expr.args)))
    _print_method.__doc__ = "Prints code for %s" % k
    setattr(cls, meth_name, _print_method)


def _attach_print_methods(cls, cont):
    for sympy_name, cxx_name in cont[cls.standard].items():
        _attach_print_method(cls, sympy_name, cxx_name)


class _CXXCodePrinterBase(object):
    printmethod = "_cxxcode"
    language = 'C++'
    _ns = 'std::'  # namespace

    def __init__(self, settings=None):
        super(_CXXCodePrinterBase, self).__init__(settings or {})

    def _print_Max(self, expr):
        from sympy import Max
        if len(expr.args) == 1:
            return self._print(expr.args[0])
        return "%smax(%s, %s)" % (self._ns, expr.args[0], self._print(Max(*expr.args[1:])))

    def _print_Min(self, expr):
        from sympy import Min
        if len(expr.args) == 1:
            return self._print(expr.args[0])
        return "%smin(%s, %s)" % (self._ns, expr.args[0], self._print(Min(*expr.args[1:])))


class CXX98CodePrinter(_CXXCodePrinterBase, C89CodePrinter):
    standard = 'C++98'
    reserved_words = set(reserved['C++98'])


_attach_print_methods(CXX98CodePrinter, _math_functions)


class CXX11CodePrinter(_CXXCodePrinterBase, C99CodePrinter):
    standard = 'C++11'
    reserved_words = set(reserved['C++11'])


_attach_print_methods(CXX11CodePrinter, _math_functions)


class CXX17CodePrinter(_CXXCodePrinterBase, C99CodePrinter):
    standard = 'C++17'
    reserved_words = set(reserved['C++17'])


_attach_print_methods(CXX17CodePrinter, _math_functions)

cxx_code_printers = {
    'c++98': CXX98CodePrinter,
    'c++11': CXX11CodePrinter,
    'c++17': CXX17CodePrinter
}

def cxxcode(expr, assign_to=None, standard='c++11', **settings):
    """ C++ equivalent of :func:`sympy.ccode`. """
    return cxx_code_printers[standard.lower()](settings).doprint(expr, assign_to)
