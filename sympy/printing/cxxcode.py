from __future__ import (absolute_import, division, print_function)
"""
C++ code printer
"""


from .ccode import (
    CCodePrinter, C99CodePrinter, known_functions as _c_known_functions,
    reserved_words as _c_reserved_words
)

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

reserved['C++11'] = reserved['C++98'][:]
reserved['C++11'].extend([
    'alignas', 'alignof', 'char16_t', 'char32_t', 'constexpr', 'decltype',
    'noexcept', 'nullptr', 'static_assert', 'thread_local'
])
reserved['C++17'] = []
# TM TS: atomic_cancel, atomic_commit, atomic_noexcept, synchronized
# concepts TS: concept, requires
# module TS: import, module


_math_functions = {
    'C++98': {
        'Mod': 'std::fmod',
        'ceiling': 'std::ceil',
    }
}

# from http://en.cppreference.com/w/cpp/header/cmath
for k in ('Abs', 'exp', 'log', 'log10', 'sqrt', 'sin', 'cos', 'tan',  # 'Pow'
          'asin', 'acos', 'atan', 'atan2', 'sinh', 'cosh', 'tanh', 'floor'):
    _math_functions['C++98'][k] = 'std::' + k.lower()

_math_functions['C++11'] = {
    'gamma': 'tgamma',
}

for k in ('asinh', 'acosh', 'atanh', 'erf', 'erfc'):
    _math_functions['C++11'][k] = 'std::' + k.lower()


def _attach_print_method(cls, k, v):
    method_name = '_print_' + k
    if hasattr(cls, method_name):
        raise ValueError("Edit method (or subclass) instead of overwriting.")
    setattr(cls, method_name, lambda self, expr:
            v+'(' + ', '.join(map(self._print, expr.args)) + ')')

def _attach_methods(cls, cont):
    for k, v in cont[cls.standard].items():
        _attach_print_method(cls, k, v)


class _CXXCodePrinterBase(object):
    language = 'C++'
    _ns = 'std::'

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


class CXX98CodePrinter(_CXXCodePrinterBase, CCodePrinter):
    standard = 'C++98'
    reserved_words = set(reserved['C++98'])


_attach_methods(CXX98CodePrinter, _math_functions)


class CXX11CodePrinter(_CXXCodePrinterBase, C99CodePrinter):
    standard = 'C++11'
    reserved_words = set(reserved['C++11'])


_attach_methods(CXX11CodePrinter, _math_functions)


# C++17 standard not finalized, below is just an example (unoffical API)
class _CXX17CodePrinter(_CXXCodePrinterBase, C99CodePrinter):
    standard = 'C++17'

    def __init__(self, settings=None):
        super(_CXX17CodePrinter, self).__init__(settings or {})
        self.reserved_words = set(reserved[self.standard])

    def _print_beta(self, expr):
        return '{0}beta({1}, {2})'.format(self._ns, *map(self._print, expr.args))

    def _print_Ei(self, expr):
        return '{0}expint({1})'.format(self._ns, self._print(expr.args[0]))

    def _print_zeta(self, expr):
        return '{0}riemann_zeta({1})'.format(self._ns, self._print(expr.args[0]))

cxx_code_printers = {
    'c++98': CXX98CodePrinter,
    'c++11': CXX11CodePrinter,
    'c++17': _CXX17CodePrinter
}

def cxxcode(expr, assign_to=None, standard='c++98', **settings):
    """ C++ equivalent of :func:`sympy.ccode`. """
    return cxx_code_printers[standard.lower()](settings).doprint(expr, assign_to)
