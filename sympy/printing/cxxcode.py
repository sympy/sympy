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
        'typedef', 'typeid', 'typename', 'union', 'unsigned', 'using(1)',
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
    }
}

# from http://en.cppreference.com/w/cpp/header/cmath
for k in ('Abs', 'Pow', 'exp', 'log', 'log10', 'sqrt', 'sin', 'cos', 'tan',
          'asin', 'acos', 'atan', 'atan2', 'sinh', 'cosh', 'tanh'):
    _math_functions['C++98'][k] = 'std::' + k.lower()

_math_functions['C++11'] = {
    'gamma': 'tgamma',
    'loggamma': 'lgamma'
}

for k in ('fma', 'fmax', 'fmin', 'exp2', 'expm1', 'log2', 'log1p', 'Cbrt',
          'hypot', 'asinh', 'acosh', 'atanh', 'erf', 'erfc'):
    _math_functions['C++11'][k] = 'std::' + k.lower()


def _attach_print_method(cls, k, v):
    setattr(cls, '_print_' + k, lambda self, expr:
            v+'(' + ', '.join(map(self._print, expr.args)) + ')')

def _attach_methods(cls, cont):
    for k, v in cont[cls.language].items():
        _attach_print_method(cls, k, v)



class CXX98CodePrinter(CCodePrinter):
    language = 'C++98'

    def __init__(self, settings=None):
        super(CXX98CodePrinter, self).__init__(settings or {})
        self.reserved_words = set(reserved[self.language])


_attach_methods(CXX98CodePrinter, _math_functions)


class CXX11CodePrinter(C99CodePrinter):
    language = 'C++11'

    def __init__(self, settings=None):
        super(CXX11CodePrinter, self).__init__(settings or {})
        self.reserved_words = set(reserved[self.language])



_attach_methods(CXX11CodePrinter, _math_functions)


# C++17 standard not finalized, below is just an example (unoffical API)
class _CXX17CodePrinter(C99CodePrinter):
    language = 'C++17'

    def __init__(self, settings=None):
        super(_CXX17CodePrinter, self).__init__(settings or {})
        self.reserved_words = set(reserved[self.language])

    def _print_beta(self, expr):
        return 'std::beta({0}, {1})'.format(*map(self._print, expr.args))

    def _print_Ei(self, expr):
        return 'std::expint({0})'.format(self._print(expr.args[0]))

    def _print_zeta(self, expr):
        return 'std::riemann_zeta({0})'.format(self._print(expr.args[0]))
