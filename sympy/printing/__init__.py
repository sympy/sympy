"""Printing subsystem"""

from sympy.core.cache import lazy_function, lazy_functions

from .pretty import pager_print, pretty, pretty_print, pprint, pprint_use_unicode, pprint_try_use_unicode
from .pycode import pycode # not lazy becaues pycode.pycode generates a name clash
from .str import StrPrinter, sstr, sstrrepr
from .tableform import TableForm

# lazy import of functions
latex, print_latex, multiline_latex = lazy_functions('sympy.printing.latex', ['latex', 'print_latex', 'multiline_latex'])
mathml, print_mathml= lazy_functions('sympy.printing.mathml', ['mathml', 'print_mathml'])
python, print_python = lazy_functions('sympy.printing.python', ['python', 'print_python'])
print_ccode, print_fcode= lazy_functions('sympy.printing.codeprinter', ['print_ccode', 'print_fcode'])
ccode, fcode, cxxcode = lazy_functions('sympy.printing.codeprinter', ['ccode', 'fcode', 'cxxcode'])
glsl_code = lazy_function('sympy.printing.glsl', 'glsl_code')
print_glsl = lazy_function('sympy.printing.glsl', 'print_glsl')
rcode, print_rcode = lazy_functions('sympy.printing.rcode', ['rcode', 'print_rcode'])
jscode, print_jscode = lazy_functions('sympy.printing.jscode', ['jscode', 'print_jscode'])
julia_code = lazy_function('sympy.printing.julia', 'julia_code')
mathematica_code = lazy_function('sympy.printing.mathematica', 'mathematica_code')
octave_code = lazy_function('sympy.printing.octave', 'octave_code')
rust_code = lazy_function('sympy.printing.rust', 'rust_code')
print_gtk = lazy_function('sympy.printing.gtk', 'print_gtk')
preview = lazy_function('sympy.printing.preview', 'preview')
srepr = lazy_function('sympy.printing.repr', 'srepr')
print_tree = lazy_function('sympy.printing.tree', 'print_tree')
dotprint = lazy_function('sympy.printing.dot', 'dotprint')
maple_code, print_maple_code= lazy_functions('sympy.printing.maple', ['maple_code', 'print_maple_code'])


__all__ = [
    # sympy.printing.pretty
    'pager_print', 'pretty', 'pretty_print', 'pprint', 'pprint_use_unicode',
    'pprint_try_use_unicode',

    # sympy.printing.latex
    'latex', 'print_latex', 'multiline_latex',

    # sympy.printing.mathml
    'mathml', 'print_mathml',

    # sympy.printing.python
    'python', 'print_python',

    # sympy.printing.pycode
    'pycode',

    # sympy.printing.codeprinter
    'ccode', 'print_ccode', 'cxxcode', 'fcode', 'print_fcode',

    # sympy.printing.smtlib
    'smtlib_code',

    # sympy.printing.glsl
    'glsl_code', 'print_glsl',

    # sympy.printing.rcode
    'rcode', 'print_rcode',

    # sympy.printing.jscode
    'jscode', 'print_jscode',

    # sympy.printing.julia
    'julia_code',

    # sympy.printing.mathematica
    'mathematica_code',

    # sympy.printing.octave
    'octave_code',

    # sympy.printing.rust
    'rust_code',

    # sympy.printing.gtk
    'print_gtk',

    # sympy.printing.preview
    'preview',

    # sympy.printing.repr
    'srepr',

    # sympy.printing.tree
    'print_tree',

    # sympy.printing.str
    'StrPrinter', 'sstr', 'sstrrepr',

    # sympy.printing.tableform
    'TableForm',

    # sympy.printing.dot
    'dotprint',

    # sympy.printing.maple
    'maple_code', 'print_maple_code',
]
