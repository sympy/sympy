"""
SymPy is a Python library for symbolic mathematics. It aims to become a
full-featured computer algebra system (CAS) while keeping the code as simple
as possible in order to be comprehensible and easily extensible.  SymPy is
written entirely in Python. It depends on mpmath, and other external libraries
may be optionally for things like plotting support.

See the webpage for more information and documentation:

    http://sympy.org

"""

from __future__ import absolute_import, print_function

__all__ = []

try:
    import mpmath
except ImportError:
    raise ImportError("SymPy now depends on mpmath as an external library. "
    "See http://docs.sympy.org/latest/install.html#mpmath for more information.")


from sympy.release import __version__

if 'dev' in __version__:
    def enable_warnings():
        import warnings
        warnings.filterwarnings('default',   '.*',   DeprecationWarning, module='sympy.*')
    enable_warnings()


import sys
if ((sys.version_info[0] == 2 and sys.version_info[1] < 7) or
    (sys.version_info[0] == 3 and sys.version_info[1] < 4)):
    raise ImportError("Python version 2.7 or 3.4 or above "
                      "is required for SymPy.")


def __sympy_debug():
    # helper function so we don't import os globally
    import os
    debug_str = os.getenv('SYMPY_DEBUG', 'False')
    if debug_str in ('True', 'False'):
        return eval(debug_str)
    else:
        raise RuntimeError("unrecognized value for SYMPY_DEBUG: %s" %
                           debug_str)
SYMPY_DEBUG = __sympy_debug()
__all__ += ["SYMPY_DEBUG"]

from . import core
__all__ += core.__all__

from . import logic
__all__ += logic.__all__

from . import assumptions
__all__ += assumptions.__all__

from . import polys
__all__ += polys.__all__

from . import series
__all__ += series.__all__

from . import functions
__all__ += functions.__all__

from . import ntheory
__all__ += ntheory.__all__

from . import concrete
__all__ += concrete.__all__

from . import simplify
__all__ += simplify.__all__

from . import sets
__all__ += sets.__all__

from . import solvers
__all__ += solvers.__all__

from . import matrices
__all__ += matrices.__all__

from . import geometry
__all__ += geometry.__all__

from . import utilities
__all__ += utilities.__all__

from . import integrals
__all__ += integrals.__all__

from . import tensor
__all__ += tensor.__all__

from . import parsing
__all__ += parsing.__all__

from . import calculus
__all__ += calculus.__all__

from . import algebras
__all__ += algebras.__all__

# Adds about .04-.05 seconds of import time
# from . import combinatorics
# __all__ += combinatorics.__all__

# This module is slow to import:
# from .physics import units
# __all__ += ["units"]

from .plotting import plot, textplot, plot_backends, plot_implicit
__all__ += ["plot", "textplot", "plot_backends", "plot_implicit"]

from .printing import (
    pretty, pretty_print, pprint, pprint_use_unicode,
    pprint_try_use_unicode, print_gtk, print_tree, pager_print, TableForm
)
__all__ += [
    "pretty", "pretty_print", "pprint", "pprint_use_unicode",
    "pprint_try_use_unicode", "print_gtk", "print_tree", "pager_print", "TableForm"
]

from .printing import (
    rcode, ccode, fcode, jscode, julia_code, mathematica_code,
    octave_code, latex, preview, rust_code, mathml, glsl_code, cxxcode
)
__all__ += [
    "rcode", "ccode", "fcode", "jscode", "julia_code", "mathematica_code",
    "octave_code", "latex", "preview", "rust_code", "mathml", "glsl_code", "cxxcode"
]

from .printing import python, print_python, srepr, sstr, sstrrepr
__all__ += ["python", "print_python", "srepr", "sstr", "sstrrepr"]

from .interactive import init_session, init_printing
__all__ += ["init_session", "init_printing"]

from .core.evalf import _create_evalf_table
_create_evalf_table()

# This is slow to import:
#import abc

from . import deprecated
__all__ += deprecated.__all__
