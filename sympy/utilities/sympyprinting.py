"""
A print function that pretty prints sympy Basic objects.

:moduleauthor: Brian Granger

Usage
=====

Once the extension is loaded, Sympy Basic objects are automatically
pretty-printed.

"""
#-----------------------------------------------------------------------------
#  Copyright (C) 2008  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from sympy import pretty, latex
from sympy.interactive.printing import init_printing

#-----------------------------------------------------------------------------
# Definitions of special display functions for use with IPython
#-----------------------------------------------------------------------------

_loaded = False

def load_ipython_extension(ip):
    """Load the extension in IPython."""
    global _loaded
    # Use extension manager to track loaded status if available
    # This is currently in IPython 0.14.dev
    if hasattr(ip.extension_manager, 'loaded'):
        loaded = 'sympy.utilities.sympyprinting' not in ip.extension_manager.loaded
    else:
        loaded = _loaded

    if not loaded:
        init_printing(ip=ip)

        _loaded = True
