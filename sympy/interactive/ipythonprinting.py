"""
A print function that pretty prints SymPy objects.

:moduleauthor: Brian Granger

Usage
=====

To use this extension, execute:

    %load_ext sympy.interactive.ipythonprinting

Once the extension is loaded, SymPy Basic objects are automatically
pretty-printed in the terminal and rendered in LaTeX in the Qt console and
notebook.

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

from sympy.interactive.printing import init_printing

#-----------------------------------------------------------------------------
# Definitions of special display functions for use with IPython
#-----------------------------------------------------------------------------

def load_ipython_extension(ip):
    """Load the extension in IPython."""
    init_printing(ip=ip)
