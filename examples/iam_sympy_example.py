"""Module to make life easier for SymPy examples"""

# hook in-tree sympy into python path if running example
import sys
import os.path
my_dir      = os.path.dirname(__file__)         # examples/
sympy_top   = os.path.split(my_dir)[0]          # ../
sympy_dir   = os.path.join(sympy_top, 'sympy')  # ../sympy/
if os.path.isdir(sympy_dir):
    #print 'working in-tree...'
    sys.path.insert(0, sympy_top)

