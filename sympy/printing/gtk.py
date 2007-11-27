from sympy import Basic
from sympy.printing.mathml import mathml
import tempfile
import os

def print_gtk(x, start_viewer=True):
    """Print to Gtkmathview, a gtk widget capable of rendering MathML.
    Needs libgtkmathview-bin"""
    from sympy.utilities.mathml import c2p

    tmp = tempfile.mktemp() # create a temp file to store the result
    file = open(tmp, 'wb')
    file.write( c2p(mathml(x), simple=True) )
    file.close()

    if start_viewer:
        os.system("mathmlviewer " + tmp)
