from sympy import Basic
from sympy.utilities.mathml import c2p
from sympy.printing.mathml import mathml
import tempfile
import os

def print_gtk(x):
    """Print to Gtkmathview, a gtk widget capable of rendering MathML.
    Needs libgtkmathview-bin"""

    tmp = tempfile.mktemp() # create a temp file to store the result
    file = open(tmp, 'wb')

    file.write( c2p(mathml(x), simple=True) )
    file.close()

    os.system("mathmlviewer " + tmp)
