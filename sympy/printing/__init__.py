"""Printing subsystem"""

from pretty import *
from latex import latex, print_latex
from mathml import mathml, print_mathml
from python import python, print_python
from gtk import *

from preview import preview, view, pngview, pdfview, dviview

from str import StrPrinter
StrPrinter = StrPrinter()

from repr import srepr

# /cyclic/
from sympy.core import basic
basic.StrPrinter = StrPrinter
basic.ReprPrinter = ReprPrinter
del basic

from tree import print_tree
