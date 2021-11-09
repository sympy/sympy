from sympy.external.importtools import version_tuple
from collections.abc import Iterable

from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.codegen.cfunctions import Sqrt
from sympy.external import import_module
from sympy.printing.precedence import PRECEDENCE
from sympy.printing.pycode import AbstractPythonCodePrinter
import sympy

tensorflow = import_module('paddle')

class PaddlePrinter(AbstractPythonCodePrinter):
  '''
  for using paddle in sympy
  '''
  
