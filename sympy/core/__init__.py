"""Core module. Provides the basic operations needed in sympy.
"""

from basic import Basic
from symbol import Symbol, Wild, symbols
from numbers import Number, Real, Rational, Integer
from power import Pow
from mul import Mul
from add import Add
from relational import Equality, Inequality, Unequality, StrictInequality
from function import Lambda, Function, Apply, FApply, Composition, FPow, WildFunction, Derivative, DefinedFunction, diff
from interval import Interval

import defined_functions
import orthogonal_polynomials
import integer_sequences

# set repr output to pretty output:
Basic.set_interactive(True)

# expose singletons like exp, log, oo, I, etc.
for _n,_cls in Basic.singleton.items():
    exec '%s = _cls()' % (_n)

sympify = Basic.sympify
