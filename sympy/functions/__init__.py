
from sympy.core.basic import Basic

import combinatorial
import elementary
import special

from elementary.trigonometric import acot, cot

for _n, _cls in Basic.singleton.items():
    exec '%s = _cls()' % (_n)
