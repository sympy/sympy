#!/usr/bin/env python
import iam_sympy_example

from sympy import Symbol, cos, sin, Plot, log, tan
from sympy.abc import x, y

Plot(cos(x)*sin(y), sin(x)*sin(y), cos(y)+log(tan(y/2))+0.2*x, [x, -0.00,
    12.4, 40], [y, 0.1, 2, 40])
