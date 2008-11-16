#!/usr/bin/env python
import iam_sympy_example

import sympy
a=sympy.Symbol('a')
b=sympy.Symbol('b')
e=(a+b)**5
print e
print e.expand()
