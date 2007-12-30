#!/usr/bin/env python
import iam_sympy_example

import sympy
a=sympy.Symbol('a')
b=sympy.Symbol('b')
e=(a+2*b)**5
print e
print e.diff(a)
print e.diff(b)
print e.diff(b).diff(a,3)
print e.expand().diff(b).diff(a,3)
