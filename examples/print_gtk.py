"""examples for print_gtk. It prints in gtkmathview using mathml"""

import sys
sys.path.append("..")

from sympy import *
from sympy.printing import print_gtk

x = Symbol('x')

#l1 = limit(sin(x)/x, x, 0, evaluate=False)
#print_gtk(l1)

l2 = integrate(exp(x), (x,0,1), evaluate=False)
print_gtk(l2)
