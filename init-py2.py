from __future__ import division
from sympy import *

from sympy.plotting import plot_parametric, plot3d

# Output of last expression
_ = None

x, y, z, t = symbols('x,y,z,t')
k, m, n = symbols('k,m,n', integer=True)
f, g, h = symbols('f,g,h', cls=Function)

