from sympy import *

x, y, z, t = symbols('x,y,z,t')
k, m, n = symbols('k,m,n', integer=True)
f, g, h = map(Function, 'fgh')

from printing import init_printing
from session import init_session

