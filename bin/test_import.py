from __future__ import division, print_function

from timeit import default_timer as clock
from get_sympy import path_hack
path_hack()
t = clock()
import sympy
t = clock() - t
print(t)
