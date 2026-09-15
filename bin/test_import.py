from __future__ import print_function

from timeit import default_timer as clock
from get_sympy import path_hack
path_hack()
t = clock()
import sympy # noqa: F401
t = clock() - t
print(t)
