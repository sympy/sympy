from sympy.abc import x, y
from sympy import ask, Q
from timeit import timeit
 
print(timeit(lambda : ask(Q.eq(x, 1) | Q.eq(y, -1), Q.integer(x) & Q.integer(y)), number=1))
print(timeit(lambda : ask(Q.eq(x, 1), Q.integer(x) & Q.integer(y)), number=1))
print(timeit(lambda : ask(Q.eq(x, 1) | Q.eq(y, -1)), number=1))
