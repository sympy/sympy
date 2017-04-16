from __future__ import division, print_function

from sympy import *

n = Symbol('n')

expr = cos(1/n) - 1
expr = (n - 1)/(n*log(n)**3) 

# s = Sum(expr, (n, 1, oo))

# s1 = s.is_convergent()

s = Sum(expr, (n, 1, oo))

s1 = s.is_convergent()




# expr4 = 1/(n + log(n))
# s1 = Sum(expr4, (n, 1, oo)).is_convergent()







# Sum(expr3, (n, 3, oo)).is_convergent()



# expr2 = expr4.subs(n, 1/n)
# s = expr2.series(n, 0, 10)
#         if type(s) is Add:
#             return s.subs(n, 1/n)
#         else:
#             return expr


# expr3 = (n - 1)/(n*log(n)**3) 
# Sum(expr3, (n, 3, oo)).is_convergent()

# print(Sum(expr3, (n, 3, oo)).is_convergent())









# expr2 = expr.subs(n, 1/n)

# t = expr2.series(n, 0, 4)
# t1 = expr2.series(n, 0, 6)
# t2 = expr2.series(n, 0, 2)


# O(expr, (n, S.Infinity))


# print(Sum(expr, (n, 1, oo)).is_convergent())

# expr2 = 1 - cos(1/n)

# print(Sum(expr2, (n, 1, oo)).is_convergent())







