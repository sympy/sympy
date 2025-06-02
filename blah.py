import sympy
from sympy import *
from sympy.assumptions.satask import satask

#     ``Q.negative``, ``Q.zero``, ``Q.positive``,
#     ``Q.nonnegative``, ``Q.nonpositive``, ``Q.nonzero``,
#     ``Q.integer``, ``Q.rational``, and ``Q.irrational`` all imply
#     ``Q.real``, as do all facts that imply those facts.



# print("buggy - should return 'None'")
# x,y,z = symbols("x y z")
# # y COULD be infinity
# print(ask(Q.zero(x*y*z), Q.zero(x*z)))
# print(ask(Q.zero(x*y*z), Q.zero(x*z)))

# assumptions during debugging
'''

Q.zero(x*z) -> x and z are not infinity

'''

x,y = symbols('x y')

# fixes inconsistent assumptions AND other thing!


# print(satask(Q.finite(x*y), ~Q.finite(x) & Q.zero(y)))
print("Tests")
print(ask(Q.finite(x*y), Q.infinite(x)))
print(ask(Q.finite(x*y), ~Q.finite(x)))
# should return true
print(satask(Q.zero(x*y),   ((Q.zero(x)) & Q.finite(y)))  )
print(satask(Q.zero(x*y),   (((Q.zero(x) | Q.zero(y)) & Q.finite(x)) & Q.finite(y)))    )

# all false
print(satask(Q.zero(x) | Q.zero(y), Q.nonzero(x*y)))
print(satask((Q.zero(x) | Q.zero(y)), Q.nonzero(x*y)))
print(satask(Q.zero(x) | Q.zero(y), Q.nonzero(x*y)))

# print(satask(Q.zero(x) | Q.zero(y), Q.nonzero(x*y)))
# print(satask(Q.zero(x) | Q.zero(y), Q.finite(x*y)))
# print(satask(Q.zero(x) | Q.zero(y), Q.positive(x*y)))
# print(satask(Q.finite(x), Q.positive(x*y)))




# print(satask(Q.zero(x) | Q.zero(y), Q.nonzero(x*y) & Q.finite(x) & Q.finite(y)))
# print(satask(Q.zero(x*y), Q.zero(x) | Q.zero(y)))
# print(satask(Q.zero(x*y), (Q.finite(x)) & (Q.finite(y)) & Q.zero(x)))

# print(satask(Q.finite(x*y), ~Q.finite(x) & Q.zero(y)))
# print(satask(~Q.finite(x) & Q.zero(y))) # gives None as expected 
# print(satask(~Q.finite(x) & Q.zero(y) & Q.finite(x*y))) # wrongly gives false
# print(satask(~Q.finite(x) & Q.zero(y) & ~Q.finite(x*y))) # also wrongly gives false

# print(satask(~Q.finite(x) & Q.zero(y))) # gives None as expected 
# print(satask(~Q.finite(x) & Q.zero(y) & Q.finite(x*y))) # wrongly gives false
# print(satask(~Q.finite(x) & Q.zero(y) & ~Q.finite(x*y)))

# print(x.assumptions0)
# print(y.assumptions0)
print("\n")


# print("normal - returns 'True'")
# x,y,z = symbols("x y z", finite = true)
# print(ask(Q.zero(x*y*z), Q.zero(x*z)))
# print(x.assumptions0)
# print(y.assumptions0)
# print(z.assumptions0)
# print("\n")


# You might be wondering why this behavior is wrong. 
# It's because variables aren't assumed to be finite by default and if y=oo, then xy=0*oo = NaN.