from sympy.combinatorics.schur_number import schur_partition
import random

"""
Checks if every subset is sum-free(There are no x,y,z in the subset such that
x + y = z)
"""
def _sum_free_test(subset):
    for i in subset:
        for j in subset:
            if i + j in subset:
                raise AssertionError("This subset is not sum free")
    return "Everything is correct"


n = random.randint(1,1000)
result = schur_partition(n)

"""
Checks if the occurance of all numbers  is exactly one
"""
t = 0
numbers = []
for item in result:
    _sum_free_test(item)
    t += len(item)
    for l in item:
        if l in numbers:
            raise AssertionError("Every number must exist once")
        numbers.append(l)
assert n == t
