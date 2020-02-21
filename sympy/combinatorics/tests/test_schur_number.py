from sympy.combinatorics.schur_number import schur_partition
from sympy.testing.randtest import _randint
import random

def _sum_free_test(subset):
    """
    Checks if subset is sum-free(There are no x,y,z in the subset such that
    x + y = z)
    """
    for i in subset:
        for j in subset:
            assert (i + j in subset) is False

def test_schur_number():
    random_number_generator = _randint(1000)
    for _ in range(5):
        n = random_number_generator(1, 1000)
        result = schur_partition(n)
        t = 0
        numbers = []
        for item in result:
            _sum_free_test(item)
            """
            Checks if the occurance of all numbers  is exactly one
            """
            t += len(item)
            for l in item:
                assert (l in numbers) is False
                numbers.append(l)
        assert n == t
