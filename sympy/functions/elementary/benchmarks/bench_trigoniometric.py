from random import random
from sympy.functions.elementary.trigonometric import atan, cos
from sympy.functions.elementary.miscellaneous import sqrt

sqrt3 = sqrt(3)

def timeit_atan_evaluate_sqrt_3():
    atan(sqrt3)


def timeit_atan_evaluate_float():
    atan(random())


def timeit_cos_evaluate_float():
    cos(random())
