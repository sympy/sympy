from random import random
from sympy.functions.elementary.trigonometric import atan, cos


def timeit_atan_evaluate_float():
    atan(random())

def timeit_cos_evaluate_float():
    cos(random())
