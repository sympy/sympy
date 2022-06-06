from sympy.core.symbol import symbols
from sympy.functions.special.spherical_harmonics import Ynm
from sympy.functions.special.beta_functions import beta

x, y = symbols('x,y')


def timeit_Ynm_xy():
    Ynm(1, 1, x, y)

def timeit_beta_evalulate_float():
    beta(3.1)
    beta(3.1, 2.1)
