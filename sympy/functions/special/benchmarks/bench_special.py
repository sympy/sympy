from sympy import symbols
from sympy.functions.special.spherical_harmonics import Ylm

x,y = symbols('x,y')

def timeit_Ylm_xy():
    Ylm(1,1, x,y)
