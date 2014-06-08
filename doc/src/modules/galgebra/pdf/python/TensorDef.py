
import sys
from sympy import symbols,sin,cos
from printer import Format,xpdf,Get_Program,Print_Function
from ga import Ga
from lt import Mlt

coords = symbols('t x y z',real=True)
(st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)

A = st4d.mv('T','bivector')

def TA(a1,a2):
    global A
    return A | (a1 ^ a2)

T = Mlt(TA,st4d) # Define multi-linear function




