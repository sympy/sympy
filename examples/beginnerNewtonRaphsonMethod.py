# Implementing Newton Raphson method in python
# Author: Haseeb

from sympy import diff
# This modules make calculations more precise
from decimal import Decimal
from math import sin, cos, exp, sqrt

def NewtonRaphson(func, a, ERROR=10**-15):
    ''' Finds root from the point a by Newton-Raphson method '''
    while True:
        x = a
        c = Decimal(a) - ( Decimal(eval(func)) / Decimal(eval(str(diff(func)))) )
        
        x = c
        a = c
        if  abs(eval(func)) < ERROR:
            return  c

# Let's Execute
if __name__ == '__main__':
    # Find root of trignometric fucntion
    # Find value of  pi
    print ('sin(x) = 0', NewtonRaphson('sin(x)', 2))
    
    # Find root of polynomial
    print ('x**2 - 5*x +2 = 0', NewtonRaphson('x**2 - 5*x +2', 0.4))
    
    # Find Square Root of 5
    print ('x**2 - 5 = 0', NewtonRaphson('x**2 - 5', 0.1))

    #  Exponential Roots
    print ('exp(x) - 1 = 0', NewtonRaphson('exp(x) - 1', 0))
    
    # Find Square Root of 25
    # Find nth root of the number by replacing 0.5 with 1/2, 1/3, ...
    print ('exp(x) - 1 = 0', NewtonRaphson('25**Decimal(0.5) - x', 0))
        

    

