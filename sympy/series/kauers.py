from sympy.abc import x
from sympy import expand

def Delta(expr, d = 1):
     
    """Takes as input function expression and returns the diffrence between final value(function's expression  	      incremented to 1) and initial value (function's expression).

       >>> from sympy import Symbol
       >>> from sympy import sin, cos
       >>> from sympy.series.integrals import Delta
       >>> x = Symbol('x')
       >>> Delta(x**2 + 3*x -2)
       2*x + 4
       >>> Delta(sin(x))
       - sin(x) + sin(x + 1)			 	  			
    
       If you want an increment different than 1, you can give it the increment value as an argument to Delta
       as shown below
       >>>Delta(x**2 - 2*x + 3, 2)
       4*x
       >>>Delta(x**3 + 3*x**2 + 4*x + 5, 3)
       W9*x**2 + 45*x + 66	 
    """ 
    expr = expand(expr)
    a = expr.subs(x, x + d)
    return expand(a) - expr
      
  
    
