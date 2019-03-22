from sympy import *
init_printing()
from functools import reduce 
from sympy.abc import x
def by_parts_integration(eq, args = []):

    diffr = [ ]
    integrates = [ ] 

    #Improve this
    if len(factor_list(eq)[1])< 2:
        if args:
            return integrate(eq, args)
        else:
            return integrate(eq)

    for i,j in factor_list(eq)[1]:
        if type(i) in [type(sec(x)), type(cos(x)), type(sin(x)), type(tan(x)), 
                       type(csc(x)), type((cot(x))), type(exp(x))]:        
            integrates.append((i, j))

        else:    
#         if type(i) in ["<class 'sympy.core.power.Pow'>", 
#             "<class 'sympy.core.add.Add'>", "<class 'sympy.core.symbol.Symbol'>"]:
            diffr.append((i, j))
  ######
    funcs = lambda l:[i**j for i,j in l]        
    diffr = funcs(diffr)
    integrates = funcs(integrates)
    j = 1        
    for i in diffr:
        j *= i
    diffr = j    
    
    j = 1        
    for i in integrates:
        j *= i
    integrates = j    
    if args:
        temp = by_parts_integration(integrates, args)
        return factor_list(eq)[0]*(  diffr * temp - by_parts_integration(diffr.diff(args[0]) * temp, args)  )
    else:
        temp = by_parts_integration(integrates)
        return factor_list(eq)[0]*(  diffr * temp - by_parts_integration(diffr.diff(list(eq.atoms(Symbol))[0]) * simplify(temp) ))
