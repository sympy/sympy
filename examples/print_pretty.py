#!/usr/bin/env python
import iam_sympy_example

from sympy import *

x=Symbol("x") 
y=Symbol("y") 


# try to enable unicode pprint
pprint_try_use_unicode()

pprint( x**x )
print '\n'# separate with two blank likes

pprint(x**2+y+x)
print '\n'

pprint(sin(x)**x)
print '\n'

pprint( sin(x)**cos(x) )
print '\n'

pprint( sin(x)/(cos(x)**2 * x**x +(2*y)) )
print '\n'

pprint( sin(x**2+exp(x)) )
print '\n'

pprint( sqrt(exp(x)) )
print '\n'

pprint( sqrt(sqrt(exp(x))) )
print '\n'

pprint( (1/cos(x)).series(x,10) )
print '\n'

