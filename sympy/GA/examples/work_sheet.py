#!/usr/bin/python
import sys

from sympy import symbols,sin,cos,Rational
from GA import *
from GAPrint import *

enhance_print()

(a,b,Qa,Qb) = MV.setup('a b Qa Qb')

Ep = (a^b)+(Qa^Qb)
Em = (a^b)-(Qa^Qb)
F  = (a^Qb)+(Qa^b)

print 'Ep x F =',Com(Ep,F)
print 'Em x F =',Com(Em,F)
