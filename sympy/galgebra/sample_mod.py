from sympy.galgebra.ga import Ga

(g3d,a,b,c,d,e) = Ga.build('a b c d e')

A =  (a|(b^c))|(b^c)
B =  (a|(d^e))|(d^e)
print (A+B) - (a|((b^c)+(d^e)))|((b^c)+(d^e))







