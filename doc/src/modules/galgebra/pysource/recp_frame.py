from sympy import expand, simplify
from ga import Ga
from printer import Format, xpdf

Format()

g = '1 # #,'+ \
    '# 1 #,'+ \
    '# # 1'

ng3d = Ga('e1 e2 e3',g=g)
(e1,e2,e3) = ng3d.mv()

print 'g_{ij} =',ng3d.g

E = e1^e2^e3
Esq = (E*E).scalar()
print 'E =',E
print '%E^{2} =',Esq
Esq_inv = 1/Esq

E1 = (e2^e3)*E
E2 = (-1)*(e1^e3)*E
E3 = (e1^e2)*E

print 'E1 = (e2^e3)*E =',E1
print 'E2 =-(e1^e3)*E =',E2
print 'E3 = (e1^e2)*E =',E3

w = (E1|e2)
w = w.expand()
print 'E1|e2 =',w

w = (E1|e3)
w = w.expand()
print 'E1|e3 =',w

w = (E2|e1)
w = w.expand()
print 'E2|e1 =',w

w = (E2|e3)
w = w.expand()
print 'E2|e3 =',w

w = (E3|e1)
w = w.expand()
print 'E3|e1 =',w

w = (E3|e2)
w = w.expand()
print 'E3|e2 =',w

w = (E1|e1)
w = (w.expand()).scalar()
Esq = expand(Esq)
print '%(E1\\cdot e1)/E^{2} =',simplify(w/Esq)

w = (E2|e2)
w = (w.expand()).scalar()
print '%(E2\\cdot e2)/E^{2} =',simplify(w/Esq)

w = (E3|e3)
w = (w.expand()).scalar()
print '%(E3\\cdot e3)/E^{2} =',simplify(w/Esq)

xpdf(paper='letter')
