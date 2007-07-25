import sys
sys.path.append("..")
from sympy import sqrt, Symbol, eye

w, x, y, z = Symbol('w'), Symbol('x'), Symbol('y'), Symbol('z')
L = [x,y,z]
V = eye(len(L))
for i in range(len(L)):
    for j in range(len(L)):
        V[i,j] = L[i]**j
det = 1
for i in range(len(L)):
    det *= L[i]-L[i-1]

print "matrix"
print V
print "det:"
print V.det().expand()
print "correct result"
print det
print det.expand()
