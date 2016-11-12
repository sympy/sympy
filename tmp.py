from sympy import *

A = MatrixSymbol("A", 1, 3)
B = MatrixSymbol("B", 1, 3)
C = MatrixSymbol("C", 1, 3)
M = MatrixSymbol("M", 1, 3)    

print( ccode(A[0,0]) )

E = A-B
F = C[0, 0]
F = F.subs(C, E)
print(ccode(F))# == "(-1)*B[0] + A[0]")

E = A - B + M
F = C[0, 0]
F = F.subs(C, E)
print(ccode(F))# == "(-1)*B[0] + A[0] + M[0]")

E = A + M
F = C[0, 1]
F = F.subs(C, E)
print(ccode(F))# == "A[1] + M[1]")
