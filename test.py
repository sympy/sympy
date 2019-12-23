from sympy import *

M = Matrix([[1, 1, -1], [1, 3, 1], [-1, 1, 3]])
M_evs = sorted(list(M.eigenvals().keys()))

M_minors = [M.minor_submatrix(i, i) for i in range(M.rows)]
M_minors_evs = [sorted(list(m.eigenvals().keys())) for m in M_minors]

def entry_LHS(i, j):
    return Abs(Symbol("v_{%s, %s}" % (i, j)))**2

def entry_RHS(i, j):
    numer = S.One
    for k in range(M.rows-1):
        numer *= M_evs[i] - M_minors_evs[j][k]

    denom = S.One
    for k in range(M.rows):
        if k == i:
            continue
        denom *= M_evs[i] - M_evs[k]

    return simplify(numer/denom)

display_mat_LHS = Matrix(M.rows, M.cols, entry_LHS)
display_mat_RHS = Matrix(M.rows, M.cols, entry_RHS)
print('Result computed from minors:')
print(Eq(display_mat_LHS.T, display_mat_RHS.T))

l = list()
for _, _, vects in M.eigenvects():
    for vect in vects:
        vect = vect.normalized()
        vect = vect.applyfunc(lambda x: simplify(x**2))
        l.append(vect)

display_mat_RHS = Matrix.hstack(*l)
print('The actual result:')
print(Eq(display_mat_LHS.T, display_mat_RHS))
