def det_gauss(self):
    ref, p = self.gauss()
    det = 1
    for i in xrange(ref.rows):
        det *= ref[i, i]
    if len(p) % 2 == 1:
        det *= -1
    return det

def det_LU(self):
    L, U, p = self.LU_sparse()
    det = 1
    for i in xrange(U.rows):
        det *= U[i, i]
    if len(p) % 2 == 1:
        det *= -1
    return det

def gauss(self):
    "Gaussian elimnation, currently tested only on square matrices"
    A = self.copy()
    row_swaps = []
    for j in xrange(A.rows):
        rlist = nz_col_lower(A, j)
        if A[j, j] == 0:
            if rlist:
                A.row_swap(j, rlist[0])
                row_swaps.append((j, rlist[0]))
                rlist.pop(0)
            else:
                continue
        for i in rlist:
            A.row_add(i, j, - A[i, j] / A[j, j])
    return A, row_swaps

def _upper_triangular_solve(self, rhs):
    X = zeros(rhs.rows, 1)
    for i in reversed(xrange(self.rows)):
        X[i, 0] = (rhs[i, 0] - sum(value * X[j, 0] for j, value in self.mat[i])) / self[i, i]
    return X

def _lower_triangular_solve(self, rhs):
    X = zeros(rhs.rows, 1)
    for i in xrange(self.rows):
        X[i, 0] = (rhs[i, 0] - sum(value * X[j, 0] for j, value in self.mat[i])) / self[i, i]
    return X

def LUsolve(self, rhs):
    L, U, p = LU_sparse(self)
    b = rhs.copy()
    b.permute(p)
    Y = _lower_triangular_solve(L, b)
    return _upper_triangular_solve(U, Y)

def solve_gauss(self, rhs):
    U, p = gauss(self)
    b = rhs.copy()
    b.permute(p)
    return _upper_triangular_solve(U, b)

def solve_rref(self, rhs):
    big = self.join_rows(rhs)
    rref = big.rref()
    return rref[:, self.cols]

def rref(self):
    "rref"
    A = self.copy()
    for j in xrange(A.rows):
        rlist = nz_col_lower(A, j)
        if A[j, j] == 0:
            if rlist:
                A.row_swap(j, rlist[0])
                rlist.pop(0)
            else:
                continue
        A.row_scale(j, 1 /  A[j, j])
        for i in A.nz_col(j):
            if i != j:
                A.row_add(i, j, - A[i, j])
    return A
                
def nz_col(self, j):
    li = []
    for i in xrange(self.rows):
        if self[i, j] != 0:
            li.append(i)
    return li

def nz_col_lower(self, j):
    " Returns the row indices of non-zero elements in column j, below the diagonal"
    li = []
    for i in xrange(j + 1, self.rows):
        if self[i, j] != 0:
            li.append(i)
    return li

def inv_rref(self):
    aug = self.join_rows(eye(self.rows, one = 1, zero = 0))
    reduced = aug.rref()
    return reduced[:,self.rows:]

def LU_sparse(self):
    row_swaps = []
    A = self.copy()
    n = self.rows
    for k in xrange(n):

        rlist = nz_col_lower(A, k)

        # Pivoting
        if A[k, k] == 0:
            if not rlist:
                raise Exception('Singular')
            A.row_swap(k, rlist[0])
            row_swaps.append((k, rlist[0]))
            rlist.pop(0)
        assert A[k, k] != 0

        # Algorithm
        for i in rlist:
            A[i, k] /= A[k, k]

        for j, val in A.mat[k]:
            if j <= k:
                continue
            # j > k
            for i in rlist:
                A[i, j] -= A[i, k] * val # i > k, j > k

    L = eye(self.rows)
    for i in xrange(L.rows):
        for j in xrange(i):
            L[i, j] = A[i, j]

    U = zeros(self.rows)
    for i in xrange(U.rows):
        for j in xrange(i, U.rows):
            U[i, j] = A[i, j]

    return L, U, row_swaps

def join_rows(A, B):
    A = A.copy()
    B = B.copy()
    assert A.rows == B.rows
    for i in xrange(B.rows):
        for ind, (j, val) in enumerate(B.mat[i]):
            B.mat[i][ind] = (j + A.cols, val)
        A.mat[i].extend(B.mat[i])
    A.cols += B.cols
    return A

def row_add(self, r1, r2, alpha):
    if r1 == r2:
        return
    row1 = self.mat[r1]
    row2 = self.mat[r2]
    self.mat[r1] = _row_add(row1, row2, alpha)

def row_scale(self, r, alpha):
    for ind, (j, val) in enumerate(self.mat[r]):
        self.mat[r][ind] = (j, alpha * val)

def _row_add(row1, row2, alpha):
    "li = row1 + alpha * row2 "
    li = []
    i1 = i2 = 0
    n1 = len(row1)
    n2 = len(row2)
    while i1 < n1 or i2 < n2:
        # print i1, i2, len(row1), len(row2)
        if i1 < n1 and (i2 >= n2 or row1[i1][0] < row2[i2][0]):
            li.append(row1[i1])
            i1 += 1
        elif i1 >= n1 or row1[i1][0] > row2[i2][0]:
            li.append((row2[i2][0], alpha * row2[i2][1]))
            i2 += 1
        else:
            val = row1[i1][1] + alpha * row2[i2][1]
            if val != 0:
                li.append((row1[i1][0], val))
            i1 += 1
            i2 += 1
    return li   

def row_swap(self, r1, r2):
    self.mat[r1], self.mat[r2] = self.mat[r2], self.mat[r1]

def zeros(*args):
    from lilmatrix import LILMatrix
    if len(args) == 1:
        m = n = args[0]
    elif len(args) == 2:
        m, n = args
    return LILMatrix(m, n, lambda i, j: 0)

def eye(n):
    matrix = zeros(n)
    for i in xrange(n):
        matrix[i, i] = 1
    return matrix

