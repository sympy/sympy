def rowdecomp(self, num):
    nmax = len(self)
    if not (0 <= num < nmax) and not (0 <= -num < nmax):
        raise ValueError("`num` must satisfy 0 <= `num` < `self.rows*" +
            "*self.cols` (%d) and 0 <= -num < " % nmax +
            "`self.rows*self.cols` (%d) to apply redecomp()." % nmax)
    i, j = 0, num
    while j >= self.cols:
        j -= self.cols
        i += 1
    return i,j

def row_del(self, k):
    newD = {}
    for (i,j) in self.mat.keys():
        if i==k:
            pass
        elif i > k:
            newD[i-1,j] = self.mat[i,j]
        else:
            newD[i,j] = self.mat[i,j]
    self.mat = newD
    self.rows -= 1

def col_del(self, k):
    newD = {}
    for (i,j) in self.mat.keys():
        if j==k:
            pass
        elif j > k:
            newD[i,j-1] = self.mat[i,j]
        else:
            newD[i,j] = self.mat[i,j]
    self.mat = newD
    self.cols -= 1

def _nonzero_elements_in_row(self, m):
    keys = sorted(self.mat.keys())
    start = 0
    end = len(keys)
    for i in xrange(len(keys)):
        if keys[i][0] == m - 1:
            start = i + 1
        if keys[i][0] == m + 1:
            end = i
            return keys[start:end]
    return keys[start:end]

def _nonzero_elements_in_col(self, n):
    keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())]
    start = 0
    end = len(keys)
    for i in xrange(len(keys)):
        if keys[i][1] == n - 1:
            start = i + 1
        if keys[i][1] == n + 1:
            end = i
            return keys[start:end]
    return keys[start:end]

def _lil_row_major(self):
    n = self.rows
    keys = sorted(self.mat.keys())
    lil = [[]] * n
    if not keys:
        return lil
    k = 0
    start = 0
    for i in xrange(len(keys)):
        if keys[i][0] > k:
            lil[k] = keys[start:i]
            start = i
            k = keys[i][0]
    lil[keys[-1][0]] = keys[start:]
    return lil

def _lil_col_major(self):
    n = self.cols
    keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())]
    lil = [[]] * n
    if not keys:
        return lil
    k = 0
    start = 0
    for i in xrange(len(keys)):
        if keys[i][1] > k:
            lil[k] = keys[start:i]
            start = i
            k = keys[i][0]
    lil[keys[-1][1]] = keys[start:]
    return lil
        
def elem_ops_row(self, row1, row2, alpha):
    ''' row1 = row1 + alpha * row2 '''

    self._lil_row_major()
    rowlist1 = lil[row1]
    rowlist2 = lil[row2]
    for key in rowlist2:
        key1 = (row1, key[1])
        if key1 not in rowlist1:
            self.mat[key1] = alpha * self.mat[key]
        else:
            self.mat[key1] += alpha * self.mat[key]

def elem_ops_col(self, col1, col2, alpha):
    "col1 = col1 + alpha * col2"

    self._lil_col_major()
    collist1 = lil[col1]
    collist2 = lil[col2]
    for key in collist2:
        key1 = (key[0], col1)
        if key1 not in collist1:
            self.mat[key1] = alpha * self.mat[key]
        else:
            self.mat[key1] += alpha * self.mat[key]

def liu(self):
    R = self._lower_row_nonzero_structure()
    parent = [None] * self.rows
    for j in xrange(self.rows):
        parent[j] = None
        for i in xrange(len(R[j])):
            r = R[j][i][1]
            while parent[r] != None:
                r = parent[r]
            if r != j:
                parent[r] = j
    return parent

def liupc(self):
    R = self._lower_row_nonzero_structure()
    parent = [None] * self.rows
    virtual = [None] * self.rows
    for j in xrange(self.rows):
        for i, (_,r) in enumerate(R[j][:-1]):
            while virtual[r] and virtual[r] < j:
                t = virtual[r]
                virtual[r] = j
                r = t
            if not virtual[r]:
                virtual[r] = j
                parent[r] = j
    return R, parent

def row_structure_symbolic_cholesky(self):
    from copy import deepcopy
    R, parent = self.liupc()
    Lrow = deepcopy(R)
    for k in xrange(self.rows):
        for _, j in R[k]:
            while j != None and j != k:
                if (k, j) not in Lrow[k]:
                    Lrow[k].append((k, j))
                    Lrow[k].sort()
                j = parent[j]
    return Lrow

def elementary_symbolic_cholesky(self):
    C = self._lower_columnar_nonzero_structure()
    for col in C:
        for iup, up in enumerate(col[1:]):
            for down in col[iup+1:]:
                if (down[0], up[0]) not in C[up[0]]:
                    C[up[0]].append((down[0], up[0]))
                    C[up[0]].sort()
    return C

def fast_symbolic_cholesky(self): 
    """ implement algorithm 1.3 from 10.1.1.163.7506 """
    C = self._lower_columnar_nonzero_structure()
    p=[]
    for k, col in enumerate(C):
        if len(C[k]) > 1:
            parent = C[k][1][0]
            p.append(parent)
        else:
            continue
        for i, _ in col[1:]: # add (parent[j], j) to C
            if (i, parent) not in C[parent]:
                C[parent].append((i, parent))
                C[parent].sort()
    return C

def _cholesky(self):

    cached_result = self._get_cache('CH')
    if cached_result and self._cached:
        return cached_result

    CHstruc = self.fast_symbolic_cholesky()
    C = DOKMatrix(self.rows, self.cols, {})
    for col in CHstruc:
        for i, j in col:
            if i != j:
                C[i, j] = (1 / C[j, j]) * (self[i, j] - sum(C[i, k] * C[j, k] 
                    for k in xrange(j)))
            elif i==j:
                C[j, j] = (self[j, j] - sum(C[j, k] ** 2
                    for k in xrange(j))) ** (0.5)

    self._set_cache('CH', C)
    return C

def _cholesky_sparse(self):

    cached_result = self._get_cache('CH')
    if cached_result and self._cached:
        return cached_result

    Crowstruc = self.row_structure_symbolic_cholesky()
    C = DOKMatrix.zeros(self.rows)
    for row in Crowstruc:
        for i, j in row:
            if i != j:
                C[i, j] = self[i, j]
                summ = 0
                for p1 in Crowstruc[i]:
                    if p1[1] < j:
                        for p2 in Crowstruc[j]:
                            if p2[1] < j:
                                if p1[1] == p2[1]:
                                    summ += C[p1] * C[p2]
                            else:
                                break
                        else:
                            break
                C[i, j] -= summ
                C[i, j] /= C[j, j]
            else:
                C[j, j] = self[j, j]
                summ = 0
                for _, k in Crowstruc[j]:
                    if k < j:
                        summ += C[j, k] ** 2
                    else:
                        break
                C[j, j] -= summ
                C[j, j] = C[j, j] ** 0.5

    self._set_cache('CH', C)
    return C
                
        
def _LDL(self):

    cached_result = self._get_cache('LDL')
    if cached_result and self._cached:
        return cached_result

    CHstruc = self.fast_symbolic_cholesky()
    L = DOKMatrix.eye(self.rows)
    D = DOKMatrix(self.rows, self.cols, {})
    
    for col in CHstruc:
        for i, j in col:
            if i != j:
                L[i, j] = (self[i, j] - sum(L[i, k] * L[j, k] * D[k, k] 
                    for k in range(j))) / D[j, j]
            elif i == j:
                D[i, i] = self[i, i] - sum(L[i, k]**2 * D[k, k] 
                    for k in range(i))

    self._set_cache('LDL', (L, D))
    return L, D


def _LDL_sparse(self):

    cached_result = self._get_cache('LDL')
    if cached_result and self._cached:
        return cached_result

    Lrowstruc = self.row_structure_symbolic_cholesky()
    L = DOKMatrix.eye(self.rows)
    D = DOKMatrix(self.rows, self.cols, {})
    
    for row in Lrowstruc:
        for i, j in row:
            if i != j:
                L[i, j] = self[i, j]
                summ = 0
                for p1 in Lrowstruc[i]:
                    if p1[1] < j:
                        for p2 in Lrowstruc[j]:
                            if p2[1] < j:
                                if p1[1] == p2[1]:
                                    summ += L[p1] * L[p2] * D[p1[1], p1[1]]
                            else:
                                break
                    else:
                        break
                L[i, j] -= summ
                L[i, j] /= D[j, j]
            elif i == j:
                D[i, i] = self[i, i]
                summ = 0
                for _, k in Lrowstruc[i]:
                    if k < i:
                        summ += L[i, k]**2 * D[k, k]
                    else:
                        break
                D[i, i] -= summ

    self._set_cache('LDL', (L, D)) 
    return L, D


def _lower_triangular_solve(self, rhs):
    rows = self._lil_row_major()
    X = DOKMatrix(rhs.rows, 1, rhs.mat)
    for i in xrange(self.rows):
        for key in rows[i]:
            if key[1] < i:
                X[i, 0] -= self[key] * X[key[1], 0]
            else:
                break
        X[i, 0] /= self[i, i]
    return X

def _upper_triangular_solve(self, rhs):
    """ Helper function of function upper_triangular_solve, without the error checks, to be used privately. """
    rows = self._lil_row_major()
    X = DOKMatrix(rhs.rows, 1, rhs.mat)
    for i in reversed(xrange(self.rows)):
        for key in reversed(rows[i]):
            if key[1] > i:
                X[i, 0] -= self[key] * X[key[1], 0]
            else:
                break
        X[i, 0] /= self[i, i]
    return X

def _diagonal_solve(self, rhs):
    return DOKMatrix(self.rows, 1, lambda i, j: rhs[i, 0] / self[i, i])

def _cholesky_solve(self, rhs):
    L = self._cholesky()
    Y = L._lower_triangular_solve(rhs)
    X = L.T._upper_triangular_solve(Y)
    return X

def LDL_solve(self, rhs):
    if not self.is_square():
        raise Exception("Make a square matrix exception")
    if not self.rows == rhs.rows:
        raise Exception
    if not self.is_symmetric():
        rhs = self.T * rhs
        self = self.T * self
    L, D = self._LDL_sparse()
    z = L._lower_triangular_solve(rhs)
    y = D._diagonal_solve(z)
    x = L.T._upper_triangular_solve(y)
    if x.has(nan) or x.has(oo): # import this
        raise Exception
    return x

def _LDLsolve(self, rhs):
    L, D = self._LDL_sparse()
    z = L._lower_triangular_solve(rhs)
    y = D._diagonal_solve(z)
    x = L.T._upper_triangular_solve(y)
    return x

def _lower_columnar_nonzero_structure(self):
    n = self.cols
    NZlist = [[]] * n
    keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())]
    k = 0
    startset = False
    for i in xrange(len(keys)):
        if startset and keys[i][1] > k:
            NZlist[k] = keys[start:i]
            startset = False
        k = keys[i][1]
        if not startset and keys[i][1] == k and keys[i][1] <= keys[i][0]:
            start = i
            startset = True
    if (n-1, n-1) == keys[-1]:
        NZlist[-1].append((n-1,n-1))
    return NZlist

def _lower_row_nonzero_structure(self):
    n = self.rows
    NZlist = [[]] * n
    keys = sorted(self.mat.keys())
    k = 0
    startset = False
    for i in xrange(len(keys)):
        if startset and (keys[i][0] > k or keys[i][1] > keys[i][0]):
            NZlist[k] = keys[start:i]
            startset = False
        k = keys[i][0]
        if not startset and keys[i][1] <= keys[i][0] and keys[i][0] == k and not NZlist[k]:
            start = i
            startset = True
    NZlist[keys[start][0]] = keys[start:]
    return NZlist

def reshape(self, _rows, _cols):
    if len(self) != _rows*_cols:
        print "Invalid reshape parameters %d %d" % (_rows, _cols)
    newD = {}
    for i in range(_rows):
        for j in range(_cols):
            m,n = self.rowdecomp(i*_cols + j)
            if (m,n) in self.mat:
                newD[(i,j)] = self.mat[(m,n)]
    return DOKMatrix(_rows, _cols, newD)

def copyin_list(self, key, value):
    self.copyin_matrix(key, DOKMatrix(value))

def inverse_solver(self, solver=None):
    from matrixutils import vecs2matrix
    if not solver or solver == 'CH':
        solver = self._cholesky_solve
    elif solver == 'LDL':
        solver = self.LDL_solve
    else:
        raise Exception('solver method not recognized')
    I = Matrix.eye(self.rows)
    return vecs2matrix([solver(I[:, i])
        for i in xrange(self.cols)])

