from sympy import SparseMatrix
from collections import defaultdict


def _doktocsr(dok):
    """Converts a sparse matrix to Compressed Sparse Row (CSR) format.

    Parameters
    ==========

    A : contains non-zero elements sorted by key (row, column)
    JA : JA[i] is the column corresponding to A[i]
    IA : IA[i] contains the index in A for the first non-zero element
        of row[i]. Thus IA[i+1] - IA[i] gives number of non-zero
        elements row[i]. The length of IA is always 1 more than the
        number of rows in the matrix.
    """
    row, JA, A = [list(i) for i in zip(*dok.row_list())]
    IA = [0]*((row[0] if row else 0) + 1)
    for i, r in enumerate(row):
        IA.extend([i]*(r - row[i - 1]))  # if i = 0 nothing is extended
    IA.extend([len(A)]*(dok.rows - len(IA) + 1))
    shape = [dok.rows, dok.cols]
    return [A, JA, IA, shape]


def _csrtodok(csr):
    """Converts a CSR representation to DOK representation"""
    smat = {}
    A, JA, IA, shape = csr
    for i in range(len(IA) - 1):
        indices = slice(IA[i], IA[i + 1])
        for l, m in zip(A[indices], JA[indices]):
            smat[i, m] = l
    return SparseMatrix(*(shape + [smat]))


def _add(dok1, dok2, K):
    if dok1.nnz() > dok2.nnz():
        dok2, dok1 = dok1, dok2
    dict1 = dok1._smat
    dict2 = dok2._smat
    residue = dict2
    smat = {}
    for k in dict1:
        smat[k] = K(dok1[k]) + K(dok2[k])
    for k in dict2:
        if k not in dict1:
            smat[k] = K(dict2[k])
    return SparseMatrix(dok1.rows, dok1.cols, smat)


def _sub(dok1, dok2, K):
    return _add(dok1, -dok2, K)


def _isIdentity(csr, K):
    a, ja, ia, shape = csr
    for i in a:
        if i == K(1):
            return True


def _symmmetric(sp, K):
    transpose = _transpose(sp)
    if transpose == sp:
        return transpose == sp
    else:
        return False


def _mulscsp(v, csr, K):
    a, ja, ia, shape = csr
    return [[K(v)*K(i) for i in a], ja, ia, shape]


def _mulspsp(dok, other, K):
    A = dok
    B = other
    # sort B's row_list into list of rows
    Blist = [[] for i in range(B.rows)]
    for i, j, v in B.row_list():
        Blist[i].append((j, K(v))) #domain added here.
    Cdict = defaultdict(int)
    for k, j, Akj in A.row_list():
        for n, Bjn in Blist[j]:
            temp = K(Akj)*K(Bjn) #domain added here
            Cdict[k, n] += temp
    rv = dok.zeros(A.rows, B.cols)
    rv._smat = dict([(k, v) for k, v in Cdict.iteritems() if v])
    return rv


def _mulspvec(csr, vec, K):
    a, ja, ia, shape = csr
    smat = {}
    for i in range(len(ia) - 1):
        stripe = slice(ia[i], ia[i + 1])
        for m, n in zip(a[stripe], ja[stripe]):
            try:
                smat[i, 0] = K(smat[i, 0]) +  K(m)*K(vec[n]) #domain added here
            except KeyError:
                smat[i, 0] = K.zero +  K(m)*K(vec[n]) #domain added here
                return SparseMatrix(shape[0], 1, smat)


def _applyfunc(csr, f, K):
    if not callable(f):
        return TypeError("f must be callable")
    else:
        a, ja, ia, shape = csr
        return [[f(i) for i in a], ja, ia, shape]


def _trace(dok, K):
    if dok.rows != dok.cols:
        raise ShapeError()
    smat = K.zero
    for k in dok._smat.keys():
        if k[0] == k[1]:
            smat += K(dok[k]) #domain added here
    return smat


def _transpose(dok, K):
    smat = {}
    for k in dok._smat.keys():
        smat[k[1], k[0]] = dok._smat[k]
    return SparseMatrix(dok.cols, dok.rows, smat)


def liupc(dok, K):
    R = [[] for r in range(dok.rows)]
    for r, c, _ in dok.row_list():
        if c <= r:
            R[r].append(c)
        inf = len(R)  # nothing will be this large
        parent = [inf]*dok.rows
        virtual = [inf]*dok.rows
        for r in range(dok.rows):
            for c in R[r][:-1]:
                while virtual[c] < r:
                    t = virtual[c]
                    virtual[c] = r
                    c = t
                if virtual[c] == inf:
                    parent[c] = virtual[c] = r
        return R, parent


def _cholesky_sparse(dok, K):
    Crowstruc = dok.row_structure_symbolic_cholesky()
    C = dok.zeros(dok.rows)
    for i in range(len(Crowstruc)):
        for j in Crowstruc[i]:
            if i != j:
                C[i, j] = K(sparse[i, j])
                summ = K.zero
                for p1 in Crowstruc[i]:
                    if p1 < j:
                        for p2 in Crowstruc[j]:
                            if p2 < j:
                                if p1 == p2:
                                    summ += C[i, p1]*C[j, p1]
                                else:
                                    break
                            else:
                                break
                    C[i, j] -= summ
                    C[i, j] /= C[j, j]
                else:
                    C[j, j] = K(dok[j, j])
                    summ = 0
                    for k in Crowstruc[j]:
                        if k < j:
                            summ += C[j, k]**2
                        else:
                            break
                    C[j, j] -= summ
                    C[j, j] = sqrt(C[j, j])
        return C

def _LDL_sparse(dok):
    Lrowstruc = dok.row_structure_symbolic_cholesky()
    L = dok.eye(dok.rows)
    D = dok.zeros(dok.rows, dok.cols)

    for i in range(len(Lrowstruc)):
        for j in Lrowstruc[i]:
            if i != j:
                L[i, j] = K(dok[i, j])
                summ = K.zero
                for p1 in Lrowstruc[i]:
                    if p1 < j:
                        for p2 in Lrowstruc[j]:
                            if p2 < j:
                                if p1 == p2:
                                    summ += K(L[i, p1]*L[j, p1]*D[p1, p1])#need a check.
                                else:
                                    break
                        else:
                            break
                    L[i, j] -= summ
                    L[i, j] /= D[j, j]
            elif i == j:
                D[i, i] = dok[i, i]
                summ = 0
                for k in Lrowstruc[i]:
                    if k < i:
                        summ += L[i, k]**2*D[k, k]
                    else:
                        break
                    D[i, i] -= summ

        return L, D

def _lower_triangular_solve(dok, rhs, K):
    rows = [[] for i in range(dok.rows)]
    for i, j, v in dok.row_list():
        if i > j:
            rows[i].append((j, K(v))) #need a check
        X = rhs.copy()
        for i in range(dok.rows):
            for j, v in rows[i]: # same as above note in this function
                X[i, 0] -= K(v)*X[j, 0] #added domains here, needs
            X[i, 0] /= dok[i, i]
        return dok._new(X)


def _upper_triangular_solve(dok, rhs, K):
    rows = [[] for i in range(dok.rows)]
    for i, j, v in dok.row_list():
        if i < j:
            rows[i].append((j, v))
        X = rhs.copy()
        for i in range(dok.rows - 1, -1, -1):
            rows[i].reverse()
            for j, v in rows[i]:
                X[i, 0] -= K(v)*K(X[j, 0])
            X[i, 0] /= K(dok[i, i])
        return dok._new(X)


def _binsearch(i, v, beg, end):
    if beg > end:
        return False
    else:
        mid = (beg + end)/2
        if v[mid] == i:
            return mid
        elif v[mid] > i:
            return _binsearch(i, v, beg, mid - 1)
        elif v[mid] < i:
            return _binsearch(i, v, mid + 1, end)
