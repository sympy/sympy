from sympy.core.singleton import S

def _slice_to_bounds(key, defmax): 
    """
        Takes slice or number and returns (min,max) for iteration
        Takes a default maxval to deal with the slice ':' which is (none, none)
    """
    lo, hi = 0, defmax
    if key.start is not None:
        if key.start >= 0:
            lo = key.start
        else:
            lo = defmax+key.start
    if key.stop is not None:
        if key.stop >= 0:
            hi = key.stop
        else:
            hi = defmax+key.stop
    return lo, hi

def _denserepr_from_dict(rows, cols, di):
    mat = [0] * (rows * cols)
    for i in xrange(rows):
        for j in xrange(cols):
            mat[i * cols + j] = di[i, j]
    return mat

def _denserepr_from_list(rows, cols, li):
    mat = [0] * (rows * cols)
    for i, val in enumerate(li):
        mat[i] = val
    return mat

def _denserepr_from_lil(rows, cols, lil):
    mat = []
    for row in lil:
        mat.extend(row)
    return mat

def _denserepr_from_callable(rows, cols, f):
    mat = [0] * (rows * cols)
    for i in xrange(rows):
        for j in xrange(cols):
            mat[i * cols + j] = f(i, j)
    return mat

def _lilrepr_from_callable(rows, cols, f):
    assert type(rows) == int, rows
    mat = [[] for i in xrange(rows)]
    for i in xrange(rows):
        for j in xrange(cols):
            val = S(f(i,j))
            if val != 0:
                mat[i].append((j, val))  
    return mat

def _lilrepr_from_list(rows, cols, li):
    mat = [[] for i in xrange(rows)]
    for i in xrange(rows):
        for j in xrange(cols):
            val = S(li[i*cols + j])
            if val != 0:
                mat[i].append((j, val))
    return mat

def _lilrepr_from_dict(rows, cols, di):
    mat = [[] for i in xrange(rows)]
    for i in xrange(rows):
        for j in xrange(cols):
            val = S(di[i, j])
            if val != 0:
                mat[i].append((j, val))
    return mat

def _lilrepr_from_lil(rows, cols, lil):
    mat = [[] for i in xrange(rows)]
    for i in xrange(rows):
        for j in xrange(cols):
            val = S(lil[i][j])
            if val != 0:
                mat[i].append((j, val))
    return mat

def _dokrepr_from_lil(rows, cols, lil):
    mat = {}
    for i in xrange(rows):
        for j in xrange(cols):
            val = lil[i][j]
            if val != 0:
                mat[i, j] = val
    return mat

def _dokrepr_from_dict(rows, cols, di):
    mat = {}
    for i in xrange(rows):
        for j in xrange(cols):
            if (i, j) in di:
                mat[i, j] = di[i, j]
    return mat

def _dokrepr_from_callable(rows, cols, func):
    mat = {}
    for i in xrange(rows):
        for j in xrange(cols):
            val = func(i, j)
            if val != 0:
                mat[i, j] = val
    return mat

def _from_callable(*args):
    rows = args[0]
    cols = args[1]
    op = args[2]
    mat = {}
    for i in xrange(rows):
        for j in xrange(cols):
            val = op(i, j)
            if val != 0:
                mat[i, j] = val
    return rows, cols, mat

def _from_list(*args):
    rows = args[0]
    cols = args[1]
    inp = args[2]
    mat = {}
    for i in xrange(rows):
        for j in xrange(cols):
            val = mat[i * rows + j]
            if val != 0:
                mat[i, j] = val
    return rows, cols, mat

def _dict_to_densematrix(rows, cols, mat):
    return DenseMatrix(rows, cols, mat)

def _dict_to_dokmatrix(rows, cols, mat):
    return DOKMatrix(rows, cols, mat)

def _dict_to_lilmatrix(rows, cols, mat):
    return LILMatrix(rows, cols, mat)

