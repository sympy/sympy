# -*- coding: utf-8 -*-

from sympy.matrices import *
from sympy.tensor.arraypy import *
from sympy import Add, diff, count, symbols, simplify


# ---------------- scal_prod g(X,Y)--------------------------------

def scal_prod(X, Y, g):
    """ Returns scalar product of vectors
        g(X,Y)=sum_{i,j}([g[i,j]*X[i]*Y[j])
        Example:
        ========
        >>> x1, x2 = symbols('x1, x2')
        >>> X = [1,2]
        >>> Y = [4,5]
        >>> A = Arraypy((2,2))
        >>> g = Tensor(A,(-1,-1))
        >>> g[0,0] = cos(x2)**2
        >>> g[0,1] = 0
        >>> g[1,0] = 0
        >>> g[1,1] = 1
        >>> print g
        >>> skal = scal_prod(X, Y, g)
    """

    # Handling of a input argument - metric tensor g
    if not isinstance(g, (Matrix, Tensor, Arraypy)):
        raise TypeError(
            'The type of metric tensor must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(g, (Arraypy, Tensor)):
            if isinstance(g, Tensor):
                if not g.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of metric tensor must be (-1,-1)')
            g = g.to_matrix()
    if not g.is_symmetric():
        raise ValueError('The metric tensor must be symmetric.')

    # Handling of a input arguments - vector or vector fields X
    if not isinstance(X, (list, Arraypy, Tensor)):
        raise TypeError('The type of vector must be list, Arraypy or Tensor')
    if isinstance(X, (Tensor, Arraypy)):
        if len(X.shape) != 1:
            raise ValueError("The dimension of X must be 1!")
        if isinstance(X, Tensor):
            if not X.type_pq == (1, 0):
                raise ValueError('The valence or ind_char of X must be (+1)')
    if isinstance(X, (Tensor, Arraypy)):
        X = X.To_list()

    # Handling of a input arguments - vector or vector fields Y
    if not isinstance(Y, (list, Arraypy, Tensor)):
        raise TypeError('The type of vector must be list, Arraypy or Tensor')
    if isinstance(Y, (Tensor, Arraypy)):
        if len(Y.shape) != 1:
            raise ValueError("The dimension of Y must be 1!")
        if isinstance(Y, Tensor):
            if not Y.type_pq == (1, 0):
                raise ValueError('The valence or ind_char of Y must be (+1)')
    if isinstance(Y, (Tensor, Arraypy)):
        Y = Y.To_list()

    if not len(X) == len(Y):
        raise ValueError('The vectors must be identical length')
    elif len(X) != g.rows:
        raise ValueError(
            'The vector fields and dimension of metric tensor must be identical length')

    # Calculation
    indices = range(len(X))
    scal = sum([g[i, j] * X[i] * Y[j] for i in indices
                for j in indices])
    # Output
    return scal


# ---------------- Christoffel_1 --------------------------------

def christoffel_1(g, var, type_output='t'):
    """Return the (-1,-1,-1)-tensor of Christoffel symbols for the given metric.
       This returns the Christoffel symbol of first kind that represents the
       Levi-Civita connection for the given metric.
       christoffel_1[i,j,k] = (diff(g[j,k],x[i])+diff(g[i,k],x[j])-diff(g[i,j],x[k]))/2.
       Example:
       ========
       >>> x1, x2 = symbols('x1, x2')
       >>> arg = [x1, x2]
       >>> A = Arraypy((2,2))
       >>> g = Tensor(A,(-1,-1))
       >>> g[0,0] = cos(x2)**2
       >>> g[0,1] = 0
       >>> g[1,0] = 0
       >>> g[1,1] = 1
       >>> print g
       >>> christoffel1=christoffel_1(g, arg, 't')

    """
    # Handling of input vector arguments var
    if not isinstance(var, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector arguments(var) must be a list, Arraypy or Tensor')
    if isinstance(var, (Tensor, Arraypy)):
        if len(var.shape) != 1:
            raise ValueError("The dimension of variables must be 1!")
        if isinstance(var, Tensor):
            if not var.type_pq == (1, 0):
                raise ValueError(
                    'The valence or ind_char of vector variables must be (+1)')
    if isinstance(var, (Tensor, Arraypy)):
        var = var.To_list()
    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    if not isinstance(g, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of metric tensor must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(g, (Arraypy, Tensor)):
            if isinstance(g, Tensor):
                if not g.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of metric tensor must be (-1,-1)')
            if not (g.to_matrix()).is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            if not (g.start_index[0] == g.start_index[1]):
                raise ValueError(
                    'The starting indices of metric tensor must be identical')
            idx_start = g.start_index[0]
        elif isinstance(g, Matrix):
            if not g.is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            idx_start = 0

    # The definition of diapason changes in an index
    [n1, n2] = g.shape
    if not n == n1:
        raise ValueError(
            'The rank of the metric tensor does not coincide with the number of variables.')

    indices = range(idx_start, idx_start + n)

    # Creating of output array with new indices
    Ch = Arraypy([3, n, idx_start])

    # Calculation
    for i in indices:
        for j in indices:
            for k in indices:
                Ch[i,
                   j,
                   k] = (diff(g[j,
                                k],
                              var[i - idx_start]) + diff(g[i,
                                                           k],
                                                         var[j - idx_start]) - diff(g[i,
                                                                                      j],
                                                                                    var[k - idx_start])) / 2

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        christoffel_1 = Ch.to_tensor((-1, -1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        christoffel_1 = Ch
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - Tensor.")

    # Output
    return christoffel_1


# ---------------- Christoffel_2 --------------------------------

def christoffel_2(g, var, type_output='t'):
    """Return the (-1, -1, +1) - tensor of Christoffel symbols for the given metric.
       This returns the Christoffel symbol of second kind that represents the
       Levi-Civita connection for the given metric.
       christoffel_2[i,j,k] =
       = Sum_{l}(g^{-1}[k,l]/2*(diff(g[j, l],x[i])+diff(g[i,l],x[j])-diff(g[i,j],x[l]))/2
       Example:
       ========
       >>> x1, x2 = symbols('x1, x2')
       >>> arg = [x1, x2]
       >>> A = Arraypy((2,2))
       >>> g = tensor(A,(-1,-1))
       >>> g[0,0] = cos(x2)**2
       >>> g[0,1] = 0
       >>> g[1,0] = 0
       >>> g[1,1] = 1
       >>> print g
       >>> christoffel_2(g,arg,'a')

    """
    # Handling of input vector arguments var
    if not isinstance(var, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector arguments(var) must be a list, Arraypy or Tensor')
    if isinstance(var, (Tensor, Arraypy)):
        if len(var.shape) != 1:
            raise ValueError("The dimension of variables must be 1!")
        if isinstance(var, Tensor):
            if not var.type_pq == (1, 0):
                raise ValueError(
                    'The valence or ind_char of vector variables must be (+1)')
    if isinstance(var, (Tensor, Arraypy)):
        var = var.To_list()
    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    if not isinstance(g, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of metric tensor must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(g, (Arraypy, Tensor)):
            if isinstance(g, Tensor):
                if not g.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of metric tensor must be (-1,-1)')
            if not (g.to_matrix()).is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            if not (g.start_index[0] == g.start_index[1]):
                raise ValueError(
                    'The starting indices of metric tensor must be identical')
            idx_start = g.start_index[0]
            g_inv = (g.to_matrix()).inv()
        elif isinstance(g, Matrix):
            if not g.is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            idx_start = 0
            g_inv = g.inv()

    # The definition of diapason changes in an index
    [n1, n2] = g.shape
    if not n == n1:
        raise ValueError(
            'The rank of the metric tensor does not coincide with the number of variables.')

    indices = range(idx_start, idx_start + n)

    # Creating of output array with new indices
    Ch = Arraypy([3, n, idx_start])

    # Calculation
    for i in indices:
        for j in indices:
            for k in indices:
                Ch[i,
                   j,
                   k] = Add(*[g_inv[k - idx_start,
                                    l - idx_start] * (diff(g[j,
                                                             l],
                                                           var[i - idx_start]) + diff(g[i,
                                                                                        l],
                                                                                      var[j - idx_start]) - diff(g[i,
                                                                                                                   j],
                                                                                                                 var[l - idx_start])) / 2 for l in indices])

    # Other variant calculation
    """christ_1 = christoffel_1(g,arg)
    for i in indices:
	for j in indices:
	    for k in indices:
		Ch[i,j,k] = Add(*[g_inv[k, l] *christ_1[i,j,l] for l in indices])"""

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        christoffel_2 = Ch.to_tensor((-1, -1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        christoffel_2 = Ch
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - Tensor.")

    # Output
    return christoffel_2


# ---------------- Covar_der --------------------------------

def covar_der(X, g, var, type_output='t'):
    """Return the covariant derivative the vector field.
       nabla X[i,j] = diff(X[j],x[i])+Sum_{k}(Gamma2[k,i,j]*X[k])
       Example:
       ========
       >>> x1, x2 = symbols('x1, x2')
       >>> arg = [x1, x2]
       >>> A = Arraypy((2,2))
       >>> g = Tensor(A,(-1,-1))
       >>> g[0,0] = cos(x2)**2
       >>> g[0,1] = 0
       >>> g[1,0] = 0
       >>> g[1,1] = 1
       >>> print g
       >>> X = [x1*x2**3,x1-cos(x2)]
       >>> covar_der(X, g, arg, 't')
    """
    # Handling of input vector arguments var
    if not isinstance(var, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector arguments(var) must be a list, Arraypy or Tensor')
    if isinstance(var, (Tensor, Arraypy)):
        if len(var.shape) != 1:
            raise ValueError("The dimension of variables must be 1!")
        if isinstance(var, Tensor):
            if not var.type_pq == (1, 0):
                raise ValueError(
                    'The valence or ind_char of vector variables must be (+1)')
    if isinstance(var, (Tensor, Arraypy)):
        var = var.To_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    if not isinstance(g, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of metric tensor must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(g, (Arraypy, Tensor)):
            if isinstance(g, Tensor):
                if not g.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of metric tensor must be (-1,-1)')
            if not (g.to_matrix()).is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            if not (g.start_index[0] == g.start_index[1]):
                raise ValueError(
                    'The starting indices of metric tensor must be identical')
            idx_g = g.start_index[0]
        elif isinstance(g, Matrix):
            if not g.is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            idx_g = 0

    # Handling of a input argument - vector field X
    if not isinstance(X, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector field must be list, Tensor or Arraypy')
    else:
        if isinstance(X, (Arraypy, Tensor)):
            if len(X.shape) != 1:
                raise ValueError("The dimension of X must be 1!")
            if isinstance(X, Tensor):
                if not X.type_pq == (1, 0):
                    raise ValueError(
                        'The valence or ind_char of vector field must be (+1)')
            idx_X = X.start_index[0]
        elif isinstance(X, list):
            idx_X = 0

    # The definition of diapason changes in an index
    [n1, n2] = g.shape
    if not n == n1:
        raise ValueError(
            'The rank of the metric Tensor does not coincide with the number of variables.')

    if (idx_g != idx_X):
        raise ValueError(
            'The start index of the metric tensor and vector field must be equal')
    else:
        idx_start = idx_g

    indices = range(idx_start, idx_start + n)

    # Creating of output array with new indices
    cov = Arraypy([2, n, idx_start])
    ch_2 = christoffel_2(g, var)
    # Calculation
    for i in indices:
        for j in indices:
            cov[i, j] = diff(X[j], var[i - idx_start]) + \
                Add(*[ch_2[k, i, j] * X[k] for k in indices])

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        cov_der = cov.to_tensor((-1, 1))
    elif type_output == str('a') or type_output == Symbol('a'):
        cov_der = cov
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - Tensor.")

    # Output
    return cov_der


# ---------------- Covar_der_XY --------------------------------

def covar_der_XY(X, Y, g, var, type_output='t'):
    """Return the covariant derivative the vector field along another field.
       nabla_XY[i] = Sum_{i}(nabla X[i,j]*Y[i])
       Example:
       >>> x1, x2 = symbols('x1, x2')
       >>> arg = [x1, x2]
       >>> A = Arraypy((2,2))
       >>> g = Tensor(A,(-1,-1))
       >>> g[0,0] = cos(x2)**2
       >>> g[0,1] = 0
       >>> g[1,0] = 0
       >>> g[1,1] = 1
       >>> print g
       >>> X = [x1*x2**3, x1-cos(x2)]
       >>> Y = [1, 2, 3]
       >>> covar_der_XY(X, Y, g, arg, 't')
    """
    # Handling of input vector arguments var
    if not isinstance(var, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector arguments(var) must be a list, Arraypy or Tensor')
    if isinstance(var, (Tensor, Arraypy)):
        if len(var.shape) != 1:
            raise ValueError("The dimension of variables must be 1!")
        if isinstance(var, Tensor):
            if not var.type_pq == (1, 0):
                raise ValueError(
                    'The valence or ind_char of vector variables must be (+1)')
    if isinstance(var, (Tensor, Arraypy)):
        var = var.To_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    if not isinstance(g, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of metric tensor must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(g, (Arraypy, Tensor)):
            if isinstance(g, Tensor):
                if not g.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of metric tensor must be (-1,-1)')
            if not (g.to_matrix()).is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            if not (g.start_index[0] == g.start_index[1]):
                raise ValueError(
                    'The starting indices of metric tensor must be identical')
            idx_g = g.start_index[0]
        elif isinstance(g, Matrix):
            if not g.is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            idx_g = 0

    # Handling of a input argument - vector field X
    if not isinstance(X, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector field must be list, Tensor or Arraypy')
    else:
        if isinstance(X, (Arraypy, Tensor)):
            if len(X.shape) != 1:
                raise ValueError("The dimension of X must be 1!")
            if isinstance(X, Tensor):
                if not X.type_pq == (1, 0):
                    raise ValueError(
                        'The valence or ind_char of vector field must be (+1)')
            idx_X = X.start_index[0]
        elif isinstance(X, list):
            idx_X = 0

    # Handling of a input argument - vector field Y
    if not isinstance(Y, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector field must be list, Tensor or Arraypy')
    else:
        if isinstance(Y, (Arraypy, Tensor)):
            if len(Y.shape) != 1:
                raise ValueError("The dimension of Y must be 1!")
            if isinstance(Y, Tensor):
                if not Y.type_pq == (1, 0):
                    raise ValueError(
                        'The valence or ind_char of vector field must be (+1)')
            idx_Y = Y.start_index[0]
        elif isinstance(Y, list):
            idx_Y = 0

    [n1, n2] = g.shape
    if not len(X) == len(Y):
        raise ValueError('The vectors must be identical length')
    elif not idx_X == idx_Y:
        raise ValueError('The start index of vector fields must be equal')
    elif not(idx_g == idx_X):
        raise ValueError(
            'The start index of the metric tensor and vector field must be equal')
    else:
        idx_start = idx_g
    if len(X) != n1:
        raise ValueError(
            'The vector fields and dimension of metric tensor must be identical length')

    # The definition of diapason changes in an index
    if not n == n1:
        raise ValueError(
            'The rank of the metric tensor does not concide with the number of variables.')
    indices = range(idx_start, idx_start + n)

    # Creating of output array with new indices
    nabla_XY = Arraypy([1, n, idx_start])
    nabla_X = covar_der(X, g, var)

    # Calculation
    for j in indices:
        nabla_XY[j] = sum([nabla_X[i, j] * Y[i] for i in indices])

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        cov_der_XY = nabla_XY.to_tensor((1))
    elif type_output == str('a') or type_output == Symbol('a'):
        cov_der_XY = nabla_XY
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - Tensor.")

    # Output
    return cov_der_XY


# ---------------- Riemann --------------------------------

def riemann(g, var, type_output='t'):
    """Return the Riemann curvature tensor of type (-1,-1,-1,+1)  for the given metric tensor.
       Riemann[i,j,k,l] = diff(Gamma_2[j,k,l],x[i])-diff(Gamma_2[i,k,l],x[j]) +
       + Sum_{p}( Gamma_2[i,p,l]*Gamma_2[j,k,p] -Gamma_2[j,p,l]*Gamma_2[i,k,p]
       Example:
       ========
       >>> x1, x2 = symbols('x1, x2')
       >>> arg = [x1, x2]
       >>> A = Arraypy((2,2))
       >>> g = Tensor(A,(-1,-1))
       >>> g[0,0] = cos(x2)**2
       >>> g[0,1] = 0
       >>> g[1,0] = 0
       >>> g[1,1] = 1
       >>> print g
       >>> riemann(g, arg, 'a')
    """
    # Handling of input vector arguments var
    if not isinstance(var, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector arguments(var) must be a list, Arraypy or Tensor')
    if isinstance(var, (Tensor, Arraypy)):
        if len(var.shape) != 1:
            raise ValueError("The dimension of variables must be 1!")
        if isinstance(var, Tensor):
            if not var.type_pq == (1, 0):
                raise ValueError(
                    'The valence or ind_char of vector variables must be (+1)')
    if isinstance(var, (Tensor, Arraypy)):
        var = var.To_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    if not isinstance(g, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of metric tensor must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(g, (Arraypy, Tensor)):
            if isinstance(g, Tensor):
                if not g.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of metric tensor must be (-1,-1)')
            if not (g.to_matrix()).is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            if not (g.start_index[0] == g.start_index[1]):
                raise ValueError(
                    'The starting indices of metric tensor must be identical')
            idx_start = g.start_index[0]
        elif isinstance(g, Matrix):
            if not g.is_symmetric():
                raise ValueError('The metric tensor must be symmetric.')
            idx_start = 0

    # The definition of diapason changes in an index
    [n1, n2] = g.shape
    if not n == n1:
        raise ValueError(
            'The rank of the metric tensor does not coincide with the number of variables.')

    indices = range(idx_start, idx_start + n)

    # Creating of output array with new indices
    R = Arraypy([4, n, idx_start])
    ch_2 = christoffel_2(g, var)

    # Calculation
    for i in indices:
        for j in indices:
            for k in indices:
                for l in indices:
                    R[i,
                      j,
                      k,
                      l] = diff(ch_2[j,
                                     k,
                                     l],
                                var[i - idx_start]) - diff(ch_2[i,
                                                                k,
                                                                l],
                                                           var[j - idx_start]) + sum([ch_2[i,
                                                                                           p,
                                                                                           l] * ch_2[j,
                                                                                                     k,
                                                                                                     p] - ch_2[j,
                                                                                                               p,
                                                                                                               l] * ch_2[i,
                                                                                                                         k,
                                                                                                                         p] for p in indices])

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        riemann = R.to_tensor((-1, -1, -1, 1))
    elif type_output == str('a') or type_output == Symbol('a'):
        riemann = R
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - Tensor.")

    # Output
    return riemann


# ---------------- Ricci --------------------------------

def ricci(riemann, var, type_output='t'):
    """Return the tensor Ricci of type (-1,-1), is symmetric tensor
       for given Riemann curvature tensor.
       Ricci[j,k] = Sum_{i}(Riemann[i,j,k,i])
       Example:
       ========
       >>> x1, x2 = symbols('x1, x2')
       >>> arg = [x1, x2]
       >>> A = Arraypy((2,2))
       >>> g = Tensor(A,(-1,-1))
       >>> g[0,0] = cos(x2)**2
       >>> g[0,1] = 0
       >>> g[1,0] = 0
       >>> g[1,1] = 1
       >>> print g
       >>> cur = riemann(g, arg, 't')
       >>> r = ricci(cur, arg, 't')
    """
    # Handling of input vector arguments var
    if not isinstance(var, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector arguments(var) must be a list, Arraypy or Tensor')
    if isinstance(var, (Tensor, Arraypy)):
        if len(var.shape) != 1:
            raise ValueError("The dimension of variables must be 1!")
        if isinstance(var, Tensor):
            if not var.type_pq == (1, 0):
                raise ValueError(
                    'The valence or ind_char of vector variables must be (+1)')
    if isinstance(var, (Tensor, Arraypy)):
        var = var.To_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument Riemann curvature tensor - riemann
    if not isinstance(riemann, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of Riemann curvature tensor must be Matrix, Arraypy or Tensor')
    else:
        if isinstance(riemann, (Arraypy, Tensor)):
            if isinstance(riemann, Tensor):
                if not riemann.type_pq == (1, 3):
                    raise ValueError(
                        'The valence or ind_char of Riemann curvature tensor must be (-1,-1,-1,+1)')
                if not (
                    count(
                        riemann.start_index.count(
                            riemann.start_index[0]) == 4)):
                    raise ValueError(
                        'The starting indices of Riemann curvature tensor must be identical')
            idx_start = riemann.start_index[0]
        else:
            idx_start = 0

    # The definition of diapason changes in an index
    [n1, n2, n3, n4] = riemann.shape
    if not n == n1:
        raise ValueError(
            'The rank of the Riemann curvature tensor does not coincide with the number of variables.')

    indices = range(idx_start, idx_start + n)

    # Creating of output array with new indices
    Ri = Arraypy([2, n, idx_start])

    # Calculation
    for j in indices:
        for k in indices:
            Ri[j, k] = sum([riemann[i, j, k, i] for i in indices])

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        ricci = Ri.to_tensor((-1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        ricci = Ri
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - Tensor.")

    # Output
    return ricci


# ---------------- Scal_curv --------------------------------

def scal_curv(g, ricci, var):
    """ The scalar curvature (or the Ricci scalar)
        is the simplest curvature invariant of a Riemannian manifold.
        S=Ricci[j,k]*g_inv[j,k]
        Example:
        ========
        >>> x1, x2 = symbols('x1, x2')
        >>> arg = [x1, x2]
        >>> A = Arraypy((2,2))
        >>> g = Tensor(A,(-1,-1))
        >>> g[0,0] = cos(x2)**2
        >>> g[0,1] = 0
        >>> g[1,0] = 0
        >>> g[1,1] = 1
        >>> print g
        >>> scal_curv(g, r, arg)
    """
    # Handling of input vector arguments var
    if not isinstance(var, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector arguments(var) must be a list, Arraypy or Tensor')
    if isinstance(var, (Tensor, Arraypy)):
        if len(var.shape) != 1:
            raise ValueError("The dimension of variables must be 1!")
        if isinstance(var, Tensor):
            if not var.type_pq == (1, 0):
                raise ValueError(
                    'The valence or ind_char of vector variables must be (+1)')
    if isinstance(var, (Tensor, Arraypy)):
        var = var.To_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    if not isinstance(g, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of metric Tensor must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(g, (Arraypy, Tensor)):
            if isinstance(g, Tensor):
                if not g.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of metric tensor must be (-1,-1)')
            g = g.to_matrix()
    if not g.is_symmetric():
        raise ValueError('The metric tensor must be symmetric.')
    # The definition of inverse matrix of the metric tensor
    g_inv = g.inv()

    # Handling of a input argument tensor Ricci - ricci
    if not isinstance(ricci, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of tensor Ricci must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(ricci, (Arraypy, Tensor)):
            if isinstance(ricci, Tensor):
                if not ricci.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of tensor Ricci must be (-1,-1)')
            ricci = ricci.to_matrix()
    if not ricci.is_symmetric():
        raise ValueError('The Ricci tensor must be symmetric.')

    if not (g.shape == ricci.shape):
        raise ValueError(
            'The rank of the metric tensor does not coincide with the rank of tensor Ricci.')

    # The definition of diapason changes in an index
    [n1, n2] = g.shape
    if not n == n1:
        raise ValueError(
            'The rank of the metric tensor does not coincide with the number of variables.')

    # Calculation
    indices = range(n)
    for i in indices:
        for j in indices:
            scal_curv = g_inv[i, j] * ricci[i, j]
    # Output
    return scal_curv


#-----------------K_sigma----------------------------

def k_sigma(X, Y, R, g, var):
    """Return Sectional curvature of thу Riemannian space
       in the direction за two-dimensional area formed by
       vectors X, Y  for the given metric tensor.
       K_sigma = Sum_{i,j,k,r,s}( g[r,s]*Riemann[i,j,k,r] *X[i]*Y[j]*Y[k]X[s])/ (scal_prod(X,X,g)*scal_prod(Y,Y,g) - scal_prod(X,Y,g)^2
       Example:
       ========
       >>> x1, x2 = symbols('x1, x2')
       >>> arg = [x1, x2]
       >>> X = [1,2]
       >>> Y = [3,4]
       >>> A = Arraypy((2,2))
       >>> g = Tensor(A,(-1,-1))
       >>> g[0,0] = cos(x2)**2
       >>> g[0,1] = 0
       >>> g[1,0] = 0
       >>> g[1,1] = 1
       >>> print g
       >>> R = riemann(g, arg)
       >>> K_sig = k_sigma(X, Y, R, g, var)
    """
    # Handling of input vector arguments var
    if not isinstance(var, (list, Arraypy, Tensor)):
        raise TypeError(
            'The type of vector arguments(var) must be a list, Arraypy or Tensor')
    if isinstance(var, (Tensor, Arraypy)):
        if len(var.shape) != 1:
            raise ValueError("The dimension of variables must be 1!")
        if isinstance(var, Tensor):
            if not var.type_pq == (1, 0):
                raise ValueError(
                    'The valence or ind_char of vector variables must be (+1)')
    if isinstance(var, (Tensor, Arraypy)):
        var = var.To_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    if not isinstance(g, (Matrix, Tensor, Arraypy)):
        raise TypeError(
            'The type of metric tensor must be Matrix, Tensor or Arraypy')
    else:
        if isinstance(g, (Arraypy, Tensor)):
            if isinstance(g, Tensor):
                if not g.type_pq == (0, 2):
                    raise ValueError(
                        'The valence or ind_char of metric tensor must be (-1,-1)')
            g = g.to_matrix()
    if not g.is_symmetric():
        raise ValueError('The metric tensor must be symmetric.')

    # Handling of a input arguments - vector or vector fields X
    if not isinstance(X, (list, Arraypy, Tensor)):
        raise TypeError('The type of vector must be list, Arraypy or Tensor')
    if isinstance(X, (Tensor, Arraypy)):
        if len(X.shape) != 1:
            raise ValueError("The dimension of X must be 1!")
        if isinstance(X, Tensor):
            if not X.type_pq == (1, 0):
                raise ValueError('The valence or ind_char of X must be (+1)')
    if isinstance(X, (Tensor, Arraypy)):
        X = X.To_list()

    # Handling of a input arguments - vector or vector fields Y
    if not isinstance(Y, (list, Arraypy, Tensor)):
        raise TypeError('The type of vector must be list, Arraypy or Tensor')
    if isinstance(Y, (Tensor, Arraypy)):
        if len(Y.shape) != 1:
            raise ValueError("The dimension of Y must be 1!")
        if isinstance(Y, Tensor):
            if not Y.type_pq == (1, 0):
                raise ValueError('The valence or ind_char of Y must be (+1)')
    if isinstance(Y, (Tensor, Arraypy)):
        Y = Y.To_list()

    if not len(X) == len(Y):
        raise ValueError('The vectors must be identical length')
    elif len(X) != g.rows:
        raise ValueError(
            'The vector fields and dimension of metric tensor must be identical length')

    # Handling of a input argument Riemann curvature tensor - R
    if not isinstance(R, (Matrix, Arraypy, Tensor)):
        raise TypeError(
            'The type of Riemann curvature tensor must be Matrix, Arraypy or Tensor')
    else:
        if isinstance(R, (Arraypy, Tensor)):
            if isinstance(R, Tensor):
                if not R.type_pq == (1, 3):
                    raise ValueError(
                        'The valence or ind_char of Riemann curvature tensor must be (-1,-1,-1,+1)')
                if not (R.start_index[0] == R.start_index[1]):
                    raise ValueError(
                        'The starting indices of Riemann curtivate tensor must be identical')
            idx_R = R.start_index[0]

    # The definition of diapason changes in an index
    [n1, n2] = g.shape
    if not n == n1:
        raise ValueError(
            'The rank of the metric tensor does not coincide with the number of variables.')
    [n1, n2, n3, n4] = R.shape
    if not n == n1:
        raise ValueError(
            'The rank of the Riemann curvature tensor does not concide with the number of variables.')

    indices = range(len(X))

    # Calculation
    Sc_pr = scal_prod(X, X, g) * scal_prod(Y, Y, g) - scal_prod(X, Y, g)**2
    if (Sc_pr == 0):
        raise ValueError('The two-dimensional area is a degenerate!')

    numerator = sum([g[r, s] * R[i + idx_R, j + idx_R, k + idx_R, r + idx_R] * X[i] * Y[j] * Y[k] * X[s] for i in indices
                     for j in indices
                     for k in indices
                     for r in indices
                     for s in indices])

    k_sigma = simplify(numerator / Sc_pr)

    # Output
    return k_sigma
