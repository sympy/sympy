# -*- coding: utf-8 -*-

from sympy.matrices import Matrix
from sympy.core import Add, diff, Symbol
from sympy.simplify import simplify
from sympy.tensor.arraypy import Arraypy, TensorArray, matrix2arraypy, \
    matrix2tensor, list2arraypy, list2tensor
from sympy.tensor.tensor_methods import is_symmetric
from sympy.tensor.helper_functions import check_vector_of_arguments, \
    check_metric_tensor, check_the_vector_field, replace_index_to_k, \
    check_the_christoffel_symbols_2

"""Module riemannian_geometry contains functions for work with tensor fields:
- the calculation of the scalar product;
- the Christoffel symbols of the first and second kind;
- the covariant derivative of the curvature tensor;
- the Ricci tensor;
- scalar and sectional curvature;
- the covariant derivative the tensor field;
- the covariant divergence of a tensor field;
- the Riemann curvature tensor and sectional curvature for left-invariant metric;
- the product of Kulkarni-Nomizu;
- the Gaussian curvature;
- the second quadratic form.


To implement the functions used modules: matrices and tensor
(with classes arraypy and tensor). All functions take arguments,
the types of which may be such as list, matrix, or array Arraypy tensor.
Some functions have optional parameter indicating the type of the function result.
Starting index of arguments with type Arraypy or TensorArray is not necessarily
and by default equal to 0. The function determines the range of the index
in array to return the object with the same range of index.

Functions are work with multidimensional arrays Arraypy and tensors,
classes and methods are contained in the module Arraypy.

"""


def scal_prod(X, Y, g):
    """Returns scalar product of vectors g(X,Y).

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import scal_prod
    >>> from sympy import symbols, cos
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> x1, x2 = symbols('x1, x2')

    X, Y it's a vector or a vector field. They can be a list,
    one-dimensional arraypy or TensorArray with valence of indices (+1):

    >>> X = [1, 2]
    >>> Y = [3, 4]

    g it's a metric tensor must be symmetric matrix, array of arraypy or
    covariant tensor with valence of indices (-1, -1):

    >>> A = Arraypy((2, 2))
    >>> g = TensorArray(A,(-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    The scalar product:
    >>> sc = scal_prod(X, Y, g)
    >>> print(sc)
    3*cos(x2)**2 + 8

    """
    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
        g = g.to_matrix()
    if not g.is_symmetric():
        raise ValueError('The metric tensor must be symmetric.')

    # Handling of a input arguments - vector or vector fields X
    check_the_vector_field(X)

    if isinstance(X, (TensorArray, Arraypy)):
        X = X.to_list()

    # Handling of a input arguments - vector or vector fields Y
    check_the_vector_field(Y)
    if isinstance(Y, (TensorArray, Arraypy)):
        Y = Y.to_list()

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


def christoffel_1(g, var, type_output='t'):
    """Return the (-1,-1,-1) - tensor of Christoffel symbols for the given metric.
    This returns the Christoffel symbol of first kind that represents the
    Levi-Civita connection for the given metric.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import christoffel_1
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var is a list of symbolic arguments. May be a list, one-dimensional
    arraypy or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]
    >>> A = Arraypy((2, 2))
    >>> g = TensorArray(A,(-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The Christoffel symbols of the first kind:
    >>> ch_1 = christoffel_1(g, var, 't')
    >>> print(ch_1)
    0  sin(x2)*cos(x2)
    -sin(x2)*cos(x2)  0
    -sin(x2)*cos(x2)  0
    0  0
    >>> ch_1.type_pq
    (0, 3)

    """
    # Handling of input vector of arguments - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()
    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
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
            "The parameter of type output result must 'a' - Arraypy or 't' and None - TensorArray.")

    # Output
    return christoffel_1


def christoffel_2(g, var, type_output='t'):
    """Return the (1, -1, -1) - tensor of Christoffel symbols for the given metric.
    This returns the Christoffel symbol of second kind that represents the
    Levi-Civita connection for the given metric.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import christoffel_2
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional
    arraypy or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valence indices (-1, -1):

    >>> A = Arraypy((2, 2))
    >>> g = TensorArray(A,(-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The Christoffel symbols of the second kind:
    >>> ch_2 = christoffel_2(g, var, 'a')
    >>> print(ch_2)
    0  sin(x2)*cos(x2)
    -sin(x2)/cos(x2)  0
    -sin(x2)/cos(x2)  0
    0  0

    """
    # Handling of input vector of arguments - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()
    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
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
    # christ_1 = christoffel_1(g, var)
    # for i in indices:
        # for j in indices:
            # for k in indices:
                # Ch[i,
                # j,
                # k] = Add(*[g_inv[k,
                # l] *christ_1[i,
                # j,
                # l] for l in indices])

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        christoffel_2 = Ch.to_tensor((1, -1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        christoffel_2 = Ch
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - TensorArray.")

    # Output
    return christoffel_2


def covar_der(X, g, var, type_output='t'):
    """Return the covariant derivative the vector field.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import covar_der
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional
    arraypy or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valence indices (-1, -1):

    >>> A = Arraypy((2, 2))
    >>> g = TensorArray(A,(-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    X it's vector field can be a list, one-dimensional arraypy, or one-dimensional
    tensor with valences of indices (+1):

    >>> X = [x1 * x2**3, x1 - cos(x2)]

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The covariant derivative:
    >>> c_v = covar_der(X, g, var, 't')
    >>> print(c_v)
    x2**3 - (x1 - cos(x2))*sin(x2)/cos(x2)  x1*x2**3*sin(x2)*cos(x2) + 1
    -x1*x2**3*sin(x2)/cos(x2) + 3*x1*x2**2  sin(x2)
    >>> c_v.type_pq
    (1, 1)

    """
    # Handling of input vector of arguments - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
        if not (g.start_index[0] == g.start_index[1]):
            raise ValueError(
                'The starting indices of metric tensor must be identical')
        idx_g = g.start_index[0]
    elif isinstance(g, Matrix):
        if not g.is_symmetric():
            raise ValueError('The metric tensor must be symmetric.')
        idx_g = 0

    # Handling of a input argument - vector field X
    check_the_vector_field(X)

    if isinstance(X, (Arraypy, TensorArray)):
        idx_X = X.start_index[0]
    elif isinstance(X, list):
        idx_X = 0

    # The definition of diapason changes in an index
    [n1, n2] = g.shape
    if not n == n1:
        raise ValueError(
            'The rank of the metric tensor does not coincide with the number of variables.')

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
        cov_der = cov.to_tensor((1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        cov_der = cov
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - TensorArray.")

    # Output
    return cov_der


def covar_der_xy(X, Y, g, var, type_output='t'):
    """Return the covariant derivative the vector field along another field.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import covar_der_xy
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional arraypy
    or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valences indices (-1, -1):

    >>> A = Arraypy((2, 2))
    >>> g = TensorArray(A,(-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    X, Y it's vector fields may be lists, one-dimensional arraypy,
    or one-dimensional tensor indices with valences (+ 1):

    >>> X = [x1 * x2**3, x1 - cos(x2)]
    >>> Y = [1, 2]

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The covariant derivative along another vector field:
    >>> c_v_XY = covar_der_xy(X, Y, g, var, 't')
    >>> print(c_v_XY)
    -2*x1*x2**3*sin(x2)/cos(x2) + 6*x1*x2**2 + x2**3 - (x1 - cos(x2))*sin(x2)/cos(x2) \
    x1*x2**3*sin(x2)*cos(x2) + 2*sin(x2) + 1

    """
    # Handling of input vector of arguments - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
        if not (g.start_index[0] == g.start_index[1]):
            raise ValueError(
                'The starting indices of metric tensor must be identical')
        idx_g = g.start_index[0]
    elif isinstance(g, Matrix):
        if not g.is_symmetric():
            raise ValueError('The metric tensor must be symmetric.')
        idx_g = 0

    # Handling of a input argument - vector field X
    check_the_vector_field(X)

    if isinstance(X, (Arraypy, TensorArray)):
        idx_X = X.start_index[0]
    elif isinstance(X, list):
        idx_X = 0

    # Handling of a input argument - vector field Y
    check_the_vector_field(Y)
    if isinstance(Y, (Arraypy, TensorArray)):
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
            "The parameter of type output result must 'a' - Arraypy or 't' and None - TensorArray.")

    # Output
    return cov_der_XY


def riemann(g, var, type_output='t'):
    """Return the Riemann curvature tensor of type (1, -1, -1, -1)
    for the given metric tensor.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import riemann
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional
    arraypy or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valence indices (-1, -1):

    >>> A = Arraypy((2, 2))
    >>> g = TensorArray(A,(-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The curvature tensor:
    >>> r = riemann(g, var, 'a')
    >>> print(r)
    0  0
    0  0
    0  -cos(x2)**2
    1  0
    0  cos(x2)**2
    -1  0
    0  0
    0  0

    """
    # Handling of input vector of arguments - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
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
        riemann = R.to_tensor((1, -1, -1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        riemann = R
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - TensorArray.")

    # Output
    return riemann


def ricci(riemann, var, type_output='t'):
    """Return the tensor Ricci of type (-1, -1), is symmetric tensor
    for given Riemann curvature tensor.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import ricci, riemann
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional
    arraypy or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valence indices (-1, -1):

    >>> A = Arraypy((2,2))
    >>> g = TensorArray(A,(-1,-1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    riemann it's a Riemann curvature tensor must be symmetric matrix,
    arraypy or tensor with valences indices (-1, -1, -1, 1):

    >>> cur = riemann(g, var, 't')

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The Ricci tensor:
    >>> r = ricci(cur, var, 't')
    >>> print(r)
    cos(x2)**2  0
    0  1
    >>> r.type_pq
    (0, 2)
    """
    # Handling of input vector of arguments - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument Riemann curvature tensor - riemann
    if not isinstance(riemann, (Matrix, Arraypy, TensorArray)):
        raise TypeError(
            'The type of Riemann curvature tensor must be Matrix, Arraypy or TensorArray')
    else:
        if isinstance(riemann, (Arraypy, TensorArray)):
            if isinstance(riemann, TensorArray):
                if not riemann.type_pq == (1, 3):
                    raise ValueError(
                        'The valence of Riemann curvature tensor must be (1, -1, -1, -1)')
                if not (
                    riemann.start_index.count(
                        riemann.start_index[0]) == 4):
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
            "The parameter of type output result must 'a' - Arraypy or 't' and None - TensorArray.")

    # Output
    return ricci


def scal_curv(g, ricci, var):
    """The scalar curvature (or the Ricci scalar) is the simplest curvature
    invariant of a Riemannian manifold.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import scal_curv, ricci, riemann
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional
    arraypy or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valence indices (-1, -1):

    >>> A = Arraypy((2,2))
    >>> g = TensorArray(A,(-1,-1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    riemann it's a Riemann curvature tensor must be symmetric matrix,
    arraypy or tensor with valences indices (-1, -1, -1, 1):

    >>> cur = riemann(g, var, 't')

    ricci it's Ricci tensor must be a matrix, arraypy or valences with
    tensor indices (-1, -1):

    >>> r = ricci(cur, var, 't')

    The Ricci tensor for the Riemann curvature tensor:
    >>> sc_c = scal_curv(g, r, var)
    >>> print(sc_c)
    1

    """
    # Handling of input vector of arguments - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
        g = g.to_matrix()
    if not g.is_symmetric():
        raise ValueError('The metric tensor must be symmetric.')
    # The definition of inverse matrix of the metric tensor
    g_inv = g.inv()

    # Handling of a input argument tensor Ricci - ricci
    if not isinstance(ricci, (Matrix, Arraypy, TensorArray)):
        raise TypeError(
            'The type of tensor Ricci must be Matrix, TensorArray or Arraypy')
    else:
        if isinstance(ricci, (Arraypy, TensorArray)):
            if isinstance(ricci, TensorArray):
                if not ricci.type_pq == (0, 2):
                    raise ValueError(
                        'The valence of tensor Ricci must be (-1,-1)')
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


def k_sigma(X, Y, R, g, var):
    """Return Sectional curvature of thу Riemannian space
    in the direction за two-dimensional area formed by
    vectors X, Y  for the given metric tensor.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import k_sigma, riemann
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional
    arraypy or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    X, Y it's a vector or a vector field. They can be a list, one-dimensional
    arraypy or tensor with valence of indices (+1):

    >>> X = [1, 2]
    >>> Y = [3, 4]

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valence indices (-1, -1):

    >>> A = Arraypy((2, 2))
    >>> g = TensorArray(A,(-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    R it's a Riemann curvature tensor must be symmetric matrix, arraypy or tensor
    with valences indices (1, -1, -1, -1):

    >>> R = riemann(g, var)

    The sectional curvature:
    >>> k_sig = k_sigma(X, Y, R, g, var)
    >>> print(k_sig)
    1
    """
    # Handling of input vector of arguments - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
        g = g.to_matrix()
    if not g.is_symmetric():
        raise ValueError('The metric tensor must be symmetric.')

    # Handling of a input arguments - vector or vector fields X
    check_the_vector_field(X)

    if isinstance(X, (TensorArray, Arraypy)):
        X = X.to_list()

    # Handling of a input arguments - vector or vector fields Y
    check_the_vector_field(Y)

    if isinstance(Y, (TensorArray, Arraypy)):
        Y = Y.to_list()

    if not len(X) == len(Y):
        raise ValueError('The vectors must be identical length')
    elif len(X) != g.rows:
        raise ValueError(
            'The vector fields and dimension of metric tensor must be identical length')

    # Handling of a input argument Riemann curvature tensor - R
    if not isinstance(R, (Matrix, Arraypy, TensorArray)):
        raise TypeError(
            'The type of Riemann curvature tensor must be Matrix, Arraypy or TensorArray')
    else:
        if isinstance(R, (Arraypy, TensorArray)):
            if isinstance(R, TensorArray):
                if not R.type_pq == (1, 3):
                    raise ValueError(
                        'The valence of Riemann curvature tensor must be (1, -1,- 1, -1)')
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


def nabla(T, ch_2, var):
    """Return the covariant derivative the tensor field.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import nabla
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy import symbols, cos, sin
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional arraypy
    or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    T it's a tensor field must be tensor:

    >>> T = Arraypy([2, 2, 0]).to_tensor((1, -1))
    >>> T[0,0] = x2
    >>> T[0,1] = -x2
    >>> T[1,0] = -x1
    >>> T[1,1] = x1

    ch_2 it's a Christoffel symbol of second kind must be arraypy or tensor
    with valence indices (1, -1, -1):

    >>> ch_2 = Arraypy([3, 2, 0]).to_tensor((1, -1, -1))
    >>> ch_2[0,0,0] = 0
    >>> ch_2[0,0,1] = sin(x2)*cos(x2)
    >>> ch_2[0,1,1] = 0
    >>> ch_2[1,1,1] = 0
    >>> ch_2[1,0,1] = 0
    >>> ch_2[1,1,0] = 0
    >>> ch_2[1,0,0] = -sin(x2)*cos(x2)
    >>> ch_2[0,1,0] = -sin(x2)*cos(x2)

    The covariant derivative of tensor field:
    >>> nabla_t = nabla(T, ch_2, var)
    >>> print(nabla_t)
    -x1*sin(x2)*cos(x2) + x2*sin(x2)*cos(x2)  0  
    x1*sin(x2)*cos(x2) + x2*sin(x2)*cos(x2)  x2*sin(x2)*cos(x2) - 1  
    -x1*sin(x2)*cos(x2) - x2*sin(x2)*cos(x2)  -x1*sin(x2)*cos(x2) - 1  
    -x1*sin(x2)*cos(x2) + x2*sin(x2)*cos(x2)  0
    """
    # Handling of a input argument - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Handling of a input argument - Christoffel symbol of second kind
    check_the_christoffel_symbols_2(ch_2)

    idx_ch = ch_2.start_index[0]

    # Handling of a input argument - tensor field T
    if not isinstance(T, TensorArray):
        raise TypeError(
            'The type of tensor field must be TensorArray')

    idx_start_T = T.start_index[0]

    if (idx_start_T != idx_ch):
        raise ValueError(
            'The start index of the tensor field and Christoffel symbol \
            of second kind must be equal')

    # The definition of diapason changes in an index
    # The number of upper indices
    p = T.type_pq[0]

    # The dimension of the input array
    n = T.shape[0]

    # The rank of the input array
    rank_T = len(T.shape)

    # The definition of the start index
    idx_char_T = T.ind_char

    idx_char_nabla_T = list(idx_char_T) + [-1]

    # upper_idx_numbers it is a list with the positions on which are the upper
    # indices
    upper_idx_numbers = [
        k for k in range(len(idx_char_T)) if idx_char_T[k] == 1]

    # low_idx_numbers it is a list with the positions on which are the lower
    # indices
    low_idx_numbers = [
        k for k in range(len(idx_char_T)) if idx_char_T[k] == -1]

    # Creating the output array in accordance with the start index
    nabla_T = Arraypy([rank_T + 1, n, idx_start_T]).to_tensor(idx_char_nabla_T)

    index_nabla_T = nabla_T.index_list

    # Calculation
    for index in index_nabla_T:
        index_T = list(index)
        del index_T[n]
        index_T = tuple(index_T)
        s = index[rank_T]
        dt = diff(T[index_T], var[index[s]])
        k = idx_start_T
        nabla_T_up = 0
        nabla_T_lo = 0
        while k < n + idx_start_T:
            for i in upper_idx_numbers:
                index_T_ik = replace_index_to_k(index_T, i, k)
                nabla_T_up += T[index_T_ik] * ch_2[index_T[i], s, k]
            for j in low_idx_numbers:
                index_T_jk = replace_index_to_k(index_T, j, k)
                nabla_T_lo += T[index_T_jk] * ch_2[index_T[j], s, k]
            k = k + 1
        nabla_T[index] = dt + nabla_T_up - nabla_T_lo

    # Output
    return nabla_T


def nabla_x(T, ch_2, X, var):
    """Return the covariant derivative the tensor field along another vector field.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import nabla_x
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy import symbols, cos, sin
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional arraypy
    or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    T it's a tensor field must be tensor:

    >>> T = Arraypy([2, 2, 0]).to_tensor((1, -1))
    >>> T[0,0] = x2
    >>> T[0,1] = -x2
    >>> T[1,0] = -x1
    >>> T[1,1] = x1

    ch_2 it's a Christoffel symbol of second kind must be arraypy or tensor
    with valence indices (1, -1, -1):

    >>> ch_2 = Arraypy([3, 2, 0]).to_tensor((1, -1, -1))
    >>> ch_2[0,0,0] = 0
    >>> ch_2[0,0,1] = sin(x2)*cos(x2)
    >>> ch_2[0,1,1] = 0
    >>> ch_2[1,1,1] = 0
    >>> ch_2[1,0,1] = 0
    >>> ch_2[1,1,0] = 0
    >>> ch_2[1,0,0] = -sin(x2)*cos(x2)
    >>> ch_2[0,1,0] = -sin(x2)*cos(x2)

    X it's vector field can be a list, one-dimensional arraypy, or one-dimensional
    tensor with valences of indices (+1):

    >>> X = [x1 * x2**3, x1 - cos(x2)]

    The covariant derivative of tensor field along another vector field:
    >>> nabla_xt = nabla_x(T, ch_2, X, var)
    >>> print(nabla_xt)
    x1*x2**3*(-x1*sin(x2)*cos(x2) + x2*sin(x2)*cos(x2))  x1*x2**3*(x1*sin(x2)*cos(x2) + \
    x2*sin(x2)*cos(x2)) + (x1 - cos(x2))*(x2*sin(x2)*cos(x2) - 1)  
    x1*x2**3*(-x1*sin(x2)*cos(x2) - x2*sin(x2)*cos(x2)) + \
    (x1 - cos(x2))*(-x1*sin(x2)*cos(x2) - 1)  x1*x2**3*(-x1*sin(x2)*cos(x2) + x2*sin(x2)*cos(x2)) 
    """
    # Handling of a input argument - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Handling of a input argument - Christoffel symbol of second kind
    check_the_christoffel_symbols_2(ch_2)

    idx_ch = ch_2.start_index[0]

    # Handling of a input argument - vector field X
    check_the_vector_field(X)

    if isinstance(X, (Arraypy, TensorArray)):
        idx_X = X.start_index[0]
    elif isinstance(X, list):
        idx_X = 0

    # Handling of a input argument - tensor field T
    if not isinstance(T, TensorArray):
        raise TypeError(
            'The type of tensor field must be TensorArray')

    idx_start_T = T.start_index[0]

    if (idx_start_T != idx_ch != idx_X):
        raise ValueError(
            'The start index of the tensor field and Christoffel symbol \
            of second kind and vector field must be equal')

    # The definition of diapason changes in an index
    # The number of upper indices
    p = T.type_pq[0]

    # The dimension of the input array
    n = T.shape[0]

    # The rank of the input array
    rank_T = len(T.shape)

    # The definition of the start index
    idx_char_T = T.ind_char

    # Creating the output array in accordance with the start index
    nabla_TX = Arraypy([rank_T, n, idx_start_T]).to_tensor(idx_char_T)

    index_nabla_TX = nabla_TX.index_list

    nabla_T = nabla(T, ch_2, var)
    # Calculation
    for index in index_nabla_TX:
        k = idx_start_T
        while k < n + idx_start_T:
            idx_nabla_T = tuple(list(index) + [k])
            nabla_TX[index] += nabla_T[idx_nabla_T] * X[k]
            k = k + 1

    # Output
    return nabla_TX


def delta(T, g, ch_2, var):
    """Return the covariant divergence of a tensor field T.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import delta
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy import symbols, cos, sin
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional arraypy
    or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    T it's a tensor field must be tensor:

    >>> T = Arraypy([2, 2, 0]).to_tensor((1, -1))
    >>> T[0,0] = x2
    >>> T[0,1] = -x2
    >>> T[1,0] = -x1
    >>> T[1,1] = x1

    ch_2 it's a Christoffel symbol of second kind must be arraypy or tensor
    with valence indices (1, -1, -1):

    >>> ch_2 = Arraypy([3, 2, 0]).to_tensor((1, -1, -1))
    >>> ch_2[0,0,0] = 0
    >>> ch_2[0,0,1] = sin(x2)*cos(x2)
    >>> ch_2[0,1,1] = 0
    >>> ch_2[1,1,1] = 0
    >>> ch_2[1,0,1] = 0
    >>> ch_2[1,1,0] = 0
    >>> ch_2[1,0,0] = -sin(x2)*cos(x2)
    >>> ch_2[0,1,0] = -sin(x2)*cos(x2)

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valence indices (-1, -1):

    >>> g = Arraypy((2, 2)).to_tensor((-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    The covariant divergence of a tensor field:
    >>> delta_T = delta(T, g, ch_2, var)
    >>> print(delta_T)
    x1*sin(x2)*cos(x2) + 1  0
    """
    # Handling of a input argument - var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
        g = g.to_matrix()
    if not g.is_symmetric():
        raise ValueError('The metric tensor must be symmetric.')

    # Handling of a input argument - Christoffel symbol of second kind
    check_the_christoffel_symbols_2(ch_2)

    idx_ch = ch_2.start_index[0]

    # Handling of a input argument - tensor field T
    if not isinstance(T, TensorArray):
        raise TypeError(
            'The type of vector field must be TensorArray')

    idx_start_T = T.start_index[0]

    # The definition of inverse matrix of the metric tensor
    g_inv = g.inv()

    # The definition of diapason changes in an index
    # The dimension of the input array
    n = T.shape[0]

    # The rank of the input array
    rank_T = len(T.shape)

    index_T = T.index_list

    idx_char_delta_T = [(-1) for i in range(rank_T - 1)]

    nabla_T = nabla(T, ch_2, var)

    # Creating the output array in accordance with the start index
    delta_T = Arraypy([rank_T - 1, n, idx_start_T]).to_tensor(idx_char_delta_T)

    # Calculation
    for index in index_T:
        k = idx_start_T
        while k < n + idx_start_T:
            for j in range(n):
                idx_nabla_T = tuple(list(index) + [k])
                idx_delta_T = list(index)
                del idx_delta_T[0]
                idx_delta_T = tuple(idx_delta_T)
                delta_T[idx_delta_T] = (-1) * \
                    nabla_T[idx_nabla_T] * g_inv[k, j]
            k = k + 1

    # Output
    return delta_T


def riemann_li(C, g, var, type_output='t'):
    """Return the Riemann curvature tensor of type (1, -1, -1, -1)
    for the given left-invariant metric tensor.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import riemann_li
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy import symbols, cos, sin
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional arraypy
    or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    C it's a structural constant must be tensor with valence indices (1,-1,-1):

    >>> C = Arraypy([3, 2, 0]).to_tensor((1, -1, -1))
    >>> C[0,0,0] = 0
    >>> C[0,0,1] = sin(x2)*cos(x2)
    >>> C[0,1,1] = 0
    >>> C[1,1,1] = 0
    >>> C[1,0,1] = 0
    >>> C[1,1,0] = 0
    >>> C[1,0,0] = -sin(x2)*cos(x2)
    >>> C[0,1,0] = -sin(x2)*cos(x2)

    g it's a left-invariant metric tensor must be symmetric matrix, arraypy or
    tensor with valence indices (-1, -1):

    >>> g = Arraypy((2, 2)).to_tensor((-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The curvature tensor:
    >>> r_li = riemann_li(C, g, var, 'a')
    >>> print(r_li)
    -0.25*sin(x2)**2*cos(x2)**2  0
    0  0
    0  0
    0  0
    0  0
    0  0
    0  0
    0  0
    """
    # Handling of input vector arguments var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    check_metric_tensor(g)

    if isinstance(g, (Arraypy, TensorArray)):
        if not (g.start_index[0] == g.start_index[1]):
            raise ValueError(
                'The starting indices of metric tensor must be identical')
        idx_g = g.start_index[0]
        g_inv = (g.to_matrix()).inv()
    elif isinstance(g, Matrix):
        if not g.is_symmetric():
            raise ValueError('The metric tensor must be symmetric.')
        idx_g = 0
        g_inv = g.inv()

    # Handling of a input argument - structure constant
    if not isinstance(C, TensorArray):
        raise TypeError(
            'The type of must be TensorArray')
    else:
        if isinstance(C, TensorArray):
            if not C.type_pq == (1, 2):
                raise ValueError(
                    'The valence or ind_char of must be (1,-1,-1)')
            idx_c = C.start_index[0]
    # The definition of diapason changes in an index
    [n1, n2] = g.shape
    if not n == n1:
        raise ValueError(
            'The rank of the metric tensor does not coincide with the number of variables.')

    if (idx_g != idx_c):
        raise ValueError(
            'The start index of the tensor field and Christoffel symbol \
            of second kind must be equal')
    else:
        idx_start = idx_g

    indices = range(idx_start, idx_start + n)

    gamma = Arraypy([3, n, idx_start])
    for p in indices:
        for i in indices:
            for j in indices:
                for s in indices:
                    for k in indices:
                        gamma[p, i, j] = 0.5 * (C[p, i, j] + g[j, s] * C[s, k, i] * g_inv[
                                                k, p] + g[i, s] * C[s, k, j] * g_inv[k, p])
    # Creating the output array in accordance with the start index
    R = Arraypy([4, n, idx_start])

    # Calculation
    for s in indices:
        for i in indices:
            for j in indices:
                for k in indices:
                    for p in indices:
                        R[i, j, k, s] = gamma[s, i, p] * gamma[p, j, k] - gamma[s, j, p] * gamma[p, i, k] - \
                            gamma[s, p, k] * gamma[p, i, j]

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        riemann = R.to_tensor((1, -1, -1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        riemann = R
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - TensorArray.")

    # Output
    return riemann


def k_sigma_li(R, g, var):
    """Return Sectional curvature  in the direction of coordinate areas.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import k_sigma_li, riemann_li
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy import symbols, cos, sin
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional
    arraypy or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    g it's a metric tensor must be symmetric matrix, arraypy or tensor
    with valence indices (-1, -1):

    >>> g = Arraypy((2, 2)).to_tensor((-1, -1))
    >>> g[0,0] = cos(x2)**2
    >>> g[0,1] = 0
    >>> g[1,0] = 0
    >>> g[1,1] = 1

    C it's a structural constant must be tensor with valence indices (1,-1,-1):
    
    >>> C = Arraypy([3, 2, 0]).to_tensor((1, -1, -1))
    >>> C[0,0,0] = 0
    >>> C[0,0,1] = sin(x2)
    >>> C[0,1,1] = cos(x2)
    >>> C[1,1,1] = cos(x2)
    >>> C[1,0,1] = cos(x2)
    >>> C[1,1,0] = 0
    >>> C[1,0,0] = -sin(x2)
    >>> C[0,1,0] = -sin(x2)

    R it's a Riemann curvature tensor must be symmetric matrix, arraypy or tensor
    with valences indices (1, -1, -1, -1):
    
    >>> R = riemann_li(C, g, var, 't')

    The sectional curvature:
    >>> k_sig_li = k_sigma_li(R, g, var)
    >>> print(k_sig_li)
    Division by zero!
    """
    # Handling of input vector arguments var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Definition of number of variables
    n = len(var)

    # Handling of a input argument - metric tensor g
    if isinstance(g, (Arraypy, TensorArray)):
        if not (g.start_index[0] == g.start_index[1]):
            raise ValueError(
                'The starting indices of metric tensor must be identical')
        idx_start = g.start_index[0]
    elif isinstance(g, Matrix):
        if not g.is_symmetric():
            raise ValueError('The metric tensor must be symmetric.')
        idx_start = 0

    # Handling of a input argument Riemann curvature tensor - R
    if not isinstance(R, (Matrix, Arraypy, TensorArray)):
        raise TypeError(
            'The type of Riemann curvature tensor must be Matrix, Arraypy or TensorArray')
    else:
        if isinstance(R, (Arraypy, TensorArray)):
            if isinstance(R, TensorArray):
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

    indices = range(n)
    k_sig_li = Arraypy([2, n, idx_start])

    # Calculation
    for i in indices:
        for j in indices:
            for k in indices:
                if (g[i, j] * g[j, j] - g[i, j]**2) == 0:
                    raise ValueError('Division by zero!')
                else:
                    k_sig_li = sum(
                        (g[k, i] * R[k, i, j, j]) / (g[i, i] * g[j, j] - g[i, j]**2))

    # Output
    return k_sig_li


def kulkarni_nomizu(h, k, var, type_output='t'):
    """Return the product of Kulkarni-Nomizu of type (-1, -1, -1, -1)
    for the given two symmetric tensor.

    Examples:
    =========

    >>> from sympy.tensor.riemannian_geometry import kulkarni_nomizu
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy import symbols, cos
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional arraypy
    or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    h,k it's a tensor must be symmetric arraypy or tensor
    with valence indices (-1, -1):

    >>> h = Arraypy((2, 2)).to_tensor((-1, -1))
    >>> h[0,0] = x1
    >>> h[0,1] = 0
    >>> h[1,0] = 0
    >>> h[1,1] = x2
    >>> k = Arraypy((2, 2)).to_tensor((-1, -1))
    >>> k[0,0] = x2
    >>> k[0,1] = 0
    >>> k[1,0] = 0
    >>> k[1,1] = x1

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The curvature tensor:
    >>> k_n = kulkarni_nomizu(h, k, var, 'a')
    >>> print(k_n)
    0  0
    0  0
    0  x1**2 + x2**2
    -x1**2 - x2**2  0
    0  -x1**2 - x2**2
    x1**2 + x2**2  0
    0  0
    0  0
    """
    # Handling of input vector arguments var
    check_vector_of_arguments(var)

    if isinstance(var, (TensorArray, Arraypy)):
        var = var.to_list()

    # Handling of input symmetric tensor h
    if not isinstance(h, TensorArray):
        raise TypeError(
            'The type of input tensor must be a TensorArray')
    if isinstance(h, TensorArray):
        if not h.type_pq == (0, 2):
            raise ValueError(
                'The valence or ind_char of tensor must be (-1,-1)')
        if not (h.to_matrix()).is_symmetric():
            raise ValueError('The tensor must be symmetric.')

    # Handling of input symmetric tensor k
    if not isinstance(k, TensorArray):
        raise TypeError(
            'The type of input tensor must be a TensorArray')
    if isinstance(k, TensorArray):
        if not k.type_pq == (0, 2):
            raise ValueError(
                'The valence or ind_char of tensor must be (-1,-1)')
        if not (k.to_matrix()).is_symmetric():
            raise ValueError('The tensor must be symmetric.')

    if (h.start_index[0] != k.start_index[0]):
        raise ValueError(
            'The start index of the tensors must be equal')
    else:
        idx_start = h.start_index[0]

    # Definition of number of variables
    n = len(var)
    kul_nom = Arraypy([4, n, idx_start])
    indices = range(idx_start, idx_start + n)

    # Calculation
    for i in indices:
        for j in indices:
            for t in indices:
                for l in indices:
                    kul_nom[i, j, t, l] = (
                        h[i, t] * k[j, l] - h[i, l] * k[j, t]) - (h[j, t] * k[i, l] - h[j, l] * k[i, t])

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        K = kul_nom.to_tensor((-1, -1, -1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        K = kul_nom
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 't' and None - TensorArray.")
    # Output
    return K


def second_surf(surf, var, type_output='t'):
    """Return the second quadratic form.

    Examples:
    =========

    >>> from sympy import symbols
    >>> from sympy.tensor.riemannian_geometry import second_surf
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional arraypy
    or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    surf it's list of functions, must be consist of one or three functions.

    type_output it's optional parameter function, indicating the type of calculation
    result and receiving the character or string value:
    - symbol 't' means that the type of the result will match TensorArray;
    - symbol 'a' means that the type of the result will be Arraypy;
    - default function takes a parameter 't', so that the result will be a TensorArray.

    The the second quadratic form.
    >>> surf3 = [x1+x2, 2*x1**2-3*x2, (1+x2)*x1+x2-4]
    >>> print(second_surf(surf3, var, 't'))
    (-x1 + x2)/(3*x1)  -(4*x1 + 3)/((x1 + 1)*(x2 + 1))
    -(4*x1 + 3)/((x1 + 1)*(x2 + 1))  0
    >>> surf1 = [x1 + 4*x2**2]
    >>> print(second_surf(surf1, var, 't'))
    0  0
    0  8
    """
    # The definition symbols i, j, k
    i = Symbol('i')
    j = Symbol('j')
    k = Symbol('k')

    b = Arraypy((2, 2))

    # Calculation
    if (len(surf) == 1):
        b[0, 0] = diff(diff(surf[0], var[0]), var[0])
        b[0, 1] = b[1, 0] = diff((diff(surf[0], var[0])), var[1])
        b[1, 1] = diff((diff(surf[0], var[1])), var[1])
    elif (len(surf) == 3):
        # The first partial derivatives
        r_u = diff(surf[0], var[0]) * i + diff(surf[1], var[0]) * j +\
            diff(surf[2], var[0]) * k
        r_v = diff(surf[0], var[1]) * i + diff(surf[1], var[1]) * j +\
            diff(surf[2], var[1]) * k

        # The vector product
        vect_prod = (r_u.coeff(j) * r_v.coeff(k) - r_v.coeff(j) * r_u.coeff(k)) * i - \
            (r_u.coeff(k) * r_v.coeff(i) - r_v.coeff(k) * r_u.coeff(i)) * j + \
            (r_u.coeff(i) * r_v.coeff(j) - r_v.coeff(i) * r_u.coeff(j)) * k

        # The length of vector product
        len_r_uv = r_u.coeff(i) * r_v.coeff(i) * i + r_u.coeff(j) * r_v.coeff(j) * j + \
            r_u.coeff(k) * r_v.coeff(k) * k

        if (len_r_uv == 0):
            raise ValueError('The two-dimensional area is a degenerate!')

        # The components of the normal vector
        n = (simplify(vect_prod.coeff(i) / len_r_uv.coeff(i)) * i +
             simplify(vect_prod.coeff(j) / len_r_uv.coeff(j)) * j +
             simplify(vect_prod.coeff(k) / len_r_uv.coeff(k)) * k)

        # The second partial derivatives
        r_uu = diff(r_u.coeff(i), var[0]) * i + diff(r_u.coeff(j), var[0]) * j + \
            diff(r_u.coeff(k), var[0]) * k
        r_uv = diff(r_u.coeff(i), var[1]) * i + diff(r_u.coeff(j), var[1]) * j + \
            diff(r_u.coeff(k), var[1]) * k
        r_vv = diff(r_v.coeff(i), var[1]) * i + diff(r_v.coeff(j), var[1]) * j + \
            diff(r_v.coeff(k), var[1]) * k

        b[0, 0] = r_uu.coeff(i) * n.coeff(i) + r_uu.coeff(j) * n.coeff(j) + \
            r_uu.coeff(k) * n.coeff(k)
        b[0, 1] = b[1, 0] = r_uv.coeff(i) * n.coeff(i) + r_uv.coeff(j) * n.coeff(j) + \
            r_uv.coeff(k) * n.coeff(k)
        b[1, 1] = r_vv.coeff(i) * n.coeff(i) + r_vv.coeff(j) * n.coeff(j) + \
            r_vv.coeff(k) * n.coeff(k)
    else:
        raise ValueError(
            "The argument surf must be consist one function or three functions")

    # Handling of an output array
    if type_output == str('t') or type_output == Symbol('t'):
        b = b.to_tensor((-1, -1))
    elif type_output == str('a') or type_output == Symbol('a'):
        b = b
    elif type_output == str('m') or type_output == Symbol('m'):
        b = b.to_matrix()
    else:
        raise ValueError(
            "The parameter of type output result must 'a' - Arraypy or 'm' - Matrix\
            't' and None - TensorArray.")

    # Output
    return b


def k_surf(surf, var):
    """Return the Gaussian curvature.

    Examples:
    =========

    >>> from sympy import symbols
    >>> from sympy.tensor.riemannian_geometry import k_surf
    >>> x1, x2 = symbols('x1, x2')

    var it's a list of symbolic arguments. May be a list, one-dimensional arraypy
    or one-dimensional tensor with valence of indices (+1):

    >>> var = [x1, x2]

    surf it's list of functions, must be consist of one or three functions.

    The Gaussian curvature:
    >>> surf3 = [x1+x2, 2*x1**2-3*x2, (1+x2)*x1+x2-4]
    >>> print(k_surf(surf3, var))
    -(4*x1 + 3)**2/((x1 + 1)**2*(x2 + 1)**2*(((x1 + 1)**2 + 10)* \
    (16*x1**2 + (x2 + 1)**2 + 1) - (-12*x1 + (x1 + 1)*(x2 + 1) + 1)**2))
    >>> surf1 = [x1 + 4*x2**2]
    >>> print(k_surf(surf1, var))
    0
    """
    # Calculation
    if (len(surf) == 1):
        K = diff(diff(surf[0], var[0]), var[0]) * diff(diff(surf[0], var[1]), var[1]) -\
            (diff(diff(surf[0], var[0]), var[1]))**2 / \
            (1 + diff(surf[0], var[0])**2 + diff(surf[0], var[1])**2)**2
    elif (len(surf) == 3):
        g = Arraypy((2, 2))
        g[0, 0] = diff(surf[0], var[0])**2 + \
            diff(surf[1], var[0])**2 + diff(surf[2], var[0])**2
        g[0, 1] = g[1, 0] = diff(surf[0], var[0]) * diff(surf[0], var[1]) + diff(
            surf[1], var[0]) * diff(surf[1], var[1]) + diff(surf[2], var[0]) * diff(surf[2], var[1])
        g[1, 1] = diff(surf[0], var[1])**2 + \
            diff(surf[1], var[1])**2 + diff(surf[2], var[1])**2

        b = second_surf(surf3, var, 't')

        K = simplify(
            (b[0, 0] * b[1, 1] - b[0, 1]**2) / (g[0, 0] * g[1, 1] - g[0, 1]**2))

    else:
        raise ValueError(
            "The argument surf must be consist one function or three functions")
    # Output
    return K
