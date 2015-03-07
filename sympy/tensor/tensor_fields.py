# -*- coding: utf-8 -*-

from sympy.matrices import Matrix, eye
from sympy.core import diff, Add, Symbol
from sympy.simplify import simplify
from sympy import sqrt
from sympy.tensor.arraypy import Arraypy, Tensor, matrix2arraypy, \
    matrix2tensor, list2arraypy, list2tensor

"""Module tensor_fields contains functions for working with the tensor fields:
-the calculation of the differential and the gradient of the function;
-curl and divergence of a vector field;
-the calculation of the derivative of Li and the external differentiation
of differential forms.

Functions are work with the multidimensional arrays arraypy and tensors,
classes and methods which are contained in the module arraypy.

"""


# ---------------- df --------------------------------


def df(f, args, output_type='l'):
    """Return an the 1-form df, differential of function f(x).

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import df
    >>> from sympy import symbols, sin
    >>> x1, x2, x3= symbols('x1 x2 x3')
    >>> f=x1**2*x2 + sin(x2*x3 - x2)

    args it’s a list of symbol arguments of the function f. It can be in list,
    array of arraypy or contravariant tensor.
    Indexing the returned object will not begin from scratch and depends on
    the initial index of the argument array args:

    >>> args_t=Arraypy([1,3,1]).to_tensor(1)
    >>> args_t[1]=x1
    >>> args_t[2]=x2
    >>> args_t[3]=x3

    output _ type  it is an optional parameter accepting  symbol value of
    'l', 'a' or  't' and indicative on the type of result of calculations:
    'l' it is  a result as a list(list);
    'a' it is a result as an unidimensional array of arraypy;
    't' it is a result as an unidimensional covariant tensor.

    Differentials:

    >>> d = df(f, args_t, 't')
    >>> d
    2*x1*x2 x1**2 + (x3 - 1)*cos(x2*x3 - x2) x2*cos(x2*x3 - x2)

    The valence of the returned tensor:

    >>> d.type_pq
    (0, 1)

    """
    # Handling of a vector of arguments
    if not isinstance(args, (list, Tensor, Arraypy)):
        raise TypeError(
            "The type of vector of arguments must be list, Tensor or Arraypy")
    if isinstance(args, (Tensor, Arraypy)):
        if len(args.shape) != 1:
            raise ValueError("The dimension of argument must be 1")
        if isinstance(args, Tensor):
            if args.type_pq != (1, 0):
                raise ValueError(
                    "The valency(ind_char) of tensor must be (+1)")
        idx = args.start_index[0]
    if isinstance(args, list):
        idx = 0

    # Creating the output array in accordance with start indexes
    n = len(args)
    array = Arraypy([1, n, idx])
    indices = range(idx, idx + n)

    # Calculation
    for k in indices:
        array[k] = diff(f, args[k])

    # Handling of an output array
    if output_type == 't' or output_type == Symbol('t'):
        differential = Arraypy.to_tensor(array, -1)
    elif output_type == 'a' or output_type == Symbol('a'):
        differential = array
    elif output_type == 'l' or output_type == Symbol('l'):
        differential = Arraypy.to_list(array)
    else:
        raise TypeError(
            "The third arguments must be 't'-tensor,'a'-massiv Arraypy, \
            'l'-list")
    # Output
    return differential

# ---------------- grad --------------------------------


def grad(f, args, g=None, output_type=None):
    """Return the vector field Gradient(f(x)) of a function f(x).

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import grad
    >>> from sympy import symbols, sin
    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> f=x1**2*x2 + sin(x2*x3 - x2)

    args it’s a list of symbol arguments of function of f.
    It can be in list, array of arraypy or contravariant tensor.
    Indexing of the returned object may be not begining from scratch and
    depends on the initial index of the metric tensor:

    >>> args=[x1,x2,x3]

    g - optional parameter, metric tensor, which can be a matrix "Matrix",
    array of arraypy or covariant tensor:

    >>> g=Arraypy([2,3,1])
    >>> g_t=g.to_tensor((-1,-1))
    >>> g_t[1,1]=2
    >>> g_t[1,2]=1
    >>> g_t[1,3]=0
    >>> g_t[2,1]=1
    >>> g_t[2,2]=3
    >>> g_t[2,3]=0
    >>> g_t[3,1]=0
    >>> g_t[3,2]=0
    >>> g_t[3,3]=1

    output _ type  it is an optional parameter accepting  symbol value of
    'l', 'a' or  't' and indicative on the type of result of calculations:
    'l' it is  a result as a list(list);
    'a' it is a result as an unidimensional array of arraypy;
    't' it is a result as an unidimensional covariant tensor.

    Gradient:
    >>> gr=grad(f,args,g_t,'a')
    >>> gr
    -x1**2/5 + 6*x1*x2/5 - (x3 - 1)*cos(x2*x3 - x2)/5 2*x1**2/5 - 2*x1*x2/5 + \
    2*(x3 - 1)*cos(x2*x3 - x2)/5 x2*cos(x2*x3 - x2)

    """
    # Handling of a vector of arguments
    if not isinstance(args, (list, Tensor, Arraypy)):
        raise TypeError(
            "The type of vector of arguments must be list, Tensor or Arraypy")
    if isinstance(args, (Tensor, Arraypy)):
        if len(args.shape) != 1:
            raise ValueError("The dimension of argument must be 1")
        if isinstance(args, Tensor):
            if args.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_args = args.start_index[0]
    if isinstance(args, list):
        idx_args = 0

    # Handling of the metric tensor
    # 1. if g is not NULL
    if g is not None:
        if output_type is None:
            output_type = 't'
        if not isinstance(g, (Tensor, Matrix, Arraypy)):
            raise ValueError("Type must be Matrix or Tensor or Arraypy")
        if isinstance(g, Tensor):
            if g.type_pq != (0, 2):
                raise ValueError("The indices of tensor must be (-1,-1)")

        # The definition of the start index
        if isinstance(g, Matrix):
            idx_st = 0
        else:
            idx_st = g.start_index[0]
        if isinstance(g, type(args)) and idx_st != idx_args:
            raise ValueError(
                "The start index of the metric tensor and vector of arguments \
                must be equal")

        if isinstance(g, (Tensor, Arraypy)):
            g = g.to_matrix()
        if not g.is_symmetric():
            raise ValueError("The metric is not symmetric")
    # 2.if g is NULL
    else:
        # g - the identity matrix
        g = eye(len(args))
        idx_st = 0

    # Creating the output array in accordance with start indexes
    n = len(args)
    array = Arraypy([1, n, idx_st])
    indices = range(idx_st, idx_st + n)

    # Calculating
    g_inv = g.inv()
    if isinstance(args, (Tensor, Arraypy)):
        args = args.to_list()
    for i in indices:
        for j in indices:
            array[i] += (g_inv[i - idx_st, j - idx_st] * diff(f, args[j -
                                                                      idx_st]))

    # Handling of an output array
    if output_type == 't' or output_type == Symbol('t'):
        gradient = Arraypy.to_tensor(array, 1)
    elif output_type == 'a' or output_type == Symbol('a'):
        gradient = array
    elif output_type == 'l' or output_type == Symbol('l') or output_type is\
            None:
        gradient = Arraypy.to_list(array)
    else:
        raise TypeError(
            "The third arguments must be 't'-tensor,'a'-massiv Arraypy, \
            'l'-list")
# Output
    return gradient

# ---------------- rot --------------------------------


def rot(X, args, output_type=None):
    """Return the vorticity vector field rot(X) of a vector field X in R^3
    (curl, rotation, rotor, vorticity).
    A rotor can be calculated for only in three-dimensional Euclidean space.

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import rot
    >>> from sympy import symbols
    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> X=Arraypy(3)
    >>> X_t=Tensor(X,(1))
    >>> X_t[0]=x1*x2**3
    >>> X_t[1]=x2-cos(x3)
    >>> X_t[2]=x3**3-x1
    >>> arg=[x1,x2,x3]
    >>> r=rot(X_t,arg,'t')
    >>> r
    -sin(x3) 1 -3*x1*x2**2

    """
    # Handling of a vector of arguments
    if not isinstance(args, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of arguments vector must be list, Tensor or Arraypy")
    if len(args) != 3:
        raise ValueError("ERROW:three variables are required")
    if isinstance(args, (Tensor, Arraypy)):
        if len(args.shape) != 1:
            raise ValueError("The lenght of argument must be 1")
        if isinstance(args, Tensor):
            if args.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_args = args.start_index[0]
    if isinstance(args, list):
        idx_args = 0

    # Handling of a vector field
    if not isinstance(X, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of vector fields must be list, Tensor or Arraypy")
    if len(X) != 3:
        raise ValueError("ERROW:a three-dimensional vector is necessary")

    if isinstance(X, (Tensor, Arraypy)):
        if len(X.shape) != 1:
            raise ValueError("The dim of argument must be 1")
        if isinstance(X, Tensor):
            out_t = 't'
            if X.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_X = X.start_index[0]
    elif isinstance(X, list):
        idx_X = 0
        out_t = 'l'

    if output_type is None:
        if out_t is not None:
            output_type = out_t
        else:
            output_type = 'a'

    # The definition of the start index
    if isinstance(X, type(args)) and (idx_X != idx_args):
        raise ValueError(
            "The start index of vector field and vector of arguments must be \
            equal")
    idx_st = idx_X

    # Creating the output array in accordance with start indexes
    array = Arraypy([1, 3, idx_st])

    # Calculation
    if isinstance(X, (Tensor, Arraypy)):
        X = X.to_list()
    if isinstance(args, (Tensor, Arraypy)):
        args = args.to_list()

    array[idx_st] = (diff(X[2], args[1]) - diff(X[1], args[2]))
    array[idx_st + 1] = diff(X[0], args[2]) - diff(X[2], args[0])
    array[idx_st + 2] = diff(X[1], args[0]) - diff(X[0], args[1])

    # Handling of an output array
    if output_type == 't' or output_type == Symbol('t'):
        rotor = Arraypy.to_tensor(array, 1)
    elif output_type == 'a' or output_type == Symbol('a'):
        rotor = array
    elif output_type == 'l' or output_type == Symbol('l'):
        rotor = Arraypy.to_list(array)
    else:
        raise TypeError(
            "The third arguments must be 't'-tensor,'a'-massiv Arraypy, \
            'l'-list")
    # Output
    return rotor


# ---------------- div --------------------------------


def div(X, args, g=None):
    """Return the divergence of a vector field X. Compute divergence of vector
    field consisting of N elements.

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import div
    >>> from sympy import symbols, cos
    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> X = [x1*x2**3,x2-cos(x3),x3**3-x1]
    >>> g = Matrix([[2,1,0],[1,3,0],[0,0,1]])
    >>> arg = [x1, x2, x3]
    >>> dv = div(X,arg,g)
    x2**3 + 3*x3**2 + 1

    """
    # Handling of a vector of arguments
    if not isinstance(args, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of arguments vector must be list, Tensor or Arraypy")
    if isinstance(args, (Tensor, Arraypy)):
        if len(args.shape) != 1:
            raise ValueError("The lenght of argument must be 1")
        if isinstance(args, Tensor):
            if args.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        args = args.to_list()

    # Handling of a vector field
    if not isinstance(X, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of vector fields must be list, Tensor or Arraypy")
    if isinstance(X, (Tensor, Arraypy)):
        if len(X.shape) != 1:
            raise ValueError("The dim of argument must be 1")
        if isinstance(X, Tensor):
            if X.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        X = X.to_list()

    # Handling of the metric tensor
    if g is not None:
        if not isinstance(g, (Tensor, Matrix, Arraypy)):
            raise ValueError("Type must be Matrix or Tensor or Arraypy")
        else:
            if isinstance(g, (Tensor, Arraypy)):
                if isinstance(g, Tensor):
                    if g.type_pq != (0, 2):
                        raise ValueError(
                            "The indices of tensor must be (-1,-1)")
                g = g.to_matrix()
            if not g.is_symmetric():
                raise ValueError("The metric is not symmetric")
    else:
        g = eye(len(args))

# Calculation
    sq = sqrt(abs(Matrix.det(g)))
    divergenc = 0
    for k in range(len(args)):
        divergenc += simplify(1 / sq * sum([diff(X[k] * sq, args[k])]))
# Output
    return divergenc


# ------------------ lie_xy -------------------------------

def lie_xy(X, Y, args, output_type=None):
    """Return the vector field [X,Y], Lie bracket (commutator) of a vector
    fields X and Y.

    Examples:
    =========
    >>> from sympy.tensor.tensor_fields import lie_xy
    >>> from sympy import symbols, sin
    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> X=[x1*x2**3,x2-cos(x3),x3**3-x1]
    >>> Y = [x1**3*x2**3, x2*x3 - sin(x1*x3), x3**3 - x1**2]
    >>> arg = [x1, x2, x3]
    >>> lie = lie_xy(X, Y, arg,'a')
    >>> lie
    2*x1**3*x2**6 + 3*x1**3*x2**2*(x2 - cos(x3)) - 3*x1*x2**2*(x2*x3 - \
    sin(x1*x3))
    -x1*x2**3*x3*cos(x1*x3) - x2*x3 + x3*(x2 - cos(x3)) + \
    (-x1 + x3**3)*(-x1*cos(x1*x3) + x2) - (-x1**2 + x3**3)*sin(x3) + sin(x1*x3)
    x1**3*x2**3 - 2*x1**2*x2**3 + 3*x3**2*(-x1 + x3**3) - 3*x3**2*(-x1**2 + x3**3)

    """
    # Handling of a vector of arguments
    if not isinstance(args, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of arguments vector must be list, Tensor or Arraypy")
    if isinstance(args, (Tensor, Arraypy)):
        if len(args.shape) != 1:
            raise ValueError("The lenght of argument must be 1")
        if isinstance(args, Tensor):
            if args.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_args = args.start_index[0]
    if isinstance(args, list):
        idx_args = 0

    # Handling of the first vector field
    if not isinstance(X, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of vector fields must be list, Tensor or Arraypy")
    if isinstance(X, (Tensor, Arraypy)):
        if len(X.shape) != 1:
            raise ValueError("The dim of argument must be 1")
        if isinstance(X, Tensor):
            out_t = 't'
            if X.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_X = X.start_index[0]
    if isinstance(X, list):
        idx_X = 0
        out_t = 'l'

    # Handling of the second vector field
    if not isinstance(Y, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of vector fields must be list, Tensor or Arraypy")
    if isinstance(Y, (Tensor, Arraypy)):
        if len(Y.shape) != 1:
            raise ValueError("The dim of argument must be 1")
        if isinstance(Y, Tensor):
            if Y.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_Y = Y.start_index[0]
    if isinstance(Y, list):
        idx_Y = 0

    if len(Y) != len(X):
        raise ValueError(
            "The different number of arguments of the vector fields")
    elif len(args) != len(X) or len(args) != len(Y):
        raise ValueError(
            "The different number of components at the vector field and \
            vector of variables")

    # Define the start index in the output tensor
    if type(Y) == type(X) == type(args):
        if idx_Y != idx_X or idx_Y != idx_args or idx_X != idx_args:
            raise ValueError(
                "The start index of vector fields and vetcor of argements \
                must be equal")
    if idx_Y != idx_X:
        raise ValueError(
            "The start index of tensor and the vector field must be equal")
    idx_st = idx_Y

    if output_type is None:
        if out_t is not None:
            output_type = out_t
        else:
            output_type = 'a'

    # Creating the output array in accordance with start indexes
    Li = Arraypy([1, len(X), idx_st])

    # Calculating
    if isinstance(Y, (Tensor, Arraypy)):
        Y = Y.to_list()
    if isinstance(X, (Tensor, Arraypy)):
        X = X.to_list()
    if isinstance(args, (Tensor, Arraypy)):
        args = args.to_list()

    if X == Y:
        return 0
    else:
        indices = range(len(args))
        for i in indices:
            for k in indices:
                Li[i +
                   idx_st] += Add(diff(Y[i], args[k]) *
                                  X[k] -
                                  diff(X[i], args[k]) *
                                  Y[k])
    # Handling of an output array
    if output_type == 't' or output_type == Symbol('t'):
        Lie = Arraypy.to_tensor(Li, 1)
    elif output_type == 'a' or output_type == Symbol('a'):
        Lie = Li
    elif output_type == 'l' or output_type == Symbol('l'):
        Lie = Arraypy.to_list(Li)
    else:
        raise TypeError(
            "The third arguments must be 't'-tensor,'a'-massiv Arraypy, \
            'l'-list")

    # Output
    return Lie


# ---------------- NotNeedElement --------------------------------

def NotNeedElement(_list, index):
    """The function returns a tuple containing the remainder of the input list
    after you remove the element at the specified index."""
    res = []
    for i in range(len(_list)):
        if i != index:
            res.append((_list[i]))
    return (tuple(res))


# ---------------- dw --------------------------------

def dw(omega, args):
    """Return a skew-symmetric tensor of type (0,p+1).
    Indexes the output tensor will start as well as the input tensor (array).
    If the input parameters of the same type, they must be equal to the initial
    indexes.

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import dw
    >>> from sympy import symbols
    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> omega=Arraypy([2,3,1]).to_tensor((-1,-1))
    >>> omega[1,2]=x3
    >>> omega[1,3]=-x2
    >>> omega[2,1]=-x3
    >>> omega[2,3]=x1
    >>> omega[3,1]=x2
    >>> omega[3,2]=-x1
    >>> domega=dw(omega, [x1,x2,x3])
    >>> domega
    0 0 0
    0 0 3
    0 -3 0
    0 0 -3
    0 0 0
    3 0 0
    0 3 0
    -3 0 0
    0 0 0

    """
# Handling of a vector of arguments
    if not isinstance(args, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of arguments vector must be list, Tensor or Arraypy")
    if isinstance(args, (Tensor, Arraypy)):
        if len(args.shape) != 1:
            raise ValueError("The lenght of argument must be 1")
        if isinstance(args, Tensor):
            if args.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_args = args.start_index[0]
    if isinstance(args, list):
        idx_args = 0

# Handling of a differential form
    if not isinstance(omega, (Tensor, Arraypy)):
        raise ValueError("Type must be Tensor or Arraypy")
    idx_omega = omega.start_index[0]

    # Define the start index in the output tensor
    if isinstance(omega, type(args)) and idx_omega != idx_args:
        raise ValueError("Raznie indeksi!!!")
    idx_st = idx_omega

    # Creating the output array in accordance with start indexes
    n = omega.shape[0]  # the dimensionality of the input array
    p = len(omega.shape)  # the rank of the input array
    a = Arraypy([p + 1, n, idx_st])
    valence_ind = [(-1) for k in range(p + 1)]
    d_omega = a.to_tensor(valence_ind)

    # Calculation
    idx = d_omega.start_index
    if isinstance(args, (Tensor, Arraypy)):
        args = args.to_list()

    for i in range(len(d_omega)):
        # list of tuple. example:[(0, 1), (0, 1), (0, 0)]
        tuple_list_indx = [NotNeedElement(idx, f) for f in range(len(idx))]
        for k in range(p + 1):
            d_omega[idx] += Add(((-1)**k) * diff(omega[tuple_list_indx[k]],
                                                 args[idx[k] - idx_st]))
        idx = d_omega.next_index(idx)

# Output
    return d_omega

# ---------------- NeedElementK --------------------------------


def NeedElementK(_list, index, k):
    """The function replaces the item "index" on the element "k".

    The result is a tuple.

    """
    output = []
    for i in range(len(_list)):
        if i != index:
            output.append((_list[i]))
        else:
            output.append(k)
    return (tuple(output))


# ---------------- lie_w --------------------------------

def lie_w(omega, X, args):
    """Return a skew-symmetric tensor of type (0,p).
    Indexes the output tensor will start as well as the input tensor (array).
    If the tensor and vector field of the same type, they must be equal to the
    initial indexes.
    If all the input parameters of the same type, they must be equal to the
    initial indexes.

    Function lie_w calculates all the components of the Lie
    derivative differential forms in a symbolic form.

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import lie_w
    >>> from sympy import symbols, cos
    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> omega=Arraypy([2,3,1]).to_tensor((-1,-1))
    >>> omega[1,2]=x3
    >>> omega[1,3]=-x2
    >>> omega[2,1]=-x3
    >>> omega[2,3]=x1
    >>> omega[3,1]=x2
    >>> omega[3,2]=-x1
    >>> arg = [x1, x2, x3]
    >>> X = [x1*x2**3,x2-cos(x3),x3**3-x1]
    >>> li = lie_w(omega,X,arg)
    >>> li
    0 x2**3*x3 + x3**3 + x3 -x2**4 - 3*x2*x3**2 - x2 + x3*sin(x3) + cos(x3)
    -x2**3*x3 - x3**3 - x3 0 -2*x1*x2**3 + 3*x1*x3**2 + x1
    x2**4 + 3*x2*x3**2 + x2 - x3*sin(x3) - cos(x3) 2*x1*x2**3 - 3*x1*x3**2 - x1 0

    """
# Handling of a vector of arguments
    if not isinstance(args, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of arguments vector must be list, Tensor or Arraypy")
    if isinstance(args, (Tensor, Arraypy)):
        if len(args.shape) != 1:
            raise ValueError("The lenght of argument must be 1")
        if isinstance(args, Tensor):
            if args.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_args = args.start_index[0]
    if isinstance(args, list):
        idx_args = 0

# Handling of a vector field
    if not isinstance(X, (list, Tensor, Arraypy)):
        raise ValueError(
            "The type of vector fields must be list, Tensor or Arraypy")
    if isinstance(X, (Tensor, Arraypy)):
        if len(X.shape) != 1:
            raise ValueError("The dim of argument must be 1")
        if isinstance(X, Tensor):
            if X.type_pq != (1, 0):
                raise ValueError("The valency of tensor must be (+1)")
        idx_X = X.start_index[0]
    if isinstance(X, list):
        idx_X = 0

# Handling of a differential form
    if not isinstance(omega, (Tensor, Arraypy)):
        raise ValueError("Type must be Tensor or Arraypy")
    idx_omega = omega.start_index[0]


# Define the start index in the output tensor
    if type(omega) == type(X) == type(args):
        if idx_omega != idx_X or idx_omega != idx_args or idx_X != idx_args:
            raise ValueError(
                "The start index of tensor,vector field and \
                vetcor of argements  must be equal")
    if isinstance(omega, type(X)) and idx_omega != idx_X:
        raise ValueError(
            "The start index of tensor and vector field must be equal")
    idx_st = idx_omega

# Creating the output array in accordance with start indexes
    n = omega.shape[0]  # the dimensionality of the input array
    r = len(omega.shape)  # the rank of the input array
    a = Arraypy([r, n, idx_st])
    valence_list = [(-1) for k in range(r)]
    diff_Lie = a.to_tensor(valence_list)

    # Calculation
    idx = diff_Lie.start_index
    if isinstance(args, (Tensor, Arraypy)):
        args = args.to_list()
    if isinstance(X, (Tensor, Arraypy)):
        X = X.to_list()

    for p in range(len(diff_Lie)):
        for k in range(len(idx) + 1):
            tuple_list_indx = [
                NeedElementK(idx, f, k + idx_st) for f in range(len(idx))]
            diff_omega = diff(omega[idx], args[k]) * X[k]
            for j in range(len(idx)):
                diff_Lie[idx] += diff(X[k], args[idx[j] - idx_st]) *\
                    omega[tuple_list_indx[j]]
            diff_Lie[idx] = diff_Lie[idx] + diff_omega
        idx = diff_Lie.next_index(idx)
# Output
    return diff_Lie
