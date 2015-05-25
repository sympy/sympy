# -*- coding: utf-8 -*-

from sympy.matrices import Matrix, eye, det
from sympy.core import diff, Add, Symbol
from sympy.simplify import simplify
from sympy.tensor.arraypy import Arraypy, TensorArray, matrix2arraypy, \
    matrix2tensor, list2arraypy, list2tensor
from sympy.tensor.tensor_methods import is_asymmetric, perm_parity, \
    lower_index, raise_index, is_symmetric, tensor_product
from sympy import sqrt
from itertools import permutations
from sympy.functions.combinatorial.factorials import factorial
from sympy.tensor.helper_functions import check_vector_of_arguments, \
    check_metric_tensor, check_the_vector_field, sign_permutations, \
    delete_index_from_list, replace_index_to_k

"""Module tensor_fields contains functions for working with the tensor fields:
-the calculation of the differential and the gradient of the function;
-curl and divergence of a vector field;
-the calculation of the derivative of Li;
-the calculation the Lie brackets of two vector fields;
-the calculation the external differentiation of differential forms;
-


All functions can be operated with input arguments specified as a tensor.
Some of the input arguments of functions can be a list, matrix or array of
arraypy. The starting index of tensor may not be equal to 0. In such instance,
the object, which returns the function will have the starting index
corresponding to the conditions of the call.
Also lie_w function and dw involves working with multidimensional arrays.

Functions work with the multidimensional arrays arraypy and tensors,
classes and methods which are contained in the module arraypy.

"""


def df(f, args, output_type='l'):
    """Return an the 1-form df, differential of function f(x).

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import df
    >>> from sympy import symbols, sin
    >>> from sympy.tensor.arraypy import Arraypy
    >>> x1, x2, x3= symbols('x1 x2 x3')

    f it's a function the differential of that is calculated:

    >>> f=x1**2*x2 + sin(x2*x3 - x2)

    args it's a list of symbol arguments of the function f. It can be in list,
    array of arraypy or contravariant tensor:

    >>> args_t=Arraypy([1,3,1]).to_tensor(1)
    >>> args_t[1]=x1
    >>> args_t[2]=x2
    >>> args_t[3]=x3

    output_type  it is an optional parameter accepting  symbol value of
    'l', 'a' or  't' and indicative on the type of result of calculations:
    - 'l' it is  a result as a list(list);
    - 'a' it is a result as an unidimensional array of arraypy;
    - 't' it is a result as an unidimensional covariant tensor.

    Differential:

    >>> d = df(f, args_t, 't')
    >>> print(d)
    2*x1*x2 x1**2 + (x3 - 1)*cos(x2*x3 - x2) x2*cos(x2*x3 - x2)

    The valence of the returned tensor:

    >>> d.type_pq
    (0, 1)

    """
    # Handling of a vector of arguments
    check_vector_of_arguments(args)

    # The definition of the start index
    if isinstance(args, list):
        idx_start = 0
    else:
        idx_start = args.start_index[0]

    # Creating the output array in accordance with the start index
    n = len(args)
    array = Arraypy([1, n, idx_start])

    # Calculation
    for k in range(idx_start, idx_start + n):
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
            "The third argument must be 't' - TensorArray,'a' - Arraypy, \
            'l' - list")
    # Output
    return differential


def grad(f, args, g=None, output_type=None):
    """Return the vector field gradient of a function f(x).

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import grad
    >>> from sympy import symbols, sin
    >>> from sympy.tensor.arraypy import Arraypy
    >>> x1, x2, x3 = symbols('x1 x2 x3')

    f it's a function the differential of that is calculated:

    >>> f=x1**2*x2 + sin(x2*x3 - x2)

    args it's a list of symbol arguments of function of f.
    It can be in list, array of arraypy or contravariant tensor:

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
    - 'l' it is  a result as a list(list);
    - 'a' it is a result as an unidimensional array of arraypy;
    - 't' it is a result as an unidimensional covariant tensor.

    Gradient:
    >>> gr=grad(f,args,g_t,'a')
    >>> print(gr)
    -x1**2/5 + 6*x1*x2/5 - (x3 - 1)*cos(x2*x3 - x2)/5 2*x1**2/5 - 2*x1*x2/5 + \
    2*(x3 - 1)*cos(x2*x3 - x2)/5 x2*cos(x2*x3 - x2)

    """
    # Handling of a vector of arguments
    check_vector_of_arguments(args)

    # The definition of the start index
    if isinstance(args, list):
        idx_args = 0
    else:
        idx_args = args.start_index[0]

    # Handling of the metric tensor
    # 1. if g is not NULL
    if g is not None:
        if output_type is None:
            output_type = 't'
        check_metric_tensor(g)

        # The definition of the start index
        if isinstance(g, Matrix):
            idx_st = 0
        else:
            idx_st = g.start_index[0]
        # The start index is the same
        if isinstance(g, type(args)) and idx_st != idx_args:
            raise ValueError(
                "The start index of the metric tensor and vector of arguments \
                must be equal")
        if isinstance(g, (TensorArray, Arraypy)):
            g = g.to_matrix()
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
    if isinstance(args, (TensorArray, Arraypy)):
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
    elif output_type == 'l' or output_type == Symbol('l') or output_type is \
            None:
        gradient = Arraypy.to_list(array)
    else:
        raise TypeError(
            "The third argument must be 't' - tensor,'a' - Arraypy, \
            'l' - list")
# Output
    return gradient


def curl(X, args, output_type=None):
    """Return the rotor vector field curl(X) of a vector field X in R^3
    (curl, rotation, rotor, vorticity).
    A rotor can be calculated for only in three-dimensional Euclidean space.

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import curl
    >>> from sympy import symbols, cos
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> x1, x2, x3 = symbols('x1 x2 x3')

    X is a vector field, args it's a list of symbol arguments of a vector field
    X. They can be a list, array arraypy or contravariant tensor:

    >>> X=Arraypy(3)
    >>> X_t=TensorArray(X,(1))
    >>> X_t[0]=x1*x2**3
    >>> X_t[1]=x2-cos(x3)
    >>> X_t[2]=x3**3-x1
    >>> arg=[x1,x2,x3]
    >>> r=curl(X_t,arg,'t')
    >>> print(r)
    -sin(x3) 1 -3*x1*x2**2
    >>> r.type_pq
    (1, 0)

    """
    # Handling of a vector of arguments
    check_vector_of_arguments(args)
    if len(args) != 3:
        raise ValueError("Three variables are required")
    # The definition of the start index
    if isinstance(args, list):
        idx_args = 0
    else:
        idx_args = args.start_index[0]

    # Handling of a vector field
    check_the_vector_field(X)
    if len(X) != 3:
        raise ValueError("A three-dimensional vector field is necessary")
    # The definition of the start index and type of output
    if isinstance(X, (TensorArray, Arraypy)):
        if isinstance(X, TensorArray):
            out_t = 't'
        idx_X = X.start_index[0]
    else:
        idx_X = 0
        out_t = 'l'

    # The definition type of output of an output array
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

    # Creating the output array in accordance with the start index
    array = Arraypy([1, 3, idx_st])

    # Calculation
    if isinstance(X, (TensorArray, Arraypy)):
        X = X.to_list()
    if isinstance(args, (TensorArray, Arraypy)):
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
            "The third argument must be 't'-tensor,'a'-Arraypy, \
            'l'-list")
    # Output
    return rotor


def diverg(X, args, g=None):
    """Return the divergence of a vector field X. Compute divergence of vector
    field consisting of N elements.

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import diverg
    >>> from sympy import symbols, cos
    >>> from sympy.matrices import Matrix
    >>> x1, x2, x3 = symbols('x1 x2 x3')

    X is a vector field, args it's a list of symbol arguments of the vector
    field X. It's can be in list, array of arraypy or contravariant tensor:

    >>> X = [x1*x2**3,x2-cos(x3),x3**3-x1]
    >>> arg = [x1, x2, x3]

    g - optional parameter, metric tensor, which can be a matrix "Matrix",
    array of arraypy or covariant tensor:

    >>> g = Matrix([[2,1,0],[1,3,0],[0,0,1]])

    >>> dv = diverg(X,arg,g)
    >>> print(dv)
    x2**3 + 3*x3**2 + 1

    """
    # Handling of a vector of arguments
    check_vector_of_arguments(args)
    if isinstance(args, list):
        idx_args = 0
    else:
        idx_args = args.start_index[0]

    # Handling of the first vector field
    check_the_vector_field(X)
    if isinstance(X, (TensorArray, Arraypy)):
        idx_X = X.start_index[0]
    else:
        idx_X = 0

    if idx_args != idx_X:
        raise ValueError(
            "The start index of vector field and vector of arguments must be \
                    equal")

    # Handling of the metric tensor
    if g is not None:
        if isinstance(g, (TensorArray, Arraypy)):
            g = g.to_matrix()
    else:
        g = eye(len(args))

    # Calculation
    sq = sqrt(abs(det(g)))
    diver = 0
    for k in range(len(args)):
        diver += simplify(1 / sq *
                          sum([diff(X[k + idx_X] * sq, args[k + idx_X])]))
    # Output
    return diver


def lie_xy(X, Y, args, output_type=None):
    """Return the vector field [X,Y], Lie bracket (commutator) of a vector
    fields X and Y.

    Examples:
    =========
    >>> from sympy.tensor.tensor_fields import lie_xy
    >>> from sympy import symbols, cos, sin
    >>> x1, x2, x3 = symbols('x1 x2 x3')

    X, Y is a vector field, args it's a list of symbol arguments.
    It's can be a list, array arraypy or contravariant tensor:

    >>> X=[x1*x2**3,x2-cos(x3),x3**3-x1]
    >>> Y = [x1**3*x2**3, x2*x3 - sin(x1*x3), x3**3 - x1**2]
    >>> arg = [x1, x2, x3]

    The Lie brackets of two vector fields:

    >>> lie = lie_xy(X, Y, arg,'a')
    >>> print(lie)
    2*x1**3*x2**6 + 3*x1**3*x2**2*(x2 - cos(x3)) - 3*x1*x2**2*(x2*x3 - \
    sin(x1*x3))
    -x1*x2**3*x3*cos(x1*x3) - x2*x3 + x3*(x2 - cos(x3)) + \
    (-x1 + x3**3)*(-x1*cos(x1*x3) + x2) - (-x1**2 + x3**3)*sin(x3) + sin(x1*x3)
    x1**3*x2**3 - 2*x1**2*x2**3 + 3*x3**2*(-x1 + x3**3) - \
    3*x3**2*(-x1**2 + x3**3)

    """
    # Handling of a vector of arguments
    check_vector_of_arguments(args)
    # The definition of the start index args
    if isinstance(args, list):
        idx_args = 0
    else:
        idx_args = args.start_index[0]

    # Handling of the first vector field
    check_the_vector_field(X)
    # The definition of the start index X and type of output
    if isinstance(X, (TensorArray, Arraypy)):
        if isinstance(X, TensorArray):
            out_t = 't'
        idx_X = X.start_index[0]
    else:
        idx_X = 0
        out_t = 'l'

    # Handling of the second vector field
    check_the_vector_field(Y)
    # The definition of the start index Y
    if isinstance(Y, (TensorArray, Arraypy)):
        idx_Y = Y.start_index[0]
    else:
        idx_Y = 0

    if len(Y) != len(X):
        raise ValueError(
            "The different number of arguments in the vector fields")
    elif len(args) != len(X) or len(args) != len(Y):
        raise ValueError(
            "The different number of components in the vector field and \
            vector of arguments")

    # Define the start index in the output tensor
    if type(Y) == type(X) == type(args):
        if idx_Y != idx_X or idx_Y != idx_args or idx_X != idx_args:
            raise ValueError(
                "The start index of vector fields and vector of argements \
                must be equal")
    if idx_Y != idx_X:
        raise ValueError(
            "The start index of vector fields must be equal")
    idx_st = idx_Y

    if output_type is None:
        if out_t is not None:
            output_type = out_t
        else:
            output_type = 'a'

    # Creating the output array in accordance with start indexes
    Li = Arraypy([1, len(X), idx_st])

    # Calculating
    if isinstance(Y, (TensorArray, Arraypy)):
        Y = Y.to_list()
    if isinstance(X, (TensorArray, Arraypy)):
        X = X.to_list()
    if isinstance(args, (TensorArray, Arraypy)):
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
            "The third argument must be 't'-TensorArray,'a'-Arraypy, \
            'l'-list")

    # Output
    return Lie


def dw(omega, args):
    """Return a skew-symmetric tensor of type (0, p+1).
    Indexes the output tensor will start as well as the input tensor (array).
    If the input parameters of the same type, they must be equal to the initial
    indexes.

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import dw
    >>> from sympy import symbols
    >>> from sympy.tensor.arraypy import Arraypy
    >>> x1, x2, x3 = symbols('x1 x2 x3')

    omega - differential form, differential which is calculated.
    It's can be a skew-symmetric tensor of type (0,p) or an array arraypy:

    >>> omega=Arraypy([2,3,1]).to_tensor((-1,-1))
    >>> omega[1,2]=x3
    >>> omega[1,3]=-x2
    >>> omega[2,1]=-x3
    >>> omega[2,3]=x1
    >>> omega[3,1]=x2
    >>> omega[3,2]=-x1

    args it's a list of symbol arguments of differential form.
    It's can be in list, array of arraypy or contravariant tensor.

    External differential of a differential forms:

    >>> domega=dw(omega, [x1,x2,x3])
    >>> print(domega)
    0 0 0
    0 0 3
    0 -3 0
    0 0 -3
    0 0 0
    3 0 0
    0 3 0
    -3 0 0
    0 0 0
    >>> domega.type_pq
    (0, 3)

    """
    # Handling of a vector of arguments
    check_vector_of_arguments(args)
    # The definition of the start index of args
    if isinstance(args, list):
        idx_args = 0
    else:
        idx_args = args.start_index[0]

    # Handling of a differential form
    if not isinstance(omega, (TensorArray, Arraypy)):
        raise ValueError(
            "The type of differential form must be TensorArray or Arraypy")
    if omega.rank > 1:
        if not is_asymmetric(omega):
            raise ValueError("The differential form must be a skew-symmetric")
    idx_omega = omega.start_index[0]

    # Define the start index in the output tensor
    if isinstance(omega, type(args)) and idx_omega != idx_args:
        raise ValueError("The start index of the differential form and \
        vector of arguments must be equal")
    idx_st = idx_omega

    # Creating the output array in accordance with start indexes
    n = omega.shape[0]  # the dimensionality of the input array
    p = len(omega.shape)  # the rank of the input array
    valence_ind = [(-1) for k in range(p + 1)]
    d_omega = Arraypy([p + 1, n, idx_st]).to_tensor(valence_ind)

    # Calculation
    idx = d_omega.start_index
    if isinstance(args, (TensorArray, Arraypy)):
        args = args.to_list()

    for i in range(len(d_omega)):
        # tuple_list_indx it's list of tuple. Example:[(0, 1), (0, 1), (0, 0)]
        tuple_list_indx = [
            delete_index_from_list(
                idx,
                f) for f in range(
                len(idx))]
        for k in range(p + 1):
            d_omega[idx] += Add(((-1)**k) * diff(omega[tuple_list_indx[k]],
                                                 args[idx[k] - idx_st]))
        idx = d_omega.next_index(idx)

    # Output
    return d_omega


def lie_w(omega, X, args):
    """Return a skew-symmetric tensor of type (0,p).
    Function lie_w calculates all the components of the Lie derivative
    differential forms in a symbolic form.

    Indexes the output tensor will start as well as the input tensor (array)
    "omega". If all the input parameters of the same type, they must be equal
    to the initial indexes.

    Examples:
    =========

    >>> from sympy.tensor.tensor_fields import lie_w
    >>> from sympy import symbols, cos
    >>> from sympy.tensor.arraypy import Arraypy
    >>> x1, x2, x3 = symbols('x1 x2 x3')

    omega - skew-symmetric tensor. Can be a tensor of type (0,p) or an array
    arraypy:

    >>> omega=Arraypy([2,3,1]).to_tensor((-1,-1))
    >>> omega[1,2]=x3
    >>> omega[1,3]=-x2
    >>> omega[2,1]=-x3
    >>> omega[2,3]=x1
    >>> omega[3,1]=x2
    >>> omega[3,2]=-x1

    X - the vector field along which the derivative is calculated:

    >>> X = [x1*x2**3,x2-cos(x3),x3**3-x1]

    args it's a list of symbol arguments.
    It's can be in list, array of arraypy or contravariant tensor:

    >>> arg = [x1, x2, x3]

    Lie derivative of a differential form:

    >>> li = lie_w(omega,X,arg)
    >>> print(li)
    0 x2**3*x3 + x3**3 + x3 -x2**4 - 3*x2*x3**2 - x2 + x3*sin(x3) + cos(x3)
    -x2**3*x3 - x3**3 - x3 0 -2*x1*x2**3 + 3*x1*x3**2 + x1
    x2**4 + 3*x2*x3**2 + x2 - x3*sin(x3) - cos(x3) 2*x1*x2**3 - \
    3*x1*x3**2 - x1 0
    >>> li.type_pq
    (0, 2)

    """
    # Handling of a vector of arguments
    check_vector_of_arguments(args)
    if isinstance(args, list):
        idx_args = 0
    else:
        idx_args = args.start_index[0]

    # Handling of a vector field
    check_the_vector_field(X)
    if isinstance(X, (TensorArray, Arraypy)):
        idx_X = X.start_index[0]
    else:
        idx_X = 0

    # Handling of a differential form
    if not isinstance(omega, (TensorArray, Arraypy)):
        raise ValueError(
            "The type of differential form must be TensorArray or Arraypy")
    if not is_asymmetric(omega):
        raise ValueError("The differential form must be a skew-symmetric")
    idx_omega = omega.start_index[0]


# Define the start index in the output tensor
    if type(omega) == type(X) == type(args):
        if idx_omega != idx_X or idx_omega != idx_args or idx_X != idx_args:
            raise ValueError(
                "The start index of differential form, vector field and \
                vector of argements must be equal")
    if isinstance(omega, type(X)) and idx_omega != idx_X:
        raise ValueError(
            "The start index of differential form and vector field must be \
            equal")
    idx_st = idx_omega

# Creating the output array in accordance with start indexes
    n = omega.shape[0]  # the dimensionality of the input array
    r = len(omega.shape)  # the rank of the input array
    a = Arraypy([r, n, idx_st])
    valence_list = [(-1) for k in range(r)]
    diff_Lie = a.to_tensor(valence_list)

    # Calculation
    idx = diff_Lie.start_index
    if isinstance(args, (TensorArray, Arraypy)):
        args = args.to_list()
    if isinstance(X, (TensorArray, Arraypy)):
        X = X.to_list()
    l_idx = len(idx)
    for p in range(len(diff_Lie)):
        for k in range(l_idx + 1):
            tuple_list_indx = [
                replace_index_to_k(idx, f, k + idx_st) for f in range(l_idx)]
            # the intermediate value
            diff_omega = diff(omega[idx], args[k]) * X[k]
            for j in range(l_idx):
                diff_Lie[idx] += diff(X[k], args[idx[j] - idx_st]) *\
                    omega[tuple_list_indx[j]]
            diff_Lie[idx] = diff_Lie[idx] + diff_omega
        idx = diff_Lie.next_index(idx)
    # Output
    return diff_Lie


class Wedge_array():

    """This class contains the elements of the array that has the numbers \ in
    the index are in ascending order.

    Examples:
    =========

    index=[1,2,3]
    index=[0,1,2]
    index=[10,40,60]

    The index and the element of this index is stored in the dictionary.

    """

    def __init__(self, *arg):
        # - dictionary. Dictionary key is an element index.
        # Dictionary value - is array value.
        self._output = dict()

    def print_index_value(self):
        for key, value in self._output.items():
            return str(key) + ':' + str(value)


def tensor2wedgearray(A):
    """
    Function takes a skew-symmetric array and will give array in
    which the indices are in ascending order. Return array of type
    Wedge_array.

    Examples:
    =========
    >>> from sympy import symbols
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy.tensor.tensor_fields import tensor2wedgearray

    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> y2 = TensorArray(Arraypy((3, 3)), (-1, -1))
    >>> y2[0, 1] = x3
    >>> y2[0, 2] = -x2
    >>> y2[1, 0] = -x3
    >>> y2[1, 2] = x1
    >>> y2[2, 0] = x2
    >>> y2[2, 1] = -x1
    >>> A=tensor2wedgearray(y2)
    >>> for key in sorted(A._output):
    ...     print(str(key) +'=>'+ str(A._output[key]))
    (0, 1)=>x3
    (0, 2)=>-x2
    (1, 2)=>x1

    """
    # Handling of a tensor
    if not isinstance(A, (Arraypy, TensorArray)):
        raise ValueError("The type of form must be Arraypy")
    if not is_asymmetric(A):
        raise ValueError("The form must be a skew-symmetric")

    # The valency (0,q)
    if isinstance(A, (TensorArray)) and A.type_pq[0] != 0:
        raise ValueError("The valency of tensor must be (0,q)")

    # The definition of the start index
    start_idx_A = A.start_index[0]

    lst_of_new_index = []

    # index it'a a one index of A(tuple)
    for index in A.index_list:
        index_list_sort = sorted(index)
        if(index_list_sort == list(index)) and \
                len(set(list(index))) == len(index):
            lst_of_new_index.append(index)

    # creating class instance
    wedge = Wedge_array()

    # wedge._output - it's a dictionary Wedge_array
    for index in lst_of_new_index:
        wedge._output[index] = A[index]
    return wedge


def wedgearray2tensor(B):
    """Function takes the Wedge_array(or dictionary) and creates a full skew
    array.

    Examples:
    =========
    >>> from sympy import symbols
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy.tensor.tensor_fields import tensor2wedgearray, \
    wedgearray2tensor

    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> y2 = TensorArray(Arraypy((3, 3)), (-1, -1))
    >>> y2[0, 1] = x3
    >>> y2[0, 2] = -x2
    >>> y2[1, 0] = -x3
    >>> y2[1, 2] = x1
    >>> y2[2, 0] = x2
    >>> y2[2, 1] = -x1
    >>> A=tensor2wedgearray(y2)
    >>> wdg=wedgearray2tensor(A)
    >>> print(wdg)
    0  x3  -x2
    -x3  0  x1
    x2  -x1  0

    """
    if isinstance(B, Wedge_array):
        B = B._output
    else:
        raise TypeError('Input tensor must be of Wedge_array or dict type')

    list_indices = sorted(B.keys())
    # (min_index, max_index) - the range of indexes
    max_index = max(max(list_indices))
    # for the start index of arraypy
    min_index = min(min(list_indices))
    # shape output array
    rank = len(list_indices[0])

    # Creating the output array
    arr = Arraypy([rank, max_index - min_index + 1, min_index])
    valence_list = [(-1) for k in range(max_index - min_index)]
    tensor = arr.to_tensor(valence_list)

    for key, value in B.items():
        tensor[key] = value
        # multiply by the sign of the permutation
        for p in permutations(list(key)):
            tensor[p] = sign_permutations(list(p)) * value

    return tensor


def int_product(w, X):
    """Calculation of the inner product of the form and the vector field.
    Return.

    Examples:
    =========
    >>> from sympy import symbols
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy.tensor.tensor_fields import int_product

    >>> x1, x2, x3, l, m, n = symbols('x1 x2 x3 l m n')
    >>> omega2=Arraypy([2,3,1]).to_tensor((-1,-1))
    >>> omega2[1,2]=x3
    >>> omega2[1,3]=-x2
    >>> omega2[2,1]=-x3
    >>> omega2[2,3]=x1
    >>> omega2[3,1]=x2
    >>> omega2[3,2]=-x1

    >>> X_t=Arraypy([1,3,1]).to_tensor((1))
    >>> X_t[1]=l
    >>> X_t[2]=m
    >>> X_t[3]=n

    >>> print(int_product(omega2, X_t))
    -m*x3 + n*x2  l*x3 - n*x1  -l*x2 + m*x1

    """
    # Handling of a differential form
    if not isinstance(w, (TensorArray, Arraypy)):
        raise ValueError(
            "The type of differential form must be TensorArray or Arraypy")

    # The valency (0,q)
    if (len(w.type_pq) == 1 and w.type_pq[0] != 1) or (
            len(w.type_pq) > 1 and w.type_pq[0] != 0):
        raise ValueError("The valency of tensor must be (0,q)")
    if not is_asymmetric(w):
        raise ValueError("The form must be a skew-symmetric")

    # Handling of a vector field
    check_the_vector_field(X)
    if isinstance(X, (TensorArray, Arraypy)):
        idx_start_X = X.start_index[0]
    else:
        idx_start_X = 0

    # Define the start index in the output form
    idx_start_w = w.start_index[0]
    if idx_start_w != idx_start_X:
        raise ValueError(
            "The start index of tensor and vector field must be equal")

    # Creating the output array in accordance with start indexes
    n = w.shape[0]  # the dimension of the input array
    q = len(w.shape)  # the rank of the input array
    if (q != 1):
        result_tensor = Arraypy([q - 1, n, idx_start_w])
        # all indices are -1
        valence_list = [(-1) for k in range(q - 1)]
        # result tensor
        int_prod = result_tensor.to_tensor(valence_list)

        for index in int_prod.index_list:
            # the sum on k.k is [0,n-1], n is len(X)
            for k in range(len(X)):
                idx = list(index)
                idx.insert(0, k + idx_start_w)
                int_prod[index] += w[tuple(idx)] * X[k + idx_start_w]
    if q == 1:
        # int_prod is number(symbolic expression)
        int_prod = 0
        for index in w.index_list:
            for k in range(len(X)):
                int_prod += w[index] * X[k + idx_start_w]
    # Output
    return int_prod


def g_tensor(T, S, g):
    """
    The scalar product in the space of tensors.
    The result is a skew-symmetric form of (0, p-1).

    Examples:
    =========

    >>> from sympy import symbols, Matrix
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy.tensor.tensor_fields import g_tensor

    >>> x, y, z, w = symbols('x, y, z, w')
    >>> omega=Arraypy([2,3,0]).to_tensor((-1,1))
    >>> omega[0,0]=w
    >>> omega[0,1]=x
    >>> omega[1,0]=y
    >>> omega[1,1]=z
    >>> omega[2,1]=y*y
    >>> omega[2,2]=x*y*w
    >>> g = Matrix([[2,1,0],[1,3,0],[0,0,1]])
    >>> print(g_tensor(omega,omega,g))
    w**2*x**2*y**2 + 3*y**4 + (-w/5 + 2*y/5)*(2*y + z) + (3*w/5 - y/5)*(2*w + x) + (w + 3*x)*(3*x/5 - z/5) + (-x/5 + 2*z/5)*(y + 3*z)

    """
    # Handling of a input tensors
    if not isinstance(T, TensorArray) and not isinstance(S, TensorArray):
        raise ValueError(
            "The type of first and second arguments must be TensorArray")

    if sum(T.type_pq) != sum(S.type_pq):
        raise ValueError("The valency of tensor must be equal")

    # Handling of the metric tensor
    check_metric_tensor(g)
    if isinstance(g, Matrix):
        g = matrix2tensor(g, (-1, -1))

    # The definition of the start index
    idx_start_T = T.start_index[0]
    idx_start_S = S.start_index[0]

    if idx_start_S != idx_start_T:
        raise ValueError("The start index of the tensors must be equal")
    if isinstance(g, Matrix):
        idx_start_g = 0
    else:
        idx_start_g = g.start_index[0]
    if idx_start_S != idx_start_g:
        raise ValueError(
            "The start index of the tensor and metric must be equal")

    # 1.Lowering of the top indices of a tensor T
    # upper_idx_numbers it is a list with the positions on which are the upper
    # indices
    upper_idx_numbers = [k + 1 for k in range(len(T.ind_char))
                         if T.ind_char[k] == 1]
    if len(upper_idx_numbers) != 0:
        for position in upper_idx_numbers:
            T = lower_index(T, g, position)

    # 2. Upping of the lower indices of a tensor S
    # low_idx_numbers it is a list with the positions on which are the lower
    # indices
    low_idx_numbers = [k + 1 for k in range(len(S.ind_char))
                       if S.ind_char[k] == -1]

    if len(low_idx_numbers) != 0:
        for position in low_idx_numbers:
            S = raise_index(S, g, position)

    # 3.Convolution of the first index
    result = 0
    for index in S.index_list:
        result += T[index] * S[index]

    return result


def g_wedge(T, S, g):
    """
    The scalar product in the space of differential forms.
    The result is a skew-symmetric form of (0, p-1).


    Examples:
    =========

    >>> from sympy import symbols
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy.tensor.tensor_fields import g_wedge

    >>> l, m, n = symbols('l, m, n')
    >>> X_w=Arraypy([1,3,1]).to_tensor((-1))
    >>> X_w[1]=l
    >>> X_w[2]=m
    >>> X_w[3]=n

    >>> b1=Arraypy([2,3,1]).to_tensor((-1,-1))
    >>> b1[1,1]=2
    >>> b1[1,2]=1
    >>> b1[1,3]=0
    >>> b1[2,1]=1
    >>> b1[2,2]=3
    >>> b1[2,3]=0
    >>> b1[3,1]=0
    >>> b1[3,2]=0
    >>> b1[3,3]=1
    >>> g_wedge(X_w,X_w,b1)
    l*(3*l/5 - m/5) + m*(-l/5 + 2*m/5) + n**2

    """
    # Handling of the input tensor
    if isinstance(S, (Wedge_array)):
        S = wedgearray2tensor(S)
    if isinstance(T, (Wedge_array)):
        T = wedgearray2tensor(T)

    g_t = g_tensor(T, S, g)

    if (len(T.type_pq) == 1 and T.type_pq[0] != 1) or (
            len(T.type_pq) > 1 and T.type_pq[0] != 0):
        raise ValueError("The valency of tensor must be (0,q)")

    p = len(T.start_index)
    result = g_t / factorial(p)

    return result


def hodge_star(T, g):
    """The calculation actions on the forms of the Hodge operator's.

    Examples:
    =========
    >>> from sympy import symbols, Matrix
    >>> from sympy.tensor.arraypy import Arraypy, TensorArray
    >>> from sympy.tensor.tensor_fields import hodge_star

    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> y3 = TensorArray(Arraypy((3, 3, 3)), (-1, -1, -1))
    >>> y3[0, 1, 2] = 3
    >>> y3[0, 2, 1] = -3
    >>> y3[1, 0, 2] = -3
    >>> y3[1, 2, 0] = 3
    >>> y3[2, 0, 1] = 3
    >>> y3[2, 1, 0] = -3
    >>> g = Matrix([[2,1,0],[1,3,0],[0,0,1]])
    >>> print(hodge_star(y3,g))
    96*sqrt(5)/5

    """

    if not isinstance(T, (TensorArray)):
        raise ValueError(
            "The type of tensor must be TensorArray")

    if (len(T.type_pq) == 1 and T.type_pq[0] != 1) or (len(T.type_pq) > 1 and
                                                       T.type_pq[0] != 0):
        raise ValueError("The valency of tensor must be (0,q)")
    if T.rank > 1:
        if not is_asymmetric(T):
            raise ValueError("The tensor must be a skew-symmetric")

    # Handling of the metric tensor
    check_metric_tensor(g)
    if isinstance(g, (TensorArray, Arraypy)):
        idx_start_g = g.start_index[0]
        det_g = det(g.to_matrix())
    else:
        idx_start_g = 0
        det_g = det(g)
        g = matrix2tensor(g)

    # The definition of the start index
    idx_start_T = T.start_index[0]

    if idx_start_T != idx_start_g:
        raise ValueError(
            "The start index of the tensor and metric must be equal")

    # 1. Calculating of tensor mu
    n = T.shape[0]  # the dimension of the input array
    k = T.rank
    sqrt_det_g = simplify(sqrt(abs(det_g)))

    valence_list_mu = [(-1) for i in range(n)]
    mu = Arraypy([n, n, idx_start_g]).to_tensor(valence_list_mu)

    for idx in mu.index_list:
        mu[idx] = simplify(sqrt_det_g * sign_permutations(list(idx)))

    # 2. Tensor product mu and T
    uT = tensor_product(mu, T)

    # 3.Convolution by the first k-index and the last k-index
    # low_idx_numbers it is a list with the positions on which are the lower
    # indices
    low_idx_numbers = [i + 1 for i in range(k)]
    # Upping of the first k-lower indices of a tensor
    for position in low_idx_numbers:
        uT = raise_index(uT, g, position)

    # Convolution
    for i in range(k):
        uT = uT.contract(1, k + 1)
        k = k - 1

    return uT


def codiff(w, g, args, eta=0):
    """Calculation of codifferential form.

    Examples:
    =========
    >>> from sympy import symbols, Matrix
    >>> from sympy.tensor.arraypy import Arraypy
    >>> from sympy.tensor.tensor_fields import codiff

    >>> x1, x2, x3 = symbols('x1 x2 x3')
    >>> omega2=Arraypy([2,3,0]).to_tensor((-1,-1))
    >>> omega2[0,1]=x3*x2
    >>> omega2[0,2]=-x2
    >>> omega2[1,0]=-x3*x2
    >>> omega2[1,2]=x1
    >>> omega2[2,0]=x2
    >>> omega2[2,1]=-x1
    >>> g = Matrix([[2,1,0],[1,3,0],[0,0,1]])
    >>> print(codiff(omega2, g, [x1, x2, x3]))
    0 0 0

    """
    n = w.shape[0]  # the dimension of the input array
    k = w.rank
    # hodge star
    hodge_star1 = hodge_star(w, g)
    hodge_dw = dw(hodge_star1, args)

    hodge_dw_hodge = hodge_star(hodge_dw, g)
    sign = (-1) ^ (n * (k + 1) + 1 + eta)
    for i in hodge_dw_hodge.index_list:
        hodge_dw_hodge[i] = simplify(sign * hodge_dw_hodge[i])

    return hodge_dw_hodge
