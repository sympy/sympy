
#external

from __future__ import print_function, division

import itertools
import collections
import functools

#internal

from sympy import Tuple, S #, diff
from sympy.core.containers import Basic

from sympy import symbols

from sympy.core.sympify import _sympify
from sympy.core.compatibility import SYMPY_INTS
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify
from sympy.core.function import Derivative
from sympy.core.expr import Expr

from sympy.utilities.iterables import flatten

from sympy.combinatorics import Permutation

from sympy.tensor import Indexed

from sympy.matrices.matrices import MatrixBase
from sympy.matrices import SparseMatrix
from sympy.matrices import Matrix

#OLD N DIM ARRAY


def deprecation_warning(ob):
    print(str(ob)+" is a deprecated class, please check the documentation and update your code")


class NDimArray:
    """
    An n-dimensional Array.
    """
    
    default_assumptions="Jack"
    _mhash=None
    def _hashable_content(self):
        return ("A",)
        
    def __setitem__(self, index, value):
        if self.immutable:
            raise TypeError('immutable N-dim array')
        else:
            """
            Allows to set items to MutableDenseNDimArray.
            """
            index = self._parse_index(index)
            if not isinstance(value, MutableNDimArray):
                value = _sympify(value)

            if isinstance(value, NDimArray):
                return NotImplementedError

            if value == 0 and index in self._sparse_array:
                self._sparse_array.pop(index)
            else:
                self._sparse_array[index] = value
                
    #this might be a cleaned up version of the above?
    #no.
    
    #def __setitem__(self, index, value):
        #"""
        #Allows to set items to MutableDenseNDimArray.
        #"""
        #index = self._parse_index(index)
        #self._setter_iterable_check(value)
        #value = _sympify(value)

        #self._array[index] = value
    
    def slice_expand(self,s, dim):
        if not isinstance(s, slice):
            return (s,)
        start, stop, step = s.indices(dim)
        return [start + i*step for i in range((stop-start)//step)]
    
    
    def __getitem__(self, index):
        """
        Allows to get items from N-dim array.
        """
        
        if self.sparse:
            syindex = self._check_symbolic_index(index)
            if syindex is not None:
                return syindex

            # `index` is a tuple with one or more slices:
            if isinstance(index, tuple) and any([isinstance(i, slice) for i in index]):

                sl_factors = [self.slice_expand(i, dim) for (i, dim) in zip(index, self.shape)]
                eindices = itertools.product(*sl_factors)
                array = [self._sparse_array.get(self._parse_index(i), S.Zero) for i in eindices]
                nshape = [len(el) for i, el in enumerate(sl_factors) if isinstance(index[i], slice)]
                return type(self)(array, nshape)
            else:
                # `index` is a single slice:
                if isinstance(index, slice):
                    start, stop, step = index.indices(self._loop_size)
                    retvec = [self._sparse_array.get(ind, S.Zero) for ind in range(start, stop, step)]
                    return retvec
                # `index` is a number or a tuple without any slice:
                else:
                    index = self._parse_index(index)
                    return self._sparse_array.get(index, S.Zero)

        else:
            
            syindex = self._check_symbolic_index(index)
            if syindex is not None:
                return syindex
            if isinstance(index, tuple) and any([isinstance(i, slice) for i in index]):
                sl_factors = [self.slice_expand(i, dim) for (i, dim) in zip(index, self.shape)]
                eindices = itertools.product(*sl_factors)
                array = [self._array[self._parse_index(i)] for i in eindices]
                nshape = [len(el) for i, el in enumerate(sl_factors) if isinstance(index[i], slice)]
                return type(self)(array, nshape)
            else:
                
                if isinstance(index, slice):
                    return self._array[index]
                    
                else:
                    index = self._parse_index(index)
                    return self._array[index]

    
    def __init__(self,iterable,shape=(),sparse=False,immutable=False, **kwargs):
        #assume dense and mutable
        self.sparse=sparse
        self.immutable=immutable
        print(sparse,immutable)
        #ImmutableDenseNDimArray
        shape, flat_list = self.handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        
        
        shape = Tuple(*map(_sympify, shape))
        flat_list = flatten(flat_list)
        flat_list = Tuple(*flat_list)
        
        
        #this is for dense
        if not self.sparse==True:
            #wtf???
            #self = Basic.__new__(type(self), flat_list, shape, **kwargs)
            self._array = list(flat_list)
            #this can't be appropriate?, surely?
        
        
        if self.sparse:
            if isinstance(flat_list, dict):
                sparse_array = dict(flat_list)
                #hmmmm
            else:
                sparse_array = {}
                for i, el in enumerate(flatten(flat_list)):
                    if el != 0:
                        sparse_array[i] = _sympify(el)

            sparse_array = dict(sparse_array)

            self._sparse_array = sparse_array
            
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0

        

    def _parse_index(self, index):

        if isinstance(index, (SYMPY_INTS, Integer)):
            if index >= self._loop_size:
                raise ValueError("index out of range")
            return index

        if len(index) != self._rank:
            raise ValueError('Wrong number of array axes')

        real_index = 0
        # check if input index can exist in current indexing
        for i in range(self._rank):
            if index[i] >= self.shape[i]:
                raise ValueError('Index ' + str(index) + ' out of border')
            real_index = real_index*self.shape[i] + index[i]

        return real_index

    def _get_tuple_index(self, integer_index):
        index = []
        for i, sh in enumerate(reversed(self.shape)):
            index.append(integer_index % sh)
            integer_index //= sh
        index.reverse()
        return tuple(index)

    def _check_symbolic_index(self, index):
        # Check if any index is symbolic:
        tuple_index = (index if isinstance(index, tuple) else (index,))
        if any([(isinstance(i, Expr) and (not i.is_number)) for i in tuple_index]):
            for i, nth_dim in zip(tuple_index, self.shape):
                if ((i < 0) == True) or ((i >= nth_dim) == True):
                    raise ValueError("index out of range")
            
            return Indexed(self, *tuple_index)
        return None

    def _setter_iterable_check(self, value):
        
        if isinstance(value, (collections.Iterable, MatrixBase, NDimArray)):
            raise NotImplementedError


    def scan_iterable_shape(self,iterable):
        #replace pointer with iterable
        if not isinstance(iterable, collections.Iterable):
            return [iterable], ()

        result = []
        
        elems, shapes = zip(*[self.scan_iterable_shape(i) for i in iterable])
        if len(set(shapes)) != 1:
            raise ValueError("could not determine shape unambiguously")
        for i in elems:
            result.extend(i)
        return result, (len(shapes),)+shapes[0]


    def handle_ndarray_creation_inputs(self,iterable=None, shape=None, **kwargs):
        
        if shape is None and iterable is None:
            shape = ()
            iterable = ()
        # Construction from another `NDimArray`:
        elif shape is None and isinstance(iterable, NDimArray):
            shape = iterable.shape
            iterable = list(iterable)
        # Construct N-dim array from an iterable (numpy arrays included):
        elif shape is None and isinstance(iterable, collections.Iterable):
            print(iterable)
            print(type(iterable))
            iterable, shape = self.scan_iterable_shape(iterable)

        # Construct N-dim array from a Matrix:
        elif shape is None and isinstance(iterable, MatrixBase):
            shape = iterable.shape

        # Construct N-dim array from another N-dim array:
        elif shape is None and isinstance(iterable, NDimArray):
            shape = iterable.shape

        # Construct NDimArray(iterable, shape)
        elif shape is not None:
            pass

        #else:
        #    shape = ()
        #    iterable = (iterable,)
        
        if isinstance(shape, (SYMPY_INTS, Integer)):
            shape = (shape,)
        if any([not isinstance(dim, (SYMPY_INTS, Integer)) for dim in shape]):
            raise TypeError("Shape should contain integers only.")

        return tuple(shape), iterable

    
    def __hash__(self):
        return Basic.__hash__(self)
    
    def __iter__(self):
        #if self.sparse:
        #    return self._array.__iter__()
        #else:
        #    return self._array.__iter__()
            
        for i in range(self._loop_size):
            yield self[i]

    def __len__(self):
        """
        Overload common function len(). Returns number of elements in array.
        """
        return self._loop_size

    
    def reshape(self, *newshape):
        """
        Returns [type] instance with new shape. Elements number
        must be suitable to new shape. The only argument of method sets
        new shape.
        """
        new_total_size = functools.reduce(lambda x,y: x*y, newshape)
        if new_total_size != self._loop_size:
            raise ValueError("Invalid reshape parameters " + newshape)

        # there is no `.func` as this class does not subtype `Basic`:
        return type(self)(self._array, newshape)

    @property
    def shape(self):
        """
        Returns array shape (dimension).
        """
        return self._shape

    def rank(self):
        """
        Returns rank of array.
        """
        return self._rank
        
    
    def diff(self, *args):
        """
        Calculate the derivative of each element in the array.
        """
        return type(self)(map(lambda x: x.diff(*args), self), self.shape)

    def applyfunc(self, f):
        """
        Apply a function to each element of the N-dim array.
        """
        return type(self)(map(f, self), self.shape)
    def recursive_string(self,sh,shape_left,i,j):
        if len(shape_left) == 1:
            return "["+", ".join([str(self[e]) for e in range(i, j)])+"]"

        sh //= shape_left[0]
        return "[" + ", ".join([self.recursive_string(sh, shape_left[1:], i+e*sh, i+(e+1)*sh) for e in range(shape_left[0])]) + "]" # + "\n"*len(shape_left)

    def __str__(self):
        """
        Returns string, allows to use standard functions print() and str().
        """
        
        r=self.recursive_string(self._loop_size, self.shape, 0, self._loop_size)
        return r

    def __repr__(self):
        return self.__str__()

    
    def tomatrix(self):
        """
        Converts [type] to Matrix. Can convert only 2-dim array, else will raise error.
        """
        if self.rank() != 2:
            raise ValueError('Dimensions must be of size of 2')

        return Matrix(self.shape[0], self.shape[1], self._array)

    def tolist(self):
        """
        Conveting MutableDenseNDimArray to one-dim list
        """
        l=[]
        for i in self:
            l.append(i)
        return l
        
        
        #if len(self.shape) == 1:
            #return [self[e] for e in range(i, j)]
        #else:
            #return self.recursive_list_conversion(self._loop_size, self.shape, 0, self._loop_size)

    def recursive_list_conversion(self,sh,shape_left,i,j):
        result = []
        sh //= shape_left[0]
        for e in range(shape_left[0]):
            result.append(self.recursive_list_conversion(sh, shape_left[1:], i+e*sh, i+(e+1)*sh))
        return result
        
    def __add__(self, other):
        if not isinstance(other, NDimArray):
            raise TypeError(str(other))

        if self.shape != other.shape:
            raise ValueError("array shape mismatch")
        result_list = [i+j for i,j in zip(self, other)]

        return type(self)(result_list, self.shape)

    def __sub__(self, other):
        if not isinstance(other, NDimArray):
            raise TypeError(str(other))

        if self.shape != other.shape:
            raise ValueError("array shape mismatch")
        result_list = [i-j for i,j in zip(self, other)]

        return type(self)(result_list, self.shape)

    def __mul__(self, other):

        if isinstance(other, (collections.Iterable,NDimArray, MatrixBase)):
            raise ValueError("scalar expected, use tensorproduct(...) for tensorial product")
        other = sympify(other)
        result_list = [i*other for i in self]
        return type(self)(result_list, self.shape)

    def __rmul__(self, other):

        if isinstance(other, (collections.Iterable,NDimArray, MatrixBase)):
            raise ValueError("scalar expected, use tensorproduct(...) for tensorial product")
        other = sympify(other)
        result_list = [other*i for i in self]
        return type(self)(result_list, self.shape)

    def __div__(self, other):

        if isinstance(other, (collections.Iterable,NDimArray, MatrixBase)):
            raise ValueError("scalar expected")
        other = sympify(other)
        result_list = [i/other for i in self]
        return type(self)(result_list, self.shape)

    def __rdiv__(self, other):
        raise NotImplementedError('unsupported operation on NDimArray')

    def __neg__(self):
        result_list = [-i for i in self]
        return type(self)(result_list, self.shape)

    def __eq__(self, other):
        """
        NDimArray instances can be compared to each other.
        Instances equal if they have same shape and data.
        """
        if not isinstance(other, NDimArray):
            return False
        return (self.shape == other.shape) and (list(self) == list(other))

    def __ne__(self, other):
        return not self == other

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def _eval_diff(self, *args, **kwargs):
        if kwargs.pop("evaluate", True):
            return self.diff(*args)
        else:
            return Derivative(self, *args, **kwargs)

    def _eval_transpose(self):
        if self.rank() != 2:
            raise ValueError("array rank not 2")
        
        return permutedims(self, (1, 0))

    def transpose(self):
        return self._eval_transpose()

    def _eval_conjugate(self):
        return type(self)([i.conjugate() for i in self], self.shape)


    def common_new(self,iterable,shape,*args,**kwargs):
        a=1
        
        
    def conjugate(self):
        return self._eval_conjugate()

    def _eval_adjoint(self):
        return self.transpose().conjugate()

    def adjoint(self):
        return self._eval_adjoint()
        
    #@classmethod
    def zeros(self, *shape):
        list_length = functools.reduce(lambda x, y: x*y, shape)
        the_zeros=[0]*list_length
        r=NDimArray(the_zeros, shape)
        return r

class ImmutableNDimArray(NDimArray):
    _op_priority = 11.0

#OLD ARRAYOP

def _arrayfy(a):

    if isinstance(a, NDimArray):
        return a
    if isinstance(a, (MatrixBase, list, tuple, Tuple)):
        return ImmutableDenseNDimArray(a)
    return a

def tensorproduct(a,b):
    """
    Tensor product among scalars or array-like objects.
    """
    
    
    #what the fuck...
    
    #if len(args) == 0:
        #return S.One
    #if len(args) == 1:
        #return NDimArray(args[0])
    #if len(args) == 2:
        #print(args)
        #NDimArray(*args)
    #if len(args) > 2:
        #return tensorproduct(tensorproduct(args[0], args[1]), *args[2:])
    #if not isinstance(a, NDimArray) or not isinstance(b, NDimArray):
        #return a*b

    #al = list(a)
    #bl = list(b)
    
    product_list = [i*j for i in al for j in bl]
    
    #erm. no. what.
    
    return ImmutableDenseNDimArray(product_list, a.shape + b.shape)


def tensorcontraction(array, *contraction_axes):
    """
    Contraction of an array-like object on the specified axes.
    """
    array = _arrayfy(array)
    #print(array)
    # Verify contraction_axes:
    taken_dims = set([])
    for axes_group in contraction_axes:
        if not isinstance(axes_group, collections.Iterable):
            raise ValueError("collections of contraction axes expected")

        dim = array.shape[axes_group[0]]

        for d in axes_group:
            if d in taken_dims:
                raise ValueError("dimension specified more than once")
            if dim != array.shape[d]:
                raise ValueError("cannot contract between axes of different dimension")
            taken_dims.add(d)

    rank = array.rank()

    remaining_shape = [dim for i, dim in enumerate(array.shape) if i not in taken_dims]
    cum_shape = [0]*rank
    _cumul = 1
    for i in range(rank):
        cum_shape[rank - i - 1] = _cumul
        _cumul *= int(array.shape[rank - i - 1])

    # DEFINITION: by absolute position it is meant the position along the one
    # dimensional array containing all the tensor components.

    # Possible future work on this module: move computation of absolute
    # positions to a class method.

    # Determine absolute positions of the uncontracted indices:
    remaining_indices = [[cum_shape[i]*j for j in range(array.shape[i])]
                         for i in range(rank) if i not in taken_dims]

    # Determine absolute positions of the contracted indices:
    summed_deltas = []
    for axes_group in contraction_axes:
        lidx = []
        for js in range(array.shape[axes_group[0]]):
            lidx.append(sum([cum_shape[ig] * js for ig in axes_group]))
        summed_deltas.append(lidx)

    # Compute the contracted array:
    #
    # 1. external for loops on all uncontracted indices.
    #    Uncontracted indices are determined by the combinatorial product of
    #    the absolute positions of the remaining indices.
    # 2. internal loop on all contracted indices.
    #    It sum the values of the absolute contracted index and the absolute
    #    uncontracted index for the external loop.
    contracted_array = []
    for icontrib in itertools.product(*remaining_indices):
        index_base_position = sum(icontrib)
        isum = S.Zero
        for sum_to_index in itertools.product(*summed_deltas):
            isum += array[index_base_position + sum(sum_to_index)]

        contracted_array.append(isum)

    if len(remaining_indices) == 0:
        assert len(contracted_array) == 1
        return contracted_array[0]

    return type(array)(contracted_array, remaining_shape)


def derive_by_array(expr, dx):
    """
    Derivative by arrays. Supports both arrays and scalars.

    Given the array `A_{i_1, \ldots, i_N}` and the array `X_{j_1, \ldots, j_M}`
    this function will return a new array `B` defined by

    `B_{j_1,\ldots,j_M,i_1,\ldots,i_N} := \frac{\partial A_{i_1,\ldots,i_N}}{\partial X_{j_1,\ldots,j_M}}`

    """
    array_types = (collections.Iterable, MatrixBase, NDimArray)

    if isinstance(dx, array_types):
        dx = ImmutableDenseNDimArray(dx)
        for i in dx:
            if not i._diff_wrt:
                raise ValueError("cannot derive by this array")

    if isinstance(expr, array_types):
        expr = ImmutableDenseNDimArray(expr)
        if isinstance(dx, array_types):
            new_array = [[y.diff(x) for y in expr] for x in dx]
            return type(expr)(new_array, dx.shape + expr.shape)
        else:
            return expr.diff(dx)
    else:
        if isinstance(dx, array_types):
            return ImmutableDenseNDimArray([expr.diff(i) for i in dx], dx.shape)
        else:
            return expr.diff(dx)


def permutedims(expr, perm):
    """
    Permutes the indices of an array.

    Parameter specifies the permutation of the indices.

    """
    #print("expr",expr)
    #print("perm",perm)
    if not isinstance(expr, NDimArray):
        raise TypeError("expression has to be an N-dim array")
        
    if not isinstance(perm, Permutation):
        perm = Permutation(list(perm))
    #print(perm.size,expr.rank())
    if perm.size != expr.rank():
        raise ValueError("wrong permutation size")

    # Get the inverse permutation:
    iperm = ~perm

    indices_span = perm([range(i) for i in expr.shape])
    
    new_array = [None]*len(expr)
    for i, idx in enumerate(itertools.product(*indices_span)):
        t = iperm(idx)
        new_array[i] = expr[t]
        
    new_shape = perm(expr.shape)
    
    
    rt=type(expr)(new_array, new_shape)
    
    return rt

class MutableNDimArray(NDimArray):
    def __new__(self, *args, **kwargs):
        deprecation_warning(MutableNDimArray)
        return NDimArray(*args, **kwargs)
    
    
class SparseNDimArray(NDimArray):
    def __new__(self, *args, **kwargs):
        deprecation_warning(SparseNDimArray)
        sparse=True
        return ImmutableSparseNDimArray(*args,sparse, **kwargs)
        
class ImmutableSparseNDimArray(NDimArray):

    def __new__(self, iterable=None, shape=None, **kwargs):
        
        deprecation_warning(ImmutableSparseNDimArray)
        sparse=True
        immutable=True        
        A=NDimArray(iterable,shape,sparse,immutable,**kwargs)
        return A

class MutableSparseNDimArray(NDimArray):

    def __new__(self, iterable=None, shape=None, **kwargs):
        
        deprecation_warning(MutableSparseNDimArray)
        
        sparse=True
        immutable=False
        
        A=NDimArray(iterable,shape,sparse,immutable,**kwargs)
        return A
        
class DenseNDimArray(NDimArray):
    def __new__(self, *args, **kwargs):
        
        deprecation_warning(DenseNDimArray)
        
        return NDimArray(*args, **kwargs)

class ImmutableDenseNDimArray(NDimArray):

    def __new__(self, iterable=None, shape=None, **kwargs):
        
        deprecation_warning(ImmutableDenseNDimArray)
        
        sparse=False
        immutable=True
        A=NDimArray(iterable,shape,sparse,immutable,**kwargs)
        return A

class MutableDenseNDimArray(NDimArray):#DenseNDimArray, MutableNDimArray):

    def __new__(self, iterable=None, shape=None, **kwargs):
        
        deprecation_warning(MutableDenseNDimArray)
        
        immutable=False
        sparse=False
        A=NDimArray(iterable,shape,sparse,immutable,**kwargs)
        return A
        
def testing():
    x,y,z,t=sympy.symbols("x y z t")
    
    NDimArray([x,y,z,t],[2,2])
    
    a=NDimArray([1,2,3,4],[2,2])
    b=NDimArray([1,2,3,4],[2,2])
    
    
    

if __name__=="__main__":
    testing()
    
