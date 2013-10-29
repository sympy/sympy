from __future__ import print_function, division

from collections import Counter
from sympy import diff, Matrix, trace

__all__ = ["Index", "Tensor", "Metric", "partial_diff"]

# Support functions

def map_n(f, *nested_list):
    """A generalization of "map" in functional programming
    to higher ranked lists and nested lists.

    Examples
    ========
    >>> print(map_n(lambda x: x+2, 3))
    5
    >>> print(map_n(lambda x: x+2, [3,[4,5]]))
    [5, [6, 7]]
    >>> print(map_n(lambda x,y: x+y+2, 3,4))
    9
    >>> print(map_n(lambda x,y: x+y+2, [3,[4,5]], [6,[7,8]]))
    [11, [13, 15]]
    >>> print(map_n(lambda x,y: x+y+2, [3,[4,5]], [6,7]))
    IndexError: Form of function arguments don't match
    """
    if len(set(map(lambda x: isinstance(x, list), nested_list))) > 1:
        raise IndexError("Form of function arguments don't match")
    elif not isinstance(nested_list[0], list): return f(*nested_list)
    else: return list(map(lambda *x: map_n(f, *x), *nested_list))

def map_on_matrix(f, nested_list):
    """Map on the innerest matrix structure

    Examples
    ========
    >>> map_on_matrix(trace, [[1,2],[3,4]])
    5
    >>> map_on_matrix(trace, [[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
    [5, 13]
    """
    if not isinstance(nested_list[0][0], list): return f(Matrix(nested_list))
    else: return list(map(lambda x: map_on_matrix(f, x), nested_list))

def transpose(rank2_list):
    """Transpose a rank-2 list

    Examples
    ========
    >>> transpose([[1,2],[3,4]])
    [[1, 3], [2, 4]]
    """
    return list(map(list, zip(*rank2_list)))

def transpose_to_outer_level(nested_list, level):
    """The generic form of transpose a list

    Note: the outerest list have level 0
    
    Examples
    ========
    >>> transpose_to_outer_level([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], 0)
    [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]
    >>> transpose_to_outer_level([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], 1)
    [[[1, 2], [5, 6]], [[3, 4], [7, 8]]]
    >>> transpose_to_outer_level([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], 2)
    [[[1, 3], [5, 7]], [[2, 4], [6, 8]]]    
    """
    if level == 0: return nested_list
    else: return transpose(map(lambda x: transpose_to_outer_level(x, level-1),
                               nested_list))

def nested_list_depth(nested_list):
    """Count the depth of a nested list

    Examples
    ========
    >>> nested_list_depth([[1,2],[3,4]])
    2
    >>> nested_list_depth([[1,2],3])
    ValueError: Nested list doesn't have the same depth
    """
    if not isinstance(nested_list, list): return 0
    else:
        sublist_len = [nested_list_depth(x) for x in nested_list]
        if len(set(sublist_len)) > 1: raise ValueError("Nested list doesn't have the same depth")
        else: return sublist_len[0] + 1
 
# Tensor calculation related codes

class Index:
    """Tensor index with the choice of either
    contravariant or covariant and a symbol

    Examples
    ========
    >>> i = Index("^", "i") # contravariant index
    >>> i
    ^i
    i = Index("_", "i") # covariant index
    >>> i
    _i
    """
    contravariant = "^"
    covariant = "_"
    def __init__(self, contravariant_or_covariant, symbol):
        self.cc = contravariant_or_covariant
        self.symbol = symbol
    def __eq__(self, other):
        return self.cc == other.cc and self.symbol == other.symbol
    def __repr__(self):
        return str(self.cc) + str(self.symbol)
    def raising(self):
        """Raising index

        Examples
        ========
        >>> Index("_", "i").contravariant()
        ^i
        """
        return Index(self.contravariant, self.symbol)
    def is_contravariant(self):
        return self.cc == self.contravariant
    def lowering(self):
        """Lowering index

        Examples
        ========
        >>> Index('^', 'i').covariant()
        _i
        """
        return Index(self.covariant, self.symbol)
    def is_covariant(self):
        return self.cc == self.covariant

class Tensor:
    """Tensor in a component form

    Note: indices on the left for outer level of components, while
    indices on the right for inner level of components.

    Examples
    ========
    >>> Tensor([Index("_", "i"), Index("_", "j")], [[1,2],[3,4]])
    Tensor[_i, _j] = [[1, 2], [3, 4]]
    >>> Tensor([Index("_", "i")], [[1, 2], [3, 4]])
    ValueError: Number of indices and the form of tensor componenets don't match
    """
    def __init__(self, indices, elements):
        if len(indices) != nested_list_depth(elements):
            raise ValueError("Number of indices and the form of tensor componenets don't match")
        else:
            self.idx = indices
            self.ele = elements
    def __repr__(self):
        return "Tensor" + str(self.idx) + " = " + str(self.ele)
    def __add__(self, other):
        """Operator Overloading of + for class Tensor.
        Add two tensors elment by element

        Examples
        ========
        >>> t1 = Tensor([Index("_", "i"), Index("_", "j")], [[1,2],[3,4]])
        >>> t2 = Tensor([Index("_", "i"), Index("_", "j")], [[5,6],[7,8]])
        >>> t1 + t2
        Tensor[_i, _j] = [[6, 8], [10, 12]]
        >>> t2 = Tensor([Index("_", "i"), Index("^", "j")], [[5,6],[7,8]])
        >>> t1 + t2
        ValueError: Tensor forms don't match
        """
        if self.idx != other.idx: raise ValueError("Tensor forms don't match")
        else: return Tensor(self.idx, map_n(lambda x,y: x+y, self.ele, other.ele))
    def __sub__(self, other):
        """Operator Overloading of - for class Tensor.
        Substract two tensors elment by element

        Examples
        ========
        >>> t1 = Tensor([Index("_", "i"), Index("_", "j")], [[1,2],[3,4]])
        >>> t2 = Tensor([Index("_", "i"), Index("_", "j")], [[5,6],[7,8]])
        >>> t1 - t2
        Tensor[_i, _j] = [[-4, -4], [-4, -4]]
        >>> t2 = Tensor([Index("_", "i"), Index("^", "j")], [[5,6],[7,8]])
        >>> t1 - t2
        ValueError: Tensor forms don't match
        """
        if self.idx != other.idx: raise ValueError("Tensor forms don't match")
        else: return Tensor(self.idx, map_n(lambda x,y: x-y, self.ele, other.ele))
    def scalar_multiplication(self, k):
        return Tensor(self.idx, map_n(lambda x: k * x, self.ele))
    def __mul__(self, other):
        if isinstance(other, Tensor):
            return Tensor(self.idx + other.idx,
                          map_n(lambda x: other.scalar_multiplication(x).ele, self.ele))
        else: return self.scalar_multiplication(other)
        """Operator Overloading of * for class Tensor: tensor product

        Examples
        ========
        >>> t1 = Tensor([Index("_", "i"), Index("_", "j")], [[1,2],[3,4]])
        >>> t1*2 # Note: 2*t1 is not working
        Tensor[_i, _j] = [[2, 4], [6, 8]]
        >>> t2 = Tensor([Index("^", "k")], [5,6])
        >>> t1*t2
        Tensor[_i, _j, ^k] = [[[5, 6], [10, 12]], [[15, 18], [20, 24]]]
        """
    def partial_diff_by_scalar(self, scalar):
        return Tensor(self.idx, map_n(lambda x: diff(x, scalar), self.ele))
    def partial_diff_on_scalar(self, scalar):
        return Tensor(self.idx, map_n(lambda x: diff(scalar, x), self.ele))
    def partial_diff(self, other):
        """Partial derivative for tensors

        Examples
        ========
        >>> from sympy import symbols
        >>> x, y, z = symbols("x y z")
        >>> t1 = Tensor([Index("_", "i"), Index("_", "j")], [[x,y],[z,2]])
        >>> t1.partial_diff(x)
        Tensor[_i, _j] = [[1, 0], [0, 0]]
        >>> t2 = Tensor([Index("^", "k")], [x, y])
        >>> t1.partial_diff(t2)
        Tensor[_i, _j, ^k] = [[[1, 0], [0, 1]], [[0, 0], [0, 0]]]
        """
        if isinstance(other, Tensor):
            return Tensor(self.idx + other.idx,
                          map_n(lambda y: other.partial_diff_on_scalar(y).ele, self.ele))
        else: return self.partial_diff_by_scalar(other)
    def change_indices(self, indices):
        """Change the indices of a tensor, but keep its elements unchanged

        Examples
        ========
        >>> t = Tensor([Index("_", "i"), Index("_", "j")], [[1,2],[3,4]])
        >>> t.change_indices([Index("_", "k"), Index("^", "l")])
        Tensor[_k, ^l] = [[1, 2], [3, 4]]
        """
        return Tensor(indices, self.ele)
    def switch_indices(self, indices):
        """Switch the indices of a tensor, and change the order of its
        components at the same time

        Examples
        ========
        >>> t = Tensor([Index("_", "i"), Index("_", "j"), Index("^", "k")], [[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
        t.switch_indices([Index("^", "k"), Index("_", "j"), Index("_", "i")])
        Tensor[^k, _j, _i] = [[[1, 5], [3, 7]], [[2, 6], [4, 8]]]
        """
        i = self.idx.index(indices[0])
        if len(indices) == 1: return self
        else: return Tensor(indices,
                            list(map(lambda x: Tensor(self.idx[:i]+self.idx[(i+1):], x).switch_indices(indices[1:]).ele,
                                     transpose_to_outer_level(self.ele, i))))
    def einstein_summation(self):
        """Einstein Summation of dummy indices

        Examples
        ========
        >>> t = Tensor([Index("_", "i"), Index("^", "i")], [[1, 2], [3, 4]])
        >>> t.einstein_summation()
        5
        >>> t = Tensor([Index("_", "i"), Index("^", "i"), Index("^", "k")], [[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
        >>> t.einstein_summation()
        Tensor[^k] = [8, 10]
        """
        indices_symbols = list(map(lambda x: x.symbol, self.idx))
        indices_dummy = [x for x,n in dict(Counter(indices_symbols)).items() if n == 2]
        if len(indices_dummy) == 0: return self
        else:
            idx = ([x for x in self.idx if x.symbol != indices_dummy[0]]
                       + [Index(Index.covariant, indices_dummy[0]), Index(Index.contravariant, indices_dummy[0])])
            t = Tensor(idx[:-2], map_on_matrix(trace, self.switch_indices(idx).ele))
            if len(t.idx) == 0: return t.ele
            else: return t.einstein_summation()

def partial_diff(t1, t2):
    """Generic function for partial derivative of tensors
    or scalars(symbols)

    Examples
    ========
    >>> from sympy import symbols
    >>> x, y, z = symbols("x y z")
    >>> partial_diff(2*x, x)
    2
    >>> partial_diff(Tensor([Index("^", "k")], [x, y]), x)
    Tensor[^k] = [1, 0]
    >>> partial_diff(x*y, Tensor([Index("^", "k")], [x, y]))
    Tensor[^k] = [y, x]
    >>> partial_diff(Tensor([Index("_", "i")], [x*y, y]),
                     Tensor([Index("^", "j")], [y, z]))
    Tensor[_i, ^j] = [[x, 0], [1, 0]]
    """
    if isinstance(t1, Tensor): return t1.partial_diff(t2)
    elif isinstance(t2, Tensor): return t2.partial_diff_on_scalar(t1)
    else: return diff(t1, t2)

class Metric(Tensor):
    """Define the metric as a rank-2 tensor

    Examples
    ========
    >>> Metric([Index("_", "i"), Index("_", "j")], [[1, 3], [3, 4]])
    Metric[_i, _j] = [[1, 3], [3, 4]]
    >>> Metric([Index("_", "i")], [1, 2])
    ValueError: Tensor form is not a metric
    >>> Metric([Index("_", "i"), Index("^", "j")], [[1, 3], [3, 4]])
    ValueError: Metric indices should be either both covariant or both contravariant
    """
    def __init__(self, indices, elements):
        if len(indices) != 2:
            raise ValueError("Tensor form is not a metric")
        elif indices[0].cc != indices[1].cc:
            raise ValueError("Metric indices should be either both covariant or both contravariant")
        else: Tensor.__init__(self, indices, elements)
    def __repr__(self):
        return "Metric" + str(self.idx) + " = " + str(self.ele)
    def change_indices(self, indices):
        """Change the indices of the metric, but keep its elements unchanged

        Examples
        ========
        >>> m = Metric([Index("_", "i"), Index("_", "j")], [[1, 3], [3, 4]])
        >>> m.change_indices([Index("_", "j"), Index("_", "i")])
        Metric[_j, _i] = [[1, 3], [3, 4]]
        """
        return Metric(indices, self.ele)
    def upperupper(self):
        """Raising the metric indices

        Examples
        ========
        >>> Metric([Index("^", "i"), Index("^", "j")], [[1, 3], [3, 4]]).upperupper()
        Metric[^i, ^j] = [[1, 3], [3, 4]]
        >>> Metric([Index("_", "i"), Index("_", "j")], [[1, 3], [3, 4]]).upperupper()
        Metric[^i, ^j] = [[-4/5, 3/5], [3/5, -1/5]]
        
        >>> m = Metric([Index("_", "i"), Index("_", "j")], [[1, 3], [3, 4]])
        >>> (m * m.upperupper().change_indices([Index("^", "j"), Index("^", "k")])).einstein_summation()
        Tensor[_i, ^k] = [[1, 0], [0, 1]]
        """
        if self.idx[0].cc == Index.contravariant: return self
        else:
            return Metric(list(map(lambda x: x.raising(), self.idx)), (Matrix(self.ele)**-1).tolist())
    def lowerlower(self):
        """Lowering the metric indices

        Examples
        ========
        >>> Metric([Index("_", "i"), Index("_", "j")], [[1, 3], [3, 4]]).lowerlower()
        Metric[^i, ^j] = [[1, 3], [3, 4]]
        >>> Metric([Index("^", "i"), Index("^", "j")], [[1, 3], [3, 4]]).lowerlower()
        Metric[^i, ^j] = [[-4/5, 3/5], [3/5, -1/5]]
        """
        if self.idx[0].cc == Index.covariant: return self
        else:
            return Metric(list(map(lambda x: x.lowering(), self.idx)), (Matrix(self.ele)**-1).tolist())
