from __future__ import print_function, division

from sympy import diff

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
    def __init__(self, contravariant_or_covariant, symbol):
        self.cc = contravariant_or_covariant
        self.symbol = symbol
    def __eq__(self, other):
        return self.cc == other.cc and self.symbol == other.symbol
    def __repr__(self):
        return str(self.cc) + str(self.symbol)
    def contravariant(self):
        """Raising index

        Examples
        ========
        >>> Index("_", "i").contravariant()
        ^i
        """
        return Index("^", self.symbol)
    def covariant(self):
        """Lowering index

        Examples
        ========
        >>> Index('^', 'i').covariant()
        _i
        """
        return Index("_", self.symbol)

class Tensor:
    """Tensor in a component form

    Note: indices on the left for outer level of components, while
    indices on the right for inner level of components.

    Examples
    ========
    >>> Tensor([Index("_", "i"), Index("_", "j")], [[1,2],[3,4]])
    Tensor[_i, _j] = [[1, 2], [3, 4]]
    """
    def __init__(self, indices, elements):
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
