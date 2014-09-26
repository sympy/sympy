# -*- coding:utf-8 -*-

"""
Definition of physical dimensions.

Unit systems will be constructed on top of these dimensions.

Most of the examples in the doc use MKS system and are presented from the
computer point of view: from a human point, adding length to time is not legal
in MKS but it is in natural system; for a computer in natural system there is
no time dimension (but a velocity dimension instead) - in the basis - so the
question of adding time to length has no meaning.
"""

from __future__ import division
from copy import copy
import numbers

from sympy.core.compatibility import reduce
from sympy.core.containers import Tuple, Dict
from sympy import sympify, nsimplify, Number, Integer, Matrix, Expr


class Dimension(Expr):
    """
    This class represent the dimension of a physical quantities.

    The dimensions may have a name and a symbol. All other
    arguments are dimensional powers. They represent a characteristic of a
    quantity, giving an interpretation to it: for example (in classical
    mechanics) we know that time is different from temperature, and dimensions
    make this difference (but they do not provide any measure of these
    quantites).

        >>> from sympy.physics.unitsystems.dimensions import Dimension
        >>> length = Dimension(length=1)
        >>> length
        {'length': 1}
        >>> time = Dimension(time=1)

    Dimensions behave like a dictionary where the key is the name and the value
    corresponds to the exponent.

    Dimensions can be composed using multiplication, division and
    exponentiation (by a number) to give new dimensions. Addition and
    subtraction is defined only when the two objects are the same dimension.

        >>> velocity = length.div(time)
        >>> velocity  #doctest: +SKIP
        {'length': 1, 'time': -1}
        >>> length.add(length)
        {'length': 1}
        >>> length.pow(2)
        {'length': 2}

    Defining addition-like operations will help when doing dimensional analysis.

    Note that two dimensions are equal if they have the same powers, even if
    their names and/or symbols differ.

        >>> Dimension(length=1) == Dimension(length=1, name="length")
        True
        >>> Dimension(length=1) == Dimension(length=1, symbol="L")
        True
        >>> Dimension(length=1) == Dimension(length=1, name="length",
        ...                                  symbol="L")
        True
    """

    is_commutative = True
    is_number = False
    # make sqrt(M**2) --> M
    is_positive = True

    def __new__(cls, *args, **kwargs):
        """
        Create a new dimension.

        Possibilities are (examples given with list/tuple work also with
        tuple/list):

            >>> from sympy.physics.unitsystems.dimensions import Dimension
            >>> Dimension(length=1)
            {'length': 1}
            >>> Dimension({"length": 1})
            {'length': 1}
            >>> Dimension([("length", 1), ("time", -1)])  #doctest: +SKIP
            {'length': 1, 'time': -1}

        """

        # before setting the dict, check if a name and/or a symbol are defined
        # if so, remove them from the dict
        name = kwargs.pop('name', None)
        symbol = kwargs.pop('symbol', None)

        # pairs of (dimension, power)
        pairs = []

        # add first items from args to the pairs
        for arg in args:
            # construction with {"length": 1}
            if isinstance(arg, dict):
                arg = copy(arg)
                pairs.extend(arg.items())
            elif isinstance(arg, (Tuple, tuple, list)):
                #TODO: add construction with ("length", 1); not trivial because
                #      e.g. [("length", 1), ("time", -1)] has also length = 2

                for p in arg:
                    #TODO: check that p is a tuple
                    if len(p) != 2:
                        raise ValueError("Length of iterable has to be 2; "
                                         "'%d' found" % len(p))

                # construction with [("length", 1), ...]
                pairs.extend(arg)
            else:
            # error if the arg is not of previous types
                raise TypeError("Positional arguments can only be: "
                                "dict, tuple, list; '%s' found" % type(arg))

        pairs.extend(kwargs.items())

        # check validity of dimension key and power
        for pair in pairs:
            #if not isinstance(p[0], str):
            #    raise TypeError("key %s is not a string." % p[0])
            if not isinstance(pair[1], (numbers.Real, Number)):
                raise TypeError("Power corresponding to '%s' is not a number"
                                % pair[0])

        # filter dimensions set to zero; this avoid the following odd result:
        # Dimension(length=1) == Dimension(length=1, mass=0) => False
        # also simplify to avoid powers such as 2.00000
        pairs = [(pair[0], nsimplify(pair[1])) for pair in pairs
                 if pair[1] != 0]
        pairs.sort(key=str)

        new = Expr.__new__(cls, Dict(*pairs))
        new.name = name
        new.symbol = symbol

        new._dict = dict(pairs)

        return new

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        return self._dict[key]

    def __setitem__(self, key, value):
        raise NotImplementedError("Dimension are Immutable")

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return self._dict.items()

    def keys(self):
        """D.keys() -> list of D's keys"""
        return self._dict.keys()

    def values(self):
        """D.values() -> list of D's values"""
        return self._dict.values()

    def __iter__(self):
        """x.__iter__() <==> iter(x)"""
        return iter(self._dict)

    def __len__(self):
        """x.__len__() <==> len(x)"""
        return self._dict.__len__()

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        return self._dict.get(key, default)

    def __contains__(self, key):
        """D.__contains__(k) -> True if D has a key k, else False"""
        return key in self._dict

    def __lt__(self, other):
        return self.args < other.args

    def __str__(self):
        """
        Display the string representation of the dimension.

        Usually one will always use a symbol to denote the dimension. If no
        symbol is defined then it uses the name or, if there is no name, the
        default dict representation.

        We do *not* want to use the dimension system to find the string
        representation of a dimension because it would imply some magic in
        order to guess the "best" form. It is better to do as if we do not
        have a system, and then to design a specific function to take it into
        account.
        """

        if self.symbol is not None:
            return self.symbol
        elif self.name is not None:
            return self.name
        else:
            return repr(self)

    def __repr__(self):

        return repr(self._dict)

    def __neg__(self):
        return self

    def add(self, other):
        """
        Define the addition for Dimension.

        Addition of dimension has a sense only if the second object is the same
        dimension (we don't add length to time).
        """

        if not isinstance(other, Dimension):
            raise TypeError("Only dimension can be added; '%s' is not valid"
                            % type(other))
        elif isinstance(other, Dimension) and self != other:
            raise ValueError("Only dimension which are equal can be added; "
                             "'%s' and '%s' are different" % (self, other))

        return self

    def sub(self, other):
        # there is no notion of ordering (or magnitude) among dimension,
        # subtraction is equivalent to addition when the operation is legal
        return self.add(other)

    def pow(self, other):
        #TODO: be sure that it works with rational numbers (e.g. when dealing
        #      with dimension under a fraction)

        #TODO: allow exponentiation by an abstract symbol x
        #      (if x.is_number is True)
        #      this would be a step toward using solve and absract powers

        other = sympify(other)
        if isinstance(other, (numbers.Real, Number)):
            return Dimension([(x, y*other) for x, y in self.items()])
        else:
            raise TypeError("Dimensions can be exponentiated only with "
                            "numbers; '%s' is not valid" % type(other))

    def mul(self, other):
        if not isinstance(other, Dimension):
            #TODO: improve to not raise error: 2*L could be a legal operation
            #      (the same comment apply for __div__)
            raise TypeError("Only dimension can be multiplied; '%s' is not "
                            "valid" % type(other))

        d = dict(self)
        for key in other:
            try:
                d[key] += other[key]
            except KeyError:
                d[key] = other[key]
        d = Dimension(d)

        # if all dimensions are zero, then return 1 so that there is no more
        # dimensions
        if d.is_dimensionless:
            return 1
        else:
            return d

    def div(self, other):
        if not isinstance(other, Dimension):
            raise TypeError("Only dimension can be divided; '%s' is not valid"
                            % type(other))

        d = dict(self)
        for key in other:
            try:
                d[key] -= other[key]
            except KeyError:
                d[key] = -other[key]
        d = Dimension(d)

        if d.is_dimensionless:
            return 1
        else:
            return d

    def rdiv(self, other):
        return other * pow(self, -1)

    @property
    def is_dimensionless(self):
        """
        Check if the dimension object really has a dimension.

        A dimension should have at least one component with non-zero power.
        """

        for key in self:
            if self[key] != 0:
                return False
        else:
            return True

    @property
    def has_integer_powers(self):
        """
        Check if the dimension object has only integer powers.

        All the dimension powers should be integers, but rational powers may
        appear in intermediate steps. This method may be used to check that the
        final result is well-defined.
        """

        for key in self:
            if not isinstance(self[key], Integer):
                return False
        else:
            return True


class DimensionSystem(object):
    """
    DimensionSystem represents a coherent set of dimensions.

    In a system dimensions are of three types:

    - base dimensions;
    - derived dimensions: these are defined in terms of the base dimensions
      (for example velocity is defined from the division of length by time);
    - canonical dimensions: these are used to define systems because one has
      to start somewhere: we can not build ex nihilo a system (see the
      discussion in the documentation for more details).

    All intermediate computations will use the canonical basis, but at the end
    one can choose to print result in some other basis.

    In a system dimensions can be represented as a vector, where the components
    represent the powers associated to each base dimension.
    """

    def __init__(self, base, dims=(), name="", descr=""):
        """
        Initialize the dimension system.

        It is important that base units have a name or a symbol such that
        one can sort them in a unique way to define the vector basis.
        """

        self.name = name
        self.descr = descr

        if (None, None) in [(d.name, d.symbol) for d in base]:
            raise ValueError("Base dimensions must have a symbol or a name")

        self._base_dims = self.sort_dims(base)
        # base is first such that named dimension are keeped
        self._dims = tuple(set(base) | set(dims))

        self._can_transf_matrix = None
        self._list_can_dims = None

        if self.is_consistent is False:
            raise ValueError("The system with basis '%s' is not consistent"
                             % str(self._base_dims))

    def __str__(self):
        """
        Return the name of the system.

        If it does not exist, then it makes a list of symbols (or names) of
        the base dimensions.
        """

        if self.name != "":
            return self.name
        else:
            return "(%s)" % ", ".join(str(d) for d in self._base_dims)

    def __repr__(self):
        return "<DimensionSystem: %s>" % repr(self._base_dims)

    def __getitem__(self, key):
        """
        Shortcut to the get_dim method, using key access.
        """

        d = self.get_dim(key)

        #TODO: really want to raise an error?
        if d is None:
            raise KeyError(key)

        return d

    def __call__(self, unit):
        """
        Wrapper to the method print_dim_base
        """

        return self.print_dim_base(unit)

    def get_dim(self, dim):
        """
        Find a specific dimension which is part of the system.

        dim can be a string or a dimension object. If no dimension is found,
        then return None.
        """

        #TODO: if the argument is a list, return a list of all matching dims

        found_dim = None

        #TODO: use copy instead of direct assignment for found_dim?
        if isinstance(dim, str):
            for d in self._dims:
                if dim in (d.name, d.symbol):
                    found_dim = d
                    break
        elif isinstance(dim, Dimension):
            try:
                i = self._dims.index(dim)
                found_dim = self._dims[i]
            except ValueError:
                pass

        return found_dim

    def extend(self, base, dims=(), name='', description=''):
        """
        Extend the current system into a new one.

        Take the base and normal units of the current system to merge
        them to the base and normal units given in argument.
        If not provided, name and description are overriden by empty strings.
        """

        base = self._base_dims + tuple(base)
        dims = self._dims + tuple(dims)

        return DimensionSystem(base, dims, name, description)

    @staticmethod
    def sort_dims(dims):
        """
        Sort dimensions given in argument using their str function.

        This function will ensure that we get always the same tuple for a given
        set of dimensions.
        """

        return tuple(sorted(dims, key=str))

    @property
    def list_can_dims(self):
        """
        List all canonical dimension names.
        """

        if self._list_can_dims is None:
            gen = reduce(lambda x, y: x.mul(y), self._base_dims)
            self._list_can_dims = tuple(sorted(map(str, gen.keys())))

        return self._list_can_dims

    @property
    def inv_can_transf_matrix(self):
        """
        Compute the inverse transformation matrix from the base to the
        canonical dimension basis.

        It corresponds to the matrix where columns are the vector of base
        dimensions in canonical basis.

        This matrix will almost never be used because dimensions are always
        define with respect to the canonical basis, so no work has to be done
        to get them in this basis. Nonetheless if this matrix is not square
        (or not invertible) it means that we have chosen a bad basis.
        """

        matrix = reduce(lambda x, y: x.row_join(y),
                        [self.dim_can_vector(d) for d in self._base_dims])

        return matrix

    @property
    def can_transf_matrix(self):
        """
        Compute the canonical transformation matrix from the canonical to the
        base dimension basis.

        It is the inverse of the matrix computed with inv_can_transf_matrix().
        """

        #TODO: the inversion will fail if the system is inconsistent, for
        #      example if the matrix is not a square
        if self._can_transf_matrix is None:
            self._can_transf_matrix = reduce(lambda x, y: x.row_join(y),
                                             [self.dim_can_vector(d)
                                              for d in self._base_dims]).inv()

        return self._can_transf_matrix

    def dim_can_vector(self, dim):
        """
        Vector representation in terms of the canonical base dimensions.
        """

        vec = []
        for d in self.list_can_dims:
            vec.append(dim.get(d, 0))

        return Matrix(vec)

    def dim_vector(self, dim):
        """
        Vector representation in terms of the base dimensions.
        """

        return self.can_transf_matrix * self.dim_can_vector(dim)

    def print_dim_base(self, dim):
        """
        Give the string expression of a dimension in term of the basis.

        Dimensions are displayed by decreasing power.
        """

        res = ""

        for (d, p) in sorted(zip(self._base_dims, self.dim_vector(dim)),
                             key=lambda x: x[1], reverse=True):
            if p == 0:
                continue
            elif p == 1:
                res += "%s " % str(d)
            else:
                res += "%s^%d " % (str(d), p)

        return res.strip()

    @property
    def dim(self):
        """
        Give the dimension of the system.

        That is return the number of dimensions forming the basis.
        """

        return len(self._base_dims)

    @property
    def is_consistent(self):
        """
        Check if the system is well defined.
        """

        # not enough or too many base dimensions compared to independent
        # dimensions
        # in vector language: the set of vectors do not form a basis
        if self.inv_can_transf_matrix.is_square is False:
            return False

        return True
