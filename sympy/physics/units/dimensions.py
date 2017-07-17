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

import collections

from sympy.core.compatibility import reduce, string_types
from sympy import sympify, Integer, Matrix, Symbol, S, Abs
from sympy.core.expr import Expr


class Dimension(Expr):
    """
    This class represent the dimension of a physical quantities.

    The ``Dimension`` constructor takes as parameters a name and an optional
    symbol.

    For example, in classical mechanics we know that time is different from
    temperature and dimensions make this difference (but they do not provide
    any measure of these quantites.

        >>> from sympy.physics.units import Dimension
        >>> length = Dimension('length')
        >>> length
        Dimension(length)
        >>> time = Dimension('time')
        >>> time
        Dimension(time)

    Dimensions can be composed using multiplication, division and
    exponentiation (by a number) to give new dimensions. Addition and
    subtraction is defined only when the two objects are the same dimension.

        >>> velocity = length / time
        >>> velocity
        Dimension(length/time)
        >>> velocity.get_dimensional_dependencies()
        {'length': 1, 'time': -1}
        >>> length + length
        Dimension(length)
        >>> l2 = length**2
        >>> l2
        Dimension(length**2)
        >>> l2.get_dimensional_dependencies()
        {'length': 2}

    """

    _op_priority = 13.0

    _dimensional_dependencies = dict()

    is_commutative = True
    is_number = False
    # make sqrt(M**2) --> M
    is_positive = True
    is_real = True

    def __new__(cls, name, symbol=None):

        if isinstance(name, string_types):
            name = Symbol(name)
        else:
            name = sympify(name)

        if not isinstance(name, Expr):
            raise TypeError("Dimension name needs to be a valid math expression")

        if isinstance(symbol, string_types):
            symbol = Symbol(symbol)
        elif symbol is not None:
            assert isinstance(symbol, Symbol)

        if symbol is not None:
            obj = Expr.__new__(cls, name, symbol)
        else:
            obj = Expr.__new__(cls, name)

        obj._name = name
        obj._symbol = symbol
        return obj

    @property
    def name(self):
        return self._name

    @property
    def symbol(self):
        return self._symbol

    def __hash__(self):
        return Expr.__hash__(self)

    def __eq__(self, other):
        if isinstance(other, Dimension):
            return self.get_dimensional_dependencies() == other.get_dimensional_dependencies()
        return False

    def __str__(self):
        """
        Display the string representation of the dimension.
        """
        if self.symbol is None:
            return "Dimension(%s)" % (self.name)
        else:
            return "Dimension(%s, %s)" % (self.name, self.symbol)

    def __repr__(self):
        return self.__str__()

    def __neg__(self):
        return self

    def _register_as_base_dim(self):
        if self.name in self._dimensional_dependencies:
            raise IndexError("already in dependecies dict")
        if not self.name.is_Symbol:
            raise TypeError("Base dimensions need to have symbolic name")
        self._dimensional_dependencies[self.name] = {self.name: 1}

    def __add__(self, other):
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

    def __sub__(self, other):
        # there is no notion of ordering (or magnitude) among dimension,
        # subtraction is equivalent to addition when the operation is legal
        return self + other

    def __pow__(self, other):
        return self._eval_power(other)

    def _eval_power(self, other):
        other = sympify(other)
        return Dimension(self.name**other)

    def __mul__(self, other):
        if not isinstance(other, Dimension):
            return self

        return Dimension(self.name*other.name)

    def __div__(self, other):
        if not isinstance(other, Dimension):
            return self

        return Dimension(self.name/other.name)

    def __rdiv__(self, other):
        return other * pow(self, -1)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def get_dimensional_dependencies(self, mark_dimensionless=False):
        name = self.name
        dimdep = self._get_dimensional_dependencies_for_name(name)
        if mark_dimensionless and dimdep == {}:
            return {'dimensionless': 1}
        return dimdep

    @classmethod
    def _from_dimensional_dependencies(cls, dependencies):
        return reduce(lambda x, y: x * y, (
            Dimension(d)**e for d, e in dependencies.items()
        ))

    @classmethod
    def _get_dimensional_dependencies_for_name(cls, name):

        if name.is_Symbol:
            if name.name in Dimension._dimensional_dependencies:
                return Dimension._dimensional_dependencies[name.name]
            else:
                return {}

        if name.is_Number:
            return {}

        if name.is_Mul:
            ret = collections.defaultdict(int)
            dicts = [Dimension._get_dimensional_dependencies_for_name(i) for i in name.args]
            for d in dicts:
                for k, v in d.items():
                    ret[k] += v
            return {k: v for (k, v) in ret.items() if v != 0}

        if name.is_Pow:
            if name.exp == 0:
                return {}
            dim = Dimension._get_dimensional_dependencies_for_name(name.base)
            return {k: v*name.exp for (k, v) in dim.items()}

        if name.is_Function:
            args = (Dimension._from_dimensional_dependencies(
                Dimension._get_dimensional_dependencies_for_name(arg)
            ) for arg in name.args)
            result = name.func(*args)

            if isinstance(result, cls):
                return result.get_dimensional_dependencies()
            # TODO shall we consider a result that is not a dimension?
            # return Dimension._get_dimensional_dependencies_for_name(result)

    @property
    def is_dimensionless(self):
        """
        Check if the dimension object really has a dimension.

        A dimension should have at least one component with non-zero power.
        """
        dimensional_dependencies = self.get_dimensional_dependencies()
        return dimensional_dependencies == {}

    @property
    def has_integer_powers(self):
        """
        Check if the dimension object has only integer powers.

        All the dimension powers should be integers, but rational powers may
        appear in intermediate steps. This method may be used to check that the
        final result is well-defined.
        """

        for dpow in self.get_dimensional_dependencies().values():
            if not isinstance(dpow, (int, Integer)):
                return False
        else:
            return True


# base dimensions (MKS)
length = Dimension(name="length", symbol="L")
mass = Dimension(name="mass", symbol="M")
time = Dimension(name="time", symbol="T")

# base dimensions (MKSA not in MKS)
current = Dimension(name='current', symbol='I')

# other base dimensions:
temperature = Dimension("temperature", "T")
amount_of_substance = Dimension("amount_of_substance")
luminous_intensity = Dimension("luminous_intensity")

# derived dimensions (MKS)
velocity = Dimension(name="velocity")
acceleration = Dimension(name="acceleration")
momentum = Dimension(name="momentum")
force = Dimension(name="force", symbol="F")
energy = Dimension(name="energy", symbol="E")
power = Dimension(name="power")
pressure = Dimension(name="pressure")
frequency = Dimension(name="frequency", symbol="f")
action = Dimension(name="action", symbol="A")
volume = Dimension("volume")

# derived dimensions (MKSA not in MKS)
voltage = Dimension(name='voltage', symbol='U')
impedance = Dimension(name='impedance', symbol='Z')
conductance = Dimension(name='conductance', symbol='G')
capacitance = Dimension(name='capacitance')
inductance = Dimension(name='inductance')
charge = Dimension(name='charge', symbol='Q')
magnetic_density = Dimension(name='magnetic_density', symbol='B')
magnetic_flux = Dimension(name='magnetic_flux')

# Create dimensions according the the base units in MKSA.
# For other unit systems, they can be derived by transforming the base
# dimensional dependency dictionary.

# Dimensional dependencies for MKS base dimensions
Dimension._dimensional_dependencies["length"] = dict(length=1)
Dimension._dimensional_dependencies["mass"] = dict(mass=1)
Dimension._dimensional_dependencies["time"] = dict(time=1)

# Dimensional dependencies for base dimensions (MKSA not in MKS)
Dimension._dimensional_dependencies["current"] = dict(current=1)

# Dimensional dependencies for other base dimensions:
Dimension._dimensional_dependencies["temperature"] = dict(temperature=1)
Dimension._dimensional_dependencies["amount_of_substance"] = dict(amount_of_substance=1)
Dimension._dimensional_dependencies["luminous_intensity"] = dict(luminous_intensity=1)

# Dimensional dependencies for derived dimensions
Dimension._dimensional_dependencies["velocity"] = dict(length=1, time=-1)
Dimension._dimensional_dependencies["acceleration"] = dict(length=1, time=-2)
Dimension._dimensional_dependencies["momentum"] = dict(mass=1, length=1, time=-1)
Dimension._dimensional_dependencies["force"] = dict(mass=1, length=1, time=-2)
Dimension._dimensional_dependencies["energy"] = dict(mass=1, length=2, time=-2)
Dimension._dimensional_dependencies["power"] = dict(length=2, mass=1, time=-3)
Dimension._dimensional_dependencies["pressure"] = dict(mass=1, length=-1, time=-2)
Dimension._dimensional_dependencies["frequency"] = dict(time=-1)
Dimension._dimensional_dependencies["action"] = dict(length=2, mass=1, time=-1)
Dimension._dimensional_dependencies["volume"] = dict(length=3)

# Dimensional dependencies for derived dimensions
Dimension._dimensional_dependencies["voltage"] = dict(mass=1, length=2, current=-1, time=-3)
Dimension._dimensional_dependencies["impedance"] = dict(mass=1, length=2, current=-2, time=-3)
Dimension._dimensional_dependencies["conductance"] = dict(mass=-1, length=-2, current=2, time=3)
Dimension._dimensional_dependencies["capacitance"] = dict(mass=-1, length=-2, current=2, time=4)
Dimension._dimensional_dependencies["inductance"] = dict(mass=1, length=2, current=-2, time=-2)
Dimension._dimensional_dependencies["charge"] = dict(current=1, time=1)
Dimension._dimensional_dependencies["magnetic_density"] = dict(mass=1, current=-1, time=-2)
Dimension._dimensional_dependencies["magnetic_flux"] = dict(length=2, mass=1, current=-1, time=-2)


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
            return "DimensionSystem(%s)" % ", ".join(str(d) for d in self._base_dims)

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
        if isinstance(dim, string_types):
            dim = Symbol(dim)

        if dim.is_Symbol:
            for d in self._dims:
                if dim in (d.name, d.symbol):
                    found_dim = d
                    break
        elif isinstance(dim, Dimension):
            for i, idim in enumerate(self._dims):
                if dim.get_dimensional_dependencies() == idim.get_dimensional_dependencies():
                    return idim

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
        dimset = set([])
        for i in self._base_dims:
            dimset.update(set(i.get_dimensional_dependencies().keys()))
        return tuple(sorted(dimset))

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

    #@cacheit
    @property
    def can_transf_matrix(self):
        """
        Return the canonical transformation matrix from the canonical to the
        base dimension basis.

        It is the inverse of the matrix computed with inv_can_transf_matrix().
        """

        #TODO: the inversion will fail if the system is inconsistent, for
        #      example if the matrix is not a square
        return reduce(lambda x, y: x.row_join(y),
                      [self.dim_can_vector(d) for d in self._base_dims]
                      ).inv()

    def dim_can_vector(self, dim):
        """
        Dimensional representation in terms of the canonical base dimensions.
        """

        vec = []
        for d in self.list_can_dims:
            vec.append(dim.get_dimensional_dependencies().get(d, 0))

        return Matrix(vec)

    def dim_vector(self, dim):
        """
        Vector representation in terms of the base dimensions.
        """
        return self.can_transf_matrix * Matrix(self.dim_can_vector(dim))

    def print_dim_base(self, dim):
        """
        Give the string expression of a dimension in term of the basis symbols.
        """
        dims = self.dim_vector(dim)
        symbols = [i.symbol if i.symbol is not None else i.name for i in self._base_dims]
        res = S.One
        for (s, p) in zip(symbols, dims):
            res *= s**p
        return res

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
