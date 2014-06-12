#TODOs: do stress-tests
from sympy.simplify import simplify, trigsimp
from sympy.core.assumptions import StdFactKB
from sympy.core.symbol import Dummy
from sympy.core import S, Add, Mul, sympify, Pow, Symbol, count_ops
from sympy.core.decorators import call_highest_priority, _sympifyit
from sympy.core.expr import Expr
from sympy.core.numbers import Zero
from sympy import diff, sqrt, ImmutableMatrix as Matrix


class Vector(Expr):
    """
    Super class for all Vector classes.
    Ideally, neither this class nor any of its subclasses should be
    instantiated by the user.
    """

    is_Vector = True
    _op_priority = 12.0

    @property
    def components(self):
        """
        Returns the components of this vector in the form of a
        Python dictionary mapping BaseVector instances to the
        corresponding measure numbers.

        Examples
        ========

        >>> from sympy.vector import i, j, k
        >>> v = 3*i + 4*j + 5*k
        >>> v.components
        {i: 3, j: 4, k: 5}

        """
        #The '_components' attribute is defined according to the
        #subclass of Vector the instance belongs to.
        return self._components

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return VectorAdd(self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return VectorAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return VectorAdd(self, -other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return VectorAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return VectorMul(self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return VectorMul(other, self)

    def __neg__(self):
        return VectorMul(S(-1), self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return _vect_div(self, other)
    
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        return TypeError("Cannot divide by a vector")

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    @call_highest_priority('__req__')
    def __eq__(self, other):
        if isinstance(other, Vector):
            return self.components == other.components
        return False

    __req__ = __eq__

    def evalf(self, *args):
        v = Vector.Zero
        for x in self.components:
            v += self.components[x].evalf(*args) * x
        return v

    def simplify(self, ratio=1.7, measure=count_ops):
        """ Returns a simplified version of this vector """
        return _func_over_vector(self, simplify, (ratio, measure))

    def trigsimp(self, *args):
        """ trigsimp method for vectors """
        return _func_over_vector(self, trigsimp, args)

    def _eval_simplify(self, ratio, measure):
        return self.simplify(ratio, measure)

    def _eval_derivative(self, wrt):
        return self.diff(wrt)

    def as_numer_denom(self):
        return (self, 1)

    def factor(self, *args, **kwargs):
        raise TypeError("Factoring not applicable for Vectors")

    def magnitude(self):
        """
        Returns the magnitude of this vector.
        """
        return sqrt(self & self)

    def normalize(self):
        """
        Returns the normalized version of this vector.
        """
        return self / sqrt(self & self)

    def dot(self, other):
        """
        Returns the dot product of this Vector, either with another
        Vector, or the Del operator.
        If 'other' is a Vector, returns the dot product scalar (Sympy
        expression).
        If 'other' is an instance of Del, returns the directional
        derivate operator as a Python function. If this function is
        applied to a scalar expression, it returns the directional
        derivative of the scalar field wrt this Vector.

        Parameters
        ==========

        other: Vector/Del
            The Vector we are dotting with, or a Del operator .

        Examples
        ========

        >>> from sympy.vector import i, j, k
        >>> i.dot(j)
        0
        >>> i & i
        1
        >>> v = 3*i + 4*j + 5*k
        >>> v.dot(k)
        5

        """

        #Check special cases
        from sympy.vector.deloperator import Del
        if not isinstance(other, Vector) and not isinstance(other, Del):
            raise TypeError(str(other) + " is not a vector or del operator")

        #Check if the other is a del operator
        if isinstance(other, Del):
            vect = self
            def directional_derivative(field):
                out = vect.dot(other._i) * \
                      diff(field, other._x)
                out += vect.dot(other._j) * \
                       diff(field, other._y)
                out += vect.dot(other._k) * \
                       diff(field, other._z)
                if out == 0 and isinstance(field, Vector):
                    out = Vector.Zero
                return out
            return directional_derivative

        dotproduct = S(0)
        for x in self.components:
            dotproduct += self.components.get(x, 0) * \
                          other.components.get(x, 0)

        return dotproduct

    def __and__(self, other):
        return self.dot(other)
    __and__.__doc__ = dot.__doc__
    __rand__ = __and__

    def cross(self, other):
        """
        Returns the cross product of this Vector with another.
        The cross product is a Vector.

        Parameters
        ==========

        other: Vector
            The Vector we are crossing with.

        Examples
        ========

        >>> from sympy.vector import i, j, k
        >>> i.cross(j)
        k
        >>> i ^ i
        0
        >>> v = 3*i + 4*j + 5*k
        >>> v ^ i
        5*j + (-4)*k

        """

        #Check special cases
        if not isinstance(other, Vector):
            raise TypeError(str(other) + " is not a vector")
        if self == Vector.Zero or other == Vector.Zero:
            return Vector.Zero

        #Compute cross product

        vx = self.components.get(j, 0)*other.components.get(k, 0) - \
             self.components.get(k, 0)*other.components.get(j, 0)
        vy = self.components.get(k, 0)*other.components.get(i, 0) - \
             self.components.get(i, 0)*other.components.get(k, 0)
        vz = self.components.get(i, 0)*other.components.get(j, 0) - \
             self.components.get(j, 0)*other.components.get(i, 0)

        return vx*i + vy*j + vz*k

    def __xor__(self, other):
        return self.cross(other)
    __xor__.__doc__ = cross.__doc__

    def as_coeff_Mul(self, rational=False):
        return (S(1), self)

    def as_coeff_add(self, *deps):
        l = [x * self.components[x] for x in self.components]
        return (0, tuple(l))

    def diff(self, *args, **kwargs):
        for x in args:
            if isinstance(x, Vector):
                raise TypeError("Cannot differentiate wrt a Vector")
        diff_components = {}
        for x in self.components:
            diff_components[x] = diff(self.components[x], *args, **kwargs)
        diff_components = [x * diff_components[x] for x in \
                               diff_components]
        return VectorAdd(*diff_components)

    def to_matrix(self):
        """
        Returns the matrix form of this vector.

        Examples
        ========

        >>> from sympy.vector import i, j, k
        >>> from sympy.abc import a, b, c
        >>> v = a*i + b*j + c*k
        >>> v.to_matrix()
        Matrix([
        [a],
        [b],
        [c]])

        """

        base_vectors = (i, j, k)
        return Matrix([self.dot(unit_vec) for unit_vec in
                       base_vectors]).reshape(3, 1)


class BaseVector(Vector, Dummy):
    """
    Class to denote a base vector.
    """

    is_Atom = True

    def __new__(cls, name, index):
        #Verify arguments
        if not index in range(0, 3):
            raise ValueError("index must be 0, 1 or 2")
        if not isinstance(name, str):
            raise ValueError("name must be a valid string")
        #Initialize an object
        obj = super(BaseVector, cls).__new__(cls, name)
        obj._index = index
        obj._components = {obj : S(1)}
        obj._base_vect = obj
        obj._measure_number = 1
        obj._name = name

        assumptions = {}
        assumptions['commutative'] = True
        obj._assumptions = StdFactKB(assumptions)

        return obj

    @call_highest_priority('__req__')
    def __eq__(self, other):
        if isinstance(other, BaseVector):
            return self._index == other._index
        return False

    __req__ = __eq__

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__


class VectorAdd(Vector, Add):
    """
    Class to denote sum of Vector instances.
    """

    def __new__(cls, *args, **options):
        components = {}

        #Check each arg and simultaneously learn the components
        for i, arg in enumerate(args):
            if not isinstance(arg, Vector):
                if isinstance(arg, Mul):
                    arg = VectorMul(*(arg.args))
                elif isinstance(arg, Add):
                    arg = VectorAdd(*(arg.args))
                else:
                    raise TypeError(str(arg) +
                                    " cannot be interpreted as a vector")
            #If argument is zero, ignore
            if arg == Vector.Zero:
                continue
            #Else, update components accordingly
            for x in arg.components:
                components[x] = components.get(x, 0) + arg.components[x]

        temp = components.keys()
        for x in temp:
            if components[x] == 0:
                del components[x]

        #Handle case of zero vector
        if len(components) == 0:
            return Vector.Zero

        newargs = [x*components[x] for x in components]

        #Build object
        obj = super(VectorAdd, cls).__new__(cls, *newargs, **options)
        if isinstance(obj, Mul):
            return VectorMul(*obj.args)
        assumptions = {}
        assumptions['commutative'] = True
        obj._assumptions = StdFactKB(assumptions)
        obj._components = components

        return obj

    __init__ = Add.__init__
        

class VectorMul(Vector, Mul):
    """
    Class to denote products of scalars and BaseVectors.
    """

    def __new__(cls, *args, **options):
        count = 0
        measure_number = S(1)
        zeroflag = False

        #Determine the component and check arguments
        #Also keep a count to ensure two vectors aren't
        #being multipled
        for arg in args:
            if isinstance(arg, VectorZero):
                count += 1
                zeroflag = True
            elif arg == S(0):
                zeroflag = True
            elif isinstance(arg, BaseVector) or \
                 isinstance(arg, VectorMul):
                count += 1
                try:
                    vect = arg._base_vect
                except:
                    print arg
                measure_number *= arg._measure_number
            elif isinstance(arg, VectorAdd):
                count += 1
                vect = arg
            else:
                measure_number *= arg
        #Make sure incompatible types weren't multipled
        if count > 1:
            raise ValueError("Cannot multiply one vector with another")
        elif count == 0:
            return Mul(*args, **options)
        #Handle zero vector case
        if zeroflag:
            return Vector.Zero

        #If one of the args was a VectorAdd, return an
        #appropriate VectorAdd instance
        if isinstance(vect, VectorAdd):
            newargs = [VectorMul(measure_number, x) for x in vect.args]
            return VectorAdd(*newargs)

        obj = super(VectorMul, cls).__new__(cls, *args, **options)
        if isinstance(obj, Add):
            return VectorAdd(*obj.args)
        obj._base_vect = vect._base_vect
        obj._measure_number = measure_number
        assumptions = {}
        assumptions['commutative'] = True
        obj._assumptions = StdFactKB(assumptions)
        
        obj._components = {vect._base_vect : measure_number}

        return obj
        
    __init__ = Mul.__init__

    def __str__(self, printer=None):
        measure_str = self._measure_number.__str__()
        if '(' in measure_str or '-' in measure_str or \
           '+' in measure_str:
            measure_str = '(' + measure_str + ')'
        return measure_str + '*' + self._base_vect.__str__()

    __repr__ = __str__
    _sympystr = __str__


class VectorZero(Vector, Zero):
    """
    Class to denote a zero vector
    """

    _op_priority = 12.1
    components = {}

    @call_highest_priority('__req__')
    def __eq__(self, other):
        return isinstance(other, VectorZero)

    __req__ = __eq__

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        if isinstance(other, Vector):
            return other
        else:
            raise TypeError(str(other) + " is not a vector")

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        if isinstance(other, Vector):
            return other
        else:
            raise TypeError(str(other) + " is not a vector")

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        if isinstance(other, Vector):
            return -other
        else:
            raise TypeError(str(other) + " is not a vector")

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        if isinstance(other, Vector):
            return other
        else:
            raise TypeError(str(other) + " is not a vector")

    def __neg__(self):
        return self

    def normalize(self):
        """
        Returns the normalized version of this vector.
        """
        return self


def _vect_div(one, other):
    """ Helper for division involving vectors. """
    if isinstance(one, Vector) and isinstance(other, Vector):
        raise TypeError("Cannot divide two vectors")
    elif isinstance(one, Vector):
        if other == S.Zero:
            raise ValueError("Cannot divide a vector by zero")
        return VectorMul(one, Pow(other, S.NegativeOne))
    else:
        raise TypeError("Invalid division involving a vector")


def _func_over_vector(vect, func, args):
    """ Applies a function over all components of a vector. """
    func_components = {}
    for x in vect.components:
        func_components[x] = func(vect.components[x], *args)
    func_components = [x * func_components[x] for x in \
                       func_components]
    return VectorAdd(*func_components)


Vector.Zero = VectorZero()

#Just some hacks for now
i = BaseVector('i', 0)
j = BaseVector('j', 1)
k = BaseVector('k', 2)
