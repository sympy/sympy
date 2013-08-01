#Note - This is a WIP
#For the time being, Symbols have been used in place of BaseVectors
#to test whether the methods behave as expected.
#Some example usage at - http://pastebin.com/WxJ15ewk
from sympy import diff
from sympy.core import Expr, Mul, Add, Pow, S, sympify
from sympy.core.decorators import call_highest_priority
from sympy.vector import BaseVector, VectAdd, VectMul, Vector
from sympy.physics.mechanics import _check_vector, _check_frame

#TODO - Add helper function to calculate outer product of two vectors
#TODO - Add case for time-differentiation of dyadics to frame.dt
#TODO - change docstrings and update doc examples as per new API


class Dyadic(Expr):
    """
    Dyadic super class
    
    See:
    http://en.wikipedia.org/wiki/Dyadic_tensor
    Kane, T., Levinson, D. Dynamics Theory and Applications. 1985 McGraw-Hill

    A more powerful way to represent a rigid body's inertia. While it is more
    complex, by choosing Dyadic components to be in body fixed basis vectors,
    the resulting matrix is equivalent to the inertia tensor.
    """
    
    _op_priority = 12.0

    @property
    def components(self):
        """
        The components of this dyadic in the form of a dict of
        BaseDyadic : measure number pairs
        """
        pass

    def __neg__(self):
        return DyadicMul(S(-1), self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return DyadicAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return DyadicAdd(other, self)

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return DyadicAdd(self, -other)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return DyadicAdd(other, -self)
    
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return DyadicMul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return DyadicMul(other, self)

    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return _dyad_div(self, other)

    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise TypeError("Cannot divide by dyadic")

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def evalf(self, *args):
        return self

    def simplify(self):
        """ Returns simplified version of the dyadic """
        simplify_components = {}
        for x in self.components:
            simplify_components[x] = simplify(self.components[x])
        simplify_components = [x * simplify_components[x] \
                               for x in simplify_components]
        return DyadicAdd(*simplify_components)

    def _eval_simplify(self, ratio, measure):
        return self.simplify()
    
    def factor(self, *args, **kwargs):
        raise TypeError("Factoring not supported for dyadics")

    def dot(self, other):
        """The inner product operator for a Dyadic and a Dyadic or Vector.

        Parameters
        ==========

        other : Dyadic or Vector
            The other Dyadic or Vector to take the inner product with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer
        >>> N = ReferenceFrame('N')
        >>> D1 = outer(N.x, N.y)
        >>> D2 = outer(N.y, N.y)
        >>> D1.dot(D2)
        (N.x|N.y)
        >>> D1.dot(N.y)
        N.x

        """
        
        if other == 0:
            return S(0)
        elif isinstance(other, Symbol):
            outvec = 0
            for x in self.components:
                vect_dot = Symbol("(" + str(x.args[1]) + "." + str(other) + ")")
                outvec += vect_dot * self.components[x] * x.args[0]
            return outvec
        elif isinstance(other, Dyadic):
            outdyad = 0
            for x in self.components:
                for y in other.components:
                    vect_dot = Symbol("(" + str(x.args[1]) + "." + str(y.args[0]) + ")")
                    outer = BaseDyadic(x.args[0], y.args[1])
                    outdyad += vect_dot * self.components[x] * \
                               other.components[y] * outer
            return outdyad
        else:
            raise TypeError(str(type(other)) + " not supported for " + \
                            "dot with dyadics")

    def rdot(self, other):
        """The inner product operator for a Vector or Dyadic, and a Dyadic

        This is for: Vector dot Dyadic

        Parameters
        ==========

        other : Vector
            The vector we are dotting with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, dot, outer
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> dot(N.x, d)
        N.x

        """
        
        if other == 0:
            return S(0)
        elif isinstance(other, Symbol):
            outvec = 0
            for x in self.components:
                vect_dot = Symbol("(" + str(x.args[0]) + "." + str(other) + ")")
                outvec += vect_dot * self.components[x] * x.args[1]
            return outvec
        else:
            raise TypeError(str(type(other)) + " not supported for " + \
                            "r-dot with dyadics")

    def cross(self, other):
        """For a cross product in the form: Dyadic x Vector.

        Parameters
        ==========

        other : Vector
            The Vector that we are crossing this Dyadic with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer, cross
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> cross(d, N.y)
        (N.x|N.z)

        """
        
        if other == 0:
            return S(0)
        elif isinstance(other, Symbol):
            outdyad = S(0)
            for x in self.components:
                cross = Symbol(str(x.args[1])+ "^" + str(other))
                outer = BaseDyadic(x.args[0], cross)
                outdyad += self.components[x] * outer
            return outdyad
        else:
            raise TypeError(str(type(other)) + " not supported for " + \
                            "cross with dyadics")

    def rcross(self, other):
        """For a cross product in the form: Vector x Dyadic

        Parameters
        ==========

        other : Vector
            The Vector that we are crossing this Dyadic with

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer, cross
        >>> N = ReferenceFrame('N')
        >>> d = outer(N.x, N.x)
        >>> cross(N.y, d)
        - (N.z|N.x)

        """
        
        if other == 0:
            return S(0)
        elif isinstance(other, Symbol):
            outdyad = S(0)
            for x in self.components:
                cross = Symbol(str(other)+ "^" + str(x.args[0]))
                outer = BaseDyadic(cross, x.args[1])
                outdyad += self.components[x] * outer
            return outdyad
        else:
            raise TypeError(str(type(other)) + " not supported for " + \
                            "r-cross with dyadics")

    def dt(self, frame):
        """Take the time derivative of this Dyadic in a frame.

        Parameters
        ==========

        frame : ReferenceFrame
            The frame to take the time derivative in

        Examples
        ========

        >>> from sympy.physics.mechanics import ReferenceFrame, outer, dynamicsymbols
        >>> N = ReferenceFrame('N')
        >>> q = dynamicsymbols('q')
        >>> B = N.orientnew('B', 'Axis', [q, N.z])
        >>> d = outer(N.x, N.x)
        >>> d.dt(B)
        - q'*(N.y|N.x) - q'*(N.x|N.y)

        """
        
        _check_frame(frame)
        outdyad = S(0)
        for x in self.components:
            measure = self.components[x]
            outer1 = BaseDyadic(frame.dt(x.args[0]), x.args[1])
            outer2 = BaseDyadic(x.args[0], frame.dt(x.args[1]))
            outdyad += diff(measure, frame.time) * x
            outdyad += measure * outer1
            outdyad += measure * outer2
        return outdyad

    def doit(self, **hints):
        """Calls .doit() on each term in the Dyadic"""
        return sum([self.components[x].doit(**hints) * x for
                    x in self.components])


class BaseDyadic(Dyadic):
    """ Class to denote a basic dyadic tensor component """
    def __new__(cls, vector1, vector2):
        if not isinstance(vector1, Symbol) or not isinstance(vector2, Symbol):
            raise TypeError("BaseDyadic cannot be composed of non-base "+ \
                            "vectors")
        elif vector1 == 0 or vector2 == 0:
            return S(0)
        obj = super(BaseDyadic, cls).__new__(cls, vector1, vector2)
        obj._base_dyad = obj
        obj._measure_number = 1
        return obj

    @property
    def components(self):
        return {self : 1}

    def __str__(self, printer=None):
        return "(" + str(self.args[0]) + "|" + str(self.args[1]) + ")"

    _sympystr = __str__
    _sympyrepr = _sympystr


class DyadicMul(Mul, Dyadic):
    """ Products of scalars and BaseDyadics """

    def __new__(cls, *args, **options):
        count = 0
        measure_number = 1
        for x in args:
            if x == 0:
                return S(0)
            elif isinstance(x, BaseDyadic) or isinstance(x, DyadicMul):
                count += 1
                dyad = x
                measure_number *= x._measure_number
            elif isinstance(x, DyadicAdd):
                count += 1
                dyad = x
            else:
                measure_number *= x
        if count > 1:
            raise ValueError("Cannot multiply two dyadics")
        elif count == 0:
            raise ValueError("No dyadics supplied")
        if isinstance(dyad, DyadicAdd):
            newargs = [DyadicMul(measure_number, x) for x in dyad.args]
            return DyadicAdd(*newargs)
        obj = super(DyadicMul, cls).__new__(cls, *args, **options)
        obj._base_dyad = dyad._base_dyad
        obj._measure_number = measure_number
        return obj

    def as_coeff_Mul(self, rational=False):
        return (S(1), self)

    __init__ = Mul.__init__

    @property
    def components(self):
        return {self._base_dyad : self._measure_number}


class DyadicAdd(Dyadic, Add):
    """ Class to hold dyadic sums """

    def __new__(cls, *args, **options):
        lookup = {}
        for x in args:
            if x == 0:
                continue
            for y in x.components:
                if y not in lookup:
                    lookup[y] = 0
                lookup[y] += x.components[y]
        newargs = [DyadicMul(x, lookup[x]) for x in lookup]
        obj = super(DyadicAdd, cls).__new__(cls, *newargs, **options)
        if isinstance(obj, Mul):
            return DyadicMul(*obj.args)
        obj._components = lookup
        return obj

    __init__ = Add.__init__

    @property
    def components(self):
        return self._components


def _dyad_div(one, other):
    """ Helper for division involving dyadics """
    if isinstance(one, Dyadic) and isinstance(other, Dyadic):
        raise TypeError("Cannot divide two dyadics")
    elif isinstance(one, Dyadic):
        return DyadicMul(one, Pow(other, S.NegativeOne))
    else:
        raise TypeError("Cannot divide by a dyadic")

def _outer(vector1, vector2):
    """ Returns the outer product of two vectors """
    #Separate the two vectors into components of base vectors, measure nos
    #Iterate over above components(one inside another) and keep adding to
    #outdyad

