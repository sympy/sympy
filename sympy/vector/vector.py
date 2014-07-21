from sympy.simplify import simplify as simp, trigsimp as tsimp
from sympy.core.assumptions import StdFactKB
from sympy.core import S, Add, Mul, sympify, Pow, Symbol, count_ops
from sympy.core.decorators import call_highest_priority, _sympifyit
from sympy.core.expr import Expr, AtomicExpr
from sympy import diff as df, sqrt, ImmutableMatrix as Matrix, \
     factor as fctr
from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.vector.functions import express


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

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> v = 3*C.i + 4*C.j + 5*C.k
        >>> v.components
        {C.i: 3, C.j: 4, C.k: 5}

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

    def evalf(self, *args):
        """
        Implements the SymPy evalf routine for vectors.

        evalf's documentation
        =====================

        """
        vec = Vector.zero
        for k, v in self.components.items():
            vec += v.evalf(*args) * k
        return vec
    evalf.__doc__ += Expr.evalf.__doc__

    def simplify(self, ratio=1.7, measure=count_ops):
        """
        Implements the SymPy simplify routine for vectors.

        simplify's documentation
        ========================

        """
        simp_components = [simp(v, ratio, measure) * k for \
                           k, v in self.components.items()]
        return VectorAdd(*simp_components)
    simplify.__doc__ += simp.__doc__

    def trigsimp(self, **opts):
        """
        Implements the SymPy trigsimp routine, for vectors.

        trigsimp's documentation
        ========================

        """
        trig_components = [tsimp(v, **opts) * k for \
                           k, v in self.components.items()]
        return VectorAdd(*trig_components)
    trigsimp.__doc__ += tsimp.__doc__

    def _eval_simplify(self, ratio, measure):
        return self.simplify(ratio, measure)

    def _eval_trigsimp(self, **opts):
        return self.trigsimp(**opts)

    def _eval_derivative(self, wrt):
        return self.diff(wrt)

    def _eval_diff(self, *args, **kwargs):
        return self.diff(*args, **kwargs)

    def as_numer_denom(self):
        """
        Returns the expression as a tuple wrt the following
        transformation -

        expression -> a/b -> a, b

        """
        return (self, 1)

    def factor(self, *args, **kwargs):
        """
        Implements the SymPy factor routine, on the scalar parts
        of a vectorial expression.

        factor's documentation
        ========================

        """
        fctr_components = [fctr(v, *args, **kwargs) * k for \
                           k, v in self.components.items()]
        return VectorAdd(*fctr_components)
    factor.__doc__ += fctr.__doc__

    def magnitude(self):
        """
        Returns the magnitude of this vector.
        """
        return sqrt(self & self)

    def normalize(self):
        """
        Returns the normalized version of this vector.
        """
        return self / self.magnitude()

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

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> C.i.dot(C.j)
        0
        >>> C.i & C.i
        1
        >>> v = 3*C.i + 4*C.j + 5*C.k
        >>> v.dot(C.k)
        5
        >>> (C.i & C.delop)(C.x*C.y*C.z)
        C.y*C.z

        """

        #Check special cases
        from sympy.vector.deloperator import Del
        if not isinstance(other, Vector) and not isinstance(other, Del):
            raise TypeError(str(other)+" is not a vector or del operator")

        #Check if the other is a del operator
        if isinstance(other, Del):
            def directional_derivative(field):
                field = express(field, other.system, variables = True)
                out = self.dot(other._i) * \
                      df(field, other._x)
                out += self.dot(other._j) * \
                       df(field, other._y)
                out += self.dot(other._k) * \
                       df(field, other._z)
                if out == 0 and isinstance(field, Vector):
                    out = Vector.zero
                return out
            return directional_derivative

        if isinstance(self, VectorZero) or isinstance(other, VectorZero):
            return S(0)

        v1 = express(self, other._sys)
        v2 = express(other, other._sys)
        dotproduct = S(0)
        for x in other._sys.base_vectors():
            dotproduct += v1.components.get(x, 0) * \
                          v2.components.get(x, 0)

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

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> C.i.cross(C.j)
        C.k
        >>> C.i ^ C.i
        0
        >>> v = 3*C.i + 4*C.j + 5*C.k
        >>> v ^ C.i
        5*C.j + (-4)*C.k

        """

        #Check special cases
        if not isinstance(other, Vector):
            raise TypeError(str(other) + " is not a vector")
        if self == Vector.zero or other == Vector.zero:
            return Vector.zero

        #Compute cross product
        def _det(mat):
            """This is needed as a little method for to find the determinant
            of a list in python.
            SymPy's Matrix won't take in Vector, so need a custom function.
            The user shouldn't be calling this.

            """

            return (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * \
                                 mat[2][1])
                    + mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] *
                    mat[2][2]) + mat[0][2] * (mat[1][0] * mat[2][1] -
                    mat[1][1] * mat[2][0]))

        outvec = Vector.zero
        for system, vect in other.separate().items():
            tempi = system.i
            tempj = system.j
            tempk = system.k
            tempm = [[tempi, tempj, tempk],
                     [self & tempi, self & tempj, self & tempk],
                     [vect & tempi, vect & tempj, vect & tempk]]
            outvec += _det(tempm)

        return outvec

    def __xor__(self, other):
        return self.cross(other)
    __xor__.__doc__ = cross.__doc__

    def as_coeff_Mul(self, rational=False):
        """Efficiently extract the coefficient of a product. """
        return (S(1), self)

    def as_coeff_add(self, *deps):
        """Efficiently extract the coefficient of a summation. """
        l = [x * self.components[x] for x in self.components]
        return (0, tuple(l))

    def diff(self, *args, **kwargs):
        """
        Implements the SymPy diff routine, for vectors.

        diff's documentation
        ========================

        """
        for x in args:
            if isinstance(x, Vector):
                raise TypeError("Cannot differentiate wrt a Vector")
        diff_components = [df(v, *args, **kwargs) * k for \
                           k, v in self.components.items()]
        return VectorAdd(*diff_components)
    diff.__doc__ += df.__doc__

    def to_matrix(self, system):
        """
        Returns the matrix form of this vector with respect to the
        specified coordinate system.

        Parameters
        ==========

        system : CoordSysCartesian
            The system wrt which the matrix form is to be computed

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> C = CoordSysCartesian('C')
        >>> from sympy.abc import a, b, c
        >>> v = a*C.i + b*C.j + c*C.k
        >>> v.to_matrix(C)
        Matrix([
        [a],
        [b],
        [c]])

        """

        return Matrix([self.dot(unit_vec) for unit_vec in
                       system.base_vectors()])

    def separate(self):
        """
        The constituents of this vector in different coordinate systems,
        as per its definition.

        Returns a dict mapping each CoordSysCartesian to the corresponding
        constituent Vector.

        Examples
        ========

        >>> from sympy.vector import CoordSysCartesian
        >>> R1 = CoordSysCartesian('R1')
        >>> R2 = CoordSysCartesian('R2')
        >>> v = R1.i + R2.i
        >>> v.separate() == {R1: R1.i, R2: R2.i}
        True

        """

        parts = {}
        for vect, measure in self.components.items():
            parts[vect.system] = parts.get(vect.system, Vector.zero) + \
                                 vect*measure
        return parts


class BaseVector(Vector, AtomicExpr):
    """
    Class to denote a base vector.
    """

    def __new__(cls, name, index, system):
        #Verify arguments
        if not index in range(0, 3):
            raise ValueError("index must be 0, 1 or 2")
        if not isinstance(name, str):
            raise TypeError("name must be a valid string")
        if not isinstance(system, CoordSysCartesian):
            raise TypeError("system should be a CoordSysCartesian")
        #Initialize an object
        obj = super(BaseVector, cls).__new__(cls, S(index),
                                             system)
        #Assign important attributes
        obj._base_vect = obj
        obj._components = {obj: S(1)}
        obj._measure_number = S(1)
        obj._name = name
        obj._system = system

        assumptions = {}
        assumptions['commutative'] = True
        obj._assumptions = StdFactKB(assumptions)

        #This attr is used for re-expression to one of the systems
        #involved in the definition of the Vector. Applies to
        #VectorMul and VectorAdd too.
        obj._sys = system

        return obj

    @property
    def system(self):
        return self._system

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
            if arg == Vector.zero:
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
            return Vector.zero

        #Build object
        newargs = [x*components[x] for x in components]
        obj = super(VectorAdd, cls).__new__(cls, *newargs, **options)
        if isinstance(obj, Mul):
            return VectorMul(*obj.args)
        assumptions = {}
        assumptions['commutative'] = True
        obj._assumptions = StdFactKB(assumptions)
        obj._components = components

        obj._sys = (list(components.keys()))[0]._sys

        return obj

    __init__ = Add.__init__

    def __str__(self, printer=None):
        ret_str = ''
        for system, vect in self.separate().items():
            base_vects = system.base_vectors()
            for x in base_vects:
                if x in vect.components:
                    temp_vect = self.components[x]*x
                    ret_str += temp_vect.__str__() + " + "
        return ret_str[:-3]

    __repr__ = __str__
    _sympystr = __str__


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
                vect = arg._base_vect
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
            return Vector.zero

        #If one of the args was a VectorAdd, return an
        #appropriate VectorAdd instance
        if isinstance(vect, VectorAdd):
            newargs = [VectorMul(measure_number, x) for x in vect.args]
            return VectorAdd(*newargs)

        obj = super(VectorMul, cls).__new__(cls, measure_number,
                                            vect._base_vect, **options)
        if isinstance(obj, Add):
            return VectorAdd(*obj.args)
        obj._base_vect = vect._base_vect
        obj._measure_number = measure_number
        assumptions = {}
        assumptions['commutative'] = True
        obj._assumptions = StdFactKB(assumptions)

        obj._components = {vect._base_vect : measure_number}

        obj._sys = vect._base_vect._sys

        return obj

    __init__ = Mul.__init__

    @property
    def base_vector(self):
        """ The BaseVector involved in the product. """
        return self._base_vect

    @property
    def measure_number(self):
        """ The scalar expression involved in the defition of
        this VectorMul.
        """
        return self._measure_number

    def __str__(self, printer=None):
        measure_str = self._measure_number.__str__()
        if '(' in measure_str or '-' in measure_str or \
           '+' in measure_str:
            measure_str = '(' + measure_str + ')'
        return measure_str + '*' + self._base_vect.__str__()

    __repr__ = __str__
    _sympystr = __str__


class VectorZero(Vector):
    """
    Class to denote a zero vector
    """

    _op_priority = 12.1
    components = {}

    def __new__(cls):
        obj = super(VectorZero, cls).__new__(cls)
        #Pre-compute a specific hash value for the zero vector
        #Use the same one always
        obj._hash = tuple([S(0), Symbol('global')]).__hash__()
        return obj

    def __hash__(self):
        return self._hash

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

    def __str__(self, printer=None):
        return '0'
    __repr__ = __str__
    _sympystr = __str__


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


Vector.zero = VectorZero()
