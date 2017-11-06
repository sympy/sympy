"""
Types used to represent a full function/module as an Abstract Syntax Tree.

Most types are small, and are merely used as tokens in the AST. A tree diagram
has been included below to illustrate the relationships between the AST types.


AST Type Tree
-------------
::

  *Basic*
       |--->Assignment
       |             |--->AugmentedAssignment
       |                                    |--->AddAugmentedAssignment
       |                                    |--->SubAugmentedAssignment
       |                                    |--->MulAugmentedAssignment
       |                                    |--->DivAugmentedAssignment
       |                                    |--->ModAugmentedAssignment
       |
       |--->CodeBlock
       |
       |--->For
       |
       |--->Token
       |        |--->Attribute
       |        |--->Type
       |                |--->IntBaseType
       |                |              |--->_SizedIntType
       |                |                               |--->SignedIntType
       |                |                               |--->UnsignedIntType
       |                |--->FloatType
       |                             |--->ComplexType
       |
       |--->Variable
       |           |---> Pointer
       |
       |--->Declaration



Predefined types
----------------
A number of ``Type`` instances are provided in the ``sympy.codegen.ast`` module
for convenience. Perhaps the two most common ones for code-generation (of numeric
codes) are ``float32`` and ``float64`` (known as single and double precision respectively).
There are also precision generic versions of Types (for which the codeprinters selects the
underlying data type at time of printing): ``real``, ``integer``, ``complex_``, ``bool_``.

The other ``Type`` instances defined are:

- ``intc``: Integer type used by C's "int".
- ``intp``: Integer type used by C's "unsigned".
- ``int8``, ``int16``, ``int32``, ``int64``: n-bit integers.
- ``uint8``, ``uint16``, ``uint32``, ``uint64``: n-bit unsigned integers.
- ``float80``: known as "extended precision" on modern x86/amd64 hardware.
- ``complex64``: Complex number represented by two ``float32`` numbers
- ``complex128``: Complex number represented by two ``float64`` numbers

"""



from __future__ import print_function, division


from sympy.core import Symbol, Tuple
from sympy.core.basic import Basic
from sympy.core.numbers import Float, Integer, oo
from sympy.core.relational import Relational
from sympy.core.sympify import _sympify, sympify
from sympy.sets import FiniteSet
from sympy.utilities.iterables import iterable

class Assignment(Relational):
    """
    Represents variable assignment for code generation.

    Parameters
    ----------
    lhs : Expr
        Sympy object representing the lhs of the expression. These should be
        singular objects, such as one would use in writing code. Notable types
        include Symbol, MatrixSymbol, MatrixElement, and Indexed. Types that
        subclass these types are also supported.

    rhs : Expr
        Sympy object representing the rhs of the expression. This can be any
        type, provided its shape corresponds to that of the lhs. For example,
        a Matrix type can be assigned to MatrixSymbol, but not to Symbol, as
        the dimensions will not align.

    Examples
    ========

    >>> from sympy import symbols, MatrixSymbol, Matrix
    >>> from sympy.codegen.ast import Assignment
    >>> x, y, z = symbols('x, y, z')
    >>> Assignment(x, y)
    Assignment(x, y)
    >>> Assignment(x, 0)
    Assignment(x, 0)
    >>> A = MatrixSymbol('A', 1, 3)
    >>> mat = Matrix([x, y, z]).T
    >>> Assignment(A, mat)
    Assignment(A, Matrix([[x, y, z]]))
    >>> Assignment(A[0, 1], x)
    Assignment(A[0, 1], x)
    """

    rel_op = ':='
    __slots__ = []

    def __new__(cls, lhs, rhs=0, **assumptions):
        from sympy.matrices.expressions.matexpr import (
            MatrixElement, MatrixSymbol)
        from sympy.tensor.indexed import Indexed
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
        # Tuple of things that can be on the lhs of an assignment
        assignable = (Symbol, MatrixSymbol, MatrixElement, Indexed)
        if not isinstance(lhs, assignable):
            raise TypeError("Cannot assign to lhs of type %s." % type(lhs))
        # Indexed types implement shape, but don't define it until later. This
        # causes issues in assignment validation. For now, matrices are defined
        # as anything with a shape that is not an Indexed
        lhs_is_mat = hasattr(lhs, 'shape') and not isinstance(lhs, Indexed)
        rhs_is_mat = hasattr(rhs, 'shape') and not isinstance(rhs, Indexed)
        # If lhs and rhs have same structure, then this assignment is ok
        if lhs_is_mat:
            if not rhs_is_mat:
                raise ValueError("Cannot assign a scalar to a matrix.")
            elif lhs.shape != rhs.shape:
                raise ValueError("Dimensions of lhs and rhs don't align.")
        elif rhs_is_mat and not lhs_is_mat:
            raise ValueError("Cannot assign a matrix to a scalar.")
        return Relational.__new__(cls, lhs, rhs, **assumptions)

# XXX: This should be handled better
Relational.ValidRelationOperator[':='] = Assignment


class AugmentedAssignment(Assignment):
    """
    Base class for augmented assignments
    """

    @property
    def rel_op(self):
        return self._symbol + '='

class AddAugmentedAssignment(AugmentedAssignment):
    _symbol = '+'


class SubAugmentedAssignment(AugmentedAssignment):
    _symbol = '-'


class MulAugmentedAssignment(AugmentedAssignment):
    _symbol = '*'


class DivAugmentedAssignment(AugmentedAssignment):
    _symbol = '/'


class ModAugmentedAssignment(AugmentedAssignment):
    _symbol = '%'


Relational.ValidRelationOperator['+='] = AddAugmentedAssignment
Relational.ValidRelationOperator['-='] = SubAugmentedAssignment
Relational.ValidRelationOperator['*='] = MulAugmentedAssignment
Relational.ValidRelationOperator['/='] = DivAugmentedAssignment
Relational.ValidRelationOperator['%='] = ModAugmentedAssignment

def aug_assign(lhs, op, rhs):
    """
    Create 'lhs op= rhs'.

    Represents augmented variable assignment for code generation. This is a
    convenience function. You can also use the AugmentedAssignment classes
    directly, like AddAugmentedAssignment(x, y).

    Parameters
    ----------
    lhs : Expr
        Sympy object representing the lhs of the expression. These should be
        singular objects, such as one would use in writing code. Notable types
        include Symbol, MatrixSymbol, MatrixElement, and Indexed. Types that
        subclass these types are also supported.

    op : str
        Operator (+, -, /, \\*, %).

    rhs : Expr
        Sympy object representing the rhs of the expression. This can be any
        type, provided its shape corresponds to that of the lhs. For example,
        a Matrix type can be assigned to MatrixSymbol, but not to Symbol, as
        the dimensions will not align.

    Examples
    --------

    >>> from sympy import symbols
    >>> from sympy.codegen.ast import aug_assign
    >>> x, y = symbols('x, y')
    >>> aug_assign(x, '+', y)
    AddAugmentedAssignment(x, y)
    """
    if op + '=' not in Relational.ValidRelationOperator:
        raise ValueError("Unrecognized operator %s" % op)
    return Relational.ValidRelationOperator[op + '='](lhs, rhs)


class CodeBlock(Basic):
    """
    Represents a block of code

    For now only assignments are supported. This restriction will be lifted in
    the future.

    Useful methods on this object are

    ``left_hand_sides``:
        Tuple of left-hand sides of assignments, in order.
    ``left_hand_sides``:
        Tuple of right-hand sides of assignments, in order.
    ``topological_sort``:
        Class method. Return a CodeBlock with assignments
        sorted so that variables are assigned before they
        are used.
    ``cse``:
        Return a new CodeBlock with common subexpressions eliminated and
        pulled out as assignments.

    Example
    =======

    >>> from sympy import symbols, ccode
    >>> from sympy.codegen.ast import CodeBlock, Assignment
    >>> x, y = symbols('x y')
    >>> c = CodeBlock(Assignment(x, 1), Assignment(y, x + 1))
    >>> print(ccode(c))
    x = 1;
    y = x + 1;

    """
    def __new__(cls, *args):
        left_hand_sides = []
        right_hand_sides = []
        for i in args:
            if isinstance(i, Assignment):
                lhs, rhs = i.args
                left_hand_sides.append(lhs)
                right_hand_sides.append(rhs)

        obj = Basic.__new__(cls, *args)

        obj.left_hand_sides = Tuple(*left_hand_sides)
        obj.right_hand_sides = Tuple(*right_hand_sides)

        return obj

    @classmethod
    def topological_sort(cls, assignments):
        """
        Return a CodeBlock with topologically sorted assignments so that
        variables are assigned before they are used.

        The existing order of assignments is preserved as much as possible.

        This function assumes that variables are assigned to only once.

        This is a class constructor so that the default constructor for
        CodeBlock can error when variables are used before they are assigned.

        Example
        =======

        >>> from sympy import symbols
        >>> from sympy.codegen.ast import CodeBlock, Assignment
        >>> x, y, z = symbols('x y z')

        >>> assignments = [
        ...     Assignment(x, y + z),
        ...     Assignment(y, z + 1),
        ...     Assignment(z, 2),
        ... ]
        >>> CodeBlock.topological_sort(assignments)
        CodeBlock(Assignment(z, 2), Assignment(y, z + 1), Assignment(x, y + z))

        """
        from sympy.utilities.iterables import topological_sort
        # Create a graph where the nodes are assignments and there is a directed edge
        # between nodes that use a variable and nodes that assign that
        # variable, like

        # [(x := 1, y := x + 1), (x := 1, z := y + z), (y := x + 1, z := y + z)]

        # If we then topologically sort these nodes, they will be in
        # assignment order, like

        # x := 1
        # y := x + 1
        # z := y + z

        # A = The nodes
        #
        # enumerate keeps nodes in the same order they are already in if
        # possible. It will also allow us to handle duplicate assignments to
        # the same variable when those are implemented.
        A = list(enumerate(assignments))

        # var_map = {variable: [assignments using variable]}
        # like {x: [y := x + 1, z := y + x], ...}
        var_map = {}

        # E = Edges in the graph
        E = []
        for i in A:
            if i[1].lhs in var_map:
                E.append((var_map[i[1].lhs], i))
            var_map[i[1].lhs] = i
        for i in A:
            for x in i[1].rhs.free_symbols:
                if x not in var_map:
                    # XXX: Allow this case?
                    raise ValueError("Undefined variable %s" % x)
                E.append((var_map[x], i))

        ordered_assignments = topological_sort([A, E])
        # De-enumerate the result
        return cls(*list(zip(*ordered_assignments))[1])

    def cse(self, symbols=None, optimizations=None, postprocess=None,
        order='canonical'):
        """
        Return a new code block with common subexpressions eliminated

        See the docstring of :func:`sympy.simplify.cse_main.cse` for more
        information.

        Examples
        ========

        >>> from sympy import symbols, sin
        >>> from sympy.codegen.ast import CodeBlock, Assignment
        >>> x, y, z = symbols('x y z')

        >>> c = CodeBlock(
        ...     Assignment(x, 1),
        ...     Assignment(y, sin(x) + 1),
        ...     Assignment(z, sin(x) - 1),
        ... )
        ...
        >>> c.cse()
        CodeBlock(Assignment(x, 1), Assignment(x0, sin(x)), Assignment(y, x0 + 1), Assignment(z, x0 - 1))

        """
        # TODO: Check that the symbols are new
        from sympy.simplify.cse_main import cse

        if not all(isinstance(i, Assignment) for i in self.args):
            # Will support more things later
            raise NotImplementedError("CodeBlock.cse only supports Assignments")

        if any(isinstance(i, AugmentedAssignment) for i in self.args):
            raise NotImplementedError("CodeBlock.cse doesn't yet work with AugmentedAssignments")

        for i, lhs in enumerate(self.left_hand_sides):
            if lhs in self.left_hand_sides[:i]:
                raise NotImplementedError("Duplicate assignments to the same "
                    "variable are not yet supported (%s)" % lhs)

        replacements, reduced_exprs = cse(self.right_hand_sides, symbols=symbols,
            optimizations=optimizations, postprocess=postprocess, order=order)
        assert len(reduced_exprs) == 1
        new_block = tuple(Assignment(var, expr) for var, expr in
            zip(self.left_hand_sides, reduced_exprs[0]))
        new_assignments = tuple(Assignment(*i) for i in replacements)
        return self.topological_sort(new_assignments + new_block)


class For(Basic):
    """Represents a 'for-loop' in the code.

    Expressions are of the form:
        "for target in iter:
            body..."

    Parameters
    ----------
    target : symbol
    iter : iterable
    body : sympy expr

    Examples
    --------

    >>> from sympy import symbols, Range
    >>> from sympy.codegen.ast import aug_assign, For
    >>> x, n = symbols('x n')
    >>> For(n, Range(10), [aug_assign(x, '+', n)])
    For(n, Range(0, 10, 1), CodeBlock(AddAugmentedAssignment(x, n)))

    """

    def __new__(cls, target, iter, body):
        target = _sympify(target)
        if not iterable(iter):
            raise TypeError("iter must be an iterable")
        if isinstance(iter, list):
            # _sympify errors on lists because they are mutable
            iter = tuple(iter)
        iter = _sympify(iter)
        if not isinstance(body, CodeBlock):
            if not iterable(body):
                raise TypeError("body must be an iterable or CodeBlock")
            body = CodeBlock(*(_sympify(i) for i in body))
        return Basic.__new__(cls, target, iter, body)

    @property
    def target(self):
        """
        Return the symbol (target) from the for-loop representation.
        This object changes each iteration.
        Target must be a symbol.
        """
        return self._args[0]

    @property
    def iterable(self):
        """
        Return the iterable from the for-loop representation.
        This is the object that target takes values from.
        Must be an iterable object.
        """
        return self._args[1]

    @property
    def body(self):
        """
        Return the sympy expression (body) from the for-loop representation.
        This is run for each value of target.
        Must be an iterable object or CodeBlock.
        """
        return self._args[2]


class Token(Basic):
    """ Similar to Symbol, but takes no assumptions.

    Defining fields are set in __slots__.
    """

    __slots__ = []

    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and not kwargs and isinstance(args[0], cls):
            return args[0]
        args = args + tuple([kwargs[k] for k in cls.__slots__[len(args):]])
        obj = Basic.__new__(cls)
        for attr, arg in zip(cls.__slots__, args):
            setattr(obj, attr, arg)
        return obj

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        for attr in self.__slots__:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def _hashable_content(self):
        return tuple([getattr(self, attr) for attr in self.__slots__])

    def __hash__(self):
        return super(Token, self).__hash__()

    def _sympystr(self, printer):
        return "{0}({1})".format(self.__class__.__name__, ', '.join(
            ['%s=%s' % (k, printer._print(getattr(self, k))) for k in self.__slots__]
        ))


class Type(Token):
    """ Represents a type.

    The naming is a super-set of NumPy naming, see [1]_. Type has a classmethod
    ``from_expr`` which offer type deduction. It also has a method
    ``cast_check`` which casts the argument to its type, possibly raising an
    exception if rounding error is not within tolerances, or if the value is not
    representable by the underlying data type (e.g. unsigned integers).

    Arguments
    ---------
    name : str
        Name of the type, e.g. ``object``, ``int16``, ``float16`` (where the latter two
        would use the ``Type`` sub-classes ``IntType`` and ``FloatType`` respectively).
        If a ``Type`` instance is given, the said instance is returned.

    Examples
    --------
    >>> from sympy.codegen.ast import Type
    >>> Type.from_expr(42).name
    'integer'
    >>> from sympy.codegen.ast import uint8
    >>> uint8.cast_check(-1)   # doctest: +ELLIPSIS
    Traceback (most recent call last):
      ...
    ValueError: Minimum value for data type bigger than new value.
    >>> from sympy.codegen.ast import float32
    >>> v6 = 0.123456
    >>> float32.cast_check(v6)
    0.123456
    >>> v10 = 12345.67894
    >>> float32.cast_check(v10)  # doctest: +ELLIPSIS
    Traceback (most recent call last):
      ...
    ValueError: Casting gives a significantly different value.
    >>> boost_mp50 = Type('boost::multiprecision::cpp_dec_float_50')
    >>> from sympy import Symbol
    >>> from sympy.printing.cxxcode import cxxcode
    >>> from sympy.codegen.ast import Declaration, Variable
    >>> cxxcode(Declaration(Variable(Symbol('x'), type_=boost_mp50)))
    'boost::multiprecision::cpp_dec_float_50 x'

    References
    ----------

    .. [1] Numpy types
        https://docs.scipy.org/doc/numpy/user/basics.types.html

    """
    __slots__ = ['name']

    default_precision_targets = {}

    def __str__(self):
        return self.name

    @classmethod
    def from_expr(cls, expr):
        """ Deduces type from an expression or a ``Symbol``.

        Parameters
        ----------
        expr : number or SymPy object
            The type will be deduced from type or properties.

        Examples
        --------
        >>> from sympy.codegen.ast import Type, integer, complex_
        >>> Type.from_expr(2) == integer
        True
        >>> from sympy import Symbol
        >>> Type.from_expr(Symbol('z', complex=True)) == complex_
        True
        >>> Type.from_expr(sum)  # doctest: +ELLIPSIS
        Traceback (most recent call last):
          ...
        ValueError: Could not deduce type from expr.

        Raises
        ------
        ValueError when type deduction fails.

        """
        if isinstance(expr, (float, Float)):
            return real
        if isinstance(expr, (int, Integer)) or getattr(expr, 'is_integer', False):
            return integer
        if getattr(expr, 'is_real', False):
            return real
        if isinstance(expr, complex) or getattr(expr, 'is_complex', False):
            return complex_
        if isinstance(expr, bool) or getattr(expr, 'is_Relational', False):
            return bool_
        else:
            raise ValueError("Could not deduce type from expr.")

    def _check(self, value):
        pass

    def cast_check(self, value, rtol=None, atol=0, limits=None, precision_targets=None):
        """ Casts a value to the data type of the instance.

        Parameters
        ----------
        value : number
        rtol : floating point number
            Relative tolerance. (will be deduced if not given).
        atol : floating point number
            Absolute tolerance (in addition to ``rtol``).
        limits : dict
            Values given by ``limits.h``, x86/IEEE754 defaults if not given.
            Default: :attr:`default_limits`.
        type_aliases : dict
            Maps substitutions for Type, e.g. {integer: int64, real: float32}

        Examples
        --------
        >>> from sympy.codegen.ast import Type, integer, float32, int8
        >>> integer.cast_check(3.0) == 3
        True
        >>> float32.cast_check(1e-40)  # doctest: +ELLIPSIS
        Traceback (most recent call last):
          ...
        ValueError: Minimum value for data type bigger than new value.
        >>> int8.cast_check(256)  # doctest: +ELLIPSIS
        Traceback (most recent call last):
          ...
        ValueError: Maximum value for data type smaller than new value.
        >>> v10 = 12345.67894
        >>> float32.cast_check(v10)  # doctest: +ELLIPSIS
        Traceback (most recent call last):
          ...
        ValueError: Casting gives a significantly different value.
        >>> from sympy.codegen.ast import float64
        >>> float64.cast_check(v10)
        12345.67894
        >>> from sympy import Float
        >>> v18 = Float('0.123456789012345646')
        >>> float64.cast_check(v18)
        Traceback (most recent call last):
          ...
        ValueError: Casting gives a significantly different value.
        >>> from sympy.codegen.ast import float80
        >>> float80.cast_check(v18)
        0.123456789012345649

        """
        from sympy.functions.elementary.complexes import im, re
        val = sympify(value)

        ten = Integer(10)
        exp10 = getattr(self, 'decimal_dig', None)

        if rtol is None:
            rtol = 1e-15 if exp10 is None else 2.0*ten**(-exp10)

        def tol(num):
            return atol + rtol*abs(num)

        new_val = self._cast_nocheck(value)
        self._check(new_val)

        delta = new_val - val
        if abs(delta) > tol(val):  # rounding, e.g. int(3.5) != 3.5
            raise ValueError("Casting gives a significantly different value.")

        return new_val


class IntBaseType(Type):
    """ Integer base type, contains no size information. """
    __slots__ = ['name']
    _cast_nocheck = Integer


class _SizedIntType(IntBaseType):
    __slots__ = ['name', 'nbits']

    def _check(self, value):
        if value < self.min:
            raise ValueError("Value is too small: %d < %d" % (value, self.min))
        if value > self.max:
            raise ValueError("Value is too big: %d > %d" % (value, self.max))


class SignedIntType(_SizedIntType):
    @property
    def min(self):
        return -2**(self.nbits-1)

    @property
    def max(self):
        return 2**(self.nbits-1) - 1


class UnsignedIntType(_SizedIntType):
    @property
    def min(self):
        return 0

    @property
    def max(self):
        return 2**self.nbits - 1

two = Integer(2)

class FloatType(Type):
    """ Represents a floating point value.

    Base 2 & one sign bit is assumed.

    Arguments
    ---------
    name : str
        Name of the type.
    nbits : integer
        Number of bits used (storage).
    nmant : integer
        Number of bits used to represent the mantissa.
    nexp : integer
        Number of bits used to represent the mantissa.

    Examples
    --------
    >>> from sympy import S, Float
    >>> from sympy.codegen.ast import FloatType
    >>> half_precision = FloatType('f16', nbits=16, nmant=10, nexp=5)
    >>> half_precision.max
    65504
    >>> half_precision.tiny == S(2)**-14
    True
    >>> half_precision.eps == S(2)**-10
    True
    >>> half_precision.dig == 3
    True
    >>> half_precision.decimal_dig == 5
    True
    >>> half_precision.cast_check(1.0)
    1.0
    >>> half_precision.cast_check(1e5)  # doctest: +ELLIPSIS
    Traceback (most recent call last):
      ...
    ValueError: Maximum value for data type smaller than new value.
    """

    __slots__ = ['name', 'nbits', 'nmant', 'nexp']

    @property
    def max_exponent(self):
        """ The largest positive number n, such that 2**(n - 1) is a representable finite value. """
        # cf. C++'s ``std::numeric_limits::max_exponent``
        return two**(self.nexp - 1)

    @property
    def min_exponent(self):
        """ The lowest negative number n, such that 2**(n - 1) is a valid normalized number. """
        # cf. C++'s ``std::numeric_limits::min_exponent``
        return 3 - self.max_exponent

    @property
    def max(self):
        """ Maximum value representable. """
        return (1 - two**-(self.nmant+1))*two**self.max_exponent

    @property
    def tiny(self):
        """ The minimum positive normalized value. """
        # See C macros: FLT_MIN, DBL_MIN, LDBL_MIN
        # or C++'s ``std::numeric_limits::min``
        # or numpy.finfo(dtype).tiny
        return two**(self.min_exponent - 1)


    @property
    def eps(self):
        """ Difference between 1.0 and the next representable value. """
        return two**(-self.nmant)

    @property
    def dig(self):
        """ Number of decimal digits that are guaranteed to be preserved in text.

        When converting text -> float -> text, you are guaranteed that at least ``dig``
        number of digits are preserved with respect to rounding or overflow.
        """
        from sympy.functions import floor, log
        return floor(self.nmant * log(2)/log(10))

    @property
    def decimal_dig(self):
        """ Number of digits needed to store & load without loss.

        Number of decimal digits needed to guarantee that two consecutive conversions
        (float -> text -> float) to be idempotent. This is useful when one do not want
        to loose precision due to rounding errors when storing a floating point value
        as text.
        """
        from sympy.functions import ceiling, log
        return ceiling((self.nmant + 1) * log(2)/log(10) + 1)

    def _cast_nocheck(self, value):
        return Float(str(sympify(value).evalf(self.decimal_dig)), self.decimal_dig)

    def _check(self, value):
        if value < -self.max:
            raise ValueError("Value is too small: %d < %d" % (value, -self.max))
        if value > self.max:
            raise ValueError("Value is too big: %d > %d" % (value, self.max))
        if abs(value) < self.tiny:
            raise ValueError("Smallest (absolute) value for data type bigger than new value.")


class ComplexType(FloatType):
    """ Represents a complex floating point number. """

    def _cast_nocheck(self, value):
        from sympy.functions import re, im
        return (
            super(ComplexType, self)._cast_nocheck(re(value)) +
            super(ComplexType, self)._cast_nocheck(im(value))*1j
        )

    def _check(self, value):
        from sympy.functions import re, im
        super(ComplexType, self)._check(re(value))
        super(ComplexType, self)._check(im(value))


# NumPy types:
intc = IntBaseType('intc')
intp = IntBaseType('intp')
int8 = SignedIntType('int8', 8)
int16 = SignedIntType('int16', 16)
int32 = SignedIntType('int32', 32)
int64 = SignedIntType('int64', 64)
uint8 = UnsignedIntType('uint8', 8)
uint16 = UnsignedIntType('uint16', 16)
uint32 = UnsignedIntType('uint32', 32)
uint64 = UnsignedIntType('uint64', 64)
float16 = FloatType('float16', 16, nexp=5, nmant=10)  # IEEE 754 binary16, Half precision
float32 = FloatType('float32', 32, nexp=8, nmant=23)  # IEEE 754 binary32, Single precision
float64 = FloatType('float64', 64, nexp=11, nmant=52)  # IEEE 754 binary64, Double precision
float80 = FloatType('float80', 80, nexp=15, nmant=63)  # x86 extended precision (1 integer part bit), "long double"
float128 = FloatType('float128', 128, nexp=15, nmant=112)  # IEEE 754 binary128, Quadruple precision
float256 = FloatType('float256', 256, nexp=19, nmant=236)  # IEEE 754 binary256, Octuple precision

complex64 = ComplexType('complex64', 64, **{k: getattr(float32, k) for k in FloatType.__slots__[2:]})
complex128 = ComplexType('complex128', 128, **{k: getattr(float64, k) for k in FloatType.__slots__[2:]})

# Generic types (precision may be chosen by code printers):
real = Type('real')
integer = IntBaseType('integer')
complex_ = Type('complex')
bool_ = Type('bool')


class Attribute(Token):
    """ Variable attribute """
    __slots__ = ['name']

value_const = Attribute('value_const')
pointer_const = Attribute('pointer_const')

class Variable(Basic):
    """ Represents a variable

    Parameters
    ----------
    symbol : Symbol
    attrs : iterable of Attribute instances
        Will be stored as a FiniteSet.
    type_ : Type (optional)
        Type of the variable.

    Examples
    --------
    >>> from sympy import Symbol
    >>> from sympy.codegen.ast import Variable, float32, integer
    >>> x = Symbol('x')
    >>> v = Variable(x, type_=float32)

    One may also construct a ``Variable`` instance with the type deduced from
    assumptions about the symbol using the ``deduced`` classmethod::
    >>> i = Symbol('i', integer=True)
    >>> v = Variable.deduced(i)
    >>> v.type == integer
    True

    """

    nargs = (2, 3)  # type is optional

    def __new__(cls, symbol, attrs=FiniteSet(), type_=None):
        args = (_sympify(symbol), attrs if isinstance(attrs, FiniteSet) else FiniteSet(*attrs))
        if type_ is not None:
            if not isinstance(type_, Type):
                raise TypeError("type_ argument should be an instance of Type")
            args += (type_,)
        return Basic.__new__(cls, *args)

    @classmethod
    def deduced(cls, symbol, attrs=FiniteSet()):
        """ Alt. constructor with type deduction from ``Type.from_expr``.

        Examples
        --------
        >>> from sympy import Symbol
        >>> from sympy.codegen.ast import Variable, complex_
        >>> n = Symbol('n', integer=True)
        >>> str(Variable.deduced(n).type)
        'integer'
        >>> x = Symbol('x', real=True)
        >>> v = Variable.deduced(x)
        >>> v.type
        Type(name='real')
        >>> z = Symbol('z', complex=True)
        >>> Variable.deduced(z).type == complex_
        True

        """
        return cls(symbol, attrs, Type.from_expr(symbol))

    @property
    def symbol(self):
        return self.args[0]

    @property
    def attributes(self):
        return self.args[1]

    @property
    def type(self):
        if len(self.args) == 3:
            return self.args[2]
        else:
            return None

    @property
    def value_const(self):
        """ Boolean value describing whether the value is constant. """
        return self.attributes.contains(value_const) == True


class Pointer(Variable):
    """ Represents a pointer """

    @property
    def pointer_const(self):
        """ Boolean value describing whether the pointer address is constant. """
        return self.attributes.contains(pointer_const) == True


class Declaration(Basic):
    """ Represents a variable declaration

    Parameters
    ----------
    var : Variable, Pointer or IndexedBase
    val : Value (optional)
        Value to be assigned upon declaration.
    cast : bool
        If val is not ``None`` val will be casted using
        ``var.Type.cast_check()``.

    Examples
    --------
    >>> from sympy import Symbol
    >>> from sympy.codegen.ast import Declaration, Type, Variable, integer
    >>> x = Symbol('x')
    >>> xvar = Variable(x)
    >>> decl = Declaration.deduced(xvar, 3)
    >>> decl.variable.type == integer
    True
    >>> k = Symbol('k', integer=True)
    >>> k_decl = Declaration.deduced(k, 3.0)
    >>> k_decl.variable.type == integer
    True

    """

    nargs = (1, 2)

    def __new__(cls, var, val=None, cast=False):
        if not isinstance(var, Variable):
            raise TypeError("var argument should be an instance of Variable")
        args = var,
        if val is not None:
            if cast:
                args += (var.type.cast_check(val),)
            else:
                args += (_sympify(val),)
        return Basic.__new__(cls, *args)

    @classmethod
    def deduced(cls, symbol, value=None, attrs=FiniteSet(), **kwargs):
        """ Deduces type primarily from ``symbol``, secondarily from ``value``.

        Examples
        --------
        >>> from sympy import Symbol
        >>> from sympy.codegen.ast import Declaration, real, integer
        >>> x = Symbol('x', real=True)
        >>> decl = Declaration.deduced(x)
        >>> decl.variable.type == real
        True
        >>> decl.value is None
        True
        >>> n = Symbol('n', integer=True)
        >>> Declaration.deduced(n).variable
        Variable(n, EmptySet(), IntBaseType(name='integer'))

        """
        try:
            type_ = Type.from_expr(symbol)
        except ValueError:
            type_ = Type.from_expr(value)
        var = Variable(symbol, attrs, type_)
        return cls(var, value, **kwargs)

    @property
    def variable(self):
        """ Variable of the declaration """
        return self.args[0]

    @property
    def value(self):
        """ Initialization value of the declaration """
        if len(self.args) == 2:
            return self.args[1]
        else:
            return None
