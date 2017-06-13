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
        Operator (+, -, /, \*, %).

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


class Type(Symbol):
    """ Represents a type.

    The naming is a super-set of NumPy naming, see [1]_. Type has a classmethod
    ``from_expr`` which offer type deduction. It also has a method
    ``cast_check`` which casts the argument to its type, possibly raising an
    exception if possible rounding error is not within tolerances.

    Arguments
    ---------
    name : str
        Either an explicit type: ``intc``, ``intp``, ``int8``, ``int16``,
        ``int32``, ``int64``, ``uint8``, ``uint16``, ``uint32``, ``uint64,
        float16``, ``float32``, ``float64``, ``complex64``, ``complex128``,
        ``bool``. Or a type category (precision decided by code-printer): ``integer``,
        ``real`` or ``complex`` (where the latter two are of floating point type).
        If a ``Type`` instance is given, the said instance is returned.

    Examples
    --------
    >>> from sympy.codegen.ast import Type
    >>> Type.from_expr(42).name
    'integer'
    >>> f32 = Type("float32")
    >>> v6 = 0.123456
    >>> f32.cast_check(v6)
    0.123456
    >>> v10 = 12345.67894
    >>> f32.cast_check(v10)  # doctest: +ELLIPSIS
    Traceback (most recent call last):
      ...
    ValueError: Casting gives a significantly different value.
    >>> f64 = Type('float64')
    >>> f64.cast_check(v10)
    12345.67894
    >>> from sympy import Float
    >>> v18 = Float('0.123456789012345646')
    >>> f64.cast_check(v18)
    Traceback (most recent call last):
      ...
    ValueError: Casting gives a significantly different value.
    >>> Type('float80').cast_check(v18)
    0.123456789012345649
    >>> boost_mp50 = Type('boost::multiprecision::cpp_dec_float_50')
    >>> from sympy import Symbol
    >>> from sympy.printing.cxxcode import cxxcode
    >>> from sympy.codegen.ast import Declaration, Variable
    >>> cxxcode(Declaration(Variable(Symbol('x'), None, boost_mp50)))
    'boost::multiprecision::cpp_dec_float_50 x;'

    References
    ----------

    .. [1] Numpy types
        https://docs.scipy.org/doc/numpy/user/basics.types.html

    """
    __slots__ = ['name']

    default_limits = {  # IEE754 data when applicable
        'int32': {  # usually the fastest integer typ (and hence 'int' on most systems)
            'max': 2**31 - 1,  # INT_MAX
            'min': -2**31  # INT_MIN
        },
        'float32': {
            'max': 3.40282347e+38,  # FLT_MAX
            'tiny': 1.17549435e-38,  # FLT_MIN, excluding subnormal numbers
            'eps': 1.1920929e-07,  # FLT_EPSILON
            'dig': 6,  # FLT_DIG
            'decimal_dig': 9,  # FLT_DECIMAL_DIG
        },
        'float64': {
            'max': 1.79769313486231571e+308,  # DBL_MAX
            'tiny': 2.22507385850720138e-308,  # DBL_MIN, excluding subnormal numbers
            'eps': 2.22044604925031308e-16,  # DBL_EPSILON
            'dig': 15,  # DBL_DIG
            'decimal_dig': 17,  # DBL_DECIMAL_DIG
        },
        'float80': {  # extended precision, usually "long double", even numpy.float128 (!)
            'max': 1.18973149535723176502e+4932,  # LDBL_MAX
            'tiny': 3.36210314311209350626e-4932,  # LDBL_MIN, excluding subnormal numbers
            'eps': 1.08420217248550443401e-19,  # LDBL_EPSILON
            'dig': 18,  # LDBL_DIG
            'decimal_dig': 21,  # LDBL_DECIMAL_DIG
        }
    }
    default_limits['complex64'] = default_limits['float32']
    default_limits['complex128'] = default_limits['float64']

    default_precision_targets = {}

    def __new__(cls, name):
        if isinstance(name, Type):
            return name
        return Symbol.__new__(cls, name)

    @classmethod
    def from_expr(cls, expr):
        """ Deduces type from an expression or a ``Symbol``.

        Parameters
        ----------
        expr : number or SymPy object
            The type will be deduced from type or properties.

        Examples
        --------
        >>> from sympy.codegen.ast import Type
        >>> Type.from_expr(2) == Type('integer')
        True
        >>> from sympy import Symbol
        >>> Type.from_expr(Symbol('z', complex=True)) == Type('complex')
        True

        Raises
        ------
        ValueError when type deduction fails.

        """
        if isinstance(expr, (float, Float)):
            return cls('real')
        if isinstance(expr, (int, Integer)) or getattr(expr, 'is_integer', False):
            return cls('integer')
        if getattr(expr, 'is_real', False):
            return cls('real')
        if isinstance(expr, complex) or getattr(expr, 'is_complex', False):
            return cls('complex')
        if isinstance(expr, bool) or getattr(expr, 'is_Relational', False):
            return cls('bool')
        else:
            raise ValueError("Could not deduce type from expr")

    def cast_check(self, value, rtol=None, atol=None, limits=None, precision_targets=None):
        """ Casts a value to the data type of the instance.

        Parameters
        ----------
        value : number
        rtol : floating point number
            Relative tolerance. (will be deduced if not given).
        atol : floating point number
            Absolute tolerance. (will be deduced if not given).
        limits : dict
            Values given by ``limits.h``, x86/IEEE754 defaults if not given.
            Default: :attr:`default_limits`.
        precision_targets : dict
            Maps substitutions for Type.name, e.g. {'integer': 'int64', 'real': 'float32'}

        Examples
        --------
        >>> from sympy.codegen.ast import Type
        >>> Type('integer').cast_check(3.0) == 3
        True
        >>> Type('float32').cast_check(1e-40)  # doctest: +ELLIPSIS
        Traceback (most recent call last):
          ...
        ValueError: Minimum value for data type bigger than new value.

        """
        from sympy.functions.elementary.complexes import im, re
        def lim(type_name, key):
            return (limits or self.default_limits).get(type_name, {}).get(key)
        val = _sympify(value)

        name = (precision_targets or self.default_precision_targets).get(self.name, self.name)
        ten = Integer(10)
        if rtol is None:
            exp10 = lim(name, 'decimal_dig')
            rtol = 1e-15 if exp10 is None else ten**(-exp10)

        if atol is None:
            if rtol == 0:
                exp10 = lim(name, 'decimal_dig')
                atol = 1e-15 if exp10 is None else 10**(-exp10)
            else:
                atol = 0

        def tol(num):
            return atol + rtol*abs(num)

        caster = lambda x: x  # identity
        _min, _max, _tiny = -oo, oo, 0  # undefined precision

        if name == 'integer':
            caster = int
        elif name.startswith('int'):
            nbits = int(name.split('int')[1])
            _min, _max = -2**(nbits - 1), 2**(nbits - 1) - 1
            caster = int
        elif name.startswith('uint'):
            nbits = int(name.split('uint')[1])
            _min, _max = 0, 2**nbits - 1
            caster = int
        elif name.startswith('float'):
            _max = +lim(name, 'max')
            _min = -lim(name, 'max')
            _tiny = lim(name, 'tiny')
            dec_dig = lim(name, 'decimal_dig')
            try:
                val = Float(str(val), dec_dig+3)
            except ValueError:
                val = val.evalf(dec_dig + 3)  # e.g. sympy.pi
            caster = lambda x: Float(str(x.evalf(dec_dig)), dec_dig+3)

        if name.startswith('complex'):
            if name == 'complex':
                nbits = 128  # assume double precision as underlying data
            else:
                nbits = int(name.split('complex')[1])
            corresponding_float = 'float%d' % (nbits // 2)
            _max = +lim(corresponding_float, 'max')
            _tiny = lim(corresponding_float, 'tiny')
            dec_dig = lim(corresponding_float, 'decimal_dig')
            caster = lambda x: (
                Float(str(re(x).evalf(dec_dig)), dec_dig+3) +
                Float(str(im(x).evalf(dec_dig)), dec_dig+3)*1j
            )
            if abs(re(val)) > _max or abs(im(val)) > _max:
                raise ValueError("Maximum value exceeded for data type.")
            if abs(re(val)) < _tiny or abs(im(val)) < _tiny:
                raise ValueError("Minimum (absolute) value for data type bigger than new value.")
            new_val = caster(re(val)) + 1j*caster(im(val))
            delta = new_val - val
            if abs(re(delta)) > tol(re(val)) or abs(im(delta)) > tol(im(val)):
                raise ValueError("Casting gives a significantly different value.")
        else:
            if val < _min:
                raise ValueError("Minumum value for data type bigger than new value.")
            if val > _max:
                raise ValueError("Maximum value exceeded for data type.")
            if abs(val) < _tiny:
                raise ValueError("Smallest (absolute) value for data type bigger than new value.")
            new_val = caster(val)
            delta = new_val - val
            if abs(delta) > tol(val):  # rounding, e.g. int(3.5) != 3.5
                raise ValueError("Casting gives a significantly different value.")

        return new_val

# NumPy types:
intc = Type('intc')
intp = Type('intp')
int8 = Type('int8')
int16 = Type('int16')
int32 = Type('int32')
int64 = Type('int64')
uint8 = Type('uint8')
uint16 = Type('uint16')
uint32 = Type('uint32')
uint64 = Type('uint64')
float16 = Type('float16')
float32 = Type('float32')
float64 = Type('float64')
float80 = Type('float80')
complex64 = Type('complex64')
complex128 = Type('complex128')
# Generic types (precision may be chosen by code printers):
real = Type('real')
integer = Type('integer')
complex_ = Type('complex')
bool_ = Type('bool')


class Attribute(Symbol):
    """ Variable attribute """
    __slots__ = ['name']

    def __new__(cls, name):
        if isinstance(name, Type):
            return name
        return Symbol.__new__(cls, name)

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
        Type of the variable. Inferred from ``symbol`` if not given.

    Examples
    --------
    >>> from sympy import Symbol
    >>> from sympy.codegen.ast import Variable, Type
    >>> i = Symbol('i', integer=True)
    >>> v = Variable.deduced(i)
    >>> v.type == Type('integer')
    True

    """

    nargs = (2, 3)  # type is optional

    def __new__(cls, symbol, attrs=None, type_=None):
        args = (_sympify(symbol), FiniteSet() if attrs is None else FiniteSet(*attrs))
        if type_ is not None:
            if not isinstance(type_, Type):
                raise NotImplementedError("Expected a Type as type_")
            args += (type_,)
        return Basic.__new__(cls, *args)

    @classmethod
    def deduced(cls, symbol, attrs=None):
        """ Alt. constructor with type deduction from ``Type.from_expr`` """
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
        return self.attributes.contains(value_const) == True


class Pointer(Variable):
    """ Represents a pointer """

    @property
    def pointer_const(self):
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
    >>> from sympy.codegen.ast import Declaration, Type, Variable
    >>> x = Symbol('x')
    >>> xvar = Variable(x)
    >>> decl = Declaration.deduced(xvar, 3)
    >>> decl.variable.type == Type('integer')
    True
    >>> k = Symbol('k', integer=True)
    >>> k_decl = Declaration.deduced(k, 3.0)
    >>> k_decl.variable.type == Type('integer')
    True

    """

    nargs = (1, 2)

    def __new__(cls, var, val=None, cast=False):
        if not isinstance(var, Variable):
            raise NotImplementedError("Expected a Variable instance as var")
        args = var,
        if val is not None:
            if cast:
                args += (var.type.cast_check(val),)
            else:
                args += (_sympify(val),)
        return Basic.__new__(cls, *args)

    @classmethod
    def deduced(cls, symbol, val=None, attrs=None, **kwargs):
        """ Deduces type primarily from symbol, secondarily from val """
        try:
            type_ = Type.from_expr(symbol)
        except ValueError:
            type_ = Type.from_expr(val)
        var = Variable(symbol, attrs, type_)
        return cls(var, val, **kwargs)

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
