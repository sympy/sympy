"""
Types used to represent a full function/module as an Abstract Syntax Tree.

Most types are small, and are merely used as tokens in the AST. A tree diagram
has been included below to illustrate the relationships between the AST types.


AST Type Tree
-------------

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
        Operator (+, -, /, *, %).

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

    ``left_hand_sides``: Tuple of left-hand sides of assignments, in order.
    ``left_hand_sides``: Tuple of right-hand sides of assignments, in order.
    ``topological_sort``: Class method. Return a CodeBlock with assignments
                          sorted so that variables are assigned before they
                          are used.
    ``cse``: Return a new CodeBlock with common subexpressions eliminated and
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
        CodeBlock(Assignment(x, 1), Assignment(x0, sin(x)), Assignment(y, x0 +
        1), Assignment(z, x0 - 1))
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


class Type(Basic):
    """ Represents a type.

    The naming is a super-set of NumPy naming, see [1]_.

    Arguments
    ---------
    name : str or Type
        Either an explicit type: ``intc``, ``intp``, ``int8``, ``int16``,
        ``int32``, ``int64``, ``uint8``, ``uint16``, ``uint32``, ``uint64,
        float16``, ``float32``, ``float64``, ``complex64``, ``complex128``,
        ``bool. Or only kind (precision decided by code-printer): ``integer``,
        ``real`` or ``complex`` (where the latter two are of floating point type).
        If a ``Type`` instance is given, the said instance is returned.

    References
    ==========
    [1] https://docs.scipy.org/doc/numpy/user/basics.types.html

    """
    allowed_names = tuple('intc intp int8 int16 int32 int64 uint8 uint16 uint32'.split() +
                          'uint64 float16 float32 float64 complex64 complex128'.split() +
                          'real integer complex bool'.split())
    __slots__ = []

    default_limits = {
        'INT_MAX': 2**31 - 1,
        'INT_MIN': -2**31,
        'FLT_MAX': 3.40282347e+38,
        'FLT_MIN': 1.17549435e-38,  # excluding subnormal numbers
        'FLT_DIG': 6,
        'FLT_EPSILON': 1.1920929e-07,
        'DBL_MAX': 1.79769313486231571e+308,
        'DBL_MIN': 2.22507385850720138e-308,  # excluding subnormal numbers
        'DBL_DIG': 15,
        'DBL_EPSILON': 2.22044604925031308e-16,
        'LDBL_MAX': 1.18973149535723176502e+4932,
        'LDBL_MIN': 3.36210314311209350626e-4932,
        'LDBL_DIG': 18,
        'LDBL_EPSILON': 1.08420217248550443401e-19
    }

    def __new__(cls, name):
        if isinstance(name, Type):
            return name
        if name not in cls.allowed_names:
            raise ValueError("Unknown type: %s" % name)
        return Basic.__new__(cls, name)

    @property
    def name(self):
        """ Return the name of the type. """
        return self._args[0]

    @classmethod
    def from_expr(cls, expr, symb=None):
        """ Infers type from an expression or a ``Symbol``.

        Parameters
        ----------
        expr : number, string or SymPy object
            The typename will be deduced from type or properties. Default is 'real'
            (e.g. when a string is given as expr).
        symb : Symbol (optional)
            If given, assumptions of ``symb`` has higher precedence than expr.

        Examples
        --------
        >>> from sympy.codegen.ast import Type
        >>> Type.from_expr(2) == Type('integer')
        True
        >>> Type.from_expr('i') == Type('integer')
        False
        >>> from sympy import Symbol
        >>> Type.from_expr(2, Symbol('j', complex=True)) == Type('complex')
        True

        """
        if symb is not None:
            if getattr(symb, 'is_integer', False):
                return cls('integer')
            if getattr(symb, 'is_complex', False):
                return cls('complex')

        if isinstance(expr, str):
            return cls('real')  # default
        else:
            if isinstance(expr, (float, Float)):
                return cls('real')
            if isinstance(expr, complex) or getattr(expr, 'is_complex', False):
                return cls('complex')
            if isinstance(expr, (int, Integer)) or getattr(expr, 'is_integer', False):
                return cls('integer')

            if getattr(expr, 'is_Relational', False):
                return cls('bool')
            else:
                return cls('real')

    def cast_check(self, value, rtol=1e-8, atol=1e-8, limits=None):
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

        Examples
        --------
        >>> from sympy.codegen.ast import Type
        >>> Type('integer').cast_check(3.0) is 3
        True
        >>> Type('float32').cast_check(1e-40)  # doctest: +ELLIPSIS
        Traceback (most recent call last):
          ...
        ValueError: Minimum value for data type bigger than new value.

        """
        from sympy.functions.elementary.complexes import im, re

        def lim(key):
            return (limits or self.default_limits)[key]

        def tol(val):
            return atol + rtol*abs(val)

        caster = lambda x: x  # identity
        flt_caster = lambda x: round(x, lim('FLT_DIG') + 3)
        _min, _max = -oo, oo  # undefined precision

        if self.name == 'integer':
            caster = int
        elif self.name.startswith('int'):
            nbits = int(self.name.split('int')[1])
            _min, _max = -2**(nbits - 1), 2**(nbits - 1) - 1
            caster = int

        elif self.name.startswith('uint'):
            nbits = int(self.name.split('uint')[1])
            _min, _max = 0, 2**nbits - 1
            caster = int

        elif self.name.startswith('float'):
            nbits = int(self.name.split('float')[1])
            if nbits == 32:
                _min, _max = lim('FLT_MIN'), lim('FLT_MAX')
                caster = flt_caster
            elif nbits == 64:
                _min, _max = lim('DBL_MIN'), lim('DBL_MAX')
                caster = float  # Python's float is double precision
            # long double is usually 80 bits long, but sometimes 128 bits
            # in NumPy float128 is actually often 80 bits (quite confusing)

        if self.name.startswith('complex'):
            if self.name != 'complex':
                nbits = int(self.name.split('complex')[1])
                if nbits == 64:
                    _min, _max = lim('FLT_MIN'), lim('FLT_MAX')
                    caster = lambda x: flt_caster(re(x)) + 1j*flt_caster(im(x))
                elif nbits == 128:
                    _min, _max = lim('DBL_MIN'), lim('DBL_MAX')
                    caster = complex  # Python's float is double precision
            if re(value) > _max or im(value) > _max:
                raise ValueError("Maximum value exceeded for data type.")
            if re(value) < _min or im(value) < _min:
                raise ValueError("Minimum value for data type bigger than new value.")
            new_val = caster(value)
            delta = new_val - value
            if abs(re(delta)) > tol(re(value)) or abs(im(delta)) > tol(im(value)):
                raise ValueError("Casting gives a significantly different value.")
        else:
            if value > _max:
                raise ValueError("Maximum value exceeded for data type.")
            if value < _min:
                raise ValueError("Minimum value for data type bigger than new value.")
            new_val = caster(value)
            delta = new_val - value
            if abs(delta) > tol(value):  # e.g. int(3.5) != 3.5
                raise ValueError("Casting gives a significantly different value.")

        return new_val


class Variable(Basic):
    """ Represents a variable

    Parameters
    ----------
    symbol : Symbol
        If a ``Variable`` instance is given as symbol, said instance is simply returned.
    type_ : Type (optional)
        Type of the variable. Inferred from ``symbol`` if not given.
    const : bool
        Constness of the variable.

    Examples
    --------
    >>> from sympy import Symbol
    >>> from sympy.codegen.ast import Variable, Type
    >>> i = Symbol('i', integer=True)
    >>> v = Variable(i)
    >>> v.type == Type('integer')
    True

    """
    def __new__(cls, symbol, type_=None, const=False):
        if isinstance(symbol, Variable):
            return symbol
        if type_ is None:
            type_ = Type.from_expr(None, symbol)
        return Basic.__new__(cls, _sympify(symbol), Type(type_), _sympify(const or False))

    @property
    def symbol(self):
        return self.args[0]

    @property
    def type(self):
        return self.args[1]

    @property
    def const(self):
        return self.args[2]


class Pointer(Basic):
    """ Represents a pointer

    Parameters
    ----------
    name : Symbol
    type_ : Type
        Type of the variable. Inferred from ``symbol`` if not given.
    value_const : bool
        Constness of the value pointed to by the variable.
    pointer_const : bool
        Constness of the pointer (i.e. opposite of the mutability of the address).
    restrict : bool
        Applicable to C > C99. If the pointer is guaranteed not to alias another pointer
        this may be set to ``True``. (allows the compiler to generate efficient assembly).

    Examples
    --------
    >>> from sympy import Symbol
    >>> from sympy.codegen.ast import Pointer, Type
    >>> x = Symbol('x')
    >>> p = Pointer(x, value_const=True, pointer_const=True, restrict=True)
    >>> p.type == Type('real')
    True

    """
    def __new__(cls, symbol, type_=None, value_const=False, pointer_const=False, restrict=False):
        if type_ is None:
            type_ = Type.from_expr(symbol)
        args = symbol, type_, value_const or False, pointer_const, restrict
        return Basic.__new__(cls, *map(_sympify, args))

    @property
    def symbol(self):
        return self.args[0]

    @property
    def type(self):
        return self.args[1]

    @property
    def value_const(self):
        return self.args[2]

    @property
    def pointer_const(self):
        return self.args[3]

    @property
    def restrict(self):
        return self.args[4]


class Declaration(Basic):
    """ Represents a variable declaration

    Parameters
    ----------
    var : Variable, Pointer or IndexedBase
    value : Value (optional)
        Value to be assigned upon declaration.
    const : bool
        Constness of value. Can not be given if ``var`` is Variable
        or Pointer since constness is then taken from ``var``.

    Examples
    --------
    >>> from sympy import Symbol
    >>> from sympy.codegen.ast import Declaration, Type
    >>> x = Symbol('x')
    >>> decl = Declaration(x, 3)
    >>> decl.variable.type == Type('integer')
    True
    >>> k = Symbol('k', integer=True)
    >>> k_decl = Declaration(k, 3.0)
    >>> k_decl.variable.type == Type('integer')
    True

    """
    def __new__(cls, var, value=None, const=None):
        from sympy.tensor.indexed import IndexedBase

        if isinstance(var, (Variable, Pointer)):
            if const is not None:
                raise ValueError("Cannot change constness of an existing Variable/Pointer")
        else:
            if value is not None:
                type_ = Type.from_expr(value, var)
                value = type_.cast_check(value)
            else:
                type_ = None

            if isinstance(var, IndexedBase):
                var = Pointer(var, type_, const)
            else:
                var = Variable(var, type_, const)

        return Basic.__new__(cls, var, sympify(value))

    @property
    def variable(self):
        return self.args[0]

    @property
    def value(self):
        return self.args[1]
