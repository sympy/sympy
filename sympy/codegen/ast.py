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
     |--->DataType
     |           |--------|--->NativeBool
     |                    |--->NativeInteger
     |                    |--->NativeFloat
     |                    |--->NativeDouble
     |                    |--->NativeVoid
     |
     |--->For
     |--->Variable
     |           |--->Argument
     |           |           |
     |           |           |--->InArgument
     |           |           |--->OutArgument
     |           |           |--->InOutArgument
     |           |--->Result
     |
     |--->FunctionDef
     |--->Import
     |--->Declare
     |--->Return
"""

from __future__ import print_function, division


from sympy.core import Symbol, Tuple
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.core.relational import Relational
from sympy.core.compatibility import string_types
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
    """

    def __new__(cls, target, iter, body):
        target = _sympify(target)
        if not iterable(iter):
            raise TypeError("iter must be an iterable")
        iter = _sympify(iter)
        if not iterable(body):
            raise TypeError("body must be an iterable")
        body = Tuple(*(_sympify(i) for i in body))
        return Basic.__new__(cls, target, iter, body)

    @property
    def target(self):
        return self._args[0]

    @property
    def iterable(self):
        return self._args[1]

    @property
    def body(self):
        return self._args[2]


# The following are defined to be sympy approved nodes. If there is something
# smaller that could be used, that would be preferable. We only use them as
# tokens.

class DataType(Basic):
    """Base class representing native datatypes"""
    pass


class NativeBool(DataType):
    _name = 'Bool'
    pass


class NativeInteger(DataType):
    _name = 'Int'
    pass


class NativeFloat(DataType):
    _name = 'Float'
    pass


class NativeDouble(DataType):
    _name = 'Double'
    pass


class NativeVoid(DataType):
    _name = 'Void'
    pass


Bool = NativeBool()
Int = NativeInteger()
Float = NativeFloat()
Double = NativeDouble()
Void = NativeVoid()


# Note: even though there is a "registry" here, these are not
# singeltonized. Don't use "is" comparison on them.
dtype_registry = {'bool': Bool,
                  'int': Int,
                  'float': Float,
                  'double': Double,
                  'void': Void}


def datatype(arg):
    """Returns a datatype instance for the given dtype.

    Parameters
    ----------
    arg : str or sympy expression
        If a str ('bool', 'int', 'float', 'double', or 'void'), return an
        instance of the corresponding dtype. If a sympy expression, return the
        datatype that best fits the expression. This is determined from the
        assumption system. For more control, use the `DataType` class
        directly.

    Returns
    -------
    DataType

    """
    from sympy.matrices import ImmutableDenseMatrix

    def infer_dtype(arg):
        if arg.is_integer:
            return Int
        elif arg.is_Boolean:
            return Bool
        else:
            return Double

    if isinstance(arg, string_types):
        if arg.lower() not in dtype_registry:
            raise ValueError("Unrecognized datatype " + arg)
        return dtype_registry[arg]
    else:
        arg = _sympify(arg)
        if isinstance(arg, ImmutableDenseMatrix):
            dts = [infer_dtype(i) for i in arg]
            if all([i == Bool for i in dts]):
                return Bool
            elif all([i == Int for i in dts]):
                return Int
            else:
                return Double
        else:
            return infer_dtype(arg)


class Variable(Basic):
    """Represents a typed variable.

    Parameters
    ----------
    dtype : str, DataType
        The type of the variable. Can be either a DataType, or a str (bool,
        int, float, double).
    name : Symbol, MatrixSymbol
        The sympy object the variable represents.

    """

    def __new__(cls, dtype, name):
        from sympy.matrices.expressions.matexpr import MatrixSymbol

        if isinstance(dtype, string_types):
            dtype = datatype(dtype)
        elif not isinstance(dtype, DataType):
            raise TypeError("datatype must be an instance of DataType.")
        if isinstance(name, string_types):
            name = Symbol(name)
        elif not isinstance(name, (Symbol, MatrixSymbol)):
            raise TypeError("Only Symbols and MatrixSymbols can be Variables.")
        return Basic.__new__(cls, dtype, name)

    @property
    def dtype(self):
        return self._args[0]

    @property
    def name(self):
        return self._args[1]


class Argument(Variable):
    """An abstract Argument data structure."""
    pass


class Result(Variable):
    """Represents a result directly returned from a routine.

    Parameters
    ----------
    dtype : str, DataType
        The type of the variable. Can be either a DataType, or a str (bool,
        int, float, double).
    name : Symbol or MatrixSymbol, optional
        The sympy object the variable represents.

    """

    def __new__(cls, dtype, name=None):
        if isinstance(dtype, string_types):
            dtype = datatype(dtype)
        elif not isinstance(dtype, DataType):
            raise TypeError("datatype must be an instance of DataType.")
        if not name:
            name = Symbol('')
        return Variable.__new__(cls, dtype, name)


class InArgument(Argument):
    """Argument provided as input only.

    Parameters
    ----------
    dtype : str, DataType
        The type of the variable. Can be either a DataType, or a str (bool,
        int, float, double).
    name : Symbol, MatrixSymbol
        The sympy object the variable represents.

    """
    pass


class OutArgument(Argument):
    """OutputArgument are always initialized in the routine.

    Parameters
    ----------
    dtype : str, DataType
        The type of the variable. Can be either a DataType, or a str (bool,
        int, float, double).
    name : Symbol, MatrixSymbol
        The sympy object the variable represents.

    """
    pass


class InOutArgument(Argument):
    """InOutArgument are never initialized in the routine.

    Parameters
    ----------
    dtype : str, DataType
        The type of the variable. Can be either a DataType, or a str (bool,
        int, float, double).
    name : Symbol, MatrixSymbol
        The sympy object the variable represents.

    """
    pass


class FunctionDef(Basic):
    """Represents a function definition.

    Parameters
    ----------
    name : str
        The name of the function.
    args : iterable
        The arguments to the function, of type `Argument`.
    body : iterable
        The body of the function.
    results : iterable
        The direct outputs of the function, of type `Result`.

    """

    def __new__(cls, name, args, body, results):
        # name
        if isinstance(name, string_types):
            name = Symbol(name)
        elif not isinstance(name, Symbol):
            raise TypeError("Function name must be Symbol or string")
        # args
        if not iterable(args):
            raise TypeError("args must be an iterable")
        if not all(isinstance(a, Argument) for a in args):
            raise TypeError("All args must be of type Argument")
        args = Tuple(*args)
        # body
        if not iterable(body):
            raise TypeError("body must be an iterable")
        body = Tuple(*(_sympify(i) for i in body))
        # results
        if not iterable(results):
            raise TypeError("results must be an iterable")
        if not all(isinstance(i, Result) for i in results):
            raise TypeError("All results must be of type Result")
        results = Tuple(*results)
        return Basic.__new__(cls, name, args, body, results)

    @property
    def name(self):
        return self._args[0]

    @property
    def arguments(self):
        return self._args[1]

    @property
    def body(self):
        return self._args[2]

    @property
    def results(self):
        return self._args[3]


class Import(Basic):
    """Represents inclusion of dependencies in the code.

    Parameters
    ----------
    fil : str
        The filepath of the module (i.e. header in C).
    funcs
        The name of the function (or an iterable of names) to be imported.

    """

    def __new__(cls, fil, funcs=None):
        fil = Symbol(fil)
        if not funcs:
            funcs = Tuple()
        elif iterable(funcs):
            funcs = Tuple(*[Symbol(f) for f in funcs])
        elif isinstance(funcs, string_types):
            funcs = Tuple(Symbol(funcs))
        else:
            raise TypeError("Unrecognized funcs type: ", funcs)
        return Basic.__new__(cls, fil, funcs)

    @property
    def fil(self):
        return self._args[0]

    @property
    def funcs(self):
        return self._args[1]


# TODO: Should Declare have an optional init value for each var?
class Declare(Basic):
    """Represents a variable declaration in the code.

    Parameters
    ----------
    dtype : DataType
        The type for the declaration.
    variable(s)
        A single variable or an iterable of Variables. If iterable, all
        Variables must be of the same type.

    """

    def __new__(cls, dtype, variables):
        if isinstance(dtype, string_types):
            dtype = datatype(dtype)
        elif not isinstance(dtype, DataType):
            raise TypeError("datatype must be an instance of DataType.")
        if isinstance(variables, Variable):
            variables = [variables]
        for var in variables:
            if not isinstance(var, Variable):
                raise TypeError("var must be of type Variable")
            if var.dtype != dtype:
                raise ValueError("All variables must have the same dtype")
        variables = Tuple(*variables)
        return Basic.__new__(cls, dtype, variables)

    @property
    def dtype(self):
        return self._args[0]

    @property
    def variables(self):
        return self._args[1]


class Return(Basic):
    """Represents a function return in the code.

    Parameters
    ----------
    expr : sympy expr
        The expression to return.

    """

    def __new__(cls, expr):
        expr = _sympify(expr)
        return Basic.__new__(cls, expr)

    @property
    def expr(self):
        return self._args[0]
