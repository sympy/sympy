"""
Types used to represent a full function/module as an Abstract Syntax Tree.

Most types are small, and are merely used as tokens in the AST. A tree diagram
has been included below to illustrate the relationships between the AST types.


AST Type Tree
-------------

*Basic*
     |--->Assign
     |--->AugAssign
     |--->NativeOp
     |           |--------------|
     |                          |--->AddOp
     |                          |--->SubOp
     |                          |--->MulOp
     |                          |--->DivOp
     |                          |--->ModOp
     |           *Singleton*----|
     |                    |
     |--->DataType        |
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
from sympy.core.singleton import Singleton
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.core.compatibility import with_metaclass
from sympy.tensor import Indexed
from sympy.matrices import ImmutableDenseMatrix
from sympy.matrices.expressions.matexpr import MatrixSymbol, MatrixElement
from sympy.utilities.iterables import iterable


class Assign(Basic):
    """Represents variable assignment for code generation.

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
    --------

    >>> from sympy import symbols, MatrixSymbol, Matrix
    >>> from sympy.codegen.ast import Assign
    >>> x, y, z = symbols('x, y, z')
    >>> Assign(x, y)
    x := y
    >>> Assign(x, 0)
    x := 0
    >>> A = MatrixSymbol('A', 1, 3)
    >>> mat = Matrix([x, y, z]).T
    >>> Assign(A, mat)
    A := Matrix([[x, y, z]])
    >>> Assign(A[0, 1], x)
    A[0, 1] := x

    """

    def __new__(cls, lhs, rhs):
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
        return Basic.__new__(cls, lhs, rhs)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{0} := {1}'.format(sstr(self.lhs), sstr(self.rhs))

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]


# The following are defined to be sympy approved nodes. If there is something
# smaller that could be used, that would be preferable. We only use them as
# tokens.


class NativeOp(with_metaclass(Singleton, Basic)):
    """Base type for native operands."""
    pass


class AddOp(NativeOp):
    _symbol = '+'


class SubOp(NativeOp):
    _symbol = '-'


class MulOp(NativeOp):
    _symbol = '*'


class DivOp(NativeOp):
    _symbol = '/'


class ModOp(NativeOp):
    _symbol = '%'


op_registry = {'+': AddOp(),
               '-': SubOp(),
               '*': MulOp(),
               '/': DivOp(),
               '%': ModOp()}


def operator(op):
    """Returns the operator singleton for the given operator"""

    if op.lower() not in op_registry:
        raise ValueError("Unrecognized operator " + op)
    return op_registry[op]


class AugAssign(Basic):
    """Represents augmented variable assignment for code generation.

    Parameters
    ----------
    lhs : Expr
        Sympy object representing the lhs of the expression. These should be
        singular objects, such as one would use in writing code. Notable types
        include Symbol, MatrixSymbol, MatrixElement, and Indexed. Types that
        subclass these types are also supported.

    op : NativeOp
        Operator (+, -, /, *, %).

    rhs : Expr
        Sympy object representing the rhs of the expression. This can be any
        type, provided its shape corresponds to that of the lhs. For example,
        a Matrix type can be assigned to MatrixSymbol, but not to Symbol, as
        the dimensions will not align.

    Examples
    --------

    >>> from sympy import symbols
    >>> from sympy.codegen.ast import AugAssign, AddOp
    >>> x, y = symbols('x, y')
    >>> AugAssign(x, AddOp(), y)
    x += y

    """

    def __new__(cls, lhs, op, rhs):
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
        if isinstance(op, str):
            op = operator(op)
        elif op not in op_registry.values():
            raise TypeError("Unrecognized Operator")
        return Basic.__new__(cls, lhs, op, rhs)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '{0} {1}= {2}'.format(sstr(self.lhs), self.op._symbol,
                sstr(self.rhs))

    @property
    def lhs(self):
        return self._args[0]

    @property
    def op(self):
        return self._args[1]

    @property
    def rhs(self):
        return self._args[2]


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

class DataType(with_metaclass(Singleton, Basic)):
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


dtype_registry = {'bool': Bool,
                  'int': Int,
                  'float': Float,
                  'double': Double,
                  'void': Void}


def datatype(arg):
    """Returns the datatype singleton for the given dtype.

    Parameters
    ----------
    arg : str or sympy expression
        If a str ('bool', 'int', 'float', 'double', or 'void'), return the
        singleton for the corresponding dtype. If a sympy expression, return
        the datatype that best fits the expression. This is determined from the
        assumption system. For more control, use the `DataType` class directly.

    Returns
    -------
    DataType

    """
    def infer_dtype(arg):
        if arg.is_integer:
            return Int
        elif arg.is_Boolean:
            return Bool
        else:
            return Double

    if isinstance(arg, str):
        if arg.lower() not in dtype_registry:
            raise ValueError("Unrecognized datatype " + arg)
        return dtype_registry[arg]
    else:
        arg = _sympify(arg)
        if isinstance(arg, ImmutableDenseMatrix):
            dts = [infer_dtype(i) for i in arg]
            if all([i is Bool for i in dts]):
                return Bool
            elif all([i is Int for i in dts]):
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
        if isinstance(dtype, str):
            dtype = datatype(dtype)
        elif not isinstance(dtype, DataType):
            raise TypeError("datatype must be an instance of DataType.")
        if isinstance(name, str):
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
        if isinstance(dtype, str):
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
        if isinstance(name, str):
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
        elif isinstance(funcs, str):
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
        if isinstance(dtype, str):
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
