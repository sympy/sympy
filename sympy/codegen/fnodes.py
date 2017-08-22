"""
AST nodes specific to the C family of languages
"""

from sympy.core.basic import Basic
from sympy.core.compatibility import string_types
from sympy.core.containers import Tuple
from sympy.core.sympify import sympify
from sympy.codegen.ast import (
    Attribute, CodeBlock, Declaration, FunctionCall, Node, none, Statement, String,
    Token, Type
)


pure = Attribute('pure')
elemental = Attribute('elemental')  # (all elemental procedures are also pure)

intent_in = Attribute('intent_in')
intent_out = Attribute('intent_out')
intent_inout = Attribute('intent_inout')

exit_ = Statement(String('exit'))  # "break" in C
cycle = Statement(String('continue'))  # "continue" in C

class Program(Token):
    __slots__ = ['name', 'body']
    _construct_name = String
    _construct_body = staticmethod(lambda body: Tuple(*body))


class Module(Token):
    __slots__ = ['name', 'declarations', 'definitions']
    _construct_name = String
    _construct_declarations = staticmethod(lambda arg: CodeBlock(*arg))
    _construct_definitions = staticmethod(lambda arg: CodeBlock(*arg))


class Subroutine(Token):
    __slots__ = ['name', 'body']
    _construct_name = String
    _construct_body = staticmethod(lambda body: Tuple(*body))


class GoTo(Token):
    __slots__ = ['labels', 'expr']
    defaults = {'expr': none}

    @classmethod
    def _construct_labels(cls, labels):
        return Tuple(*labels)

    _construct_expr = staticmethod(sympify)


class FortranReturn(Token):
    """ AST node explicitly mapped to a fortran "return".

    Because a return statement in fortran is different from C, and
    in order to aid reuse of our codegen ASTs the ordinary
    ``.codegen.ast.ReturnStatement`` is interpreted as assignment to
    the result variable of the function. If one for some reason needs
    to generate a fortran RETURN statement, this node should be used.
    """
    __slots__ = ['return_value']
    defaults = {'return_value': none}
    _construct_return_value = staticmethod(sympify)


class Extent(Basic):
    """ Represents a dimension extent.

    Examples
    --------
    >>> from sympy.codegen.fnodes import Extent
    >>> e = Extent(-3, 3)  # -3, -2, -1, 0, 1, 2, 3

    """
    def __new__(cls, *args):
        if len(args) == 2:
            low, high = args
            return Basic.__new__(cls, sympify(low), sympify(high))
        elif len(args) == 0 or (len(args) == 1 and args[0] in (':', None)):
            return Basic.__new__(cls)  # assumed shape
        else:
            raise ValueError("Expected 0 or 2 args (or one argument == None or ':')")

    def _sympystr(self, printer):
        if len(self.args) == 0:
            return ':'
        return '%d:%d' % self.args

assumed_extent = Extent() # or Extent(':'), Extent(None)


def dimension(*args):
    if len(args) > 7:
        raise ValueError("Fortran only supports up to 7 dimensional arrays")
    parameters = []
    for arg in args:
        if isinstance(arg, Extent):
            parameters.append(arg)
        elif isinstance(arg, string_types):
            if arg == ':':
                parameters.append(Extent())
            else:
                parameters.append(String(arg))
        elif iterable(arg):
            parameters.append(Extent(*arg))
    return Attribute('dimension', parameters)

assumed_size = dimension('*')

def bind_C(name=None):
    return Attribute('bind_C', [String(name)] if name else [])

def _printable(arg):
    return String(arg) if isinstance(arg, string_types) else arg


def allocated(array):
    """ Creates an AST node for a function call to Fortran's "allocated(...)"  """
    return FunctionCall('allocated', [_printable(array)])


def lbound(array, dim=None, kind=None):
    return FunctionCall('lbound',
                        [_printable(array)],
                        [_printable(dim)] if dim else [] +
                        [_printable(kind)] if kind else [])


def shape(source, kind=None):
    """ Creates an AST node for a function call to Fortran's "shape(...)"  """
    return FunctionCall('shape',
                        [_printable(source)] +
                        [_printable(kind)] if kind else [])


def size(array, dim=None, kind=None):
    """ Creates an AST node for a function call to Fortran's "size(...)"  """
    return FunctionCall('size',
                        [_printable(array)] +
                        [_printable(dim)] if dim else [] +
                        [_printable(kind)] if kind else [])


def reshape(source, shape, pad=None, order=None):
    """ Creates an AST node for a function call to Fortran's "reshape(...)"  """
    return FunctionCall('reshape',
                        [_printable(source), _printable(shape)] +
                        [_printable(pad)] if pad else [] +
                        [_printable(order)] if pad else []
    )

def ubound(array, dim=None, kind=None):
    return FunctionCall('ubound',
                        [_printable(array)],
                        [_printable(dim)] if dim else [] +
                        [_printable(kind)] if kind else [])
