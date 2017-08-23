"""
AST nodes specific to Fortran.

The functions defined in this module allows the user to express functions such as ``dsign``
as a SymPy function for symbolic manipulation.
"""

from sympy.core.basic import Basic
from sympy.core.compatibility import string_types
from sympy.core.containers import Tuple
from sympy.core.function import Function
from sympy.core.numbers import Float, Integer
from sympy.core.sympify import sympify
from sympy.codegen.ast import (
    Attribute, CodeBlock, Declaration, FunctionCall, Node, none, Statement, String,
    Token, Type, _mk_Tuple, Variable
)
from sympy.logic import true, false
from sympy.utilities.iterables import iterable



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
    _construct_body = staticmethod(lambda body: CodeBlock(*body))


class use_rename(Token):
    __slots__ = ['local', 'original']
    _construct_local = String
    _construct_original = String

def _name(arg):
    if hasattr(arg, 'name'):
        return arg.name
    else:
        return String(arg)

class use(Token):
    __slots__ = ['namespace', 'rename', 'only']
    defaults = {'rename': none, 'only': none}
    _construct_namespace = staticmethod(_name)
    _construct_rename = staticmethod(lambda args: Tuple(*[arg if isinstance(arg, use_rename) else use_rename(*arg) for arg in args]))
    _construct_only = staticmethod(lambda args: Tuple(*[arg if isinstance(arg, use_rename) else _name(arg) for arg in args]))


class Module(Token):
    __slots__ = ['name', 'declarations', 'definitions']
    defaults = {'declarations': Tuple()}
    _construct_name = String
    _construct_declarations = staticmethod(lambda arg: CodeBlock(*arg))
    _construct_definitions = staticmethod(lambda arg: CodeBlock(*arg))


class Subroutine(Node):
    __slots__ = ['name', 'parameters', 'body', 'attrs']
    _construct_name = String
    _construct_parameters = staticmethod(lambda params: Tuple(*map(Variable.deduced, params)))

    @classmethod
    def _construct_body(cls, itr):
        if isinstance(itr, CodeBlock):
            return itr
        else:
            return CodeBlock(*itr)

class SubroutineCall(Token):
    __slots__ = ['name', 'subroutine_args']
    _construct_name = staticmethod(_name)
    _construct_function_args = staticmethod(_mk_Tuple)


class Do(Token):
    __slots__ = ['body', 'counter', 'first', 'last', 'step', 'concurrent']
    defaults = {'step': Integer(1), 'concurrent': false}
    _construct_body = staticmethod(lambda body: CodeBlock(*body))
    _construct_var = staticmethod(sympify)
    _construct_start = staticmethod(sympify)
    _construct_stop = staticmethod(sympify)
    _construct_step = staticmethod(sympify)
    _construct_concurrent = staticmethod(lambda arg: true if arg else false)


class GoTo(Token):
    __slots__ = ['labels', 'expr']
    defaults = {'expr': none}
    _construct_labels = staticmethod(_mk_Tuple)
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
        else:
            parameters.append(sympify(arg))
    if len(args) == 0:
        raise ValueError("Need at least one dimension")
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
    return FunctionCall(
        'lbound',
        [_printable(array)] +
        ([_printable(dim)] if dim else []) +
        ([_printable(kind)] if kind else [])
    )


def shape(source, kind=None):
    """ Creates an AST node for a function call to Fortran's "shape(...)"  """
    return FunctionCall(
        'shape',
        [_printable(source)] +
        ([_printable(kind)] if kind else [])
    )


def size(array, dim=None, kind=None):
    """ Creates an AST node for a function call to Fortran's "size(...)"  """
    return FunctionCall(
        'size',
        [_printable(array)] +
        ([_printable(dim)] if dim else []) +
        ([_printable(kind)] if kind else [])
    )


def reshape(source, shape, pad=None, order=None):
    """ Creates an AST node for a function call to Fortran's "reshape(...)"  """
    return FunctionCall(
        'reshape',
        [_printable(source), _printable(shape)] +
        ([_printable(pad)] if pad else []) +
        ([_printable(order)] if pad else [])
    )

def ubound(array, dim=None, kind=None):
    return FunctionCall(
        'ubound',
        [_printable(array)] +
        ([_printable(dim)] if dim else []) +
        ([_printable(kind)] if kind else [])
    )

class FFunction(Function):
    _required_standard = 77

    def _fcode(self, printer):
        name = self.__class__.__name__
        if printer._settings['standard'] < self._required_standard:
            raise NotImplementedError("%s requires Fortran %d or newer" %
                                      (name, self._required_standard))
        return '{0}({1})'.format(name, ', '.join(map(printer._print, self.args)))

class F95Function(FFunction):
    _required_standard = 95


class isign(FFunction):
    """ Fortran sign intrinsic for integer arguments. """
    nargs = 2


class dsign(FFunction):
    """ Fortran sign intrinsic for double precision arguments. """
    nargs = 2


class cmplx(FFunction):
    """ Fortran complex conversion function. """
    nargs = 2  # may be extended to (2, 3) at a later point


class kind(FFunction):
    """ Fortran kind function. """
    nargs = 1


class merge(F95Function):
    """ Fortran merge function """
    nargs = 3


class _literal(Float):
    _token = None
    _decimals = None

    def _fcode(self, printer):
        mantissa, sgnd_ex = ('%.{0}e'.format(self._decimals) % self).split('e')
        mantissa = mantissa.strip('0').rstrip('.')
        ex_sgn, ex_num = sgnd_ex[0], sgnd_ex[1:].lstrip('0')
        ex_sgn = '' if ex_sgn == '+' else ex_sgn
        return (mantissa or '0') + self._token + ex_sgn + (ex_num or '0')


class literal_sp(_literal):
    """ Fortran single precision real literal """
    _token = 'e'
    _decimals = 9


class literal_dp(_literal):
    """ Fortran double precision real literal """
    _token = 'd'
    _decimals = 17


class sum_(Token):
    __slots__ = ['array', 'dim', 'mask']
    defaults = {'dim': none, 'mask': none}
    _construct_array = staticmethod(sympify)


class product_(Token):
    __slots__ = ['array', 'dim', 'mask']
    defaults = {'dim': none, 'mask': none}
    _construct_array = staticmethod(sympify)
