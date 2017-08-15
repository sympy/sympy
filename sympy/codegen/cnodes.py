"""
AST nodes specific to the C family of languages
"""

from sympy.core.basic import Basic
from sympy.core.compatibility import string_types
from sympy.core.containers import Tuple
from sympy.core.sympify import sympify
from sympy.codegen.ast import Attribute, Declaration, Node, Statement, String, Token, Type, none, FunctionCall


restrict = Attribute('restrict')  # guarantees no pointer aliasing
volatile = Attribute('volatile')
static = Attribute('static')

break_ = Statement(String('break'))
continue_ = Statement(String('continue'))


def alignof(arg):
    return FunctionCall('alignof', [String(arg) if isinstance(arg, string_types) else arg])


def sizeof(arg):
    return FunctionCall('sizeof', [String(arg) if isinstance(arg, string_types) else arg])


class CommaOperator(Basic):
    def __new__(cls, *args):
        return Basic.__new__(cls, *[sympify(arg) for arg in args])


class Label(String):
    """ Label for use with e.g. goto statement. """

class goto(Token):
    __slots__ = ['label']
    _construct_label = Label


class PreDecrement(Basic):
    nargs = 1


class PostDecrement(Basic):
    nargs = 1


class PreIncrement(Basic):
    nargs = 1


class PostIncrement(Basic):
    nargs = 1


class struct(Node):
    __slots__ = ['name', 'declarations']
    defaults = {'name': none}
    _construct_name = String

    @classmethod
    def _construct_declarations(cls, args):
        return Tuple(*[Declaration(arg) for arg in args])


class union(struct):
    pass
