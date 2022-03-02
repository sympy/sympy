from .abstract_nodes import List as AbstractList
from .ast import Token


class List(AbstractList):
    pass


class NumExprEvaluate(Token):
    """represents a call to `numexpr`s `evaluate`"""
    __slots__ = ('expr',)
