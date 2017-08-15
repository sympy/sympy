"""
AST nodes specific for C++
"""

from sympy.codegen.ast import String, Token, Type, none

class using(Token):
    __slots__ = ['type', 'alias']
    defaults = {'alias': none}
    _construct_type = Type
    _construct_alias = String
