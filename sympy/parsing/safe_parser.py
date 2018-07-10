# -*- coding: utf-8 -*-
"""
Parser that allows evaling expressions (mostly) safely

Implements AST whitelisting.

Warning: there are still limitations here, and this should not be used as a
substitute for proper sandboxing when parsing untrusted input.

When working with this file, following is an invaluable reference:
https://greentreesnakes.readthedocs.io/en/latest/nodes.html#expressions
"""
import sys

import ast

# Only Expr nodes are allowed (we only care about eval). Any node that isn't
# needed for sympify() shouldn't be included, even if it is nominally safe to
# minimize the attack surface. The dangerous nodes include comprehensions
# (they create their own namespace), and attribute access (see
# https://nedbatchelder.com/blog/201206/eval_really_is_dangerous.html).
WHITELISTED_NODES = [
    ast.Module,
    ast.Expr,

    # Literals
    ast.Num,
    ast.Str,
    # Not included: FormattedValue, JoinedStr (f-strings)
    ast.List,
    ast.Tuple,
    ast.Set,
    ast.Dict,
    ast.Ellipsis,

    # Variables
    ast.Name,
    ast.Load,

    # Expressions

    # Unary operations
    ast.UnaryOp,
    ast.UAdd,
    ast.USub,
    ast.Invert,
    # Not included: Not

    # Binary operations
    ast.BinOp,
    ast.Add,
    ast.Sub,
    ast.Mult,
    ast.Div,
    ast.FloorDiv,
    ast.Mod,
    ast.Pow,
    ast.LShift,
    ast.RShift,
    ast.BitOr,
    ast.BitXor,
    ast.BitAnd,

    # Not included: BoolOp, And, Or

    # Comparisons
    ast.Compare,
    ast.Eq,
    ast.NotEq,
    ast.Lt,
    ast.LtE,
    ast.Gt,
    ast.GtE,
    # Not included: Is, IsNot, In, NotIn

    # Other
    ast.Call,
    ast.keyword,
    ast.Lambda,
    ast.arguments,
    ast.Param,
    # Not included: IfExp, Attribute

    # Subscripting
    ast.Subscript,
    ast.Index,
    ast.Slice,
    ast.ExtSlice,

    # Not included: ListComp, SetComp, GeneratorExp, DictComp, comprehension
]

# Python 3-only
if sys.version_info >= (3,):
    WHITELISTED_NODES += [
        ast.Bytes,
        ast.arg,
    ]

if sys.version_info >= (3, 4):
    WHITELISTED_NODES += [
        # True, False, None
        ast.NameConstant,
    ]

if sys.version_info >= (3, 5):
    WHITELISTED_NODES += [
        # Matmul (A @ B)
        ast.MatMult,
    ]

WHITELISTED_NODES = tuple(WHITELISTED_NODES)

BLACKLISTED_NAMES = ['eval', 'exec', 'sympify', 'parse_expr', 'locals', 'globals']

def check_string_for_safety(s, whitelisted_nodes=WHITELISTED_NODES,
    blacklisted_names=BLACKLISTED_NAMES):
    """
    Checks if the input string is safe for parsing.

    Returns None if the string s is safe to parse and raises
    UnsafeSympifyError if it is not.

    whitelisted_nodes should be a tuple of AST types to be whitelisted.

    blacklisted_names should be a list of variable names that are disallowed.

    """
    from ..core.sympify import UnsafeSympifyError

    p = ast.parse(s)
    for node in ast.walk(p):
        if not isinstance(node, whitelisted_nodes):
            raise UnsafeSympifyError(s, reason='Non-whitelisted AST node %s found' % node)

        if isinstance(node, ast.Name) and node.id in blacklisted_names:
            # Note, the AST module normalizes Unicode characters
            # automatically, so we don't need to worry about things like evªl.
            raise UnsafeSympifyError(s, reason='disallowed name %r found' % node.id)
