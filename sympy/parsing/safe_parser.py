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
whitelisted_nodes = [
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
    whitelisted_nodes += [
        ast.Bytes
    ]

if sys.version_info >= (3, 4):
    whitelisted_nodes += [
        # True, False, None
        ast.NamedConstant,
    ]

if sys.version_info >= (3, 5):
    whitelisted_nodes += [
        # Matmul (A @ B)
        ast.MatMult,
    ]
