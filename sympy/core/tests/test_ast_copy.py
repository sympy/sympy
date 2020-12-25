import pickle
from copy import deepcopy

import sympy as sp
from sympy.core.ast_copy import clone, to_ast

x, y, z = sp.var("x y z")
neg_x = sp.Symbol("x", negative=True)
p = sp.Symbol("p", prime=True)


def test_clone():
    for expr in [
        sp.Add(
            sp.Pow(sp.Integer(2), 2, evaluate=False),
            sp.Pow(sp.Float(3.1), 2, evaluate=False),
            evaluate=False,
        ),
        # from Issue #5783
        sp.Mul(sp.Float(-11.1), sp.Integer(19)),
        # from PR #622
        sp.Pow(3, 2, evaluate=False),
        # symbols
        x ** 2 + sp.sin(y) + x * y * z,
        # rational
        sp.Add(sp.Rational(2, 3), sp.Float(sp.pi), evaluate=False),
        # assumptions
        neg_x ** 3 + sp.sin(p),
    ]:
        assert clone(expr) == expr


def test_pickle_unpickle_ast():
    for expr in [
        sp.Add(
            sp.Pow(sp.Integer(2), 2, evaluate=False),
            sp.Pow(sp.Float(3.1), 2, evaluate=False),
            evaluate=False,
        ),
        # from Issue #5783
        sp.Mul(sp.Float(-11.1), sp.Integer(19)),
        # from PR #622
        sp.Pow(3, 2, evaluate=False),
        # symbols
        x ** 2 + sp.sin(y) + x * y * z,
        # rational
        sp.Add(sp.Rational(2, 3), sp.Float(sp.pi), evaluate=False),
        # assumptions
        neg_x ** 3 + sp.sin(p),
    ]:
        ast_expr = to_ast(expr)
        orig = compile(ast_expr, filename="<ast>", mode="eval")
        unpickled_ast_expr = pickle.loads(pickle.dumps(ast_expr))
        new = compile(unpickled_ast_expr, filename="<ast>", mode="eval")
        assert eval(orig) == eval(new)


def test_7672():
    # deepcopy interferes with original Issue #7672
    x, xx = sp.symbols("x"), sp.symbols("x", real=True)
    assert x != xx
    xxx = deepcopy(xx)
    assert x == xx

    # clone does not
    y, yy = sp.symbols("y"), sp.symbols("y", real=True)
    assert y != yy
    # deepcopy broken (should not be equal)
    yyy = clone(yy)
    assert y != yy
