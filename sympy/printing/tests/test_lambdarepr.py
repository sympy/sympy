from sympy import symbols, sin, Matrix, Interval, Piecewise, Sum, lambdify
from sympy.utilities.pytest import raises

from sympy.printing.lambdarepr import lambdarepr

x, y, z = symbols("x,y,z")
i, a, b = symbols("i,a,b")
j, c, d = symbols("j,c,d")


def test_basic():
    assert lambdarepr(x*y) == "x*y"
    assert lambdarepr(x + y) in ["y + x", "x + y"]
    assert lambdarepr(x**y) == "x**y"


def test_matrix():
    A = Matrix([[x, y], [y*x, z**2]])
    assert lambdarepr(A) == "MutableDenseMatrix([[x, y], [x*y, z**2]])"

    # Test printing a Matrix that has an element that is printed differently
    # with the LambdaPrinter than in the StrPrinter.
    p = Piecewise((x, True), evaluate=False)
    A = Matrix([p])
    assert lambdarepr(A) == "MutableDenseMatrix([[((x) if (True) else None)]])"


def test_piecewise():
    # In each case, test eval() the lambdarepr() to make sure there are a
    # correct number of parentheses. It will give a SyntaxError if there aren't.

    h = "lambda x: "

    p = Piecewise((x, True), evaluate=False)
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((x) if (True) else None)"

    p = Piecewise((x, x < 0))
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((x) if (x < 0) else None)"

    p = Piecewise(
        (1, x < 1),
        (2, x < 2),
        (0, True)
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((1) if (x < 1) else (((2) if (x < 2) else " \
        "(((0) if (True) else None)))))"

    p = Piecewise(
        (1, x < 1),
        (2, x < 2),
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((1) if (x < 1) else (((2) if (x < 2) else None)))"

    p = Piecewise(
        (x, x < 1),
        (x**2, Interval(3, 4, True, False).contains(x)),
        (0, True),
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((x) if (x < 1) else (((x**2) if (((x <= 4) and " \
        "(x > 3))) else (((0) if (True) else None)))))"

    p = Piecewise(
        (x**2, x < 0),
        (x, Interval(0, 1, False, True).contains(x)),
        (2 - x, x >= 1),
        (0, True)
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((x**2) if (x < 0) else (((x) if (((x >= 0) and (x < 1))) " \
        "else (((-x + 2) if (x >= 1) else (((0) if (True) else None)))))))"

    p = Piecewise(
        (x**2, x < 0),
        (x, Interval(0, 1, False, True).contains(x)),
        (2 - x, x >= 1),
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((x**2) if (x < 0) else (((x) if (((x >= 0) and " \
        "(x < 1))) else (((-x + 2) if (x >= 1) else None)))))"

    p = Piecewise(
        (1, x >= 1),
        (2, x >= 2),
        (3, x >= 3),
        (4, x >= 4),
        (5, x >= 5),
        (6, True)
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == ("((1) if (x >= 1) else (((2) if (x >= 2) else (((3) if "
        "(x >= 3) else (((4) if (x >= 4) else (((5) if (x >= 5) else (((6) if "
        "(True) else None)))))))))))")

    p = Piecewise(
        (1, x <= 1),
        (2, x <= 2),
        (3, x <= 3),
        (4, x <= 4),
        (5, x <= 5),
        (6, True)
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((1) if (x <= 1) else (((2) if (x <= 2) else (((3) if " \
        "(x <= 3) else (((4) if (x <= 4) else (((5) if (x <= 5) else (((6) if " \
        "(True) else None)))))))))))"

    p = Piecewise(
        (1, x > 1),
        (2, x > 2),
        (3, x > 3),
        (4, x > 4),
        (5, x > 5),
        (6, True)
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == ("((1) if (x > 1) else (((2) if (x > 2) else (((3) if "
        "(x > 3) else (((4) if (x > 4) else (((5) if (x > 5) else (((6) if "
        "(True) else None)))))))))))")

    p = Piecewise(
        (1, x < 1),
        (2, x < 2),
        (3, x < 3),
        (4, x < 4),
        (5, x < 5),
        (6, True)
    )
    l = lambdarepr(p)
    eval(h + l)
    assert l == "((1) if (x < 1) else (((2) if (x < 2) else (((3) if " \
        "(x < 3) else (((4) if (x < 4) else (((5) if (x < 5) else (((6) if " \
        "(True) else None)))))))))))"


def test_sum():
    # In each case, test eval() the lambdarepr() to make sure that
    # it evaluates to the same results as the symbolic expression

    s = Sum(x ** i, (i, a, b))

    l = lambdarepr(s)
    assert l == "(builtins.sum(x**i for i in range(a, b+1)))"

    assert (lambdify((x, a, b), s)(2, 3, 8) ==
            s.subs([(x, 2), (a, 3), (b, 8)]).doit())

    s = Sum(i * x, (i, a, b))

    l = lambdarepr(s)
    assert l == "(builtins.sum(i*x for i in range(a, b+1)))"

    assert (lambdify((x, a, b), s)(2, 3, 8) ==
            s.subs([(x, 2), (a, 3), (b, 8)]).doit())


def test_multiple_sums():
    s = Sum(i * x + j, (i, a, b), (j, c, d))

    l = lambdarepr(s)
    assert l == "(builtins.sum(i*x + j for i in range(a, b+1) for j in range(c, d+1)))"

    assert (lambdify((x, a, b, c, d), s)(2, 3, 4, 5, 6) ==
            s.subs([(x, 2), (a, 3), (b, 4), (c, 5), (d, 6)]).doit())

def test_settings():
    raises(TypeError, lambda: lambdarepr(sin(x), method="garbage"))
