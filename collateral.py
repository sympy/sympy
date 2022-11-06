from sympy import Symbol, count_ops, Ge, Rational, Le, Lt, Piecewise, ITE, S, symbols, integrate, sin, pi, Eq, Or, And

def collateral_damage_tests():
    # Tests that have failed due to changes in the simplifier.
    # measure function is count_ops.

    q = Symbol("q")
    r = Symbol("r")
    x = Symbol("x")
    y = Symbol("y")
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    e = Symbol("e")
    f = Symbol("f")
    g = Symbol("g")
    z = Symbol("z")


    expr = ((-q + r) - (q - r) <= 0)
    actual = expr.simplify()
    print("File \"/home/runner/work/sympy/sympy/sympy/core/tests/test_relational.py\", line 727, "
          "in test_simplify_relational")
    print(f"Before: {expr}")
    print("Expected: q >= r")
    print(f"Actual: {actual}, count_ops: {count_ops(actual)}\n")

    expr = Ge(3 * x * (x + 1) + 4, 3 * x)
    actual = expr.simplify()
    expected = [Ge(x ** 2, -Rational(4, 3)), Le(-x ** 2, Rational(4, 3))]
    print("File \"/home/runner/work/sympy/sympy/sympy/core/tests/test_relational.py\", line 1143, "
          "in test_polynomial_relation_simplification")
    print(f"Before: {expr}")
    print(f"Expected: {expected}, count_ops: {[count_ops(i) for i in expected]}")
    print(f"Actual: {actual}, count_ops: {count_ops(actual)}\n")

    print("File \"/home/runner/work/sympy/sympy/sympy/core/tests/test_relational.py\", line 1156, "
          "in test_multivariate_linear_function_simplification")
    expr = Lt(a + b + c + 2 * d, 3 * d - f + g)
    expected = Lt(a, -b - c + d - f + g)
    actual = expr.simplify()
    print(f"Before: {expr}")
    print(f"Expected: {expected}, count_ops: {count_ops(expected)}")
    print(f"Actual: {actual}, count_ops: {count_ops(actual)}\n")

    from sympy.logic.boolalg import ITE
    from sympy.functions.elementary.piecewise import piecewise_fold
    from sympy import print_latex

    p = (Piecewise((0, ITE((x - y > 1) | (2 * x - 2 * y > 1), False,
                           ITE(x - y > 1, 2 * y - 2 < -1, 2 * x - 2 * y > 1))),
                   (Piecewise((0, ITE(x - y > 1, True, 2 * x - 2 * y > 1)),
                              (2 * Piecewise((0, x - y > 1), (y, True)), True)), True))
         + 2 * Piecewise((1, ITE((x - y > 1) | (2 * x - 2 * y > 1), False,
                                 ITE(x - y > 1, 2 * y - 2 < -1, 2 * x - 2 * y > 1))),
                         (Piecewise((1, ITE(x - y > 1, True, 2 * x - 2 * y > 1)),
                                    (2 * Piecewise((1, x - y > 1), (x, True)), True)), True)))
    print("File \"/home/runner/work/sympy/sympy/sympy/functions/elementary/tests/test_piecewise.py\", line 454, "
          "in test_issue_22917")
    expr = piecewise_fold(p)
    expected = Piecewise((2, (x - y > S.Half) | (x - y > 1)),
                                          (2*y + 4, x - y > 1),
                                          (4*x + 2*y, True))
    actual = expr.simplify()
    print(f"Before: {expr}")
    print(f"Expected: {expected}, count_ops: {count_ops(expected)}")
    print(f"Actual: {actual}, count_ops: {count_ops(actual)}\n")

    print("File \"/home/runner/work/sympy/sympy/sympy/integrals/tests/test_integrals.py\", line 1121, "
          "in test_issue_4527")
    k, m = symbols('k m', integer=True)
    expr = integrate(sin(k * x) * sin(m * x), (x, 0, pi))
    expected = Piecewise((0, Eq(k, 0) | Eq(m, 0)),
                     (-pi / 2, Eq(k, -m) | (Eq(k, 0) & Eq(m, 0))),
                     (pi / 2, Eq(k, m) | (Eq(k, 0) & Eq(m, 0))),
                     (0, True))
    actual = expr.simplify()
    print(f"Before: {expr}")
    print(f"Expected: {expected}, count_ops: {count_ops(expected)}")
    print(f"Actual: {actual}, count_ops: {count_ops(actual)}\n")

    print("File \"/home/runner/work/sympy/sympy/sympy/logic/tests/test_boolalg.py\", line 316, "
          "in test_simplification_boolalg")
    expr = Or(x <= y, And(y > x, z))
    expected = x <= y
    actual = expr.simplify()
    print(f"Before: {expr}")
    print(f"Expected: {expected}, count_ops: {count_ops(expected)}")
    print(f"Actual: {actual}, count_ops: {count_ops(actual)}\n")

    print("File \"/home/runner/work/sympy/sympy/sympy/logic/tests/test_boolalg.py\", line 1111, "
          "in test_relational_simplification")
    expr = And(Eq(x, y), Eq(x, -y))
    expected = And(Eq(x, 0), Eq(y, 0))
    actual = expr.simplify()
    print(f"Before: {expr}")
    print(f"Expected: {expected}, count_ops: {count_ops(expected)}")
    print(f"Actual: {actual}, count_ops: {count_ops(actual)}\n")

    expr = (0 >= -2*(1 + 2**(pi/2) + 3**pi)**(1/pi) + (2**pi + 2**(3*pi/2) + 6**pi)**(1/pi))
    expected = True
    actual = expr.simplify()
    print(f"Before: {expr}")
    print(f"Expected: {expected}, count_ops: {count_ops(expected)}")
    print(f"Actual: {actual}, count_ops: {count_ops(actual)}\n")
collateral_damage_tests()
