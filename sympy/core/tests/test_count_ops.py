from sympy import symbols, sin, exp, cos, Derivative, Integral, Basic, \
    count_ops, S, And, I, pi, Eq

x, y, z = symbols('x,y,z')


def test_count_ops_non_visual():
    def count(val):
        return count_ops(val, visual=False)
    assert count(x) == 0
    assert count(x) is not S.Zero
    assert count(x + y) == 1
    assert count(x + y) is not S.One
    assert count(x + y*x + 2*y) == 4
    assert count({x + y: x}) == 1
    assert count({x + y: S(2) + x}) is not S.One


def test_count_ops_visual():
    ADD, MUL, POW, SIN, COS, EXP, AND, D, G = symbols(
        'Add Mul Pow sin cos exp And Derivative Integral'.upper())
    DIV, SUB, NEG = symbols('DIV SUB NEG')

    def count(val):
        return count_ops(val, visual=True)

    assert count(7) is S.Zero
    assert count(S(7)) is S.Zero
    assert count(-1) == NEG
    assert count(-2) == NEG
    assert count(S(2)/3) == DIV
    assert count(pi/3) == DIV
    assert count(-pi/3) == DIV + NEG
    assert count(I - 1) == SUB
    assert count(1 - I) == SUB
    assert count(1 - 2*I) == SUB + MUL

    assert count(x) is S.Zero
    assert count(-x) == NEG
    assert count(-2*x/3) == NEG + DIV + MUL
    assert count(1/x) == DIV
    assert count(1/(x*y)) == DIV + MUL
    assert count(-1/x) == NEG + DIV
    assert count(-2/x) == NEG + DIV
    assert count(x/y) == DIV
    assert count(-x/y) == NEG + DIV

    assert count(x**2) == POW
    assert count(-x**2) == POW + NEG
    assert count(-2*x**2) == POW + MUL + NEG

    assert count(x + pi/3) == ADD + DIV
    assert count(x + S(1)/3) == ADD + DIV
    assert count(x + y) == ADD
    assert count(x - y) == SUB
    assert count(y - x) == SUB
    assert count(-1/(x - y)) == DIV + NEG + SUB
    assert count(-1/(y - x)) == DIV + NEG + SUB
    assert count(1 + x**y) == ADD + POW
    assert count(1 + x + y) == 2*ADD
    assert count(1 + x + y + z) == 3*ADD
    assert count(1 + x**y + 2*x*y + y**2) == 3*ADD + 2*POW + 2*MUL
    assert count(2*z + y + x + 1) == 3*ADD + MUL
    assert count(2*z + y**17 + x + 1) == 3*ADD + MUL + POW
    assert count(2*z + y**17 + x + sin(x)) == 3*ADD + POW + MUL + SIN
    assert count(2*z + y**17 + x + sin(x**2)) == 3*ADD + MUL + 2*POW + SIN
    assert count(2*z + y**17 + x + sin(
        x**2) + exp(cos(x))) == 4*ADD + MUL + 2*POW + EXP + COS + SIN

    assert count(Derivative(x, x)) == D
    assert count(Integral(x, x) + 2*x/(1 + x)) == G + DIV + MUL + 2*ADD
    assert count(Basic()) is S.Zero

    assert count({x + 1: sin(x)}) == ADD + SIN
    assert count([x + 1, sin(x) + y, None]) == ADD + SIN + ADD
    assert count({x + 1: sin(x), y: cos(x) + 1}) == SIN + COS + 2*ADD
    assert count({}) is S.Zero
    assert count([x + 1, sin(x)*y, None]) == SIN + ADD + MUL
    assert count([]) is S.Zero

    # XXX: These are a bit surprising, only Expr-compatible ops are counted.
    assert count(And(x, y, z)) == 0
    assert count(Basic(x, x + y)) == ADD
    assert count(Eq(x + y, S(2))) == ADD
