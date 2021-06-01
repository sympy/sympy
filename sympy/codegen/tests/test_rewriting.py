from sympy import log, exp, cos, Symbol, Pow, sin, MatrixSymbol, sinc
from sympy.assumptions import assuming, Q
from sympy.printing import ccode
from sympy.codegen.matrix_nodes import MatrixSolve
from sympy.codegen.cfunctions import log2, exp2, expm1, log1p
from sympy.codegen.numpy_nodes import logaddexp, logaddexp2
from sympy.codegen.scipy_nodes import cosm1
from sympy.codegen.rewriting import (
    optimize, cosm1_opt, log2_opt, exp2_opt, expm1_opt, log1p_opt, optims_c99,
    create_expand_pow_optimization, matinv_opt, logaddexp_opt, logaddexp2_opt,
    optims_numpy, sinc_opts
)
from sympy.testing.pytest import XFAIL


def test_log2_opt():
    x = Symbol('x')
    expr1 = 7*log(3*x + 5)/(log(2))
    opt1 = optimize(expr1, [log2_opt])
    assert opt1 == 7*log2(3*x + 5)
    assert opt1.rewrite(log) == expr1

    expr2 = 3*log(5*x + 7)/(13*log(2))
    opt2 = optimize(expr2, [log2_opt])
    assert opt2 == 3*log2(5*x + 7)/13
    assert opt2.rewrite(log) == expr2

    expr3 = log(x)/log(2)
    opt3 = optimize(expr3, [log2_opt])
    assert opt3 == log2(x)
    assert opt3.rewrite(log) == expr3

    expr4 = log(x)/log(2) + log(x+1)
    opt4 = optimize(expr4, [log2_opt])
    assert opt4 == log2(x) + log(2)*log2(x+1)
    assert opt4.rewrite(log) == expr4

    expr5 = log(17)
    opt5 = optimize(expr5, [log2_opt])
    assert opt5 == expr5

    expr6 = log(x + 3)/log(2)
    opt6 = optimize(expr6, [log2_opt])
    assert str(opt6) == 'log2(x + 3)'
    assert opt6.rewrite(log) == expr6


def test_exp2_opt():
    x = Symbol('x')
    expr1 = 1 + 2**x
    opt1 = optimize(expr1, [exp2_opt])
    assert opt1 == 1 + exp2(x)
    assert opt1.rewrite(Pow) == expr1

    expr2 = 1 + 3**x
    assert expr2 == optimize(expr2, [exp2_opt])


def test_expm1_opt():
    x = Symbol('x')

    expr1 = exp(x) - 1
    opt1 = optimize(expr1, [expm1_opt])
    assert expm1(x) - opt1 == 0
    assert opt1.rewrite(exp) == expr1

    expr2 = 3*exp(x) - 3
    opt2 = optimize(expr2, [expm1_opt])
    assert 3*expm1(x) == opt2
    assert opt2.rewrite(exp) == expr2

    expr3 = 3*exp(x) - 5
    assert expr3 == optimize(expr3, [expm1_opt])

    expr4 = 3*exp(x) + log(x) - 3
    opt4 = optimize(expr4, [expm1_opt])
    assert 3*expm1(x) + log(x) == opt4
    assert opt4.rewrite(exp) == expr4

    expr5 = 3*exp(2*x) - 3
    opt5 = optimize(expr5, [expm1_opt])
    assert 3*expm1(2*x) == opt5
    assert opt5.rewrite(exp) == expr5


@XFAIL  # ideally this test should pass: need to improve `expm1_opt`
def test_expm1_two_exp_terms():
    x, y = map(Symbol, 'x y'.split())
    expr1 = exp(x) + exp(y) - 2
    opt1 = optimize(expr1, [expm1_opt])
    assert opt1 == expm1(x) + expm1(y)


def test_cosm1_opt():
    x = Symbol('x')

    expr1 = cos(x) - 1
    opt1 = optimize(expr1, [cosm1_opt])
    assert cosm1(x) - opt1 == 0
    assert opt1.rewrite(cos) == expr1

    expr2 = 3*cos(x) - 3
    opt2 = optimize(expr2, [cosm1_opt])
    assert 3*cosm1(x) == opt2
    assert opt2.rewrite(cos) == expr2

    expr3 = 3*cos(x) - 5
    assert expr3 == optimize(expr3, [cosm1_opt])

    expr4 = 3*cos(x) + log(x) - 3
    opt4 = optimize(expr4, [cosm1_opt])
    assert 3*cosm1(x) + log(x) == opt4
    assert opt4.rewrite(cos) == expr4

    expr5 = 3*cos(2*x) - 3
    opt5 = optimize(expr5, [cosm1_opt])
    assert 3*cosm1(2*x) == opt5
    assert opt5.rewrite(cos) == expr5


@XFAIL  # ideally this test should pass: need to improve `cosm1_opt`
def test_cosm1_two_cos_terms():
    x, y = map(Symbol, 'x y'.split())
    expr1 = cos(x) + cos(y) - 2
    opt1 = optimize(expr1, [cosm1_opt])
    assert opt1 == cosm1(x) + cosm1(y)


@XFAIL  # ideally this test should pass: need to add a new combined expm1_cosm1_opt?
def test_expm1_cosm1_mixed():
    x = Symbol('x')
    expr1 = exp(x) + cos(x) - 2
    opt1 = optimize(expr1, [expm1_opt, cosm1_opt])  # need a combined opt pass?
    assert opt1 == cosm1(x) + expm1(x)


def test_log1p_opt():
    x = Symbol('x')
    expr1 = log(x + 1)
    opt1 = optimize(expr1, [log1p_opt])
    assert log1p(x) - opt1 == 0
    assert opt1.rewrite(log) == expr1

    expr2 = log(3*x + 3)
    opt2 = optimize(expr2, [log1p_opt])
    assert log1p(x) + log(3) == opt2
    assert (opt2.rewrite(log) - expr2).simplify() == 0

    expr3 = log(2*x + 1)
    opt3 = optimize(expr3, [log1p_opt])
    assert log1p(2*x) - opt3 == 0
    assert opt3.rewrite(log) == expr3

    expr4 = log(x+3)
    opt4 = optimize(expr4, [log1p_opt])
    assert str(opt4) == 'log(x + 3)'


def test_optims_c99():
    x = Symbol('x')

    expr1 = 2**x + log(x)/log(2) + log(x + 1) + exp(x) - 1
    opt1 = optimize(expr1, optims_c99).simplify()
    assert opt1 == exp2(x) + log2(x) + log1p(x) + expm1(x)
    assert opt1.rewrite(exp).rewrite(log).rewrite(Pow) == expr1

    expr2 = log(x)/log(2) + log(x + 1)
    opt2 = optimize(expr2, optims_c99)
    assert opt2 == log2(x) + log1p(x)
    assert opt2.rewrite(log) == expr2

    expr3 = log(x)/log(2) + log(17*x + 17)
    opt3 = optimize(expr3, optims_c99)
    delta3 = opt3 - (log2(x) + log(17) + log1p(x))
    assert delta3 == 0
    assert (opt3.rewrite(log) - expr3).simplify() == 0

    expr4 = 2**x + 3*log(5*x + 7)/(13*log(2)) + 11*exp(x) - 11 + log(17*x + 17)
    opt4 = optimize(expr4, optims_c99).simplify()
    delta4 = opt4 - (exp2(x) + 3*log2(5*x + 7)/13 + 11*expm1(x) + log(17) + log1p(x))
    assert delta4 == 0
    assert (opt4.rewrite(exp).rewrite(log).rewrite(Pow) - expr4).simplify() == 0

    expr5 = 3*exp(2*x) - 3
    opt5 = optimize(expr5, optims_c99)
    delta5 = opt5 - 3*expm1(2*x)
    assert delta5 == 0
    assert opt5.rewrite(exp) == expr5

    expr6 = exp(2*x) - 3
    opt6 = optimize(expr6, optims_c99)
    delta6 = opt6 - (exp(2*x) - 3)
    assert delta6 == 0

    expr7 = log(3*x + 3)
    opt7 = optimize(expr7, optims_c99)
    delta7 = opt7 - (log(3) + log1p(x))
    assert delta7 == 0
    assert (opt7.rewrite(log) - expr7).simplify() == 0

    expr8 = log(2*x + 3)
    opt8 = optimize(expr8, optims_c99)
    assert opt8 == expr8


def test_create_expand_pow_optimization():
    cc = lambda x: ccode(
        optimize(x, [create_expand_pow_optimization(4)]))
    x = Symbol('x')
    assert cc(x**4) == 'x*x*x*x'
    assert cc(x**4 + x**2) == 'x*x + x*x*x*x'
    assert cc(x**5 + x**4) == 'pow(x, 5) + x*x*x*x'
    assert cc(sin(x)**4) == 'pow(sin(x), 4)'
    # gh issue 15335
    assert cc(x**(-4)) == '1.0/(x*x*x*x)'
    assert cc(x**(-5)) == 'pow(x, -5)'
    assert cc(-x**4) == '-x*x*x*x'
    assert cc(x**4 - x**2) == '-x*x + x*x*x*x'
    i = Symbol('i', integer=True)
    assert cc(x**i - x**2) == 'pow(x, i) - x*x'
    # gh issue 20753
    cc2 = lambda x: ccode(optimize(x, [create_expand_pow_optimization(
        4, base_req=lambda b: b.is_Function)]))
    assert cc2(x**3 + sin(x)**3) == "pow(x, 3) + sin(x)*sin(x)*sin(x)"

def test_matsolve():
    n = Symbol('n', integer=True)
    A = MatrixSymbol('A', n, n)
    x = MatrixSymbol('x', n, 1)

    with assuming(Q.fullrank(A)):
        assert optimize(A**(-1) * x, [matinv_opt]) == MatrixSolve(A, x)
        assert optimize(A**(-1) * x + x, [matinv_opt]) == MatrixSolve(A, x) + x


def test_logaddexp_opt():
    x, y = map(Symbol, 'x y'.split())
    expr1 = log(exp(x) + exp(y))
    opt1 = optimize(expr1, [logaddexp_opt])
    assert logaddexp(x, y) - opt1 == 0
    assert logaddexp(y, x) - opt1 == 0
    assert opt1.rewrite(log) == expr1


def test_logaddexp2_opt():
    x, y = map(Symbol, 'x y'.split())
    expr1 = log(2**x + 2**y)/log(2)
    opt1 = optimize(expr1, [logaddexp2_opt])
    assert logaddexp2(x, y) - opt1 == 0
    assert logaddexp2(y, x) - opt1 == 0
    assert opt1.rewrite(log) == expr1


def test_sinc_opts():
    def check(d):
        for k, v in d.items():
            assert optimize(k, sinc_opts) == v

    x = Symbol('x')
    check({
        sin(x)/x       : sinc(x),
        sin(2*x)/(2*x) : sinc(2*x),
        sin(3*x)/x     : 3*sinc(3*x),
        x*sin(x)       : x*sin(x)
    })

    y = Symbol('y')
    check({
        sin(x*y)/(x*y)       : sinc(x*y),
        y*sin(x/y)/x         : sinc(x/y),
        sin(sin(x))/sin(x)   : sinc(sin(x)),
        sin(3*sin(x))/sin(x) : 3*sinc(3*sin(x)),
        sin(x)/y             : sin(x)/y
    })


def test_optims_numpy():
    def check(d):
        for k, v in d.items():
            assert optimize(k, optims_numpy) == v

    x = Symbol('x')
    check({
        sin(2*x)/(2*x) + exp(2*x) - 1: sinc(2*x) + expm1(2*x),
        log(x+3)/log(2) + log(x**2 + 1): log1p(x**2) + log2(x+3)
    })


@XFAIL  # room for improvement, ideally this test case should pass.
def test_optims_numpy_TODO():
    def check(d):
        for k, v in d.items():
            assert optimize(k, optims_numpy) == v

    x, y = map(Symbol, 'x y'.split())
    check({
        log(x*y)*sin(x*y)*log(x*y+1)/(log(2)*x*y): log2(x*y)*sinc(x*y)*log1p(x*y),
        exp(x*sin(y)/y) - 1: expm1(x*sinc(y))
    })
