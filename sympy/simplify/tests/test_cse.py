import itertools

from sympy import (Add, Pow, Symbol, exp, sqrt, symbols, sympify, cse,
                   Matrix, S, cos, sin, Eq, Function, Tuple, RootOf,
                   IndexedBase, Idx, Piecewise, O, Mul)
from sympy.simplify.cse_opts import sub_pre, sub_post
from sympy.simplify.cse_main import (replace_subsequence,
    shortest_repeated_subsequence, match_common_args, match_common_args_nc)
from sympy.functions.special.hyper import meijerg
from sympy.simplify import cse_main, cse_opts
from sympy.utilities.pytest import XFAIL, raises
from sympy.matrices import (MutableDenseMatrix, MutableSparseMatrix,
    ImmutableDenseMatrix, ImmutableSparseMatrix)
from sympy.matrices.expressions import MatrixSymbol, MatMul

from sympy.core.compatibility import range


w, x, y, z = symbols('w,x,y,z')
x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12 = symbols('x:13')


def test_numbered_symbols():
    ns = cse_main.numbered_symbols(prefix='y')
    assert list(itertools.islice(
        ns, 0, 10)) == [Symbol('y%s' % i) for i in range(0, 10)]
    ns = cse_main.numbered_symbols(prefix='y')
    assert list(itertools.islice(
        ns, 10, 20)) == [Symbol('y%s' % i) for i in range(10, 20)]
    ns = cse_main.numbered_symbols()
    assert list(itertools.islice(
        ns, 0, 10)) == [Symbol('x%s' % i) for i in range(0, 10)]

# Dummy "optimization" functions for testing.


def opt1(expr):
    return expr + y


def opt2(expr):
    return expr*z


def test_preprocess_for_cse():
    assert cse_main.preprocess_for_cse(x, [(opt1, None)]) == x + y
    assert cse_main.preprocess_for_cse(x, [(None, opt1)]) == x
    assert cse_main.preprocess_for_cse(x, [(None, None)]) == x
    assert cse_main.preprocess_for_cse(x, [(opt1, opt2)]) == x + y
    assert cse_main.preprocess_for_cse(
        x, [(opt1, None), (opt2, None)]) == (x + y)*z


def test_postprocess_for_cse():
    assert cse_main.postprocess_for_cse(x, [(opt1, None)]) == x
    assert cse_main.postprocess_for_cse(x, [(None, opt1)]) == x + y
    assert cse_main.postprocess_for_cse(x, [(None, None)]) == x
    assert cse_main.postprocess_for_cse(x, [(opt1, opt2)]) == x*z
    # Note the reverse order of application.
    assert cse_main.postprocess_for_cse(
        x, [(None, opt1), (None, opt2)]) == x*z + y


def test_cse_single():
    # Simple substitution.
    e = Add(Pow(x + y, 2), sqrt(x + y))
    substs, reduced = cse([e])
    assert substs == [(x0, x + y)]
    assert reduced == [sqrt(x0) + x0**2]


def test_cse_single2():
    # Simple substitution, test for being able to pass the expression directly
    e = Add(Pow(x + y, 2), sqrt(x + y))
    substs, reduced = cse(e)
    assert substs == [(x0, x + y)]
    assert reduced == [sqrt(x0) + x0**2]
    substs, reduced = cse(Matrix([[1]]))
    assert isinstance(reduced[0], Matrix)


def test_cse_not_possible():
    # No substitution possible.
    e = Add(x, y)
    substs, reduced = cse([e])
    assert substs == []
    assert reduced == [x + y]
    # issue 6329
    eq = (meijerg((1, 2), (y, 4), (5,), [], x) +
          meijerg((1, 3), (y, 4), (5,), [], x))
    assert cse(eq) == ([], [eq])


def test_nested_substitution():
    # Substitution within a substitution.
    e = Add(Pow(w*x + y, 2), sqrt(w*x + y))
    substs, reduced = cse([e])
    assert substs == [(x0, w*x + y)]
    assert reduced == [sqrt(x0) + x0**2]


def test_subtraction_opt():
    # Make sure subtraction is optimized.
    e = (x - y)*(z - y) + exp((x - y)*(z - y))
    substs, reduced = cse(
        [e], optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)])
    assert substs == [(x0, (x - y)*(y - z))]
    assert reduced == [-x0 + exp(-x0)]
    e = -(x - y)*(z - y) + exp(-(x - y)*(z - y))
    substs, reduced = cse(
        [e], optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)])
    assert substs == [(x0, (x - y)*(y - z))]
    assert reduced == [x0 + exp(x0)]
    # issue 4077
    n = -1 + 1/x
    e = n/x/(-n)**2 - 1/n/x
    assert cse(e, optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)]) == \
        ([], [0])


def test_multiple_expressions():
    e1 = (x + y)*z
    e2 = (x + y)*w
    substs, reduced = cse([e1, e2])
    assert substs == [(x0, x + y)]
    assert reduced == [x0*z, x0*w]
    l = [w*x*y + z, w*y]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    assert substs == rsubsts
    assert reduced == [z + x*x0, x0]
    l = [w*x*y, w*x*y + z, w*y]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    assert substs == rsubsts
    assert reduced == [x1, x1 + z, x0]
    l = [(x - z)*(y - z), x - z, y - z]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    assert substs == [(x0, -z), (x1, x + x0), (x2, x0 + y)]
    assert rsubsts == [(x0, -z), (x1, x0 + y), (x2, x + x0)]
    assert reduced == [x1*x2, x1, x2]
    l = [w*y + w + x + y + z, w*x*y]
    assert cse(l) == ([(x0, w*y)], [w + x + x0 + y + z, x*x0])
    assert cse([x + y, x + y + z]) == ([(x0, x + y)], [x0, z + x0])
    assert cse([x + y, x + z]) == ([], [x + y, x + z])
    assert cse([x*y, z + x*y, x*y*z + 3]) == \
        ([(x0, x*y)], [x0, z + x0, 3 + x0*z])

def test_non_commutative_cse():
    A, B, C = symbols('A B C', commutative=False)
    x0 = symbols('x0', commutative=False)

    l = [A*B*C, A*C]
    assert cse(l) == ([], l)
    l = [B+C, A*(B+C)]
    assert cse(l) == ([(x0, B+C)], [x0, A*x0])

    repl, exprs = cse(B*A*A + A*A)
    assert (repl, exprs) == ([(x0, A**2)], [B*x0 + x0])
    repl, exprs = cse(B*A*A - A*A*B)
    assert (repl, exprs) == ([(x0, A**2)], [B*x0 - x0*B])
    assert repl[0][0].is_commutative == False
    assert exprs[0] != 0

# cse on noncommutative Mul. When this works the test below it should be removed.
@XFAIL
def test_noncommutative_mul_cse():
    A, B, C = symbols('A B C', commutative=False)
    l = [A*B*C, A*B]
    repls, exprs = cse(l)
    assert (repls, exprs) == ([(x0, A*B)], [x0*C, x0])
    assert repls[0][0].is_commutative == False

# Test if CSE of non-commutative Mul terms is disabled
def test_bypass_non_commutatives():
    A, B, C = symbols('A B C', commutative=False)
    l = [A*B*C, A*C]
    assert cse(l) == ([], l)
    l = [A*B*C, A*B]
    assert cse(l) == ([], l)
    l = [B*C, A*B*C]
    assert cse(l) == ([], l)



@XFAIL
def test_powers():
    assert cse(x*y**2 + x*y) == ([(x0, x*y)], [x0*y + x0])


def test_issue_4498():
    assert cse(w/(x - y) + z/(y - x), optimizations='basic') == \
        ([], [(w - z)/(x - y)])


def test_issue_4020():
    assert cse(x**5 + x**4 + x**3 + x**2, optimizations='basic') \
        == ([(x0, x**2)], [x0*(x**3 + x + x0 + 1)])


def test_issue_4203():
    assert cse(sin(x**x)/x**x) == ([(x0, x**x)], [sin(x0)/x0])


def test_issue_6263():
    e = Eq(x*(-x + 1) + x*(x - 1), 0)
    assert cse(e, optimizations='basic') == ([], [True])


def test_dont_cse_tuples():
    from sympy import Subs
    f = Function("f")
    g = Function("g")

    name_val, (expr,) = cse(
        Subs(f(x, y), (x, y), (0, 1))
        + Subs(g(x, y), (x, y), (0, 1)))

    assert name_val == []
    assert expr == (Subs(f(x, y), (x, y), (0, 1))
            + Subs(g(x, y), (x, y), (0, 1)))

    name_val, (expr,) = cse(
        Subs(f(x, y), (x, y), (0, x + y))
        + Subs(g(x, y), (x, y), (0, x + y)))

    assert name_val == [(x0, x + y)]
    assert expr == Subs(f(x, y), (x, y), (0, x0)) + \
        Subs(g(x, y), (x, y), (0, x0))


def test_pow_invpow():
    assert cse(1/x**2 + x**2) == \
        ([(x0, x**2)], [x0 + 1/x0])
    assert cse(x**2 + (1 + 1/x**2)/x**2) == \
        ([(x0, x**2), (x1, 1/x0)], [x0 + x1*(x1 + 1)])
    assert cse(1/x**2 + (1 + 1/x**2)*x**2) == \
        ([(x0, x**2), (x1, 1/x0)], [x0*(x1 + 1) + x1])
    assert cse(cos(1/x**2) + sin(1/x**2)) == \
        ([(x0, x**(-2))], [sin(x0) + cos(x0)])
    assert cse(cos(x**2) + sin(x**2)) == \
        ([(x0, x**2)], [sin(x0) + cos(x0)])
    assert cse(y/(2 + x**2) + z/x**2/y) == \
        ([(x0, x**2)], [y/(x0 + 2) + z/(x0*y)])
    assert cse(exp(x**2) + x**2*cos(1/x**2)) == \
        ([(x0, x**2)], [x0*cos(1/x0) + exp(x0)])
    assert cse((1 + 1/x**2)/x**2) == \
        ([(x0, x**(-2))], [x0*(x0 + 1)])
    assert cse(x**(2*y) + x**(-2*y)) == \
        ([(x0, x**(2*y))], [x0 + 1/x0])


def test_postprocess():
    eq = (x + 1 + exp((x + 1)/(y + 1)) + cos(y + 1))
    assert cse([eq, Eq(x, z + 1), z - 2, (z + 1)*(x + 1)],
        postprocess=cse_main.cse_separate) == \
        [[(x1, y + 1), (x2, z + 1), (x, x2), (x0, x + 1)],
        [x0 + exp(x0/x1) + cos(x1), z - 2, x0*x2]]


def test_issue_4499():
    # previously, this gave 16 constants
    from sympy.abc import a, b
    B = Function('B')
    G = Function('G')
    t = Tuple(*
        (a, a + S(1)/2, 2*a, b, 2*a - b + 1, (sqrt(z)/2)**(-2*a + 1)*B(2*a -
        b, sqrt(z))*B(b - 1, sqrt(z))*G(b)*G(2*a - b + 1),
        sqrt(z)*(sqrt(z)/2)**(-2*a + 1)*B(b, sqrt(z))*B(2*a - b,
        sqrt(z))*G(b)*G(2*a - b + 1), sqrt(z)*(sqrt(z)/2)**(-2*a + 1)*B(b - 1,
        sqrt(z))*B(2*a - b + 1, sqrt(z))*G(b)*G(2*a - b + 1),
        (sqrt(z)/2)**(-2*a + 1)*B(b, sqrt(z))*B(2*a - b + 1,
        sqrt(z))*G(b)*G(2*a - b + 1), 1, 0, S(1)/2, z/2, -b + 1, -2*a + b,
        -2*a))
    c = cse(t)
    ans = (
        [(x0, 2*a), (x1, -b), (x2, x1 + 1), (x3, x0 + x2), (x4, sqrt(z)), (x5,
        B(x0 + x1, x4)), (x6, G(b)), (x7, G(x3)), (x8, -x0), (x9,
        (x4/2)**(x8 + 1)), (x10, x6*x7*x9*B(b - 1, x4)), (x11, x6*x7*x9*B(b,
        x4)), (x12, B(x3, x4))], [(a, a + S(1)/2, x0, b, x3, x10*x5,
        x11*x4*x5, x10*x12*x4, x11*x12, 1, 0, S(1)/2, z/2, x2, b + x8, x8)])
    assert ans == c


def test_issue_6169():
    r = RootOf(x**6 - 4*x**5 - 2, 1)
    assert cse(r) == ([], [r])
    # and a check that the right thing is done with the new
    # mechanism
    assert sub_post(sub_pre((-x - y)*z - x - y)) == -z*(x + y) - x - y


def test_cse_Indexed():
    len_y = 5
    y = IndexedBase('y', shape=(len_y,))
    x = IndexedBase('x', shape=(len_y,))
    Dy = IndexedBase('Dy', shape=(len_y-1,))
    i = Idx('i', len_y-1)

    expr1 = (y[i+1]-y[i])/(x[i+1]-x[i])
    expr2 = 1/(x[i+1]-x[i])
    replacements, reduced_exprs = cse([expr1, expr2])
    assert len(replacements) > 0


def test_cse_MatrixSymbol():
    # MatrixSymbols have non-Basic args, so make sure that works
    A = MatrixSymbol("A", 3, 3)
    assert cse(A) == ([], [A])

    n = symbols('n', integer=True)
    B = MatrixSymbol("B", n, n)
    assert cse(B) == ([], [B])

def test_cse_MatrixExpr():
    from sympy import MatrixSymbol
    A = MatrixSymbol('A', 3, 3)
    y = MatrixSymbol('y', 3, 1)

    expr1 = (A.T*A).I * A * y
    expr2 = (A.T*A) * A * y
    replacements, reduced_exprs = cse([expr1, expr2])
    assert len(replacements) > 0

    replacements, reduced_exprs = cse([expr1 + expr2, expr1])
    assert replacements

    replacements, reduced_exprs = cse([A**2, A + A**2])
    assert replacements

def test_Piecewise():
    f = Piecewise((-z + x*y, Eq(y, 0)), (-z - x*y, True))
    ans = cse(f)
    actual_ans = ([(x0, -z), (x1, x*y)], [Piecewise((x0+x1, Eq(y, 0)), (x0 - x1, True))])
    assert ans == actual_ans


def test_ignore_order_terms():
    eq = exp(x).series(x,0,3) + sin(y+x**3) - 1
    assert cse(eq) == ([], [sin(x**3 + y) + x + x**2/2 + O(x**3)])


def test_name_conflict():
    z1 = x0 + y
    z2 = x2 + x3
    l = [cos(z1) + z1, cos(z2) + z2, x0 + x2]
    substs, reduced = cse(l)
    assert [e.subs(reversed(substs)) for e in reduced] == l


def test_name_conflict_cust_symbols():
    z1 = x0 + y
    z2 = x2 + x3
    l = [cos(z1) + z1, cos(z2) + z2, x0 + x2]
    substs, reduced = cse(l, symbols("x:10"))
    assert [e.subs(reversed(substs)) for e in reduced] == l


def test_symbols_exhausted_error():
    l = cos(x+y)+x+y+cos(w+y)+sin(w+y)
    sym = [x, y, z]
    with raises(ValueError) as excinfo:
        cse(l, symbols=sym)


def test_issue_7840():
    # daveknippers' example
    C393 = sympify( \
        'Piecewise((C391 - 1.65, C390 < 0.5), (Piecewise((C391 - 1.65, \
        C391 > 2.35), (C392, True)), True))'
    )
    C391 = sympify( \
        'Piecewise((2.05*C390**(-1.03), C390 < 0.5), (2.5*C390**(-0.625), True))'
    )
    C393 = C393.subs('C391',C391)
    # simple substitution
    sub = {}
    sub['C390'] = 0.703451854
    sub['C392'] = 1.01417794
    ss_answer = C393.subs(sub)
    # cse
    substitutions,new_eqn = cse(C393)
    for pair in substitutions:
        sub[pair[0].name] = pair[1].subs(sub)
    cse_answer = new_eqn[0].subs(sub)
    # both methods should be the same
    assert ss_answer == cse_answer

    # GitRay's example
    expr = sympify(
        "Piecewise((Symbol('ON'), Equality(Symbol('mode'), Symbol('ON'))), \
        (Piecewise((Piecewise((Symbol('OFF'), StrictLessThan(Symbol('x'), \
        Symbol('threshold'))), (Symbol('ON'), S.true)), Equality(Symbol('mode'), \
        Symbol('AUTO'))), (Symbol('OFF'), S.true)), S.true))"
    )
    substitutions, new_eqn = cse(expr)
    # this Piecewise should be exactly the same
    assert new_eqn[0] == expr
    # there should not be any replacements
    assert len(substitutions) < 1


def test_issue_8891():
    for cls in (MutableDenseMatrix, MutableSparseMatrix,
            ImmutableDenseMatrix, ImmutableSparseMatrix):
        m = cls(2, 2, [x + y, 0, 0, 0])
        res = cse([x + y, m])
        ans = ([(x0, x + y)], [x0, cls([[x0, 0], [0, 0]])])
        assert res == ans
        assert isinstance(res[1][-1], cls)

def test_replace_subsequence():
    l = [1, 2, 3, 2, 3, 4]
    replace_subsequence(l, [2, 3], 5)
    assert l == [1, 5, 5, 4]

    l = [1, 2, 1, 2, 3]
    replace_subsequence(l, [1, 2, 3], 4)
    assert l == [1, 2, 4]

    l = [1, 2, 1, 2, 1, 2]
    replace_subsequence(l, [1, 2], 3)
    assert l == [3, 3, 3]


def test_shortest_repeated_subsequence():
    a, b, c = 'aabacde', 'aabade', 'abcde'
    abc = '$'.join([a, b, c])
    assert shortest_repeated_subsequence(abc, ignore='$') in ['ab', 'de']

    aa = '$'.join([a, a])
    assert shortest_repeated_subsequence(aa, ignore='$') == a

    ab = '$'.join([a, b])
    assert shortest_repeated_subsequence(ab, ignore='$') == 'de'

    assert shortest_repeated_subsequence('xxx') == 'xx'
    assert shortest_repeated_subsequence('xxxx') == 'xx'
    assert shortest_repeated_subsequence('xxxxx') == 'xx'

    def _cse_str(S):
        """
        A mock version of the cse algorithm using strings

        """
        S = '$'.join(S)
        n = 0
        repls = {}
        while True:
            rep = shortest_repeated_subsequence(S, ignore='$')
            if not rep:
                break
            repls[str(n)] = rep
            S = S.replace(rep, str(n))
            n += 1

        return repls, S.split('$')

    def _replace_cse(strs, replacement_dict):
        """
        Undo _cse_str.
        """
        for i in sorted(replacement_dict, reverse=True):
            strs = [x.replace(i, replacement_dict[i]) for x in strs]
        return strs

    repl, strs = _cse_str([a, b, c])
    assert (repl, strs) == ({'0': 'de', '1': 'c0', '2': 'ab', '3': 'a2a'}, ['31', '30', '21'])
    assert _replace_cse(strs, repl) == [a, b, c]

    S = ['sandollar', 'sandlot', 'handler', 'grand', 'pantry']
    repl, strs = _cse_str(S)
    assert (repl, strs) == ({'0': 'an', '1': '0d', '2': '1l'}, ['s1ollar',
        's2ot', 'h2er', 'gr1', 'p0try'])
    assert _replace_cse(strs, repl) == S

    repl, strs = _cse_str(['xxx'])
    assert (repl, strs) == ({'0': 'xx'}, ['0x'])
    assert _replace_cse(strs, repl) == ['xxx']

    repl, strs = _cse_str(['xxxx'])
    assert (repl, strs) == ({'0': 'xx'}, ['00'])
    assert _replace_cse(strs, repl) == ['xxxx']

    repl, strs = _cse_str(['xxxxx'])
    assert (repl, strs) == ({'0': 'xx'}, ['00x'])
    assert _replace_cse(strs, repl) == ['xxxxx']

    repl, strs = _cse_str(['xxxxxx'])
    assert (repl, strs) == ({'0': 'xx', '1': '00'}, ['10'])
    assert _replace_cse(strs, repl) == ['xxxxxx']

def test_match_common_args():
    # Test that evaluate=False actually works
    assert Mul(x, Mul(y, z, evaluate=False), evaluate=False).args == (x, Mul(y, z,
        evaluate=False))

    opt_subs = match_common_args(Mul, [x*y*z, x*y*w])
    assert set(opt_subs.keys()) == set([x*y*z, x*y*w])

    # The evaluate=False Muls could be in any order. cse uses a custom
    # substitution routine that handles this.
    assert set(opt_subs[x*y*z].args) in [set([Mul(x, y, evaluate=False), z]),
        set([Mul(y, x, evaluate=False), z])]
    assert set(opt_subs[x*y*w].args) in [set([Mul(x, y, evaluate=False), w]),
        set([Mul(y, x, evaluate=False), w])]

def test_match_common_args_nc():
    A, B, C, D, E = symbols('A B C D E', commutative=False)

    # We have to use evaluate=False to prevent A*A from becoming A**2. This is
    # handled in cse as a pre-transformation.
    m = lambda *x: Mul(*x, evaluate=False)

    # Test that evaluate=False actually works
    assert m(A, m(B, C)).args == (A, m(B, C))

    a = m(A, A, B, A, C, D, E)
    b = m(A, A, B, A, D, E)
    c = m(A, B, C, D, E)

    opt_subs = match_common_args_nc(Mul, [a, b, c])

    # This is the same example from test_shortest_repeated_subsequence()
    # above.

    # ({'0': 'de', '1': 'c0', '2': 'ab', '3': 'a2a'}, ['31', '30', '21'])

    assert opt_subs == {
        a: m(m(A, m(A, B), A), m(C, m(D, E))),
        b: m(m(A, m(A, B), A), m(D, E)),
        c: m(m(A, B), m(C, m(D, E))),
    }

    X = Symbol("X", commutative=False)

    expr = m(X, X, X, X, X)
    opt_subs = match_common_args_nc(Mul, [expr])

    assert opt_subs == {expr: m(m(X, X), m(X, X), X)}

    n = symbols('n', integer=True)
    M = MatrixSymbol('M', n, n)
    N = MatrixSymbol('N', n, n)
    x = MatrixSymbol('x', n, 1)

    # Make sure evaluate=False works with MatMul
    assert MatMul(M, MatMul(N, x, evaluate=False), evaluate=False).args ==\
        (M, MatMul(N, x, evaluate=False))

    expr1 = N*M*x
    expr2 = N*M
    opt_subs = match_common_args_nc(MatMul, [expr1, expr2])

    assert opt_subs == {
        expr1: MatMul(MatMul(N, M, evaluate=False), x, evaluate=False),
        expr2: MatMul(N, M, evaluate=False),
    }
