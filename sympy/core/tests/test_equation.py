from sympy.testing.pytest import raises
from sympy import functions, symbols, Equation

#####
# Testing that sympy functions work with Equations
#####

# Overridden elsewhere
_extended_ = ('sqrt', 'root')

# Either not applicable to Equations or have not yet figured out a way
# to systematically apply to an Equation.
# TODO examine these more carefully (top priority: real_root, Ynm_c).
_not_applicable_to_equations_ = ('Min', 'Max', 'Id', 'real_root',
        'unbranched_argument', 'polarify', 'unpolarify',
        'piecewise_fold', 'E1', 'Eijk', 'bspline_basis',
        'bspline_basis_set', 'interpolating_spline', 'jn_zeros',
        'jacobi_normalized', 'Ynm_c', 'piecewise_exclusive', 'Piecewise',
        'motzkin', 'hyper','meijerg', 'chebyshevu_root', 'chebyshevt_root',
        'betainc_regularized')
_skip_ = _extended_ + _not_applicable_to_equations_

import importlib
temp = importlib.import_module('sympy', package="functions")
for func in functions.__all__:
    globals()[func] = getattr(temp, func)


def test_functions_extensions():
    from inspect import signature
    failures = []
    a, b, c = symbols('a b c')
    eq = Equation(a, b / c)
    n = symbols('n', positive=True, integer=True)
    for func in functions.__all__:
        if func not in _skip_ or func in _extended_:
            obj = globals()[func]
            sig = signature(obj).parameters
            if func == 'betainc' or func == 'betainc_regularized':
                # The signature is undefined need 4 complex numbers:
                # a, b, x1, x2.
                sig = {'arg1': 'a', 'arg2': 'b', 'arg3': 'x1',
                       'arg4': 'x2'}
            keylist = [key for key in sig]
            tempargs = [eq]
            largs = [eq.lhs]
            rargs = [eq.rhs]
            for key in sig:
                if (str(sig[key]).find("=")) == -1 and (str(sig[key]).
                    find("**")) == -1 and key != keylist[0]:
                    tempargs.append(n)
                    largs.append(n)
                    rargs.append(n)
            try:
                tst = obj(*tempargs)
                if not (tst == Equation(obj(*largs), obj(*rargs))):
                    failures.append(func + ' extended but did not work.')
            except Exception as e:
                failures.append(str(func) + ': ' + str(e))
    assert (failures == [])
    pass


def test_functions_extensions_eqn_not_arg1():
    from inspect import signature
    failures = []
    a, b, c = symbols('a b c')
    eq = Equation(a, b / c)
    n = symbols('n', positive=True, integer=True)
    for func in functions.__all__:
        if func not in _skip_ or func in _extended_:
            obj = globals()[func]
            sig = signature(obj).parameters
            if func == 'betainc' or func == 'betainc_regularized':
                # The signature is undefined need 4 complex numbers:
                # a, b, x1, x2.
                sig = {'arg1': 'a', 'arg2': 'b', 'arg3': 'x1',
                       'arg4': 'x2'}
            keylist = [key for key in sig]
            for j in range(1, len(sig)):
                tempargs = [n]
                largs = [n]
                rargs = [n]
                for k in range(1, len(sig)):
                    if ((str(sig[keylist[k]]).find("=")) == -1 and
                        (str(sig[keylist[k]]).find("**")) == -1):
                        if k == j:
                            tempargs.append(eq)
                            largs.append(eq.lhs)
                            rargs.append(eq.rhs)
                        else:
                            tempargs.append(n)
                            largs.append(n)
                            rargs.append(n)
                try:
                    tst = obj(*tempargs)
                    if (isinstance(tst, Equation) and not
                    (tst == Equation(obj(*largs), obj(*rargs)))):
                        failures.append(func + '(' + str(*tempargs) + ') ' \
                                                                      'extended but did not work.')
                except Exception as e:
                    failures.append(str(func) + ': ' + str(e))
    assert (failures == [])
    pass


    def test_two_eqn():
        a, b, c = symbols('a b c')
        eq = Equation(a, b / c)
        obj = globals()['besselj']
        raises(NotImplementedError, lambda: obj(eq, eq))

#####
# Testing of Equation class
#####
from sympy import integrate, simplify, expand, factor, Add
from sympy import diff, FiniteSet, Equality, functions, Matrix, S
from sympy import sin, cos, log, exp, I, collect, Eqn
from sympy import sqrt, root, Heaviside

def test_define_equation():
    a, b, c = symbols('a b c')
    raises(TypeError, lambda: Equation(FiniteSet(a), FiniteSet(b, c)))
    assert(Equation(1, 0).check() == False)
    assert Eqn(1, 0) == Equation(1, 0)
    tsteqn = Equation(a, b/c)
    assert tsteqn.args == (a, b/c)
    assert tsteqn.lhs == a
    assert tsteqn.rhs == b/c
    assert tsteqn.free_symbols == {a, b, c}

def test_convert_equation():
    a, b, c = symbols('a b c')
    tsteqn = Equation(a, b/c)
    assert tsteqn.as_Boolean() == Equality(a, b/c)
    assert tsteqn.reversed == Equation(b/c, a)
    assert tsteqn.swap == Equation(b/c, a)


def test_binary_op():
    a, b, c = symbols('a b c')
    tsteqn = Equation(a, b/c)
    assert tsteqn + c == Equation(a + c, b/c + c)
    assert c + tsteqn == Equation(c + a, c + b/c)
    assert tsteqn*c == Equation(a*c, b)
    assert c*tsteqn == Equation(c*a, b)
    assert tsteqn - c == Equation(a - c, b/c - c)
    assert c - tsteqn == Equation(c - a, c - b/c)
    assert tsteqn/ c == Equation(a/c, b/c**2)
    assert c/tsteqn == Equation(c/a, c**2/b)
    assert tsteqn % c == Equation(a % c, (b/c) % c)
    assert c % tsteqn == Equation(c % a, c % (b/c))
    assert tsteqn**c == Equation(a**c, (b/c)**c)
    assert c**tsteqn == Equation(c**a, c**(b/c))
    assert tsteqn + tsteqn == Equation(2*a, 2*b/c)
    assert tsteqn*tsteqn == Equation(a**2, b**2/c**2)
    assert tsteqn - tsteqn == Equation(0, 0)
    assert tsteqn/tsteqn == Equation(1, 1)
    assert tsteqn % tsteqn == Equation(0, 0)
    assert tsteqn**tsteqn == Equation(a**a, (b/c)**(b/c))
    assert tsteqn**a == Equation(a**a, (b/c)**a)
    assert tsteqn._eval_power(tsteqn) == Equation(a**a, (b/c)**(b/c))
    assert tsteqn._eval_power(a) == Equation(a**a, (b/c)**a)

def test_sympy_functions():
    a, b, c = symbols('a b c')
    tsteqn = Equation(a, b/c)
    assert sin(tsteqn) == Equation(sin(a),sin(b/c))
    assert log(tsteqn) == Equation(log(a),log(b/c))
    # Check matrix exponentiation is not overridden.
    assert exp(tsteqn) == Equation(exp(tsteqn.lhs),exp(tsteqn.rhs))
    tsteqn5 = Equation(a, Matrix([[1, 1], [1, 1]]))
    assert exp(tsteqn5).lhs == exp(a)
    assert exp(tsteqn5).rhs == exp(Matrix([[1, 1], [1, 1]]))

def test_helper_functions():
    a, b, c, x= symbols('a b c x')
    tsteqn = Equation(a, b/c)
    raises(ValueError, lambda: integrate(tsteqn, c))
    raises(AttributeError, lambda: integrate(tsteqn, c, side='right'))
    assert tsteqn.evalf(4, {b: 2.0, c: 4}) == Equation(a, 0.5000)
    assert diff(tsteqn, c) == Equation(diff(a, c, evaluate=False), -b/c**2)
    tsteqn = Equation(a*c, b/c)
    assert diff(tsteqn, c) == Equation(a, -b/c**2)
    assert integrate(tsteqn, c, side='rhs') == integrate(tsteqn.rhs, c)
    assert integrate(tsteqn, c, side='lhs') == integrate(tsteqn.lhs, c)

    def adsq(eqn):
        # Arbitrary python function
        return eqn + eqn**2

    assert adsq(Equation(a*c, b/c)) == Equation(a**2*c**2 + a*c, b**2/c**2 +
                                                b/c)
    assert Equation((a - 1)*(a + 1), (2*b + c)**2).expand() == Equation(
        a**2 - 1, 4*b**2 + 4*b*c + c**2)
    assert expand(Equation((a - 1)*(a + 1), (2*b + c)**2)) == Equation(
        a**2 - 1, 4*b**2 + 4*b*c + c**2)
    assert Equation(a**2 - 1, 4*b**2 + 4*b*c + c**2).factor() == Equation(
        (a - 1)*(a + 1), (2*b + c)**2)
    assert factor(Equation(a**2 - 1, 4*b**2 + 4*b*c + c**2)) == Equation(
        (a - 1)*(a + 1), (2*b + c)**2)
    assert Equation(a**2 - 1, 4*b**2 + 4*b*c + c*a).collect(c) == Equation(
        a**2- 1, 4*b**2 + c*(a + 4*b))
    assert collect(Equation(a**2 - 1, 4*b**2 + 4*b*c + c*a), c) == Equation(
        a**2- 1, 4*b**2 + c*(a + 4*b))
    assert Equation((a + 1)**2/(a + 1), exp(log(c))).simplify() == Equation(
        a + 1, c)
    assert simplify(Equation((a + 1)**2/(a + 1), exp(log(c)))) == Equation(
        a + 1, c)
    assert root(Eqn(a,b/c),3) == Equation(a**(S(1)/S(3)), (b/c)**(S(1)/S(3)))
    assert root(b/c,3) == (b/c)**(S(1)/S(3))
    assert sqrt(Eqn(a,b/c)) == Equation(sqrt(a), sqrt(b/c))

def test_Heaviside():
    a, b, c, x = symbols('a b c x')
    tsteqn = Equation(a, b / c)
    assert (Heaviside(tsteqn) ==
            Equation(Heaviside(tsteqn.lhs), Heaviside(tsteqn.rhs)))
    assert Heaviside(0) == S(1)/S(2)

def test_apply_syntax():
    a, b, c, x = symbols('a b c x')
    tsteqn = Equation(a, b/c)
    assert tsteqn.apply(log) == Equation(log(a), log(b/c))
    assert tsteqn.applylhs(log) == Equation(log(a), b / c)
    assert tsteqn.applyrhs(log) == Equation(a, log(b / c))
    poly = Equation(a*x**2 + b*x + c*x**2, a*x**3 + b*x**3 + c*x)
    assert poly.applyrhs(collect, x) == Equation(poly.lhs, poly.rhs.collect(x))


def test_do_syntax():
    a, b, c, x = symbols('a b c x')
    tsteqn = Equation(a, b/c)
    raises(AttributeError, lambda: tsteqn.do.log())
    poly = Equation(a*x**2 + b*x + c*x**2, a*x**3 + b*x**3 + c*x)
    assert poly.dorhs.collect(x) == Eqn(poly.lhs, poly.rhs.collect(x))
    assert poly.dolhs.collect(x) == Eqn(poly.lhs.collect(x), poly.rhs)
    assert poly.do.collect(x) == Eqn(poly.lhs.collect(x), poly.rhs.collect(x))


def test_rewrite_add():
    b, x = symbols("x, b")
    eq = Equation(x + b, x - b)
    assert eq.rewrite(Add) == Equation(2 * b, 0)
    assert set(eq.rewrite(Add, evaluate=None).lhs.args) == set({b, x, b, -x})
    assert set(eq.rewrite(Add, evaluate=False).lhs.args) == set({b, x, b, -x})
    assert eq.rewrite(Add, eqn=False) == 2 * b
    assert set(eq.rewrite(Add, eqn=False, evaluate=False).args) == set({b, x,
                                                                        b, -x})


def test_rewrite():
    x = symbols("x")
    eq = Equation(exp(I*x),cos(x) + I*sin(x))

    # NOTE: Must use `sexp` otherwise the test is going to fail.
    # This reflects the fact that rewrite pulls the fuction exp internally
    # from the definitions of functions in sympy and not from the globally
    # redefined functions that are Equation aware.
    from sympy import exp as sexp
    assert eq.rewrite(exp) == Equation(exp(I*x), sexp(I*x))
    assert eq.rewrite(Add) == Equation(exp(I*x) - I*sin(x) - cos(x), 0)


def test_subs():
    a, b, c, x = symbols('a b c x')
    eq1 = Equation(x + a + b + c, x * a * b * c)
    eq2 = Equation(x + a, 4)
    assert eq1.subs(a, 2) == Equation(x + b + c + 2, 2 * x * b * c)
    assert eq1.subs([(a, 2), (b, 3)]) == Equation(x + c + 5, 6 * x * c)
    assert eq1.subs({a: 2, b: 3}) == Equation(x + c + 5, 6 * x * c)
    assert eq1.subs(eq2) == Equation(4 + b + c, x * a * b * c)

    # verify that proper errors are raised
    eq3 = Equation(b, 5)
    raises(TypeError, lambda: eq1.subs([eq2, eq3]))
    raises(ValueError, lambda: eq1.subs(eq2, {b: 5}))

    # verify that substituting an Equation into an expression is not supported
    raises(ValueError, lambda: eq1.dolhs.subs(eq2))
    raises(ValueError, lambda: eq1.dorhs.subs(eq2))
    raises(ValueError, lambda: (x + a + b + c).subs(eq2))

    # verify the effectiveness of `simultaneous`
    eq = Equation((x + a) / a, b * c)
    sd = {x + a: a, a: x + a}
    assert eq.subs(sd) == Equation(1, b * c)
    assert eq.subs(sd, simultaneous=True) == Equation(a / (x + a), b * c)
