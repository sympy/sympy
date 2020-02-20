from sympy.external import import_module
from sympy.utilities.decorator import doctest_depends_on
from sympy.core import Integer, Float
from sympy import Pow, Add, Integral, Mul, S, Function, E
from sympy.functions import exp as sym_exp
import inspect
import re
from sympy import powsimp
matchpy = import_module("matchpy")

if matchpy:
    from matchpy import ManyToOneReplacer, ManyToOneMatcher
    from sympy.integrals.rubi.utility_function import (
        rubi_exp, rubi_unevaluated_expr, process_trig
    )

    from sympy.utilities.matchpy_connector import op_iter, op_len

    @doctest_depends_on(modules=('matchpy',))
    def get_rubi_object():
        """
        Returns rubi ManyToOneReplacer by adding all rules from different modules.

        Uncomment the lines to add integration capabilities of that module.

        Currently, there are parsing issues with special_function,
        derivative and miscellaneous_integration. Hence they are commented.
        """
        from sympy.integrals.rubi.rules.integrand_simplification import integrand_simplification
        from sympy.integrals.rubi.rules.linear_products import linear_products
        from sympy.integrals.rubi.rules.quadratic_products import quadratic_products
        from sympy.integrals.rubi.rules.binomial_products import binomial_products
        from sympy.integrals.rubi.rules.trinomial_products import trinomial_products
        from sympy.integrals.rubi.rules.miscellaneous_algebraic import miscellaneous_algebraic
        from sympy.integrals.rubi.rules.exponential import exponential
        from sympy.integrals.rubi.rules.logarithms import logarithms
        from sympy.integrals.rubi.rules.sine import sine
        from sympy.integrals.rubi.rules.tangent import tangent
        from sympy.integrals.rubi.rules.secant import secant
        from sympy.integrals.rubi.rules.miscellaneous_trig import miscellaneous_trig
        from sympy.integrals.rubi.rules.inverse_trig import inverse_trig
        from sympy.integrals.rubi.rules.hyperbolic import hyperbolic
        from sympy.integrals.rubi.rules.inverse_hyperbolic import inverse_hyperbolic
        from sympy.integrals.rubi.rules.special_functions import special_functions
        #from sympy.integrals.rubi.rules.derivative import derivative
        #from sympy.integrals.rubi.rules.piecewise_linear import piecewise_linear
        from sympy.integrals.rubi.rules.miscellaneous_integration import miscellaneous_integration

        rules = []

        rules += integrand_simplification()
        rules += linear_products()
        rules += quadratic_products()
        rules += binomial_products()
        rules += trinomial_products()
        rules += miscellaneous_algebraic()
        rules += exponential()
        rules += logarithms()
        rules += special_functions()
        rules += sine()
        rules += tangent()
        rules += secant()
        rules += miscellaneous_trig()
        rules += inverse_trig()
        rules += hyperbolic()
        rules += inverse_hyperbolic()
        #rubi = piecewise_linear(rubi)
        rules += miscellaneous_integration()

        rubi = ManyToOneReplacer(*rules)
        return rubi, rules
    _E = rubi_unevaluated_expr(E)


class LoadRubiReplacer:
    """
    Class trick to load RUBI only once.
    """

    _instance = None

    def __new__(cls):
        if matchpy is None:
            print("MatchPy library not found")
            return None
        if LoadRubiReplacer._instance is not None:
            return LoadRubiReplacer._instance
        obj = object.__new__(cls)
        obj._rubi = None
        obj._rules = None
        LoadRubiReplacer._instance = obj
        return obj

    def load(self):
        if self._rubi is not None:
            return self._rubi
        rubi, rules = get_rubi_object()
        self._rubi = rubi
        self._rules = rules
        return rubi

    def to_pickle(self, filename):
        import pickle
        rubi = self.load()
        with open(filename, "wb") as fout:
            pickle.dump(rubi, fout)

    def to_dill(self, filename):
        import dill
        rubi = self.load()
        with open(filename, "wb") as fout:
            dill.dump(rubi, fout)

    def from_pickle(self, filename):
        import pickle
        with open(filename, "rb") as fin:
            self._rubi = pickle.load(fin)
        return self._rubi

    def from_dill(self, filename):
        import dill
        with open(filename, "rb") as fin:
            self._rubi = dill.load(fin)
        return self._rubi


@doctest_depends_on(modules=('matchpy',))
def process_final_integral(expr):
    """
    Rubi's `rubi_exp` need to be replaced back to SymPy's general `exp`.

    Examples
    ========
    >>> from sympy import Function, E, Integral
    >>> from sympy.integrals.rubi.rubimain import process_final_integral
    >>> from sympy.integrals.rubi.utility_function import rubi_unevaluated_expr
    >>> from sympy.abc import a, x
    >>> _E = rubi_unevaluated_expr(E)
    >>> process_final_integral(Integral(a, x))
    Integral(a, x)
    >>> process_final_integral(_E**5)
    exp(5)

    """
    if expr.has(_E):
        expr = expr.replace(_E, E)
    return expr


@doctest_depends_on(modules=('matchpy',))
def rubi_powsimp(expr):
    """
    This function is needed to preprocess an expression as done in matchpy
    `x^a*x^b` in matchpy auotmatically transforms to `x^(a+b)`

    Examples
    ========

    >>> from sympy.integrals.rubi.rubimain import rubi_powsimp
    >>> from sympy.abc import a, b, x
    >>> rubi_powsimp(x**a*x**b)
    x**(a + b)

    """
    lst_pow = []
    lst_non_pow = []
    if isinstance(expr, Mul):
        for i in expr.args:
            if isinstance(i, (Pow, rubi_exp, sym_exp)):
                lst_pow.append(i)
            else:
                lst_non_pow.append(i)
        return powsimp(Mul(*lst_pow))*Mul(*lst_non_pow)
    return expr


@doctest_depends_on(modules=('matchpy',))
def rubi_integrate(expr, var, showsteps=False):
    """
    Rule based algorithm for integration. Integrates the expression by applying
    transformation rules to the expression.

    Returns `Integrate` if an expression cannot be integrated.

    Parameters
    ==========
    expr : integrand expression
    var : variable of integration

    Returns Integral object if unable to integrate.
    """
    rubi = LoadRubiReplacer().load()
    expr = expr.replace(sym_exp, rubi_exp)
    expr = process_trig(expr)
    expr = rubi_powsimp(expr)
    if isinstance(expr, (int, Integer)) or isinstance(expr, (float, Float)):
        return S(expr)*var
    if isinstance(expr, Add):
        results = 0
        for ex in expr.args:
            results += rubi.replace(Integral(ex, var))
        return process_final_integral(results)

    results = util_rubi_integrate(Integral(expr, var))
    return process_final_integral(results)


@doctest_depends_on(modules=('matchpy',))
def util_rubi_integrate(expr, showsteps=False, max_loop=10):
    rubi = LoadRubiReplacer().load()
    expr = process_trig(expr)
    expr = expr.replace(sym_exp, rubi_exp)
    for i in range(max_loop):
        results = expr.replace(
            lambda x: isinstance(x, Integral),
            lambda x: rubi.replace(x, max_count=10)
        )
        if expr == results:
            return results
    return results


@doctest_depends_on(modules=('matchpy',))
def get_matching_rule_definition(expr, var):
    """
    Prints the list or rules which match to `expr`.

    Parameters
    ==========
    expr : integrand expression
    var : variable of integration
    """
    rubi = LoadRubiReplacer()
    matcher = rubi.matcher
    miter = matcher.match(Integral(expr, var))
    for fun, e in miter:
        print("Rule matching: ")
        print(inspect.getsourcefile(fun))
        code, lineno = inspect.getsourcelines(fun)
        print("On line: ", lineno)
        print("\n".join(code))
        print("Pattern matching: ")
        pattno = int(re.match(r"^\s*rule(\d+)", code[0]).group(1))
        print(matcher.patterns[pattno-1])
        print(e)
        print()
