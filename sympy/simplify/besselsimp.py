from __future__ import print_function, division

from sympy.core.numbers import I, pi
from sympy.core import Mul, S, Dummy
from sympy.functions import exp_polar, exp, log
from sympy.functions.special.bessel import besselj, besseli, besselk, jn, bessely
from sympy.simplify.trigsimp import trigsimp, exptrigsimp
from sympy.functions.elementary.complexes import unpolarify


def besselsimp(expr):
    """
    Simplify bessel-type functions.

    This routine tries to simplify bessel-type functions. Currently it only
    works on the Bessel J and I functions, however. It works by looking at all
    such functions in turn, and eliminating factors of "I" and "-1" (actually
    their polar equivalents) in front of the argument. Then, functions of
    half-integer order are rewritten using strigonometric functions and
    functions of integer order (> 1) are rewritten using functions
    of low order.  Finally, if the expression was changed, compute
    factorization of the result with factor().

    >>> from sympy import besselj, besseli, besselsimp, polar_lift, I, S
    >>> from sympy.abc import z, nu
    >>> besselsimp(besselj(nu, z*polar_lift(-1)))
    exp(I*pi*nu)*besselj(nu, z)
    >>> besselsimp(besseli(nu, z*polar_lift(-I)))
    exp(-I*pi*nu/2)*besselj(nu, z)
    >>> besselsimp(besseli(S(-1)/2, z))
    sqrt(2)*cosh(z)/(sqrt(pi)*sqrt(z))
    >>> besselsimp(z*besseli(0, z) + z*(besseli(2, z))/2 + besseli(1, z))
    3*z*besseli(0, z)/2
    """
    # TODO
    # - better algorithm?
    # - simplify (cos(pi*b)*besselj(b,z) - besselj(-b,z))/sin(pi*b) ...
    # - use contiguity relations?

    def replacer(fro, to, factors):
        factors = set(factors)

        def repl(nu, z):
            if factors.intersection(Mul.make_args(z)):
                return to(nu, z)
            return fro(nu, z)
        return repl

    def torewrite(fro, to):
        def tofunc(nu, z):
            return fro(nu, z).rewrite(to)
        return tofunc

    def tominus(fro):
        def tofunc(nu, z):
            return exp(I*pi*nu)*fro(nu, exp_polar(-I*pi)*z)
        return tofunc

    orig_expr = expr

    ifactors = [I, exp_polar(I*pi/2), exp_polar(-I*pi/2)]
    expr = expr.replace(
        besselj, replacer(besselj,
        torewrite(besselj, besseli), ifactors))
    expr = expr.replace(
        besseli, replacer(besseli,
        torewrite(besseli, besselj), ifactors))

    minusfactors = [-1, exp_polar(I*pi)]
    expr = expr.replace(
        besselj, replacer(besselj, tominus(besselj), minusfactors))
    expr = expr.replace(
        besseli, replacer(besseli, tominus(besseli), minusfactors))

    z0 = Dummy('z')

    def expander(fro):
        def repl(nu, z):
            if (nu % 1) == S(1)/2:
                return exptrigsimp(trigsimp(unpolarify(
                        fro(nu, z0).rewrite(besselj).rewrite(jn).expand(
                            func=True)).subs(z0, z)))
            elif nu.is_Integer and nu > 1:
                return fro(nu, z).expand(func=True)
            return fro(nu, z)
        return repl

    expr = expr.replace(besselj, expander(besselj))
    expr = expr.replace(bessely, expander(bessely))
    expr = expr.replace(besseli, expander(besseli))
    expr = expr.replace(besselk, expander(besselk))

    if expr != orig_expr:
        expr = expr.factor()

    return expr
