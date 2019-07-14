#!/usr/bin/env python

"""
A test suite that uses Hypothesis to generate and compare expressions.

See https://hypothesis.readthedocs.io/ and "stateful testing" for how
this works, and feel free to add new methods to generate or manipulate
expressions.
"""

import operator
import os
import string
import time

import hypothesis.strategies as st
from hypothesis import Verbosity, assume, reject, settings
from hypothesis.stateful import (
    Bundle,
    RuleBasedStateMachine,
    precondition,
    rule,
    run_state_machine_as_test,
)

import sympy
from sympy.core.cache import clear_cache

if os.name == "nt":
    from stopit import ThreadingTimeout as Timeout
else:
    from stopit import SignalTimeout as Timeout


@settings(deadline=None, verbosity=Verbosity.normal, max_examples=1000)
class SympyRules(RuleBasedStateMachine):
    def __init__(self, *args, **kwargs):
        RuleBasedStateMachine.__init__(self, *args, **kwargs)
        clear_cache()

    def execute_step(self, step):
        with Timeout(2, swallow_exc=False) as timeout_ctx:
            val = super(SympyRules, self).execute_step(step)
        assert timeout_ctx.state == timeout_ctx.EXECUTED
        return val

    # Bundles are like groups of variables.  We assign into them by returning
    # a value from a rule, and can later use them as rule arguments.
    c = Bundle("c")
    v = Bundle("v")
    expr = Bundle("expr")
    term = Bundle("term")
    poly = Bundle("poly")

    # For example, this asks Hypothesis for any integer 0 <= x <= 10, and adds
    # it to the `c` Bundle.  Any `c` can be an `expr` via the `c_is_expr` rule.
    @rule(target=c, x=st.integers(0, 10))
    def gen_c(self, x):
        return sympy.Integer(x)

    @rule(target=v, x=st.sampled_from(string.ascii_lowercase))
    def gen_v(self, x):
        return sympy.Symbol(x, positive="e" <= x <= "l")

    @rule(target=expr, m=c, n=c)
    def gen_rat(self, m, n):
        return sympy.Rational(m, n)

    @rule(target=expr, x=c)
    def c_is_expr(self, x):
        return x

    @rule(target=expr, x=v)
    def v_is_expr(self, x):
        return x

    @rule(target=term, x=c, y=v, z=c)
    def gen_term(self, x, y, z):
        return x * (y ** z)

    @rule(target=poly, x=term)
    def term_is_poly(self, x):
        return x

    @rule(target=poly, x=poly, y=poly)
    def add_poly(self, x, y):
        return x + y

    @rule(target=expr, x=poly)
    def poly_is_expr(self, x):
        return x

    # This precondition skips all `%`-expressions, as we trigger a known bug.
    @precondition(lambda _: False)  # TODO: fix #14563
    @rule(target=expr, x=expr, y=expr)
    def mod_expr(self, x, y):
        try:
            return x % y
        except ZeroDivisionError:
            reject()

    @rule(target=expr, x=st.sampled_from([sympy.pi, sympy.E, sympy.I]))
    def constants_are_expr(self, x):
        return x

    @rule(
        target=expr,
        op=st.sampled_from(
            [
                operator.add,
                operator.sub,
                operator.mul,
                operator.truediv,
                # Turns out that this can be pathologically slow in some cases.
                # TODO: work out why, and fix it.  (on Unix so signal timeouts work)
                # operator.pow
            ]
        ),
        x=expr,
        y=expr,
    )
    def add_expr_ops(self, op, x, y):
        return op(x, y)

    @rule(
        target=expr,
        fn=st.sampled_from(
            [
                sympy.sin,
                sympy.cos,
                sympy.tan,
                sympy.asin,
                sympy.acos,
                sympy.atan,
                sympy.sqrt,
                sympy.powsimp,
                sympy.expand_power_exp,
                sympy.expand_power_base,
                sympy.powdenest,
                sympy.expand_log,
                sympy.logcombine,
                sympy.expand_func,
                sympy.combsimp,
                sympy.trigsimp,
                sympy.expand_trig,
                sympy.factorial,
            ]
        ),
        x=expr,
    )
    def fn_expr(self, fn, x):
        return fn(x)

    @rule(target=expr, x=expr, y=expr, z=expr)
    def subs_expr(self, x, y, z):
        return x.subs(y, z)

    @rule(target=expr, x=expr)
    def simplify_expr(self, x):
        return sympy.simplify(x)

    @rule(target=expr, x=expr)
    def cancel_expr(self, x):
        return sympy.cancel(x)

    @rule(target=expr, x=expr, force=st.booleans())
    def expand_expr(self, x, force):
        return sympy.expand(x, force=force)

    @rule(target=expr, x=expr)
    def factor_expr(self, x):
        try:
            return sympy.factor(x)
        except ValueError:
            reject()

    @rule(target=expr, x=expr)
    def apart_expr(self, x):
        try:
            return sympy.apart(x)
        except (sympy.PolynomialError, sympy.NotAlgebraic, NotImplementedError):
            reject()

    @rule(target=expr, x=expr, y=expr)
    def collect(self, x, y):
        return sympy.collect(x, y)

    @rule(x=expr)
    def expr_doit(self, x):
        x.doit()

    @rule(x=expr)
    def expr_evalf(self, x):
        x.evalf()

    @rule(
        target=expr,
        x=expr,
        var=v,
        f=c,
        t=c,
        op=st.sampled_from([sympy.Sum, sympy.Product]),
    )
    def combine_many(self, x, var, f, t, op):
        return op(x, (var, f, t))

    @rule(source=expr)
    def doit_is_idempotent(self, source):
        done = source.doit()
        if done == done:
            assert done == done.doit()

    @rule(source=expr, data=st.data())
    def substitution_commutes_with_simplify(self, source, data):
        symbols = sorted(source.free_symbols, key=lambda a: a.name)
        assume(symbols and sympy.simplify(source).free_symbols)

        # TODO: this is quite noisy, but does ensure that we can tell
        # which expression is taking a very long time if it gets stuck.
        print("    checking: `{!r}`".format(source))

        variables = data.draw(st.permutations(symbols), label="variables")
        assignment = [sympy.Integer(data.draw(st.integers())) for _ in variables]

        r1 = source.subs(dict(zip(variables, assignment)))
        assume(not isinstance(r1, sympy.numbers.NaN))

        r2 = source
        for a, v in zip(variables, assignment):
            r2 = sympy.simplify(r2.subs({a: v}))
        assume(not isinstance(r2, sympy.numbers.NaN))

        assert r1.doit().equals(r2.doit())


# This line ensures that the unittest or pytest runners can execute the tests.
TestSympy = SympyRules.TestCase

if __name__ == "__main__":
    start = time.time()
    print("Running Hypothesis tests...")
    try:
        run_state_machine_as_test(SympyRules)
    except Exception:
        print(
            "\n\nHypothesis tests found some problems.\n"
            "These errors may or may not be related to your pull request.\n"
            "If not, please open an issue with the reported examples.\n"
        )
    finally:
        print("Finished Hypothesis tests in {:.2f} seconds".format(time.time() - start))
