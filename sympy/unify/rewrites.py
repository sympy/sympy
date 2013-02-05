from sympy.unify.usympy import construct, deconstruct
from sympy.unify.core import unify, reify, Variable
from sympy.rules.tools import subs
from sympy.rules.branch import exhaust, multiplex
import functools
import itertools as it

def crl(source, target, variables, condition=None, unify=unify, reify=reify):
    """ Conditional rewrite rule """
    def f(expr):
        for match in unify(source, expr, {}):
            if not condition or condition(*map(match.get, variables, variables)):
                yield reify(target, match)
    return f

def rewriterules(sources, targets, variabless, conditions,
        construct=construct,
        deconstruct=deconstruct,
        unify=unify,
        reify=reify,
        strategy=lambda rules: exhaust(multiplex(*rules))):

    sources2 = map(deconstruct, sources, variabless)
    targets2 = map(deconstruct, targets, variabless)
    variabless2 = map(lambda l: map(Variable, l), variabless)
    conditions2 = map(functools.partial(fn_deconstruct, construct=construct),
                      conditions)

    crl2 = functools.partial(crl, reify=reify)

    rules = map(crl2, sources2, targets2, variabless2, conditions2)

    return lambda expr: it.imap(construct, strategy(rules)(deconstruct(expr)))

def types(expr):
    return set([type(expr)]).union(
            reduce(set.union, map(types, expr.args), set()))

def fn_deconstruct(fn, construct=construct):
    if not fn:
        return fn
    return lambda *args: fn(*map(construct, args))
