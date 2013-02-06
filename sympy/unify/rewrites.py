from sympy.unify.usympy import construct, deconstruct
from sympy.unify.core import reify, Variable, Compound
from sympy.unify import core as ucore
from sympy.rules.tools import subs
from sympy.rules.branch import exhaust, multiplex
import functools
import itertools as it

from sympy import Add, MatAdd, MatMul, Union, Intersection, FiniteSet
from sympy.core.operations import AssocOp

comm_ops = (Add, MatAdd, Union, Intersection, FiniteSet)
def is_commutative(x):
    return isinstance(x, Compound) and x.op in comm_ops

assoc_ops = (AssocOp, MatAdd, MatMul, Union, Intersection, FiniteSet)
def is_associative(x):
    return isinstance(x, Compound) and any(issubclass(x.op, aop) for aop in assoc_ops)

unify = functools.partial(ucore.unify, is_associative=is_associative,
                                       is_commutative=is_commutative)

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
