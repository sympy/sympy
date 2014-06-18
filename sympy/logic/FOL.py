"""
First Order Logic module for SymPy
"""

from __future__ import print_function
from itertools import product

from sympy.core import Symbol
from sympy.core.compatibility import ordered
from sympy.logic.boolalg import (BooleanFunction, And, Or, Not,
    eliminate_implications)


class FOL(BooleanFunction):
    """
    Base class for all First Order Logic
    """
    is_Fol = True


class Callable(FOL):

    def __init__(self, name):
        self._name = name

    def __call__(self, *args):
        return self.apply()(self, *args)

    def _sympystr(self, *args, **kwargs):
        return self.name

    def __eq__(self, other):
        if isinstance(other, self.func):
            return self.name == other.name
        else:
            return False

    @property
    def name(self):
        return self._name


class Applied(FOL):

    def __init__(self, func, *args):
        if not args:
            raise ValueError("Use a constant instead of %s()" % func)
        self._func = func
        self._args = tuple(args)

    def _sympystr(self, *args, **kwargs):
        return "%s(%s)" % (self.name,
            ', '.join(str(arg) for arg in self.args))

    def __eq__(self, other):
        if isinstance(other, self.func.__class__):
            return (self.name, self.args) == (other.name, other.args)
        else:
            return False

    @property
    def name(self):
        return self._func.name

    @property
    def func(self):
        return self._func


class Predicate(Callable):
    """
    Creates a Predicate with the given name.
    To apply, simply call the Predicate object with arguments.

    Examples
    ========

    >>> from sympy.abc import A, B
    >>> from sympy.logic.FOL import Predicate
    >>> Knows = Predicate('Knows')
    >>> Knows('John', 'Jack')
    Knows(John, Jack)
    >>> Knows(A, B)
    Knows(A, B)
    """

    @classmethod
    def apply(cls):
        return AppliedPredicate


class AppliedPredicate(Applied):
    pass


class Function(Callable):
    """
    Creates a Function with the given name.
    To apply, simply call the Function object with arguments.

    Examples
    ========

    >>> from sympy.abc import A, B
    >>> from sympy.logic.FOL import Function
    >>> f = Function('f')
    >>> f(1, 2)
    f(1, 2)
    >>> g = Function('g')
    >>> f(A, g(A, B))
    f(A, g(A, B))
    """

    @classmethod
    def apply(cls):
        return AppliedFunction


class AppliedFunction(Applied):
    pass


class Quantifier(FOL):

    def __init__(self, var, expr):
        try:
            var = set(var)
        except TypeError:
            var = set([var])

        if isinstance(expr, self.func):
            v, e = expr.args
            var = var.union(v)
            expr = e

        for x in expr.atoms(Quantifier):
            v = var.intersection(x.args[0])
            if v:
                raise ValueError("Variable %s is already bound" % tuple(v))

        self._args = (tuple(ordered(var)), expr)

    @property
    def vars(self):
        return self.args[0]

    @property
    def expr(self):
        return self.args[1]


class ForAll(Quantifier):
    """
    Applies the Universal Quantifier on the given variable to the expr.

    Examples
    ========

    >>> from sympy.abc import X, Y
    >>> from sympy.logic.FOL import ForAll, Predicate

    >>> Man = Predicate('Man')
    >>> Mortal = Predicate('Mortal')
    >>> ForAll(X, Man(X) >> Mortal(X))
    ForAll((X,), Implies(Man(X), Mortal(X)))

    >>> Knows = Predicate('Knows')
    >>> ForAll(X, ForAll(Y, Knows(X, Y)))
    ForAll((X, Y), Knows(X, Y))
    """


class Exists(Quantifier):
    """
    Applies the Existential Quantifier on the given variable to the expr.

    Examples
    ========

    >>> from sympy.abc import X, Y
    >>> from sympy.logic.FOL import Exists, Predicate

    >>> Man = Predicate('Man')
    >>> Smart = Predicate('Smart')
    >>> Exists(X, Man(X) >> Smart(X))
    Exists((X,), Implies(Man(X), Smart(X)))

    >>> Knows = Predicate('Knows')
    >>> Exists(X, Exists(Y, Knows(X, Y)))
    Exists((X, Y), Knows(X, Y))
    """


def fol_true(expr, model={}):
    """
    Return whether given expr is satisfied by the given model.

    Parameters
    ==========

    model : dict
        Mapping of all the symbols to their corresponding values.
        Constants need not be mapped.
        Free variables should be mapped to a single value {X: 1, Y: 2}
        Bound variables should be mapped to a domain {X: [1, 2], Y: [2, 3]}
        Functions and predicates can be mapped in 2 ways:
            dict: {P: {(1, 2): True, (1, 3): True, 'default': False}}
                'default' indicates default value when key is not found.
            Callable: Any callable with the same arity as the predicate/function.


    Examples
    ========

    >>> from sympy.abc import X, T
    >>> from sympy.logic.FOL import Predicate, ForAll, Exists, fol_true
    >>> Person = Predicate('Person')
    >>> Time = Predicate('Time')
    >>> CanFool = Predicate('CanFool')
    >>> X_domain = ['John', 'Jack']
    >>> T_domain = [1, 2, 3]
    >>> def person(X): return X in X_domain
    >>> def time(T): return T in T_domain
    >>> CanFoolMap = {('John',2):False, ('John',3):False, 'default':True}
    >>> model = {X:X_domain, T:T_domain, Person:person, Time:time, CanFool:CanFoolMap}

    # You can fool some of the people all of the time
    >>> expr = Exists(X, ForAll(T, (Person(X) & Time(T)) >> CanFool(X, T)))
    >>> fol_true(expr, model)
    True

    # You can fool all of the people some of the time
    >>> expr = ForAll(X, Exists(T, (Person(X) & Time(T)) >> CanFool(X, T)))
    >>> fol_true(expr, model)
    True

    # You can fool all of the people all of the time
    >>> expr = ForAll(X, ForAll(T, (Person(X) & Time(T)) >> CanFool(X, T)))
    >>> fol_true(expr, model)
    False
    """
    return _fol_true(expr, model)


def _fol_true(expr, model={}):

    # Variables
    if isinstance(expr, Symbol):
        return model.get(expr)

    # Constants
    if not isinstance(expr, BooleanFunction):
        return expr

    # Quantifiers
    if isinstance(expr, Quantifier):
        if isinstance(expr, ForAll):
            flag = False
        elif isinstance(expr, Exists):
            flag = True
        else:
            raise ValueError()

        var = expr.args[0]
        domain = [model.get(v) for v in var]
        if None in domain:
            return None
        values = product(*domain)
        N = False
        for value in values:
            m = dict(zip(var, value))
            result = _fol_true(expr.args[1].xreplace(m), model)
            if result is None:
                N = True
                continue
            if flag == result:
                return result
        if N:
            return None
        else:
            return not flag

    args = tuple([_fol_true(arg, model) for arg in expr.args])

    # Functions / Predicates
    if isinstance(expr, Applied):
        mapping = model.get(expr.func)
        if mapping is None:
            return None
        if hasattr(mapping, '__call__'):
            try:
                return mapping(*args)
            except:
                return None
        default = mapping.get('default')
        return mapping.get(args, default)

    # PL Operators
    return expr.func(*args)


def to_pnf(expr):
    """
    Converts the given FOL expression into Prenex Normal Form

    Examples
    ========

    >>> from sympy.abc import a, X, Y
    >>> from sympy.logic.FOL import Predicate, to_pnf, ForAll, Exists
    >>> F = Predicate('F')
    >>> G = Predicate('G')
    >>> H = Predicate('H')
    >>> to_pnf((F(a) | Exists(X, G(X))) >> ForAll(Y, H(Y)))
    ForAll((X, Y), Or(And(Not(F(a)), Not(G(X))), H(Y)))


    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Prenex_normal_form
    """
    expr = eliminate_implications(expr)
    return _to_pnf(expr)


def _to_pnf(expr):

    if not isinstance(expr, BooleanFunction):
        return expr

    if isinstance(expr, Applied):
        return expr

    if isinstance(expr, Not):
        args = expr.args[0]
        if isinstance(args, Quantifier):
            if isinstance(args, ForAll):
                func = Exists
            elif isinstance(args, Exists):
                func = ForAll
            else:
                raise ValueError()
            return func(args.vars, _to_pnf(~args.expr))
        return expr

    if isinstance(expr, (And, Or)):
        prefix = []
        matrix = []
        args = [_to_pnf(arg) for arg in expr.args]

        for arg in args:
            while isinstance(arg, Quantifier):
                prefix.append((arg.func, arg.vars))
                arg = arg.expr
            matrix.append(arg)

        expr = expr.func(*matrix)
        while prefix:
            func, var = prefix.pop()
            expr = func(var, expr)
        return expr

    if isinstance(expr, Quantifier):
        return expr.func(expr.vars, _to_pnf(expr.expr))

    raise ValueError()
