"""
First Order Logic module for SymPy
"""

from __future__ import print_function
from itertools import product

from sympy.core import Symbol
from sympy.core.compatibility import ordered
from sympy.logic.boolalg import (And, BooleanFunction,
    eliminate_implications, Not, Or)
from sympy.utilities.iterables import numbered_symbols


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
        if isinstance(other, self.__class__):
            return (self.func, self.args) == (other.func, other.args)
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

    def __new__(cls, var, expr, **kwargs):
        
        if hasattr(var, '__iter__'):
            var = set(var)
        else:
            var = set([var])

        if isinstance(expr, cls):
            v, e = expr.args
            var = var.union(v)
            expr = e

        for x in expr.atoms(Quantifier):
            v = var.intersection(x.args[0])
            if v:
                raise ValueError("Variable %s is already bound" % tuple(v))

        var = var.intersection(expr.atoms())
        if not var:
            return expr

        obj = super(Quantifier, cls).__new__(cls, var, expr, **kwargs)
        obj._args = (tuple(ordered(var)), expr)
        return obj

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

def standardize(expr):
    """
    Rename variables so that each quantifier has its own unique variables

    Examples
    ========

    >>> from sympy.abc import X, Y
    >>> from sympy.logic.FOL import Predicate, ForAll, standardize
    >>> P = Predicate('P')
    >>> Q = Predicate('Q')
    >>> standardize(ForAll(X, P(X) >> Q(X)) | ForAll(X, Q(X) >> P(X)))
    """
    return _standardize(expr, {})


def _standardize(expr, var_set):

    if not isinstance(expr, BooleanFunction):
        return expr

    if isinstance(expr, Quantifier):
        d = {}
        for var in expr.vars:
            if var in var_set:
                if not var_set[var]:
                    var_set[var] = numbered_symbols(var.name)
                v = next(var_set[var])
                d[var] = v
            else:
                var_set[var] = None
                d[var] = var
        e = _standardize(expr.expr, var_set)
        return expr.func(d.values(), e.subs(d))

    return expr.func(*[_standardize(arg, var_set) for arg in expr.args])


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


def to_snf(expr):
    """
    Converts the given FOL expression into Skolem Normal Form

    Examples
    ========

    >>> from sympy.logic.FOL import to_snf, Predicate, ForAll, Exists
    >>> from sympy.abc import X, Y
    >>> P = Predicate('P')
    >>> R = Predicate('R')
    >>> to_snf(ForAll(X, P(X) | Exists(Y, R(X, Y))))
    ForAll((X,), Or(P(X), R(X, f0(X))))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Skolem_normal_form
    """

    expr = to_pnf(expr)
    var_list = []

    while isinstance(expr, Quantifier):
        skolemFunc = numbered_symbols('f', Function)

        if isinstance(expr, ForAll):
            var_list.extend(expr.vars)
            expr = expr.expr

        elif isinstance(expr, Exists):
            if var_list:
                d = {}
                for var in expr.vars:
                    d[var] = next(skolemFunc)(*var_list)
                expr = expr.expr.subs(d)
            else:
                expr = expr.expr

        else:
            raise ValueError()

    return ForAll(var_list, expr)
