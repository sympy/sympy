"""
First Order Logic module for SymPy
"""

from __future__ import print_function
from itertools import combinations, product

from sympy.core import Symbol
from sympy.core.compatibility import ordered
from sympy.logic.boolalg import (And, Boolean, BooleanFunction,
    eliminate_implications, false, Not, Or, true)
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

    def __eq__(self, other):
        if isinstance(other, self.func):
            return self.name == other.name
        else:
            return False

    def _sympystr(self, *args, **kwargs):
        return self.name

    @classmethod
    def apply(cls):
        raise NotImplementedError()

    @property
    def name(self):
        return self._name


class Applied(FOL):

    def __init__(self, func, *args):
        if not args:
            raise ValueError("Use a constant instead of %s()" % func)
        args = [arg if isinstance(arg, (Boolean, Symbol))
                    else Constant(arg) for arg in args]
        self._func = func
        self._args = tuple(args)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.func, self.args) == (other.func, other.args)
        else:
            return False

    def _sympystr(self, *args, **kwargs):
        return "%s(%s)" % (self.name,
            ', '.join(str(arg) for arg in self.args))

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


class Constant(Boolean):
    """
    Creates a constant with the given value.
    """
    is_Constant = True

    def __new__(cls, name, **kwargs):
        if isinstance(name, cls):
            return name
        return super(Boolean, cls).__new__(cls, name, **kwargs)

    def __init__(self, name):
        self._name = name

    def __eq__(self, other):
        return isinstance(other, self.func) and self.name == other.name

    def _sympystr(self, *args, **kwargs):
        return str(self.name)

    @property
    def name(self):
        return self._name


class Quantifier(FOL):
    """
    Abstract base class for ForAll and Exists.
    """

    def __new__(cls, *args, **kwargs):

        var = args[:-1]
        expr = args[-1]
        if len(var) == 1 and hasattr(var[0], '__iter__'):
            var = set(var[0])
        else:
            var = set(var)

        if isinstance(expr, cls):
            v, e = expr.vars, expr.expr
            var = var.union(v)
            expr = e

        for x in expr.atoms(Quantifier):
            v = var.intersection(x.vars)
            if v:
                raise ValueError("Variable %s is already bound" % tuple(v))

        var = var.intersection(expr.atoms())
        if not var:
            return expr

        args = tuple(ordered(var)) + (expr, )
        obj = super(Quantifier, cls).__new__(cls, *args, **kwargs)
        return obj

    def _sympystr(self, *args, **kwargs):
        return "%s((%s), %s)" % (self.func, ', '.join(
                str(v) for v in self.vars), self.expr)

    @property
    def vars(self):
        return self.args[:-1]

    @property
    def expr(self):
        return self.args[-1]


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
    ForAll((X), Implies(Man(X), Mortal(X)))

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
    Exists((X), Implies(Man(X), Smart(X)))

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
    >>> person = lambda X: X in X_domain
    >>> time = lambda T: T in T_domain
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
    model = dict(model)
    for key, val in model.items():
        if isinstance(key, Symbol):
            if hasattr(val, '__iter__'):
                model[key] = [Constant(v) for v in val]
            else:
                model[key] = Constant(val)

        elif isinstance(key, Callable):
            if hasattr(val, '__call__'):
                continue
            mapping = {}
            for k, v in val.items():
                if k != 'default':
                    k = tuple(k) if hasattr(k, '__iter__') else (k,)
                if v is None:
                    mapping[k] = None
                else:
                    if isinstance(key, Predicate):
                        mapping[k] = true if v else false
                    else:
                        mapping[k] = Constant(v)
            model[key] = mapping

        else:
            raise ValueError()

    result = _fol_true(expr, model)
    if result is None:
        return None
    else:
        return bool(result)


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

        var = expr.vars
        domains = [model.get(v) for v in var]
        if None in domains:
            return None
        values = product(*domains)
        N = False
        for value in values:
            m = dict(zip(var, value))
            result = _fol_true(expr.expr.xreplace(m), model)
            if result is None:
                N = True
                continue
            if flag == result:
                return result
        if N:
            return None
        else:
            return not flag

    args = [_fol_true(arg, model) for arg in expr.args]

    # Functions / Predicates
    if isinstance(expr, Applied):
        args = [a.name if isinstance(a, Constant) else a for a in args]
        mapping = model.get(expr.func)
        if mapping is None:
            return None
        if hasattr(mapping, '__call__'):
            return mapping(*args)
        default = mapping.get('default')
        return mapping.get(tuple(args), default)

    # PL Operators
    return expr.func(*args)


def standardize(expr):
    """
    Rename variables so that each quantifier has its own unique variables.

    Examples
    ========

    >>> from sympy.abc import X
    >>> from sympy.logic.FOL import Predicate, ForAll, standardize
    >>> P = Predicate('P')
    >>> Q = Predicate('Q')
    >>> standardize(ForAll(X, P(X) & Q(X)) | ForAll(X, Q(X) >> P(X)))
    Or(ForAll((X), And(P(X), Q(X))), ForAll((X0), Implies(Q(X0), P(X0))))
    """
    return _standardize(expr, {})


def _standardize(expr, var_set):

    def update_var_set(vars):
        """ Adds variables to var_set and returns subsitutions to be made. """
        d = {}
        for var in vars:
            if var in var_set:
                if not var_set[var]:
                    var_set[var] = numbered_symbols(var.name)
                v = next(var_set[var])
                d[var] = v
            else:
                var_set[var] = None
                d[var] = var
        return d

    if not isinstance(expr, BooleanFunction):
        return expr

    # Prevent renaming of variables based on following equivalences
    # ForAll(X, P(X)) & ForAll(X, Q(X)) == ForAll(X, P(X) & Q(X))
    # Exists(X, P(X)) | Exists(X, Q(X)) == Exists(X, P(X) | Q(X))
    if isinstance(expr, (And, Or)):
        if isinstance(expr, And):
            cls = ForAll
        else:
            cls = Exists

        v = set()
        for arg in expr.args:
            if isinstance(arg, cls):
                v.update(arg.vars)
        d = update_var_set(v)
        expr = expr.subs(d)

        args = []
        for arg in expr.args:
            if isinstance(arg, cls):
                a = _standardize(arg.expr, var_set)
                args.append(arg.func(arg.vars, a))
            else:
                args.append(_standardize(arg, var_set))
        return expr.func(*args)

    if isinstance(expr, Quantifier):
        d = update_var_set(expr.vars)
        e = _standardize(expr.expr, var_set)
        return expr.func(d.values(), e.subs(d))

    return expr.func(*[_standardize(arg, var_set) for arg in expr.args])


def to_pnf(expr):
    """
    Converts the given FOL expression into Prenex Normal Form.

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
    expr = standardize(eliminate_implications(expr))
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
    Converts the given FOL expression into Skolem Normal Form.
    The formula in SNF is only equisatisfiable to the original
    formula and not necessarily equivalent.

    Examples
    ========

    >>> from sympy.logic.FOL import to_snf, Predicate, ForAll, Exists
    >>> from sympy.abc import X, Y
    >>> P = Predicate('P')
    >>> R = Predicate('R')
    >>> to_snf(ForAll(X, P(X) | Exists(Y, R(X, Y))))
    ForAll((X), Or(P(X), R(X, f0(X))))

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Skolem_normal_form
    """

    expr = to_pnf(expr)
    var_list = []

    skolemFunc = numbered_symbols('f', Function)
    skolemConst = numbered_symbols('c', Constant)
    while isinstance(expr, Quantifier):

        if isinstance(expr, ForAll):
            var_list.extend(expr.vars)
            expr = expr.expr

        elif isinstance(expr, Exists):
            d = {}
            for var in expr.vars:
                if var_list:
                    d[var] = next(skolemFunc)(*var_list)
                else:
                    d[var] = next(skolemConst)
            expr = expr.expr.subs(d)

        else:
            raise ValueError()

    return ForAll(var_list, expr)


def to_cnf(expr):
    """
    Converts a given FOL formula into Conjunctive Normal Form.
    The given expr is first converted to an equisatisfiable formula
    in SNF followed by dropping of implicit universal quantification
    and distribution of conjuction over disjunction.

    Examples
    ========

    >>> from sympy.abc import X, Y
    >>> from sympy.logic.FOL import Predicate, ForAll, Exists, to_cnf
    >>> P = Predicate('P')
    >>> Q = Predicate('Q')
    >>> to_cnf(ForAll(X, Exists(Y, P(X, Y) >> Q(X, Y))))
    Or(Not(P(X, f0(X))), Q(X, f0(X)))
    """
    from sympy.logic.boolalg import to_cnf as to_cnf_prop
    expr = to_snf(expr)
    while isinstance(expr, ForAll):
        expr = expr.expr
    return to_cnf_prop(expr)


def to_dnf(expr):
    """
    Converts a given FOL formula into Disjunctive Normal Form.
    The given expr is first converted to an equisatisfiable formula
    in SNF followed by dropping of implicit universal quantification
    and distribution of disjunction over conjuction.

    Examples
    ========

    >>> from sympy.abc import X, Y
    >>> from sympy.logic.FOL import Predicate, ForAll, Exists, to_dnf
    >>> P = Predicate('P')
    >>> Q = Predicate('Q')
    >>> to_dnf(ForAll(X, Exists(Y, P(X, Y) >> Q(X, Y))))
    Or(Not(P(X, f0(X))), Q(X, f0(X)))
    """
    from sympy.logic.boolalg import to_dnf as to_dnf_prop
    expr = to_snf(expr)
    while isinstance(expr, ForAll):
        expr = expr.expr
    return to_dnf_prop(expr)


def mgu(expr1, expr2):
    """
    Returns the Most General Unifier of two Predicates if it exists.

    Examples
    ========

    >>> from sympy.abc import X, Y, Z
    >>> from sympy.logic.FOL import Predicate, Function, Constant, mgu
    >>> P = Predicate('P')
    >>> f = Function('f')
    >>> g = Function('g')
    >>> a = Constant('a')
    >>> mgu(P(f(X), Z), P(Y, a))
    {Y: f(X), Z: a}
    >>> mgu(P(f(a), g(X)), P(Y, Y))
    False

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Unification_(computer_science)
    """

    if any([not isinstance(expr1, Applied), not isinstance(expr2, Applied),
            expr1.func != expr2.func, len(expr1.args) != len(expr2.args)]):
        return False

    subs = {}
    args1, args2 = list(expr1.args), list(expr2.args)

    while args1 and args2:
        arg1, arg2 = args1.pop(), args2.pop()

        if arg1 == arg2:
            continue

        sub = {}
        if arg1.is_Symbol or arg2.is_Symbol:
            if arg2.is_Symbol:
                arg1, arg2 = arg2, arg1
            if isinstance(arg2, BooleanFunction) and arg1 in arg2.atoms():
                return False
            sub[arg1] = arg2

        elif isinstance(arg1, Applied) and isinstance(arg2, Applied):
            sub = mgu(arg1, arg2)
            if not sub:
                return False

        else:
            return False

        args1 = [arg.subs(sub) for arg in args1]
        args2 = [arg.subs(sub) for arg in args2]
        for v, s in sub.items():
            subs[v] = s

    if subs == {}:
        return {true: true}

    for var, sub in subs.items():
        for v in subs:
            subs[v] = subs[v].subs({var:sub})

    return subs


def resolve(*expr):
    """
    Returns the resolution of given set of FOL formulas.

    Examples
    ========

    >>> from sympy.abc import X, Y, Z
    >>> from sympy.logic.FOL import Predicate, Function, resolve
    >>> P = Predicate('P')
    >>> Q = Predicate('Q')
    >>> f = Function('f')
    >>> a = Constant('a')
    >>> resolve((P(X) | ~Q(f(Z))), ~P(f(a)) & Q(Y))
    False

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Resolution_(logic)
    """

    expr = to_cnf(And(*expr))
    clauses = []
    for clause in expr.args:
        c = {}
        if isinstance(clause, AppliedPredicate):
            c[clause.func] = clause
        elif isinstance(clause, Not):
            c[Not(clause.args[0].func)] = clause
        elif isinstance(clause, (And, Or)):
            for literal in clause.args:
                if literal.is_Not:
                    c[Not(literal.args[0].func)] = literal
                else:
                    c[literal.func] = literal
        else:
            raise ValueError()
        clauses.append(c)

    visited = set()
    while True:
        temp = []
        for c1, c2 in combinations(clauses, 2):
            key = (tuple(c1.values()), tuple(c2.values()))

            if not key in visited:
                visited.add(key)
                t = _resolve(c1, c2)
                if {} in t:
                    return False
                temp.extend(t)

        if not temp:
            return True
        clauses.extend(temp)


def _resolve(clause1, clause2):

    clauses = []
    if len(clause1) > len(clause2):
        clause1, clause2 = clause2, clause1

    for literal in clause1:
        if Not(literal) in clause2:
            pred1 = clause1[literal]
            pred2 = clause2[Not(literal)]

            if pred1.is_Not:
                pred1 = pred1.args[0]
            elif pred2.is_Not:
                pred2 = pred2.args[0]
            subs = mgu(pred1, pred2)

            if subs:
                c = dict(list(clause1.items()) + list(clause2.items()))
                c.pop(literal)
                c.pop(Not(literal))
                for pred, literal in c.items():
                    c[pred] = literal.subs(subs)
                clauses.append(c)

    return clauses


def entails(expr, formula_set=[]):
    """
    Check whether the formula_set entails the given expr.

    Examples
    ========

    >>> from sympy.abc import X
    >>> from sympy.logic.FOL import Predicate, Constant, entails
    >>> Man = Predicate('Man')
    >>> Mortal = Predicate('Mortal')
    >>> Socrates = Constant('Socrates')
    >>> entails(Mortal(Socrates), [Man(X) >> Mortal(X), Man(Socrates)])
    True
    """

    formula_set = list(formula_set)
    formula_set.append(~expr)
    return not resolve(*formula_set)
