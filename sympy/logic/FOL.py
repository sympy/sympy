"""
First Order Logic module for SymPy
"""

from __future__ import print_function
from collections import defaultdict, deque
from itertools import chain, combinations, product

from sympy.core import Symbol
from sympy.core.compatibility import iterable, ordered
from sympy.logic.boolalg import (And, Boolean, BooleanFunction,
    conjuncts, disjuncts, eliminate_implications, false, Implies,
    Not, Or, to_cnf, true)
from sympy.utilities.iterables import numbered_symbols


class FOL(BooleanFunction):
    """
    Abstract base class for all First Order Logic classes.
    Only attributes and dispatcher methods that need to be inherited by
    every other classes in the module should go here.
    """

    def to_nnf(self, simplify=True):
        return self


class Callable(FOL):
    """
    Abstract base class for 'Predicate' and 'Function'.
    This class provides the functionality for the 'Predicate' and
    'Function' objects to be called to yield its 'Applied' version.
    The classes extending 'Callable' simply need to override 'apply'
    method to return the appropriate 'Applied' class. This class is
    then called with the arguments supplied to the call to return an
    object of type 'AppliedPredicate' or 'AppliedFunction'.
    """

    def __init__(self, name):
        self._name = name

    def __call__(self, *args):
        """
        Uses internal dispatching to return the appropriate Applied
        object with the given arguments.
        """
        return self.apply()(self, *args)

    def __eq__(self, other):
        return isinstance(other, self.func) and self.name == other.name

    def __hash__(self):
        return super(Callable, self).__hash__()

    def _hashable_content(self):
        return (self.func, self.name)

    def _sympystr(self, *args, **kwargs):
        return self.name

    @classmethod
    def apply(cls):
        """
        Returns the 'Applied' version of the class.
        This method is intended to be overridden by the subclass
        returning the corresponding applied class.
        """
        raise NotImplementedError()

    @property
    def name(self):
        return self._name


class Applied(FOL):
    """
    Abstract base class for 'AppliedPredicate' and 'AppliedFunction'.
    This class provides common functionality for all subclasses to
    sanitize given arguments such that any non-Boolean argument is
    converted to a Constant. It also provides methods to return the
    original Callable object which was called to obtain this object.
    """

    def __init__(self, func, *args):
        if not args:
            raise ValueError("Use a constant instead of %s()" % func)
        args = [arg if isinstance(arg, (Boolean, Symbol))
                    else Constant(arg) for arg in args]
        self._args = tuple(args)
        self._func = func

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.func, self.args) == (other.func, other.args)
        else:
            return False

    def __hash__(self):
        return super(Applied, self).__hash__()

    def _hashable_content(self):
        return (self.__class__, self.name) + self.args

    def _sympystr(self, *args, **kwargs):
        return "%s(%s)" % (self.name,
            ', '.join(str(arg) for arg in self.args))

    @property
    def name(self):
        """
        Returns the name of the original 'Predicate' or 'Function'
        """
        return self.func.name

    @property
    def func(self):
        """
        Returns the class from which the given class was applied.
        This functionality is different from the usual SymPy convention
        of returning the __class__ of the object.
        """
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
    """
    Applied version of Predicate.
    All AppliedPredicate objects are intended to be created by
    calling the corresponding 'Predicate' object with arguments.
    """


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
    """
    Applied version of Function.
    All AppliedFunction objects are intended to be created by
    calling the corresponding 'Function' object with arguments.
    """


class Constant(Boolean):
    """
    Creates a constant with the given value.
    All non-Boolean objects in the FOL universe are Constants and are
    implicitly converted when 'Applied' or used for interpretation.
    Boolean objects include all 'BooleanFunctions', true/ false constants
    and Symbols (which also extends 'Boolean').

    Examples
    ========

    >>> from sympy.logic.FOL import Constant
    >>> Cons = Constant('Cons')
    >>> Cons
    Cons
    >>> isinstance(Cons, Constant)
    True

    Notes
    =====
    It is possible to make do without a separate class for Constants
    simply by using Symbols in its place. However during unification,
    which is critical to the inference system, it is important to be
    able to differentiate between Symbols and Constants as Symbols can be
    unified with some other object but the same is not true for Constants.
    In future if some technique can be used to differentiate between these
    without using the Constants class or using some pre-existing SymPy
    construct then this class can be safely removed.
    """

    def __new__(cls, name, **kwargs):
        return super(Boolean, cls).__new__(cls, str(name), **kwargs)

    def __init__(self, name):
        if isinstance(name, self.func):
            self._name = name.name
        else:
            self._name = name

    def __eq__(self, other):
        return isinstance(other, self.func) and self.name == other.name

    def __hash__(self):
        return super(Constant, self).__hash__()

    def _hashable_content(self):
        return (self.func, str(self.name))

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
        if len(var) == 1 and iterable(var[0]):
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
        """
        Returns the list of bound variables.
        """
        return self.args[:-1]

    @property
    def expr(self):
        """
        Returns the quantified expression.
        """
        return self.args[-1]

    def to_nnf(self, simplify=True):
        from sympy.logic.boolalg import to_nnf
        return self.func(*[to_nnf(arg, simplify=simplify)
                                        for arg in self.args])


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
    >>> X_ = ['John', 'Jack']
    >>> T_ = [1, 2, 3]
    >>> Person_ = lambda X: X in X_
    >>> Time_ = lambda T: T in T_
    >>> CanFool_ = {('John', 2): False, ('John', 3): False, 'default': True}
    >>> model = {X: X_, T: T_, Person: Person_, Time: Time_, CanFool: CanFool_}

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
            if iterable(val):
                model[key] = [Constant(v) for v in val]
            else:
                model[key] = Constant(val)

        elif isinstance(key, Callable):
            if hasattr(val, '__call__'):
                continue
            mapping = {}
            for k, v in val.items():
                if k != 'default':
                    k = tuple(k) if iterable(k) else (k,)
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
        none = False
        for value in values:
            m = dict(zip(var, value))
            result = _fol_true(expr.expr.xreplace(m), model)
            if result is None:
                none = True
                continue
            if flag == result:
                return result
        if none:
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


def standardize(expr, variables=None):
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
    return _standardize(expr, {}, variables)


def _standardize(expr, var_set, variables=None):

    def update_var_set(vars):
        """ Adds variables to var_set and returns subsitutions to be made. """
        d = {}
        for var in vars:
            if var in var_set:
                if variables is None:
                    if not var_set[var]:
                        var_set[var] = numbered_symbols(var.name)
                    v = next(var_set[var])
                    d[var] = v
                else:
                    d[var] = next(variables)
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
                a = _standardize(arg.expr, var_set, variables)
                args.append(arg.func(arg.vars, a))
            else:
                args.append(_standardize(arg, var_set, variables))
        return expr.func(*args)

    if isinstance(expr, Quantifier):
        d = update_var_set(expr.vars)
        e = _standardize(expr.expr, var_set, variables)
        return expr.func(d.values(), e.subs(d))

    return expr.func(*[_standardize(arg, var_set, variables)
                                        for arg in expr.args])


def to_pnf(expr, variables=None):
    """
    Converts the given FOL expression into Prenex Normal Form.
    A FOL formula is in Prenex Normal Form if it can be expressed as
    a collection of quantifiers (prefix) followed by a quantifier-free
    expr (matrix). The expr in PNF is equivalent to the given formula.


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
    expr = standardize(eliminate_implications(expr), variables)
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


def to_snf(expr, functions=None, variables=None, constants=None):
    """
    Converts the given FOL expression into Skolem Normal Form.
    A FOL formula is in Skolem Normal Form if it is in PNF with no
    existential quantifier and all existentially quantified variables
    replaced by Skolem functions/ constants.
    The returned formula has universal quantifiers dropped and all
    variables are assumed to be universally quantified.

    The formula in SNF is only equisatisfiable to the original
    formula (satisfiable if and only if original formula is
    satisfiable) and not necessarily equivalent (same truth table).


    Parameters
    ==========
    expr :       The formula to be converted to SNF.
    functions :  Generator/ Iterator for Skolem Functions.
    variables :  Generator/ Iterator for new variables for standardization.
    Constants :  Generator/ Iterator for Skolem Constants.


    Examples
    ========

    >>> from sympy.logic.FOL import to_snf, Predicate, ForAll, Exists
    >>> from sympy.abc import X, Y
    >>> P = Predicate('P')
    >>> R = Predicate('R')
    >>> to_snf(ForAll(X, P(X) | Exists(Y, R(X, Y))))
    Or(P(X), R(X, f0(X)))


    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Skolem_normal_form
    """

    expr = to_pnf(expr, variables)
    var_list = []

    if functions is None:
        skolemFunc = numbered_symbols('f', Function)
    else:
        skolemFunc = functions
    if constants is None:
        skolemConst = numbered_symbols('c', Constant)
    else:
        skolemConst = constants

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

    return expr


def mgu(expr1, expr2):
    """
    Returns the Most General Unifier of two Predicates if it exists.
    This function is critical for the entire inference system as it
    determines if two clauses can be resolved together, and the value
    of the resolved clause.

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
    Returns the resolution of set of given FOL formulas.
    Resolution in FOL is a refutation-complete inference system
    that gives the (un)satisfiability of a set of clauses.

    Examples
    ========

    >>> from sympy.abc import X, Y, Z
    >>> from sympy.logic.FOL import Predicate, Function, Constant, resolve
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

    expr = to_cnf(to_snf(And(*expr)))
    clauses = []
    for clause in conjuncts(expr):
        c = defaultdict(list)
        for literal in disjuncts(clause):
            func = literal.func
            if isinstance(literal, Not):
                literal = literal.args[0]
                func = ~literal.func
            if isinstance(literal, AppliedPredicate):
                c[func].append(literal)
            else:
                raise ValueError()
        clauses.append(c)

    visited = set()
    new_clauses = []
    generator = combinations(clauses, 2)
    while True:
        temp = []
        for c1, c2 in generator:
            key = frozenset(frozenset(chain.from_iterable(c.values())) for c in (c1, c2))
            if key not in visited:
                visited.add(key)
                t = _resolve(c1, c2)
                if {} in t:
                    return False
                if t:
                    temp.extend(t)

        if not temp:
            return True
        clauses.extend(new_clauses)
        new_clauses = temp
        generator = product(clauses, new_clauses)


def _resolve(clause1, clause2):
    clauses = []
    if len(clause1) > len(clause2):
        clause1, clause2 = clause2, clause1

    for literal in clause1:
        if ~literal in clause2:
            for P1, P2 in product(clause1[literal], clause2[~literal]):
                subs = mgu(P1, P2)
                if subs:
                    c = defaultdict(list)
                    for P, l in chain(clause1.items(), clause2.items()):
                        c[P].extend([lit.subs(subs) for lit in l])
                    P = P1.subs(subs)
                    for lit in (literal, ~literal):
                        c[lit].remove(P)
                        if not c[lit]:
                            c.pop(lit)
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


class FOL_KB():
    """
    First Order Logic Knowledge Base.
    This KB allows addition of pure Horn clauses only and uses
    backward chaining to provide results.
    """

    def __init__(self):
        self.vars = numbered_symbols()
        self.clauses = defaultdict(deque)
        self.visited = None
        self.max_limit = 8
        self.limit_increment_func = lambda x: x * 2

    def tell(self, clause):
        """
        Add a clause to the Knowledge Base.
        All facts must be of the form: P,
        All rules must be of the form: (P1 & P2 & ... & Pn) >> Q,
        where P and Q are non-negative (Applied)Predicates.
        """

        def _validate(pred):
            """
            Helper method to check for valid Predicates.
            """
            if isinstance(pred, Not):
                raise ValueError("Negative Predicate found: %s" % pred)
            if not isinstance(pred, AppliedPredicate):
                raise ValueError("Invalid Predicate: %s" % pred)

        variables = clause.atoms()
        subs = dict((v, next(self.vars)) for v in variables)
        clause = clause.subs(subs)
        if isinstance(clause, Implies):
            ante, cons = clause.args
            _validate(cons)
            ante = tuple(conjuncts(ante))
            for pred in ante:
                _validate(pred)
            self.clauses[cons.func].append((cons, ante))
        else:
            _validate(clause)
            self.clauses[clause.func].appendleft((clause, None))

    def ask(self, query, all_answers=False):
        """
        Ask a query.
        If query contains variables then returns a possible answer.
        If all_answers is True then returns all possible answers.
        Setting all_answers does nothing if query contains no variables.

        Examples
        ========

        >>> from sympy.abc import X
        >>> from sympy.logic.FOL import Predicate, Constant, FOL_KB
        >>> KB = FOL_KB()
        >>> Croaks = Predicate('Croaks')
        >>> EatsFlies = Predicate('EatsFlies')
        >>> Chirps = Predicate('Chirps')
        >>> Sings = Predicate('Sings')
        >>> Frog = Predicate('Frog')
        >>> Green = Predicate('Green')
        >>> Canary = Predicate('Canary')
        >>> Yellow = Predicate('Yellow')
        >>> Fritz = Constant('Fritz')
        >>> Tweety = Constant('Tweety')
        >>> KB.tell((Croaks(X) & EatsFlies(X)) >> Frog(X))
        >>> KB.tell((Chirps(X) & Sings(X)) >> Canary(X))
        >>> KB.tell(Frog(X) >> Green(X))
        >>> KB.tell(Canary(X) >> Yellow(X))
        >>> KB.tell(Croaks(Fritz))
        >>> KB.tell(EatsFlies(Fritz))
        >>> KB.tell(Yellow(Tweety))
        >>> KB.tell(EatsFlies(Tweety))
        >>> KB.ask(Green(Fritz))
        True
        >>> KB.ask(Frog(Tweety))
        False
        >>> KB.ask(EatsFlies(X), all_answers=True)
        [EatsFlies(Fritz), EatsFlies(Tweety)]
        >>> KB.ask(Frog(X))
        [Frog(Fritz)]
        """
        limit = 1
        self.query = query
        self.query_vars = query.atoms()
        while limit <= self.max_limit:
            self.models = set()
            result = self._ask(list(conjuncts(query)), limit,
                                self.query_vars and all_answers)
            limit = self.limit_increment_func(limit)
            if result:
                break
        if all_answers or self.query_vars:
            return list(ordered(self.models))
        return result

    def _ask(self, query, limit, all_answers=False, level=0, unifiers=[]):
        if not query:
            self.models.add(self.query.subs(dict(unifiers)))
            return True
        if level > limit:
            return False
        literal = query.pop()
        if literal.func in self.clauses:
            for clause in self.clauses[literal.func]:
                goal = query
                cons, ante = clause
                unifier = mgu(literal, cons)
                if unifier:
                    t = unifiers + list((key, value) for key, value in
                            unifier.items() if key in self.query_vars)
                    goal = [l.subs(unifier) for l in goal]
                    if ante:
                        u = list(unifier.items()) + [(a, next(self.vars))
                        for a in chain.from_iterable((l.args for l in ante))
                            if isinstance(a, Symbol)]
                        goal.extend(l.subs(u) for l in ante)
                    if self._ask(goal, limit, all_answers, level+1, t) \
                            and not all_answers:
                        return True
        return False
