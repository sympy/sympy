"""Boolean algebra module for SymPy"""
from collections import defaultdict
from itertools import product

from sympy.core.basic import Basic
from sympy.core.numbers import Number
from sympy.core.decorators import deprecated
from sympy.core.operations import LatticeOp
from sympy.core.function import Application, sympify
from sympy.core.compatibility import ordered


class Boolean(Basic):
    """A boolean object is an object for which logic operations make sense."""

    __slots__ = []

    def __and__(self, other):
        """Overloading for & operator"""
        return And(self, other)

    def __or__(self, other):
        """Overloading for |"""
        return Or(self, other)

    def __invert__(self):
        """Overloading for ~"""
        return Not(self)

    def __rshift__(self, other):
        """Overloading for >>"""
        return Implies(self, other)

    def __lshift__(self, other):
        """Overloading for <<"""
        return Implies(other, self)

    def __xor__(self, other):
        return Xor(self, other)


class BooleanFunction(Application, Boolean):
    """Boolean function is a function that lives in a boolean space
    It is used as base class for And, Or, Not, etc.
    """
    is_Boolean = True

    def __call__(self, *args):
        return self.func(*[arg(*args) for arg in self.args])

    def _eval_simplify(self, ratio, measure):
        return simplify_logic(self)


class And(LatticeOp, BooleanFunction):
    """
    Logical AND function.

    It evaluates its arguments in order, giving False immediately
    if any of them are False, and True if they are all True.

    Examples
    ========

    >>> from sympy.core import symbols
    >>> from sympy.abc import x, y
    >>> from sympy.logic.boolalg import And
    >>> x & y
    And(x, y)

    Notes
    =====

    The operator operator ``&`` will perform bitwise operations
    on integers and for convenience will construct an Add when
    the arguments are symbolic, but the And function does not
    perform bitwise operations and when any argument is True it
    is simply removed from the arguments:

    >>> And(x, y).subs(x, 1)
    y
    """
    zero = False
    identity = True

    @classmethod
    def _new_args_filter(cls, args):
        newargs = []
        for x in args:
            if isinstance(x, Number) or x in (0, 1):
                newargs.append(True if x else False)
            else:
                newargs.append(x)
        return LatticeOp._new_args_filter(newargs, And)


class Or(LatticeOp, BooleanFunction):
    """
    Logical OR function

    It evaluates its arguments in order, giving True immediately
    if any of them are True, and False if they are all False.

    Examples
    ========

    >>> from sympy.core import symbols
    >>> from sympy.abc import x, y
    >>> from sympy.logic.boolalg import Or
    >>> x | y
    Or(x, y)

    Notes
    =====

    The operator operator ``|`` will perform bitwise operations
    on integers and for convenience will construct an Or when
    the arguments are symbolic, but the Or function does not
    perform bitwise operations and when any argument is False it
    is simply removed from the arguments:

    >>> Or(x, y).subs(x, 0)
    y
    """
    zero = True
    identity = False

    @classmethod
    def _new_args_filter(cls, args):
        newargs = []
        for x in args:
            if isinstance(x, Number) or x in (0, 1):
                newargs.append(True if x else False)
            else:
                newargs.append(x)
        return LatticeOp._new_args_filter(newargs, Or)


class Not(BooleanFunction):
    """
    Logical Not function (negation)

    Notes
    =====

    De Morgan rules are applied automatically.
    """

    is_Not = True

    @classmethod
    def eval(cls, arg):
        """
        Logical Not function (negation)

        Returns True if the statement is False
        Returns False if the statement is True

        Examples
        ========

        >>> from sympy.logic.boolalg import Not, And, Or
        >>> from sympy.abc import x
        >>> Not(True)
        False
        >>> Not(False)
        True
        >>> Not(And(True, False))
        True
        >>> Not(Or(True, False))
        False
        >>> Not(And(And(True, x), Or(x, False)))
        Not(x)
        """
        if isinstance(arg, Number) or arg in (0, 1):
            return False if arg else True
        # apply De Morgan Rules
        if arg.func is And:
            return Or(*[Not(a) for a in arg.args])
        if arg.func is Or:
            return And(*[Not(a) for a in arg.args])
        if arg.func is Not:
            return arg.args[0]


class Xor(BooleanFunction):
    """
    Logical XOR (exclusive OR) function.
    """
    @classmethod
    def eval(cls, *args):
        """
        Logical XOR (exclusive OR) function.

        Returns True if an odd number of the arguments are True
            and the rest are False.
        Returns False if an even number of the arguments are True
            and the rest are False.

        Examples
        ========

        >>> from sympy.logic.boolalg import Xor
        >>> Xor(True, False)
        True
        >>> Xor(True, True)
        False

        >>> Xor(True, False, True, True, False)
        True
        >>> Xor(True, False, True, False)
        False
        """
        if not args:
            return False
        args = list(args)
        A = args.pop()
        while args:
            B = args.pop()
            A = Or(And(A, Not(B)), And(Not(A), B))
        return A


class Nand(BooleanFunction):
    """
    Logical NAND function.

    It evaluates its arguments in order, giving True immediately if any
    of them are False, and False if they are all True.
    """
    @classmethod
    def eval(cls, *args):
        """
        Logical NAND function.

        Returns True if any of the arguments are False
        Returns False if all arguments are True

        Examples
        ========

        >>> from sympy.logic.boolalg import Nand
        >>> Nand(False, True)
        True
        >>> Nand(True, True)
        False
        """
        return Not(And(*args))


class Nor(BooleanFunction):
    """
    Logical NOR function.

    It evaluates its arguments in order, giving False immediately if any
    of them are True, and True if they are all False.
    """
    @classmethod
    def eval(cls, *args):
        """
        Logical NOR function.

        Returns False if any argument is True
        Returns True if all arguments are False

        Examples
        ========

        >>> from sympy.logic.boolalg import Nor
        >>> Nor(True, False)
        False
        >>> Nor(True, True)
        False
        >>> Nor(False, True)
        False
        >>> Nor(False, False)
        True
        """
        return Not(Or(*args))


class Implies(BooleanFunction):
    """
    Logical implication.

    A implies B is equivalent to !A v B
    """
    @classmethod
    def eval(cls, *args):
        """
        Logical implication.

        Accepts two Boolean arguments; A and B.
        Returns False if A is True and B is False
        Returns True otherwise.

        Examples
        ========

        >>> from sympy.logic.boolalg import Implies
        >>> Implies(True, False)
        False
        >>> Implies(False, False)
        True
        >>> Implies(True, True)
        True
        >>> Implies(False, True)
        True
        """
        try:
            newargs = []
            for x in args:
                if isinstance(x, Number) or x in (0, 1):
                    newargs.append(True if x else False)
                else:
                    newargs.append(x)
            A, B = newargs
        except ValueError:
            raise ValueError(
                "%d operand(s) used for an Implies "
                "(pairs are required): %s" % (len(args), str(args)))
        if A is True or A is False or B is True or B is False:
            return Or(Not(A), B)
        else:
            return Basic.__new__(cls, *args)


class Equivalent(BooleanFunction):
    """
    Equivalence relation.

    Equivalent(A, B) is True iff A and B are both True or both False
    """
    @classmethod
    def eval(cls, *args):
        """
        Equivalence relation.

        Returns True if all of the arguments are logically equivalent.
        Returns False otherwise.

        Examples
        ========

        >>> from sympy.logic.boolalg import Equivalent, And
        >>> from sympy.abc import x
        >>> Equivalent(False, False, False)
        True
        >>> Equivalent(True, False, False)
        False
        >>> Equivalent(x, And(x, True))
        True

        """

        newargs = []
        for x in args:
            if isinstance(x, Number) or x in (0, 1):
                newargs.append(True if x else False)
            else:
                newargs.append(x)
        argset = set(newargs)
        if len(argset) <= 1:
            return True
        if True in argset:
            argset.discard(True)
            return And(*argset)
        if False in argset:
            argset.discard(False)
            return Nor(*argset)
        return Basic.__new__(cls, *set(args))


class ITE(BooleanFunction):
    """
    If then else clause.
    """
    @classmethod
    def eval(cls, *args):
        """
        If then else clause

        ITE(A, B, C) evaluates and returns the result of B if A is true
        else it returns the result of C

        Examples
        ========

        >>> from sympy.logic.boolalg import ITE, And, Xor, Or
        >>> from sympy.abc import x, y, z
        >>> x = True
        >>> y = False
        >>> z = True
        >>> ITE(x, y, z)
        False
        >>> ITE(Or(x, y), And(x, z), Xor(z, x))
        True
        """
        args = list(args)
        if len(args) == 3:
            return Or(And(args[0], args[1]), And(Not(args[0]), args[2]))
        raise ValueError("ITE expects 3 arguments, but got %d: %s" %
                         (len(args), str(args)))

### end class definitions. Some useful methods


def conjuncts(expr):
    """Return a list of the conjuncts in the expr s.

    Examples
    ========

    >>> from sympy.logic.boolalg import conjuncts
    >>> from sympy.abc import A, B
    >>> conjuncts(A & B)
    frozenset([A, B])
    >>> conjuncts(A | B)
    frozenset([Or(A, B)])

    """
    return And.make_args(expr)


def disjuncts(expr):
    """Return a list of the disjuncts in the sentence s.

    Examples
    ========

    >>> from sympy.logic.boolalg import disjuncts
    >>> from sympy.abc import A, B
    >>> disjuncts(A | B)
    frozenset([A, B])
    >>> disjuncts(A & B)
    frozenset([And(A, B)])

    """
    return Or.make_args(expr)


def distribute_and_over_or(expr):
    """
    Given a sentence s consisting of conjunctions and disjunctions
    of literals, return an equivalent sentence in CNF.

    Examples
    ========

    >>> from sympy.logic.boolalg import distribute_and_over_or, And, Or, Not
    >>> from sympy.abc import A, B, C
    >>> distribute_and_over_or(Or(A, And(Not(B), Not(C))))
    And(Or(A, Not(B)), Or(A, Not(C)))
    """
    return _distribute((expr, And, Or))


def distribute_or_over_and(expr):
    """
    Given a sentence s consisting of conjunctions and disjunctions
    of literals, return an equivalent sentence in DNF.

    Note that the output is NOT simplified.

    Examples
    ========

    >>> from sympy.logic.boolalg import distribute_or_over_and, And, Or, Not
    >>> from sympy.abc import A, B, C
    >>> distribute_or_over_and(And(Or(Not(A), B), C))
    Or(And(B, C), And(C, Not(A)))
    """
    return _distribute((expr, Or, And))


def _distribute(info):
    """
    Distributes info[1] over info[2] with respect to info[0].
    """
    if info[0].func is info[2]:
        for arg in info[0].args:
            if arg.func is info[1]:
                conj = arg
                break
        else:
            return info[0]
        rest = info[2](*[a for a in info[0].args if a is not conj])
        return info[1](*list(map(_distribute,
            [(info[2](c, rest), info[1], info[2]) for c in conj.args])))
    elif info[0].func is info[1]:
        return info[1](*list(map(_distribute,
            [(x, info[1], info[2]) for x in info[0].args])))
    else:
        return info[0]


def to_cnf(expr, simplify=False):
    """
    Convert a propositional logical sentence s to conjunctive normal form.
    That is, of the form ((A | ~B | ...) & (B | C | ...) & ...)

    Examples
    ========

    >>> from sympy.logic.boolalg import to_cnf
    >>> from sympy.abc import A, B, D
    >>> to_cnf(~(A | B) | D)
    And(Or(D, Not(A)), Or(D, Not(B)))

    """
    expr = sympify(expr)
    if not isinstance(expr, BooleanFunction):
        return expr

    if simplify:
        simplified_expr = distribute_and_over_or(simplify_logic(expr))
        if len(simplified_expr.args) < len(to_cnf(expr).args):
            return simplified_expr
        else:
            return to_cnf(expr)

    # Don't convert unless we have to
    if is_cnf(expr):
        return expr

    expr = eliminate_implications(expr)
    return distribute_and_over_or(expr)


def to_dnf(expr, simplify=False):
    """
    Convert a propositional logical sentence s to disjunctive normal form.
    That is, of the form ((A & ~B & ...) | (B & C & ...) | ...)

    Examples
    ========

    >>> from sympy.logic.boolalg import to_dnf
    >>> from sympy.abc import A, B, C, D
    >>> to_dnf(B & (A | C))
    Or(And(A, B), And(B, C))

    """
    expr = sympify(expr)
    if not isinstance(expr, BooleanFunction):
        return expr

    if simplify:
        simplified_expr = distribute_or_over_and(simplify_logic(expr))
        if len(simplified_expr.args) < len(to_dnf(expr).args):
            return simplified_expr
        else:
            return to_dnf(expr)

    # Don't convert unless we have to
    if is_dnf(expr):
        return expr

    expr = eliminate_implications(expr)
    return distribute_or_over_and(expr)


def is_cnf(expr):
    """
    Test whether or not an expression is in conjunctive normal form.

    Examples
    ========

    >>> from sympy.logic.boolalg import is_cnf
    >>> from sympy.abc import A, B, C
    >>> is_cnf(A | B | C)
    True
    >>> is_cnf(A & B & C)
    True
    >>> is_cnf((A & B) | C)
    False

    """
    return _is_form(expr, And, Or)


def is_dnf(expr):
    """
    Test whether or not an expression is in disjunctive normal form.

    Examples
    ========

    >>> from sympy.logic.boolalg import is_dnf
    >>> from sympy.abc import A, B, C
    >>> is_dnf(A | B | C)
    True
    >>> is_dnf(A & B & C)
    True
    >>> is_dnf((A & B) | C)
    True
    >>> is_dnf(A & (B | C))
    False

    """
    return _is_form(expr, Or, And)


def _is_form(expr, function1, function2):
    """
    Test whether or not an expression is of the required form.

    """
    expr = sympify(expr)

    # Special case of an Atom
    if expr.is_Atom:
        return True

    # Special case of a single expression of function2
    if expr.func is function2:
        for lit in expr.args:
            if lit.func is Not:
                if not lit.args[0].is_Atom:
                    return False
            else:
                if not lit.is_Atom:
                    return False
        return True

    # Special case of a single negation
    if expr.func is Not:
        if not expr.args[0].is_Atom:
            return False

    if expr.func is not function1:
        return False

    for cls in expr.args:
        if cls.is_Atom:
            continue
        if cls.func is Not:
            if not cls.args[0].is_Atom:
                return False
        elif cls.func is not function2:
            return False
        for lit in cls.args:
            if lit.func is Not:
                if not lit.args[0].is_Atom:
                    return False
            else:
                if not lit.is_Atom:
                    return False

    return True


def eliminate_implications(expr):
    """
    Change >>, <<, and Equivalent into &, |, and ~. That is, return an
    expression that is equivalent to s, but has only &, |, and ~ as logical
    operators.

    Examples
    ========

    >>> from sympy.logic.boolalg import Implies, Equivalent, \
         eliminate_implications
    >>> from sympy.abc import A, B, C
    >>> eliminate_implications(Implies(A, B))
    Or(B, Not(A))
    >>> eliminate_implications(Equivalent(A, B))
    And(Or(A, Not(B)), Or(B, Not(A)))
    """
    expr = sympify(expr)
    if expr.is_Atom:
        return expr  # (Atoms are unchanged.)
    args = list(map(eliminate_implications, expr.args))
    if expr.func is Implies:
        a, b = args[0], args[-1]
        return (~a) | b
    elif expr.func is Equivalent:
        a, b = args[0], args[-1]
        return (a | Not(b)) & (b | Not(a))
    else:
        return expr.func(*args)


@deprecated(
    useinstead="sympify", issue=3451, deprecated_since_version="0.7.3")
def compile_rule(s):
    """
    Transforms a rule into a SymPy expression
    A rule is a string of the form "symbol1 & symbol2 | ..."

    Note: This function is deprecated.  Use sympify() instead.

    """
    import re
    return sympify(re.sub(r'([a-zA-Z_][a-zA-Z0-9_]*)', r'Symbol("\1")', s))


def to_int_repr(clauses, symbols):
    """
    Takes clauses in CNF format and puts them into an integer representation.

    Examples
    ========

    >>> from sympy.logic.boolalg import to_int_repr
    >>> from sympy.abc import x, y
    >>> to_int_repr([x | y, y], [x, y]) == [set([1, 2]), set([2])]
    True

    """

    # Convert the symbol list into a dict
    symbols = dict(list(zip(symbols, list(range(1, len(symbols) + 1)))))

    def append_symbol(arg, symbols):
        if arg.func is Not:
            return -symbols[arg.args[0]]
        else:
            return symbols[arg]

    return [set(append_symbol(arg, symbols) for arg in Or.make_args(c))
            for c in clauses]


def _check_pair(minterm1, minterm2):
    """
    Checks if a pair of minterms differs by only one bit. If yes, returns
    index, else returns -1.
    """
    index = -1
    for x, (i, j) in enumerate(zip(minterm1, minterm2)):
        if i != j:
            if index == -1:
                index = x
            else:
                return -1
    return index


def _convert_to_varsSOP(minterm, variables):
    """
    Converts a term in the expansion of a function from binary to it's
    variable form (for SOP).
    """
    temp = []
    for i, m in enumerate(minterm):
        if m == 0:
            temp.append(Not(variables[i]))
        elif m == 1:
            temp.append(variables[i])
        else:
            pass  # ignore the 3s
    return And(*temp)


def _convert_to_varsPOS(maxterm, variables):
    """
    Converts a term in the expansion of a function from binary to it's
    variable form (for POS).
    """
    temp = []
    for i, m in enumerate(maxterm):
        if m == 1:
            temp.append(Not(variables[i]))
        elif m == 0:
            temp.append(variables[i])
        else:
            pass  # ignore the 3s
    return Or(*temp)


def _simplified_pairs(terms):
    """
    Reduces a set of minterms, if possible, to a simplified set of minterms
    with one less variable in the terms using QM method.
    """
    simplified_terms = []
    todo = list(range(len(terms)))
    for i, ti in enumerate(terms[:-1]):
        for j_i, tj in enumerate(terms[(i + 1):]):
            index = _check_pair(ti, tj)
            if index != -1:
                todo[i] = todo[j_i + i + 1] = None
                newterm = ti[:]
                newterm[index] = 3
                if newterm not in simplified_terms:
                    simplified_terms.append(newterm)
    simplified_terms.extend(
        [terms[i] for i in [_ for _ in todo if _ is not None]])
    return simplified_terms


def _compare_term(minterm, term):
    """
    Return True if a binary term is satisfied by the given term. Used
    for recognizing prime implicants.
    """
    for i, x in enumerate(term):
        if x != 3 and x != minterm[i]:
            return False
    return True


def _rem_redundancy(l1, terms):
    """
    After the truth table has been sufficiently simplified, use the prime
    implicant table method to recognize and eliminate redundant pairs,
    and return the essential arguments.
    """
    essential = []
    for x in terms:
        temporary = []
        for y in l1:
            if _compare_term(x, y):
                temporary.append(y)
        if len(temporary) == 1:
            if temporary[0] not in essential:
                essential.append(temporary[0])
    for x in terms:
        for y in essential:
            if _compare_term(x, y):
                break
        else:
            for z in l1:
                if _compare_term(x, z):
                    if z not in essential:
                        essential.append(z)
                    break

    return essential


def SOPform(variables, minterms, dontcares=None):
    """
    The SOPform function uses simplified_pairs and a redundant group-
    eliminating algorithm to convert the list of all input combos that
    generate '1' (the minterms) into the smallest Sum of Products form.

    The variables must be given as the first argument.

    Return a logical Or function (i.e., the "sum of products" or "SOP"
    form) that gives the desired outcome. If there are inputs that can
    be ignored, pass them as a list, too.

    The result will be one of the (perhaps many) functions that satisfy
    the conditions.

    Examples
    ========

    >>> from sympy.logic import SOPform
    >>> minterms = [[0, 0, 0, 1], [0, 0, 1, 1],
    ...             [0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 1, 1]]
    >>> dontcares = [[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 1]]
    >>> SOPform(['w','x','y','z'], minterms, dontcares)
    Or(And(Not(w), z), And(y, z))

    References
    ==========

    .. [1] en.wikipedia.org/wiki/Quine-McCluskey_algorithm

    """
    from sympy.core.symbol import Symbol

    variables = [Symbol(v) if not isinstance(v, Symbol) else v
                 for v in variables]
    if minterms == []:
        return False

    minterms = [list(i) for i in minterms]
    dontcares = [list(i) for i in (dontcares or [])]
    for d in dontcares:
        if d in minterms:
            raise ValueError('%s in minterms is also in dontcares' % d)

    old = None
    new = minterms + dontcares
    while new != old:
        old = new
        new = _simplified_pairs(old)
    essential = _rem_redundancy(new, minterms)
    return Or(*[_convert_to_varsSOP(x, variables) for x in essential])


def POSform(variables, minterms, dontcares=None):
    """
    The POSform function uses simplified_pairs and a redundant-group
    eliminating algorithm to convert the list of all input combinations
    that generate '1' (the minterms) into the smallest Product of Sums form.

    The variables must be given as the first argument.

    Return a logical And function (i.e., the "product of sums" or "POS"
    form) that gives the desired outcome. If there are inputs that can
    be ignored, pass them as a list, too.

    The result will be one of the (perhaps many) functions that satisfy
    the conditions.

    Examples
    ========

    >>> from sympy.logic import POSform
    >>> minterms = [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1],
    ...             [1, 0, 1, 1], [1, 1, 1, 1]]
    >>> dontcares = [[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 1]]
    >>> POSform(['w','x','y','z'], minterms, dontcares)
    And(Or(Not(w), y), z)

    References
    ==========

    .. [1] en.wikipedia.org/wiki/Quine-McCluskey_algorithm

    """
    from sympy.core.symbol import Symbol

    variables = [Symbol(v) if not isinstance(v, Symbol) else v
                 for v in variables]
    if minterms == []:
        return False

    minterms = [list(i) for i in minterms]
    dontcares = [list(i) for i in (dontcares or [])]
    for d in dontcares:
        if d in minterms:
            raise ValueError('%s in minterms is also in dontcares' % d)

    maxterms = []
    for t in product([0, 1], repeat=len(variables)):
        t = list(t)
        if (t not in minterms) and (t not in dontcares):
            maxterms.append(t)
    old = None
    new = maxterms + dontcares
    while new != old:
        old = new
        new = _simplified_pairs(old)
    essential = _rem_redundancy(new, maxterms)
    return And(*[_convert_to_varsPOS(x, variables) for x in essential])


def simplify_logic(expr):
    """
    This function simplifies a boolean function to its
    simplified version in SOP or POS form. The return type is an
    Or or And object in SymPy. The input can be a string or a boolean
    expression.

    Examples
    ========

    >>> from sympy.logic import simplify_logic
    >>> from sympy.abc import x, y, z
    >>> from sympy import S
    >>> b = '(~x & ~y & ~z) | ( ~x & ~y & z)'
    >>> simplify_logic(b)
    And(Not(x), Not(y))

    >>> S(b)
    Or(And(Not(x), Not(y), Not(z)), And(Not(x), Not(y), z))
    >>> simplify_logic(_)
    And(Not(x), Not(y))

    """
    expr = sympify(expr)
    if not isinstance(expr, BooleanFunction):
        return expr
    variables = list(expr.free_symbols)
    truthtable = []
    for t in product([0, 1], repeat=len(variables)):
        t = list(t)
        if expr.subs(list(zip(variables, t))) == True:
            truthtable.append(t)
    if (len(truthtable) >= (2 ** (len(variables) - 1))):
        return SOPform(variables, truthtable)
    else:
        return POSform(variables, truthtable)


def _finger(eq):
    """
    Assign a 5-item fingerprint to each symbol in the equation:
    [
    # of times it appeared as a Symbol,
    # of times it appeared as a Not(symbol),
    # of times it appeared as a Symbol in an And or Or,
    # of times it appeared as a Not(Symbol) in an And or Or,
    sum of the number of arguments with which it appeared,
    counting Symbol as 1 and Not(Symbol) as 2
    ]

    >>> from sympy.logic.boolalg import _finger as finger
    >>> from sympy import And, Or, Not
    >>> from sympy.abc import a, b, x, y
    >>> eq = Or(And(Not(y), a), And(Not(y), b), And(x, y))
    >>> dict(finger(eq))
    {(0, 0, 1, 0, 2): [x], (0, 0, 1, 0, 3): [a, b], (0, 0, 1, 2, 8): [y]}

    So y and x have unique fingerprints, but a and b do not.
    """
    f = eq.free_symbols
    d = dict(list(zip(f, [[0] * 5 for fi in f])))
    for a in eq.args:
        if a.is_Symbol:
            d[a][0] += 1
        elif a.is_Not:
            d[a.args[0]][1] += 1
        else:
            o = len(a.args) + sum(ai.func is Not for ai in a.args)
            for ai in a.args:
                if ai.is_Symbol:
                    d[ai][2] += 1
                    d[ai][-1] += o
                else:
                    d[ai.args[0]][3] += 1
                    d[ai.args[0]][-1] += o
    inv = defaultdict(list)
    for k, v in ordered(iter(d.items())):
        inv[tuple(v)].append(k)
    return inv


def bool_equal(bool1, bool2, info=False):
    """Return True if the two expressions represent the same logical
    behavior for some correspondence between the variables of each
    (which may be different). For example, And(x, y) is logically
    equivalent to And(a, b) for {x: a, y: b} (or vice versa). If the
    mapping is desired, then set ``info`` to True and the simplified
    form of the functions and the mapping of variables will be
    returned.

    Examples
    ========

    >>> from sympy import SOPform, bool_equal, Or, And, Not, Xor
    >>> from sympy.abc import w, x, y, z, a, b, c, d
    >>> function1 = SOPform(['x','z','y'],[[1, 0, 1], [0, 0, 1]])
    >>> function2 = SOPform(['a','b','c'],[[1, 0, 1], [1, 0, 0]])
    >>> bool_equal(function1, function2, info=True)
    (And(Not(z), y), {y: a, z: b})

    The results are not necessarily unique, but they are canonical. Here,
    ``(w, z)`` could be ``(a, d)`` or ``(d, a)``:

    >>> eq =  Or(And(Not(y), w), And(Not(y), z), And(x, y))
    >>> eq2 = Or(And(Not(c), a), And(Not(c), d), And(b, c))
    >>> bool_equal(eq, eq2)
    True
    >>> bool_equal(eq, eq2, info=True)
    (Or(And(Not(y), w), And(Not(y), z), And(x, y)), {w: a, x: b, y: c, z: d})
    >>> eq = And(Xor(a, b), c, And(c,d))
    >>> bool_equal(eq, eq.subs(c, x), info=True)
    (And(Or(Not(a), Not(b)), Or(a, b), c, d), {a: a, b: b, c: d, d: x})

    """

    def match(function1, function2):
        """Return the mapping that equates variables between two
        simplified boolean expressions if possible.

        By "simplified" we mean that a function has been denested
        and is either an And (or an Or) whose arguments are either
        symbols (x), negated symbols (Not(x)), or Or (or an And) whose
        arguments are only symbols or negated symbols. For example,
        And(x, Not(y), Or(w, Not(z))).

        Basic.match is not robust enough (see issue 1736) so this is
        a workaround that is valid for simplified boolean expressions
        """

        # do some quick checks
        if function1.__class__ != function2.__class__:
            return None
        if len(function1.args) != len(function2.args):
            return None
        if function1.is_Symbol:
            return {function1: function2}

        # get the fingerprint dictionaries
        f1 = _finger(function1)
        f2 = _finger(function2)

        # more quick checks
        if len(f1) != len(f2):
            return False

        # assemble the match dictionary if possible
        matchdict = {}
        for k in f1.keys():
            if k not in f2:
                return False
            if len(f1[k]) != len(f2[k]):
                return False
            for i, x in enumerate(f1[k]):
                matchdict[x] = f2[k][i]
        return matchdict

    a = simplify_logic(bool1)
    b = simplify_logic(bool2)
    m = match(a, b)
    if m and info:
        return a, m
    return m is not None
