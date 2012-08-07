"""Boolean algebra module for SymPy"""
from sympy.core.basic import Basic
from sympy.core.operations import LatticeOp
from sympy.core.function import Application, sympify

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

class And(LatticeOp, BooleanFunction):
    """
    Logical AND function.

    It evaluates its arguments in order, giving False immediately if any of them
    are False, and True if they are all True.

    Examples
    ========

        >>> from sympy.core import symbols
        >>> from sympy.abc import x, y
        >>> x & y
        And(x, y)
    """
    zero = False
    identity = True

class Or(LatticeOp, BooleanFunction):
    """
    Logical OR function

    It evaluates its arguments in order, giving True immediately if any of them are
    True, and False if they are all False.
    """
    zero = True
    identity = False

class Xor(BooleanFunction):
    """
    Logical XOR (exclusive OR) function.
    """
    @classmethod
    def eval(cls, *args):
        """
        Logical XOR (exclusive OR) function.

        Returns True if an odd number of the arguments are True, and the rest are False.
        Returns False if an even number of the arguments are True, and the rest are False.

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
        if not args: return False
        args = list(args)
        A = args.pop()
        while args:
            B = args.pop()
            A = Or(And(A, Not(B)), And(Not(A), B))
        return A

class Not(BooleanFunction):
    """
    Logical Not function (negation)

    Note: De Morgan rules applied automatically
    """

    is_Not = True

    @classmethod
    def eval(cls, *args):
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

        If multiple statements are given, returns an array of each result

        >>> Not(True, False)
        [False, True]
        >>> Not(True and False, True or False, True)
        [True, False, False]

        >>> Not(And(And(True, x), Or(x, False)))
        Not(x)
        """
        if len(args) > 1:
            return map(cls, args)
        arg = args[0]
        if type(arg) is bool:
            return not arg
        # apply De Morgan Rules
        if arg.func is And:
            return Or(*[Not(a) for a in arg.args])
        if arg.func is Or:
            return And(*[Not(a) for a in arg.args])
        if arg.func is Not:
            return arg.args[0]

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
            A, B = args
        except ValueError:
            raise ValueError("%d operand(s) used for an Implies (pairs are required): %s" % (len(args), str(args)))
        if A is True or A is False or B is True or B is False:
            return Or(Not(A), B)
        else:
            return Basic.__new__(cls, *args)

class Equivalent(BooleanFunction):
    """
    Equivalence relation.

    Equivalent(A, B) is True if and only if A and B are both True or both False
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

        argset = set(args)
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
        >>> from sympy.abc import x,y,z
        >>> x = True
        >>> y = False
        >>> z = True
        >>> ITE(x,y,z)
        False
        >>> ITE(Or(x, y), And(x, z), Xor(z, x))
        True
        """
        args = list(args)
        if len(args) == 3:
            return Or(And(args[0], args[1]), And(Not(args[0]), args[2]))
        raise ValueError("ITE expects 3 arguments, but got %d: %s" % (len(args), str(args)))

### end class definitions. Some useful methods

def fuzzy_not(arg):
    """
    Not in fuzzy logic

    Will return Not if arg is a boolean value, and None if argument
    is None.

    Examples:

    >>> from sympy.logic.boolalg import fuzzy_not
    >>> fuzzy_not(True)
    False
    >>> fuzzy_not(None)
    >>> fuzzy_not(False)
    True

    """
    if arg is None:
        return
    return not arg

def conjuncts(expr):
    """Return a list of the conjuncts in the expr s.

    Examples:

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

    Examples:

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
    if expr.func is Or:
        for arg in expr.args:
            if arg.func is And:
                conj = arg
                break
        else:
            return expr
        rest = Or(*[a for a in expr.args if a is not conj])
        return And(*map(distribute_and_over_or,
                   [Or(c, rest) for c in conj.args]))
    elif expr.func is And:
        return And(*map(distribute_and_over_or, expr.args))
    else:
        return expr

def to_cnf(expr):
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
    # Don't convert unless we have to
    if is_cnf(expr):
        return expr

    expr = sympify(expr)
    expr = eliminate_implications(expr)
    return distribute_and_over_or(expr)

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
    expr = sympify(expr)

    # Special case of a single disjunction
    if expr.func is Or:
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

    if expr.func is not And:
        return False

    for cls in expr.args:
        if cls.is_Atom:
            continue
        if cls.func is Not:
            if not cls.args[0].is_Atom:
                return False
        elif cls.func is not Or:
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

    >>> from sympy.logic.boolalg import Implies, Equivalent, eliminate_implications
    >>> from sympy.abc import A, B, C
    >>> eliminate_implications(Implies(A, B))
    Or(B, Not(A))
    >>> eliminate_implications(Equivalent(A, B))
    And(Or(A, Not(B)), Or(B, Not(A)))
    """
    expr = sympify(expr)
    if expr.is_Atom:
        return expr     ## (Atoms are unchanged.)
    args = map(eliminate_implications, expr.args)
    if expr.func is Implies:
        a, b = args[0], args[-1]
        return (~a) | b
    elif expr.func is Equivalent:
        a, b = args[0], args[-1]
        return (a | Not(b)) & (b | Not(a))
    else:
        return expr.func(*args)

def compile_rule(s):
    """
    Transforms a rule into a sympy expression
    A rule is a string of the form "symbol1 & symbol2 | ..."
    See sympy.assumptions.known_facts for examples of rules

    TODO: can this be replaced by sympify ?

    Examples
    ========

    >>> from sympy.logic.boolalg import compile_rule
    >>> compile_rule('A & B')
    And(A, B)
    >>> compile_rule('(~B & ~C)|A')
    Or(A, And(Not(B), Not(C)))
    """
    import re
    from sympy.core import Symbol
    return eval(re.sub(r'([a-zA-Z0-9_.]+)', r'Symbol("\1")', s), {'Symbol' : Symbol})


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
    symbols = dict(zip(symbols, xrange(1, len(symbols) + 1)))

    def append_symbol(arg, symbols):
        if arg.func is Not:
            return -symbols[arg.args[0]]
        else:
            return symbols[arg]

    return [set(append_symbol(arg, symbols) for arg in Or.make_args(c)) \
                                                            for c in clauses]
