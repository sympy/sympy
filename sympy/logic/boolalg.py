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

class And(LatticeOp, BooleanFunction):
    """
    Logical AND function.

    It evaluates its arguments in order, giving False immediately if any of them
    are False, and True if they are all True.

    Examples:
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
    """Logical XOR (exclusive OR) function.
    returns True if an odd number of the arguments are True, and the rest are False.
    returns False if an even number of the arguments are True, and the rest are False.
    """
    @classmethod
    def eval(cls, *args):
        if not args: return False
        args = list(args)
        A = args.pop()
        while args:
            B = args.pop()
            A = Or(And(A, Not(B)), And(Not(A), B))
        return A

class Not(BooleanFunction):
    """Logical Not function (negation)

    Note: De Morgan rules applied automatically"""

    is_Not = True

    @classmethod
    def eval(cls, *args):
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
    """Logical NAND function.

    It evaluates its arguments in order, giving True immediately if any
    of them are False, and False if they are all True.
    """
    @classmethod
    def eval(cls, *args):
        return Not(And(*args))

class Nor(BooleanFunction):
    """Logical NOR function.

    It evaluates its arguments in order, giving False immediately if any
    of them are True, and True if they are all False.
    """
    @classmethod
    def eval(cls, *args):
        return Not(Or(*args))

class Implies(BooleanFunction):
    """Logical implication.

    A implies B is equivalent to !A v B
    """
    @classmethod
    def eval(cls, *args):
        if len(args) != 2:
            raise ValueError, "%d operand(s) used for an Implies (pairs are required): %s" % (len(args), str(args))
        else:
            return Or(Not(args[0]), args[1])

class Equivalent(BooleanFunction):
    """Equivalence relation.

    Equivalent(A, B) is True if and only if A and B are both True or both False
    """
    @classmethod
    def eval(cls, *args):
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

### end class definitions. Some useful methods

def fuzzy_not(arg):
    """
    Not in fuzzy logic

    will return Not if arg is a boolean value, and None if argument
    is None

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
    """Convert a propositional logical sentence s to conjunctive normal form.
    That is, of the form ((A | ~B | ...) & (B | C | ...) & ...)

    Examples:

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
    """Test whether or not an expression is in conjunctive normal form.

    Examples:

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

    if not expr.func is And:
        return False

    for cls in expr.args:
        if cls.is_Atom:
            continue
        if cls.func is Not:
            if not cls.args[0].is_Atom:
                return False
        elif not cls.func is Or:
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
    """Change >>, <<, and Equivalent into &, |, and ~. That is, return an
    expression that is equivalent to s, but has only &, |, and ~ as logical
    operators.
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
    """Transforms a rule into a sympy expression
    A rule is a string of the form "symbol1 & symbol2 | ..."
    See sympy.assumptions.known_facts for examples of rules

    TODO: can this be replaced by sympify ?
    """
    import re
    from sympy.core import Symbol
    return eval(re.sub(r'([a-zA-Z0-9_.]+)', r'Symbol("\1")', s), {'Symbol' : Symbol})


def to_int_repr(clauses, symbols):
    """
    takes clauses in CNF puts them into integer representation

    Examples:
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
