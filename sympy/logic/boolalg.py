"""Boolean algebra module for SymPy"""
from sympy.core import Basic, Function, sympify, Symbol
from sympy.utilities import flatten, make_list
from sympy.core.operations import LatticeOp


class BooleanFunction(Function):
    """Boolean function is a function that lives in a boolean space
    It is used as base class for And, Or, Not, etc.
    """
    pass

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
            A = Or(A & Not(B), (Not(A) & B))
        return A

class Not(BooleanFunction):
    """Logical Not function (negation)

    Note: De Morgan rules applied automatically"""
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
        if not args:
            return False
        args = list(args)
        A = Not(args.pop())
        while args:
            B = args.pop()
            A = Or(A, Not(B))
        return A

class Nor(BooleanFunction):
    """Logical NOR function.
    It evaluates its arguments in order, giving False immediately if any
    of them are True, and True if they are all False.
    """
    @classmethod
    def eval(cls, *args):
        if not args:
            return False
        args = list(args)
        A = Not(args.pop())
        while args:
            B = args.pop()
            A = And(A, Not(B))
        return A

class Implies(BooleanFunction):
    pass

class Equivalent(BooleanFunction):
    """Equivalence relation.
    Equivalent(A, B) is True if and only if A and B are both True or both False
    """
    @classmethod
    def eval(cls, *args):
        return Basic.__new__(cls, *sorted(args))

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
    [A, B]
    >>> conjuncts(A | B)
    [Or(A, B)]

    """
    return make_list(expr, And)

def disjuncts(expr):
    """Return a list of the disjuncts in the sentence s.
    >>> from sympy.logic.boolalg import disjuncts
    >>> from sympy.abc import A, B
    >>> disjuncts(A | B)
    [A, B]
    >>> disjuncts(A & B)
    [And(A, B)]

    """
    return make_list(expr, Or)

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
    expr = sympify(expr)
    expr = eliminate_implications(expr)
    return distribute_and_over_or(expr)

def eliminate_implications(expr):
    """Change >>, <<, and Equivalent into &, |, and ~. That is, return an
    expression that is equivalent to s, but has only &, |, and ~ as logical
    operators.
    """
    expr = sympify(expr)
    if expr.is_Atom:
        return expr     ## (Atoms are unchanged.)
    args = map(eliminate_implications, expr.args)
    a, b = args[0], args[-1]
    if expr.func is Implies:
        return (~a) | b
    elif expr.func is Equivalent:
        return (a | Not(b)) & (b | Not(a))
    else:
        return type(expr)(*args)

def compile_rule(s):
    """Transforms a rule into a sympy expression
    A rule is a string of the form "symbol1 & symbol2 | ..."
    See sympy.assumptions.known_facts for examples of rules

    TODO: can this be replaced by sympify ?
    """
    import re
    return eval(re.sub(r'([a-zA-Z0-9_.]+)', r'Symbol("\1")', s), {'Symbol' : Symbol})


def to_int_repr(clauses, symbols):
    """
    takes clauses in CNF puts them into integer representation

    Examples:
        >>> from sympy.logic.boolalg import to_int_repr
        >>> from sympy.abc import x, y
        >>> to_int_repr([x | y, y], [x, y])
        [[1, 2], [2]]

    """
    def append_symbol(arg, symbols):
        if arg.func is Not:
            return -(symbols.index(arg.args[0])+1)
        else:
            return symbols.index(arg)+1

    output = []
    for c in clauses:
        if c.func is Or:
            t = []
            for arg in c.args:
                t.append(append_symbol(arg, symbols))
            output.append(t)
        else:
            output.append([append_symbol(c, symbols)])
    return output
