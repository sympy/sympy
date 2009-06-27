"""Boolean algebra module for SymPy"""
from sympy.core import Basic, Function, sympify, Symbol
from sympy.utilities import flatten


class BooleanFunction(Function):
    """Boolean function is a function that lives in a boolean space
    It is used as base class for And, Or, Not, etc.
    """
    pass

class And(BooleanFunction):
    """Logical AND function.
    It evaluates its arguments in order, giving False immediately if any of them
    are False, and True if they are all True.

    Examples:
        >>> from sympy import *
        >>> x, y = symbols('xy')
        >>> x & y
        And(x, y)
    """
    @classmethod
    def eval(cls, *args):
        out_args = []
        for arg in args: # we iterate over a copy or args
            if isinstance(arg, bool):
                if not arg: return False
                else: continue
            out_args.append(arg)
        if len(out_args) == 0: return True
        if len(out_args) == 1: return out_args[0]
        sargs = sorted(flatten(out_args, cls=cls))
        return Basic.__new__(cls, *sargs)

class Or(BooleanFunction):
    """Logical OR function
     It evaluates its arguments in order, giving True immediately if any of them are
     True, and False if they are all False.
     """
    @classmethod
    def eval(cls, *args):
        out_args = []
        for arg in args: # we iterate over a copy or args
            if isinstance(arg, bool):
                if arg: return True
                else: continue
            out_args.append(arg)
        if len(out_args) == 0: return False
        if len(out_args) == 1: return out_args[0]
        sargs = sorted(flatten(out_args, cls=cls))
        return Basic.__new__(cls, *sargs)

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
        # apply De Morgan Rules
        if isinstance(arg, And):
            return Or( *[Not(a) for a in arg.args])
        if isinstance(arg, Or):
            return And(*[Not(a) for a in arg.args])
        if isinstance(arg, bool): return not arg
        if isinstance(arg, Not):
            return arg.args[0]

class Nand(BooleanFunction):
    """Logical NAND function.
    It evaluates its arguments in order, giving True immediately if any
    of them are False, and False if they are all True.
    """
    @classmethod
    def eval(cls, *args):
        if not args: return False
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
        if not args: return False
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

def conjuncts(expr):
    """Return a list of the conjuncts in the expr s.
    >>> from sympy import symbols
    >>> A, B = symbols('AB')
    >>> conjuncts(A & B)
    [A, B]
    >>> conjuncts(A | B)
    [Or(A, B)]
    """
    if isinstance(expr, And):
        return list(expr.args)
    else:
        return [expr]

def disjuncts(expr):
    """Return a list of the disjuncts in the sentence s.
    >>> from sympy import symbols
    >>> A, B = symbols('AB')
    >>> disjuncts(A | B)
    [A, B]
    >>> disjuncts(A & B)
    [And(A, B)]
    """
    if isinstance(expr, Or):
        return list(expr.args)
    else:
        return [expr]

def distribute_and_over_or(expr):
    """Given a sentence s consisting of conjunctions and disjunctions
    of literals, return an equivalent sentence in CNF.
    """
    if isinstance(expr, Or):
        for arg in expr.args:
            if isinstance(arg, And):
                conj = arg
                break
        else: return type(expr)(*expr.args)
        rest = Or(*[a for a in expr.args if a is not conj])
        return And(*map(distribute_and_over_or,
                   [Or(c, rest) for c in conj.args]))
    elif isinstance(expr, And):
        return And(*map(distribute_and_over_or, expr.args))
    else:
        return expr

def to_cnf(expr):
    """Convert a propositional logical sentence s to conjunctive normal form.
    That is, of the form ((A | ~B | ...) & (B | C | ...) & ...)
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
    if expr.is_Atom: return expr     ## (Atoms are unchanged.)
    args = map(eliminate_implications, expr.args)
    a, b = args[0], args[-1]
    if isinstance(expr, Implies):
        return (~a) | b
    elif isinstance(expr, Equivalent):
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
