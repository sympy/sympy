"""
Congruence Closure Engine for EUF

Implements:
    Nieuwenhuis & Oliveras, "Congruence Closure with Integer Offsets"
    https://link.springer.com/chapter/10.1007/978-3-540-39813-4_5

Open source version of the paper (with more context):
    Nieuwenhuis & Oliveras, "Congruence Closure with Integer Offsets"
    https://www.cs.upc.edu/~oliveras/IC.pdf

This algorithm efficiently computes the equivalence closure of a set of
ground equalities over uninterpreted functions, using union-find,
curryfication/flattening, and congruence propagation.

Terminology, variable names, and algorithm steps are consistent with those
shown in the paper, especially Section 4.

Short introduction to terms:
    Constants: non-function terms. in e.g f(a,b) = d, d is a constant.
    Classes: Congruence classes. All the terms in the same class have
        relational and congruence properties. Shortly everything is
        knwon to be equal.
    Representative (repr): A constant term that represents the entire class.
        terms that have the same representative are in the same class.

Example usage (doctest):

>>> from sympy import symbols, Function
>>> from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
>>> from sympy.assumptions.ask import Q
>>> f = Function('f')
>>> a, b, x, y = symbols('a b x y')
>>> cc = EUFCongruenceClosure([Q.eq(a, b), Q.eq(f(a), x), Q.eq(f(b), y)])
>>> cc.are_congruent(x, y)
True
>>> cc.are_congruent(a, x)
False

Classes
-------
    EUFCongruenceClosure: Implements the congruence closure algorithm for EUF.
"""

from collections import defaultdict, deque
from sympy.core.symbol import Symbol
from sympy.core.function import Lambda
from sympy.core.symbol import Dummy
from sympy.core.numbers import Number
from sympy.core import Basic
from sympy.utilities.iterables import numbered_symbols
from sympy.assumptions.assume import AppliedPredicate


class EUFUnhandledInput(Exception):
    """
    Raised while creating an EUFCongruenceClosure if unhandled input is present.
    """


class EUFCongruenceClosure:
    """
    Congruence closure algorithm for ground Equality with Uninterpreted Functions (EUF).

    See:
        Nieuwenhuis & Oliveras, "Congruence Closure with Integer Offsets"
        https://link.springer.com/chapter/10.1007/978-3-540-39813-4_5

    Major data structures (using algorithm's variable names):
        pending_unions:   deque, list of pairs of constants yet to be merged (PENDING).
        representative_table: dict, mapping: constant -> its class representative (REPRESENTATIVE).
        classlist: defaultdict(set), rep -> set of all elements in class (CLASSLIST).
        lookup_table: dict, maps (function, tuple of args) to a constant (LOOKUP).
        use_list: defaultdict(list), rep -> list of (func, args, result) triples (USELIST).
    """

    def __init__(self, equations):
        """
        Parameters
        ----------
        equations : list of Q.eq or SymPy expressions
            The ground equalities to be saturated.
        """
        self.pending_unions = deque()
        self.representative_table = {}           # Representative[c]
        self.classlist = defaultdict(set)        # ClassList[rep]
        self.lookup_table = {}                   # Lookup_table[function, args]
        self.use_list = defaultdict(list)        # UseList[rep]

        self._dummies = numbered_symbols('c', Dummy)
        # the terms that flattened to a constant
        self._term_to_const = {}                 # _term_to_const[expr] -> const
        self._lambda_cache = {}

        # Transform every term of the input equations first, then merge.
        for eq in equations:
            left_id = self._flatten(eq.lhs)
            right_id = self._flatten(eq.rhs)
            self.pending_unions.append((left_id, right_id))
        self._process_pending_unions()

    def _register(self, const):
        """
        Initialize the constant as its own representative in one element class.
        """
        if const not in self.representative_table:
            self.representative_table[const] = const
            self.classlist[const].add(const)

    def _new_dummy(self):
        d = next(self._dummies)
        self._register(d)
        return d

    def _flatten(self, expr):
        """
        Curryfy, and flatten the expression.
        This method MUST be called before any merging. It is a required step.

        Returns
        -------
        Symbol/Dummy : the new and unique constant that replaced the term.
        """
        if expr in self._term_to_const:
            return self._term_to_const[expr]

        if isinstance(expr, (Dummy, Symbol)):
            self._register(expr)
            const = expr
        elif isinstance(expr, Number) or getattr(expr, "is_Atom", False):
            const = self._new_dummy()
        elif isinstance(expr, AppliedPredicate):
            arg_ids = tuple(self._flatten(arg) for arg in expr.arguments)
            const = self._record_func_eq(expr.function, arg_ids)
        elif isinstance(expr, Lambda):
            lam = expr if len(expr.variables) == 1 else expr.curry()
            body_id = self._flatten(lam.expr)
            lam_key = Lambda(lam.variables[0], body_id)
            if lam_key not in self._lambda_cache:
                self._lambda_cache[lam_key] = self._new_dummy()
            const = self._lambda_cache[lam_key]
        else:
            func = expr.func
            func_id = self._flatten(func) if isinstance(func, Basic) else func
            arg_ids = tuple(self._flatten(arg) for arg in expr.args)
            const = self._record_func_eq(func_id, arg_ids)

        self._term_to_const[expr] = const
        return const

    def _record_func_eq(self, func, arg_ids):
        """Record the equation func(arg_ids) = d and return d."""
        key = (func, arg_ids)
        if key in self.lookup_table:
            return self.lookup_table[key]
        d = self._new_dummy()
        self.lookup_table[key] = d
        for arg_id in set(arg_ids):
            self.use_list[arg_id].append((func, arg_ids, d))
        return d

    def _const_of(self, term):
        """
        Return the constant that replaced the transformed term.
        If the term was constant to begin with, returns itself.
        """
        return self._term_to_const[term]

    def _find_repr(self, const):
        """
        Return the unique class representative for const.
        """
        return self.representative_table[const]

    def _union(self, a, b):
        rep_a, rep_b = self._find_repr(a), self._find_repr(b)
        if rep_a == rep_b:
            return
        # Ensure |ClassList(a)| <= |ClassList(b)|
        if len(self.classlist[rep_a]) > len(self.classlist[rep_b]):
            rep_a, rep_b = rep_b, rep_a
        # Move all members of ClassList(rep_a) into ClassList(rep_b)
        for c in list(self.classlist[rep_a]):
            self.representative_table[c] = rep_b
            self.classlist[rep_b].add(c)
        del self.classlist[rep_a]
        # For each application (func, args, term) in UseList(rep_a)
        for func, arg_ids, term in list(self.use_list.pop(rep_a, [])):
            rep_args = tuple(self._find_repr(arg) for arg in arg_ids)
            rep_term = self._find_repr(term)
            key = (func, rep_args)
            if key in self.lookup_table:
                other = self._find_repr(self.lookup_table[key])
                if other != rep_term:
                    self.pending_unions.append((rep_term, other))
            else:
                self.lookup_table[key] = rep_term
                self.use_list[rep_b].append((func, arg_ids, term))

    def _process_pending_unions(self):
        """
        Saturates pending_unions queue (Main loop, Paper Section 4).
        """
        while self.pending_unions:
            self._union(*self.pending_unions.popleft())

    def merge(self, lhs, rhs):
        """
        Merge the classes of two already-transformed terms and propagate
        closure.
        Raises KeyError if a term was not transformed.

        Examples
        --------
        >>> from sympy import symbols
        >>> from sympy.assumptions.ask import Q
        >>> from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
        >>> a, b, x, y = symbols('a b x y')
        >>> cc = EUFCongruenceClosure([Q.eq(a, x), Q.eq(b, y)])
        >>> cc.merge(a, b)
        >>> cc.are_congruent(x, y)
        True
        """
        self.pending_unions.append((self._const_of(lhs), self._const_of(rhs)))
        self._process_pending_unions()

    def are_congruent(self, lhs, rhs):
        """
        Query whether two terms are in the same class under the closure.
        Terms that were never transformed always return False.

        Examples
        --------
        >>> from sympy import symbols, Function
        >>> from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
        >>> from sympy.assumptions.ask import Q
        >>> f = Function('f')
        >>> x, y = symbols('x y')
        >>> cc = EUFCongruenceClosure([Q.eq(x, y), Q.eq(f(x), f(y))])
        >>> cc.are_congruent(x, y)
        True
        >>> cc.are_congruent(f(x), f(y))
        True
        >>> cc.are_congruent(x, f(x))
        False
        """
        try:
            lhs_id = self._const_of(lhs)
            rhs_id = self._const_of(rhs)
        except KeyError:
            return False
        return self._find_repr(lhs_id) == self._find_repr(rhs_id)
