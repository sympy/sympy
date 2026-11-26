"""
Congruence Closure Engine for EUF
Implements:
    Nieuwenhuis & Oliveras, "Congruence Closure with Integer Offsets"
    https://www.cs.upc.edu/~oliveras/dpllt.pdf
This algorithm efficiently computes the equivalence closure of a set of
ground equalities over uninterpreted functions, using union-find,
curryfication/flattening, and congruence propagation.
Terminology, variable names, and algorithm steps are consistent with those
shown in the paper, especially Section 4.
Example usage (doctest):
>>> from sympy import symbols, Function, Eq
>>> from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
>>> from sympy.assumptions.ask import Q
>>> f, g = symbols('f g', cls=Function)
>>> a, b, x, c = symbols('a b x c')
>>> cc = EUFCongruenceClosure([Q.eq(a, b), Q.eq(f(a), x)])
>>> cc.are_equal(f(b), x)
True
>>> cc.add_equality(a, c)
>>> cc.are_equal(a, c)
True
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
        https://www.cs.upc.edu/~oliveras/dpllt.pdf
    Major data structures (using algorithm's variable names):
        pending_unions:   deque, list of pairs of constants yet to be merged (PENDING).
        representative_table: dict, mapping: constant -> its class representative (REPRESENTATIVE).
        classlist: defaultdict(set), rep -> set of all elements in class (CLASSLIST).
        lookup_table: dict, maps (function, tuple of args) to a constant (LOOKUP).
        use_list: defaultdict(list), rep -> list of (func, args, result) triples (USELIST).
    Curryfication and flattening (see Sec 3) are handled as in the paper.
    """

    def __init__(self, equations):
        """
        Parameters
        ----------
        equations : list of Q.eq or SymPy expressions
            The ground equalities to be saturated.
        """
        # ----- Section 4, Paper -----
        self.pending_unions = deque()
        self.representative_table = {}           # Representative[c]
        self.classlist = defaultdict(set)        # ClassList[rep]
        self.lookup_table = {}                   # Lookup_table[function, args]
        self.use_list = defaultdict(list)        # UseList[rep]

        self._dummies = numbered_symbols('c', Dummy)
        self._num_dummy_cache = {}
        self._atom_dummy_cache = {}
        self._flatten_cache = {}
        self._predicate_cache = {}
        self._lambda_pattern_cache = {}

        # Curryfication and flattening (Sec 3/4): every unique subterm assigned constant
        for eq in equations:
            left_id = self._flatten(eq.lhs)
            right_id = self._flatten(eq.rhs)
            self._register(left_id)
            self._register(right_id)
            self.pending_unions.append((left_id, right_id))
        self._process_pending_unions()

    def _register(self, const):
        """Ensure const is in class structures as its own singleton."""
        if const not in self.representative_table:
            self.representative_table[const] = const
            self.classlist[const].add(const)

    def _new_dummy(self):
        d = next(self._dummies)
        self._register(d)
        return d

    def _flatten(self, expr):
        """
        Curryfy and flatten expr, assign new dummy for each unique term as in Sec. 3.
        Returns
        -------
        Symbol/Dummy : unique id for the term subtree.
        """
        if isinstance(expr, (Dummy, Symbol)):
            self._register(expr)
            return expr
        if isinstance(expr, Number):
            if expr not in self._num_dummy_cache:
                d = self._new_dummy()
                self._num_dummy_cache[expr] = d
            return self._num_dummy_cache[expr]
        if getattr(expr, "is_Atom", False):
            if expr not in self._atom_dummy_cache:
                d = self._new_dummy()
                self._atom_dummy_cache[expr] = d
            return self._atom_dummy_cache[expr]
        if isinstance(expr, AppliedPredicate):
            pred = expr.function
            args = expr.arguments
            arg_ids = tuple(self._flatten(arg) for arg in args)
            key = (pred, arg_ids)
            if key not in self._flatten_cache:
                d = self._new_dummy()
                self._flatten_cache[key] = d
                self.lookup_table[key] = d
                for arg_id in arg_ids:
                    self.use_list[arg_id].append((pred, arg_ids, d))
            return self._flatten_cache[key]
        if isinstance(expr, Lambda):
            lam = expr if len(expr.variables) == 1 else expr.curry()
            body_id = self._flatten(lam.expr)
            lam_key = Lambda(lam.variables[0], body_id)
            if lam_key not in self._flatten_cache:
                d = self._new_dummy()
                self._flatten_cache[lam_key] = d
            return self._flatten_cache[lam_key]
        func = expr.func
        func_id = self._flatten(func) if isinstance(func, Basic) and not isinstance(expr, Lambda) else func
        arg_ids = tuple(self._find(self._flatten(arg)) for arg in expr.args)
        key = (func_id, arg_ids)
        if key not in self._flatten_cache:
            d = self._new_dummy()
            self._flatten_cache[key] = d
            self.lookup_table[key] = d
            for arg in arg_ids:
                self.use_list[arg].append((func_id, arg_ids, d))
        return self._flatten_cache[key]

    def _find(self, const):
        """
        Return the unique class representative for const (with path compression).
        """
        if const not in self.representative_table:
            self._register(const)
        root = const
        # Find root
        while root != self.representative_table[root]:
            root = self.representative_table[root]
        # Path compression
        while const != root:
            parent = self.representative_table[const]
            self.representative_table[const] = root
            const = parent
        return root

    def _union(self, a, b):
        rep_a, rep_b = self._find(a), self._find(b)
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
            rep_args = tuple(self._find(arg) for arg in arg_ids)
            rep_term = self._find(term)
            key = (func, rep_args)
            if key in self.lookup_table:
                other = self._find(self.lookup_table[key])
                if other != rep_term:
                    self.pending_unions.append((rep_term, other))
            self.lookup_table[key] = rep_term
            self.use_list[rep_b].append((func, arg_ids, term))
            self._flatten_cache[key] = rep_term

    def _process_pending_unions(self):
        """
        Saturates pending_unions queue (Main loop, Paper Section 4).
        """
        while self.pending_unions:
            self._union(*self.pending_unions.popleft())

    def add_equality(self, lhs, rhs):
        """
        Incrementally add a new equality and propagate closure.
        Examples
        --------
        >>> from sympy import symbols
        >>> from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
        >>> a, b = symbols('a b')
        >>> cc = EUFCongruenceClosure([])
        >>> cc.add_equality(a, b)
        >>> cc.are_equal(a, b)
        True
        """
        self.pending_unions.append((self._flatten(lhs), self._flatten(rhs)))
        self._process_pending_unions()


    def _replace_with_representative(self, expr):
        """
        Recursively replace all symbols and subexpressions in expr with their
        current class representatives using union-find _find.
        """
        from sympy import Symbol, Lambda

        if isinstance(expr, (Symbol, Dummy)):
            # Return the representative of the symbol (already registered)
            rep = self._find(expr)
            return rep

        # Terminal cases: Number or atom without proper class
        if getattr(expr, "is_Atom", False):
            return expr

        # Recursive case: rebuild expression with replaced args
        new_args = tuple(self._replace_with_representative(arg) for arg in expr.args)

        # For Lambda and other special classes, preserve structure
        if isinstance(expr, Lambda):
            # Rebuild the Lambda with replaced body
            return Lambda(expr.variables[0], new_args[0])

        # Rebuild with same function and replaced args
        return expr.func(*new_args)


    def are_equal(self, lhs, rhs):
        """
        Query whether two terms are in the same class under the closure.
        Examples
        --------
        >>> from sympy import symbols, Function
        >>> from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
        >>> from sympy.assumptions.ask import Q
        >>> f = Function('f')
        >>> x, y = symbols('x y')
        >>> cc = EUFCongruenceClosure([Q.eq(x, y), Q.eq(f(x), f(y))])
        >>> cc.are_equal(x, y)
        True
        >>> cc.are_equal(f(x), f(y))
        True
        """
        lhs_repl = self._replace_with_representative(lhs)
        rhs_repl = self._replace_with_representative(rhs)
        lhs_id, rhs_id = self._flatten(lhs_repl), self._flatten(rhs_repl)
        if self._find(lhs_id) == self._find(rhs_id):
            return True
