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
from __future__ import annotations

from collections import defaultdict, deque
from sympy.core.symbol import Symbol
from sympy.core.function import Lambda
from sympy.core.symbol import Dummy
from sympy.core.numbers import Number
from sympy.core import Basic
from sympy.utilities.iterables import numbered_symbols
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.ask import Q


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
        pending:   deque, list of pairs of constants yet to be merged (PENDING).
        representative: dict, mapping: constant -> its class representative (REPRESENTATIVE).
        classlist: defaultdict(set), rep -> set of all elements in class (CLASSLIST).
        lookup_table: dict, maps (function, tuple of args) to the equation
            recording that application, as a (func, args, result) triple (LOOKUP).
        use_list: defaultdict(list), rep -> list of (func, args, result) triples (USELIST).
    """

    def __init__(self, equations):
        """
        Parameters
        ----------
        equations : list of Q.eq or SymPy expressions
            The ground equalities to be saturated.
        """
        self.pending = deque()
        self.representative = {}                 # Representative[c]
        self.classlist = defaultdict(set)        # ClassList[rep]
        self.lookup_table = {}                   # Lookup_table[function, args]
        self.use_list = defaultdict(list)        # UseList[rep]

        self._dummies = numbered_symbols('c', Dummy)
        # the terms that flattened to a constant
        self._term_to_const = {}                 # _term_to_const[expr] -> const
        self._lambda_cache = {}

        """
        Proof forest data structures. This i
        """
        self.pf_parent = {}                      # proof forest: const -> parent const
        self.pf_label = {}                       # const -> label of edge to parent
        # additional Union-Find data structure for explain
        self._aux_parent = {}

        # Transform every term of the input equations first, then merge.
        for eq in equations:
            if not (isinstance(eq, AppliedPredicate) and eq.function == Q.eq):
                raise EUFUnhandledInput
            left_id = self._flatten(eq.lhs)
            right_id = self._flatten(eq.rhs)
            self.pending.append((left_id, right_id, eq))
        self._process_pending_unions()

    def _register(self, const):
        if const not in self.representative:
            self.representative[const] = const
            self.classlist[const].add(const)

    def _new_dummy(self):
        d = next(self._dummies)
        self._register(d)
        return d

    def _flatten(self, expr):
        """
        Curryfy, and flatten the expression. This method, in parallel, registers.
        This method should be called before any merging.

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
            body_id = self._find_repr(self._flatten(lam.expr))
            lam_key = Lambda(lam.variables[0], body_id)
            if lam_key not in self._lambda_cache:
                self._lambda_cache[lam_key] = self._new_dummy()
            const = self._lambda_cache[lam_key]
        else:
            func = expr.func
            func_id = self._find_repr(self._flatten(func)) if isinstance(func, Basic) else func
            arg_ids = tuple(self._flatten(arg) for arg in expr.args)
            const = self._record_func_eq(func_id, arg_ids)

        self._term_to_const[expr] = const
        return const

    def _record_func_eq(self, func, arg_ids):
        rep_args = tuple(self._find_repr(arg) for arg in arg_ids)
        key = (func, rep_args)
        d = self._new_dummy()
        eq = (func, arg_ids, d)
        if key in self.lookup_table:
            other = self.lookup_table[key]
            self.pending.append((d, other[2], (eq, other)))
            self._process_pending_unions()
        else:
            self.lookup_table[key] = eq
            for arg_id in set(rep_args):
                self.use_list[arg_id].append(eq)
        return d

    def _find_repr(self, const):
        return self.representative[const]

    def _union(self, a, b, label=None):
        rep_a, rep_b = self._find_repr(a), self._find_repr(b)
        if rep_a == rep_b:
            return
        # Ensure |ClassList(a)| <= |ClassList(b)|
        if len(self.classlist[rep_a]) > len(self.classlist[rep_b]):
            rep_a, rep_b = rep_b, rep_a
            a, b = b, a
        self._insert_edge(a, b, label)
        # Move all members of ClassList(rep_a) into ClassList(rep_b)
        for c in self.classlist[rep_a]:
            self.representative[c] = rep_b
            self.classlist[rep_b].add(c)
        del self.classlist[rep_a]
        # For each application (func, args, term) in UseList(rep_a)
        for eq in self.use_list.pop(rep_a, []):
            func, arg_ids, term = eq
            rep_args = tuple(self._find_repr(arg) for arg in arg_ids)
            key = (func, rep_args)
            if key in self.lookup_table:
                other = self.lookup_table[key]
                if self._find_repr(other[2]) != self._find_repr(term):
                    self.pending.append((term, other[2], (eq, other)))
            else:
                self.lookup_table[key] = eq
                self.use_list[rep_b].append(eq)

    def _process_pending_unions(self):
        """
        Saturates pending_unions queue (Main loop, Paper Section 4).
        """
        while self.pending:
            self._union(*self.pending.popleft())

    def merge(self, lhs, rhs):
        """
        Merge the classes of two already-transformed terms and propagate
        closure.

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
        self.pending.append((self._flatten(lhs), self._flatten(rhs), Q.eq(lhs, rhs)))
        self._process_pending_unions()

    def are_congruent(self, lhs, rhs):
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
        >>> cc.are_congruent(x, y)
        True
        >>> cc.are_congruent(f(x), f(y))
        True
        >>> cc.are_congruent(x, f(x))
        False
        """
        lhs_id = self._flatten(lhs)
        rhs_id = self._flatten(rhs)
        return self._find_repr(lhs_id) == self._find_repr(rhs_id)

    def _insert_edge(self, a, b, label):
        """
        Reverses the edges on the paath between a and the root of its tree
        Adds an directional edge a -> b with a label (i.e justification).

        This method is parallel to the addition to the Union in the paper.
        """
        path = []
        cursor = a
        # all roots have property of having no keys, so cursor will stop at root.
        while cursor in self.pf_parent:
            path.append((cursor, self.pf_parent[cursor], self.pf_label[cursor]))
            cursor = self.pf_parent[cursor]

        # reverse the path, self explanators
        for child, parent, path_label in path:
            self.pf_parent[parent] = child
            self.pf_label[parent] = path_label
        # add edge a -> b
        self.pf_parent[a] = b
        # label (justification) is symmetric so reversing path
        # does not affect it. (label is in the form of Q.eq)
        self.pf_label[a] = label

    def _highest_node(self, x):
        """
        Return the highest node among all the nodes of the proof tree in the class of x.
        I.e return the node that is closes to root.
        """
        parent = self._aux_parent
        root = x
        # find the root
        while root in parent:
            root = parent[root]
        # path compression
        while x in parent:
            parent[x], x = root, parent[x]
        return root

    def _nearest_common_ancestor(self, a, b):
        seen = set()
        x = a
        while True:
            x = self._highest_node(x)
            seen.add(x)
            if x not in self.pf_parent:
                break
            x = self.pf_parent[x]
        y = b
        while True:
            y = self._highest_node(y)
            if y in seen:
                return y
            y = self.pf_parent[y]

    def _explain_along_path(self, a, c, output, pending_proofs):
        a = self._highest_node(a)
        while a != c:
            b = self.pf_parent[a]
            label = self.pf_label[a]
            if isinstance(label, AppliedPredicate):
                output.add(label)
            elif label is not None:
                (_, args1, _), (_, args2, _) = label
                for x, y in zip(args1, args2):
                    if x != y:
                        pending_proofs.append((x, y))
            self._aux_parent[a] = self._highest_node(b)
            a = self._highest_node(a)

    def explain(self, lhs, rhs):
        a = self._flatten(lhs)
        b = self._flatten(rhs)
        self._process_pending_unions()
        if self._find_repr(a) != self._find_repr(b):
            return None
        self._aux_parent = {}
        output = set()
        pending_proofs = [(a, b)]
        while pending_proofs:
            x, y = pending_proofs.pop()
            c = self._nearest_common_ancestor(x, y)
            self._explain_along_path(x, c, output, pending_proofs)
            self._explain_along_path(y, c, output, pending_proofs)
        return output
