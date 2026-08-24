"""
Congruence Closure Engine for EUF

The engine has 3 major parts:
    1) Union-find data structue for the congruence closure engine. Ref: [1], Section 2
    2) explain() method to produce explanations for the engine's outputs. Ref: [1], Section 2-3
    3) backtrack() method

    Note that everything on explanation method side, (2)  use lazy versions of the algorithms.

References:
    [1] Fast congruence closure and extensions,
        https://www.cs.upc.edu/~oliveras/IC.pdf
    [2] Producing Shorter Congruence Closure Proofs in a State-of-the-Art SMT Solver
        https://link.springer.com/chapter/10.1007/978-3-032-15700-3_1

[2] Has a good introduction and background section, so it is recommended to read it in general.
It may be also helpful to read the important papers mentioned in [2] and/or [1]

Short introduction to terminologies:
    Constants: non-function terms. in e.g f(a,b) = d, d is a constant.
    Classes: Congruence classes. All the terms in the same class have
        relational and congruence properties. Shortly everything is
        knwon to be equal.
    Representative (repr): A constant term that represents the entire class.
        terms that have the same representative are in the same class.
        A representative of a constant a is usually written as a'.
    Application (app): A function with it args. I.e f(a) is an app, f is a function,
        and a is an argument.
    Ground term: a term with only constants and functions. f(a,g(b)) is a ground term, f(x, g(b)) is
        not as x is a variable.
    Ground equation: An equation where both sides are ground terms

Classes
-------
    EUFCongruenceClosure: Implements the congruence closure algorithm for EUF.
"""
from __future__ import annotations

from collections import defaultdict, deque
from sympy.core import Atom
from sympy.core.symbol import Symbol
from sympy.core.function import Lambda
from sympy.core.symbol import Dummy
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
    """

    def __init__(self, equations):
        """
        Parameters
        ----------
        equations : list of Q.eq or SymPy expressions
            The ground equalities to be saturated. Ground equalities are equalities with no variable args (see the header file)
        """

        """
        Part 1) of the engine

        pending:  list of pairs of constants yet to be merged.
        representative[c] -> rep:  returns the class representative of the constant
        classlist[repr] -> list of constants: set of all the constants in the class.
        lookup_table[func, repr_args] -> (func, args, d): maps an app (usually with repr args) to a eqaution.
            if it returns an app, what it means is that  the key and the output are congruent, and all the methods
            that use this data structure should  send them to pending for the merge (unless they are in the same class). So, there should never be a case that lookup_table should return more than one equation
            If it returns nothing/absent, usually the equation with non-rep args will be the output of the key application with the same equation but with repr args. (TODO: explain this better)
        use_list[rep] -> list of (func, args, d): all the apps that use the input rep as an argument.
            This is used in  two cases:
                1) whenever lookup_table returns nothing/absent, app the rep_args of the key of lookup_table store the equation.
                2) during _union, after the smaller class a is merged to class b, we check all the equation of a with _use_list[a] for rekeying

        """
        self.pending = deque()
        self.representative = {}
        self.classlist = defaultdict(set)
        self.lookup_table = {}
        self.use_list = defaultdict(list)

        # _flatten caches/stuff
        self._dummies = numbered_symbols('c', Dummy)
        self._term_to_const = {}                 # _term_to_const[expr] -> const. USED for not doing _flatten twice
        self._const_to_app = {}                  # const -> (func, arg consts). USED for Greedy algorithm

        """
        Part 2) of the engine
        These are proof forest data structures
        pf_parent: Since every edge is directional in pf, inputs a child and outputs the parent.
        pf_label (some papers call it justification): the label of the edge, why the edge was constructed.
        _aux_parent: secondary union-find data structure talked in the papers, used to not rewalk already walked
            paths in explain_classical().
        """
        self.pf_parent = {}                      # proof forest: const -> parent const
        self.pf_label = {}                       # const -> label of edge to parent
        self._aux_parent = {}                    # const ->  const

        """
        Part 3) of the engine
        _asserted: the equations currently asserted, in order. backtrack() drops
            the last n of them and rebuilds the closure from the rest.
        """
        self._asserted = []

        # Transform every term of the input equations first, then merge.
        for eq in equations:
            if not (isinstance(eq, AppliedPredicate) and eq.function == Q.eq):
                raise EUFUnhandledInput
            left_id = self._flatten(eq.lhs)
            right_id = self._flatten(eq.rhs)
            self.pending.append((left_id, right_id, eq))
            self._asserted.append(eq)
        self._process_pending_unions()

    def _register(self, const):
        """Initialize a constant as its own representative and an unique class for it."""
        if const not in self.representative:
            self.representative[const] = const
            self.classlist[const].add(const)

    def _new_dummy(self):
        d = next(self._dummies)
        self._register(d)
        return d

    def _flatten(self, expr):
        """
        flatten the expression. This method, in parallel, registers.
        This method should be called before any merging.

        Returns
        -------
        Symbol/Dummy : the new and unique constant that replaced the term.
        """
        # check if expr is flattened already
        if expr in self._term_to_const:
            return self._term_to_const[expr]

        # symbols are (should be) already proper constants so register them only
        if isinstance(expr, Symbol):
            self._register(expr)
            const = expr
        elif isinstance(expr, Atom):
            const = self._new_dummy()
        elif isinstance(expr, Lambda):
            # lambda-like, should be also currified
            lam = expr.curry()
            const = self._record_app((Lambda, lam.variables[0]),
                                     (self._flatten(lam.expr),))
        else:
            # function-like f(args) = const
            arg_ids = tuple(self._flatten(arg) for arg in expr.args)
            const = self._record_app(expr.func, arg_ids)

        self._term_to_const[expr] = const
        return const

    def _record_app(self, func, arg_ids):
        """
        Record the application f(a) in the related data structures, and
        return a new dummy d that replaced it in _flatten i.e f(a) = d
        """
        d = self._new_dummy()
        self._const_to_app[d] = (func, arg_ids)
        self._index_app(func, arg_ids, d)
        return d

    def _index_app(self, func, arg_ids, d):
        """
        Put the application f(a) = d into lookup_table and use_list under the
        representatives the classes have right now
        """
        rep_args = tuple(self._find_repr(arg) for arg in arg_ids)
        key = (func, rep_args)
        eq = (func, arg_ids, d)
        # check if there's another equation that are congruent
        if key in self.lookup_table:
            other_eq = self.lookup_table[key]
            self.pending.append((d, other_eq[2], (eq, other_eq)))
            self._process_pending_unions()
        else:
            # otherwise, add its own registry in lookup_table for future checks
            self.lookup_table[key] = eq
            for arg_id in set(rep_args):
                self.use_list[arg_id].append(eq)

    def _find_repr(self, const):
        return self.representative[const]

    def _union(self, a, b, label=None):
        """
        A method for the union of two classes (possibly the same) of terms a and b.
        """
        rep_a, rep_b = self._find_repr(a), self._find_repr(b)
        if rep_a == rep_b:
            return
        # Ensure |ClassList(a)| <= |ClassList(b)|
        if len(self.classlist[rep_a]) > len(self.classlist[rep_b]):
            rep_a, rep_b = rep_b, rep_a
            a, b = b, a
        # add c-graph and pf graph edges from a to b
        self._insert_pf_edge(a, b, label)
        # Move all members of ClassList(rep_a) into ClassList(rep_b)
        for c in self.classlist[rep_a]:
            self.representative[c] = rep_b
            self.classlist[rep_b].add(c)
        del self.classlist[rep_a]
        # For each equation (func, args, term) in UseList(rep_a)
        # i.e all the equation that includes rep_a:
        for eq in self.use_list.pop(rep_a, []):
            func, arg_ids, term = eq
            rep_args = tuple(self._find_repr(arg) for arg in arg_ids)
            key = (func, rep_args)
            # if there's a constant other = key
            if key in self.lookup_table:
                other_eq = self.lookup_table[key]
                # if they are not yet in the same class, add into pending for them to be in the
                # same class
                if self._find_repr(other_eq[2]) != self._find_repr(term):
                    self.pending.append((term, other_eq[2], (eq, other_eq)))
            else:
                # if there's no such key in lookup_table, add it to the register these for future
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
        eq = Q.eq(lhs, rhs)
        self._asserted.append(eq)
        self.pending.append((self._flatten(lhs), self._flatten(rhs), eq))
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

    # classical explain()

    def _insert_pf_edge(self, a, b, label):
        """
        Reverses the edges on the paath between a and the root of its tree
        Adds an directional edge a -> b with a label (i.e justification).

        This method was discussed in section 5. of [1]
        """
        path = []
        cursor = a
        # all roots have property of having no keys, so cursor will stop at root.
        while cursor in self.pf_parent:
            path.append((cursor, self.pf_parent[cursor], self.pf_label[cursor]))
            cursor = self.pf_parent[cursor]

        # reverse the path
        for child, parent, path_label in path:
            self.pf_parent[parent] = child
            self.pf_label[parent] = path_label
        # add edge a -> b
        self.pf_parent[a] = b
        # label (justification) is symmetric so reversing path
        # does not affect it. (label is in the form of Q.eq)
        self.pf_label[a] = label

    def _highest_node(self, a):
        """
        Return the highest node among all the nodes of the proof tree in the class of a.
        I.e return the node that is closest to root.
        """
        parent = self._aux_parent
        root = a
        # find the root
        while root in parent:
            root = parent[root]
        # path compression
        while a in parent:
            parent[a], a = root, parent[a]
        return root

    def _nearest_common_ancestor(self, a, b):
        """
        Compute the nearest common ancestor of the nodes a and b.
        """
        seen = set()
        cursor = a
        # walk to the root and store it in seen
        while True:
            cursor = self._highest_node(cursor)
            seen.add(cursor)
            # if it is the root
            if cursor not in self.pf_parent:
                break
            # move one edge
            cursor = self.pf_parent[cursor]
        cursor = b
        # when we find cursor in seen, stop
        while True:
            cursor = self._highest_node(cursor)
            if cursor in seen:
                return cursor
            # move one edge
            cursor = self.pf_parent[cursor]

    def _explain_along_path(self, a, c, output, pending_proofs):
        a = self._highest_node(a)
        while a != c:
            b = self.pf_parent[a]
            label = self.pf_label[a]
            # If label is Q.eq(a,b)
            if isinstance(label, AppliedPredicate):
                output.add(label)
            # if label is Q.eq(f(args),a) and Q.eq(g(args),b)
            elif label is not None:
                (_, args1, _), (_, args2, _) = label
                for x, y in zip(args1, args2):
                    if x != y:
                        pending_proofs.append((x, y))
            self._aux_parent[a] = self._highest_node(b)
            a = self._highest_node(a)

    def _explain_classical(self, c1, c2, output):
        """
        TODO: add more docs
        """
        self._aux_parent = {}
        pending_proofs = [(c1, c2)]
        while pending_proofs:
            x, y = pending_proofs.pop()
            c = self._nearest_common_ancestor(x, y)
            self._explain_along_path(x, c, output, pending_proofs)
            self._explain_along_path(y, c, output, pending_proofs)

    def explain(self, lhs, rhs):
        a = self._flatten(lhs)
        b = self._flatten(rhs)
        self._process_pending_unions()
        if self._find_repr(a) != self._find_repr(b):
            return None
        output = set()
        self._explain_classical(a, b, output)
        return output

    # backtracking
    # this is the simplest backtracking possible

    def _rebuild(self):
        """
        Discard everything that depends on which equations are asserted, then
        derive it again from self._asserted.
        """
        consts = list(self.representative)          # in registration order
        apps = list(self._const_to_app.items())
        asserted = self._asserted

        self.pending = deque()
        self.representative = {}
        self.classlist = defaultdict(set)
        self.lookup_table = {}
        self.use_list = defaultdict(list)
        self.pf_parent = {}
        self.pf_label = {}
        self._aux_parent = {}

        for const in consts:
            self._register(const)
        for d, (func, arg_ids) in apps:
            self._index_app(func, arg_ids, d)
        for eq in asserted:
            self.pending.append((self._flatten(eq.lhs), self._flatten(eq.rhs), eq))
            self._process_pending_unions()

    def backtrack(self, n):
        """
        Retract the last n asserted equations.

        TODO: Currently this just rebuilds the CC. The sources
        used in this file did not have proper backtracking information, and I am
        not even sure if it existed it would be faster than this kind of backtracking?

        Examples
        --------
        >>> from sympy import symbols
        >>> from sympy.assumptions.ask import Q
        >>> from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure
        >>> a, b, c = symbols('a b c')
        >>> cc = EUFCongruenceClosure([Q.eq(a, b)])
        >>> cc.merge(b, c)
        >>> cc.are_congruent(a, c)
        True
        >>> cc.backtrack(1)
        >>> cc.are_congruent(a, c)
        False
        >>> cc.are_congruent(a, b)
        True
        """
        if n == 0:
            return
        if not 0 < n <= len(self._asserted):
            raise ValueError("cannot backtrack %s of %s merge steps"
                             % (n, len(self._asserted)))
        del self._asserted[-n:]
        self._rebuild()
