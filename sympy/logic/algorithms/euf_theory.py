"""
Congruence Closure Engine for EUF

The engine has 3 major parts:
    1) Union-find data structue for the congruence closure engine. Ref: [1], Section 2
    2) explain_classical() method to produce explanations for the engine's outputs. Ref: [1], Section 2-3
    3) Greedy algorithm with its expain() method to produce shorter explanations from explan_classical(). Ref: [2] Section 2.3, 2.4, 3 ,4

    Note that everything on explanation() method side, (2) and (3) use lazy versions of the algorithms.


Each consequent part builds from previous one. For example, 2) extends _union from 1) with extra data structures, and
3) extends from 2) by adding c-graphs and levels to edges.

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

Classes
-------
    EUFCongruenceClosure: Implements the congruence closure algorithm for EUF.
"""
from __future__ import annotations

from collections import defaultdict, deque
from heapq import heappush, heappop
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
    """

    def __init__(self, equations):
        """
        Parameters
        ----------
        equations : list of Q.eq or SymPy expressions
            The ground equalities to be saturated.
        """

        """
        Part 1) of the engine

        pending:  list of pairs of constants yet to be merged.
        representative:  returns the class representative of the constant
        classlist: set of all the constants in the class.
        lookup_table: maps an app to a constant.
        use_list: all the apps that use the input rep as an argument.

        """
        self.pending = deque()
        self.representative = {}                 # Representative[c]
        self.classlist = defaultdict(set)        # ClassList[rep]
        self.lookup_table = {}                   # Lookup_table[function, args]
        self.use_list = defaultdict(list)        # UseList[rep]

        # _flatten caches/stuff
        self._dummies = numbered_symbols('c', Dummy)
        self._term_to_const = {}                 # _term_to_const[expr] -> const. USED for not doing _flatten twice
        self._const_to_app = {}                  # const -> (func, arg consts). USED for Greedy algorithm
        self._lambda_cache = {}                  # Lambda -> const. USED for mapping equivalent lambdas to the same const

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
        Greedy algorithm and C-graph data structures
        c-graph is more extended version of the proof tree. It is an undirected graph s.t:
            1) it mirrors all the edges on the proof tree in union time with _insert_cgraph_edge
            2) it does not ignore redundant edges during union time with _insert_cgraph_edge
                E.g, take a->b, b->c, a->c. Proof tree would ingore a->c as we already have
                a path from a to c, in c-graph we also add a-c edge.
            3) it also contains lazily created extra possible edges with _compute_extra_edges


        adjacency: cgraph data structure
        _level: level of the edge is an integer that describes the date of the edge during _union,
            e.g lower level edge was created before high level edge. If an edge was created without
            _union (e.g _compute_extra_edges during explain()), the level is max level of its arguments
            plus 1, i.e the possible level it could have if it was created in _union.
        _level_counter: global counter that counts per edge.
        _n_edges_during_union: number of edges added to c-graph during the union
        _n_edges_extra: number of edges added to c-graph with _compute_extra_edges.
            From the paper's words, _compute_extra_edges could produce O(n^2) algorithm. We instead use _n_edges_extra
            and _n_edges_during_union to limit the number of operations. See section 4.2 of [2]
        _extra_edges_seen: the set of pairs _compute_extra_edges already checked. Used to not make duplicate edges
        greedy_fuel: huerustic discussed in the paper, paper default is 10



        """
        self.adjacency = defaultdict(list)       # const -> list of [(neighbor, label, level)]
        self._level = {}
        self._level_counter = 0
        self._n_edges_during_union = 0
        self._n_edges_extra = 0
        self._extra_edges_seen = set()
        self.greedy_fuel = 10

        # Transform every term of the input equations first, then merge.
        for eq in equations:
            if not (isinstance(eq, AppliedPredicate) and eq.function == Q.eq):
                raise EUFUnhandledInput
            left_id = self._flatten(eq.lhs)
            right_id = self._flatten(eq.rhs)
            self.pending.append((left_id, right_id, eq))
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
            const = self._record_app(expr.function, arg_ids)
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
            const = self._record_app(func_id, arg_ids)

        self._term_to_const[expr] = const
        return const

    def _record_app(self, func, arg_ids):
        """
        Record the application f(a) in the related data structures, and
        return a new dummy d that replaced it in _flatten i.e f(a) = d
        """
        rep_args = tuple(self._find_repr(arg) for arg in arg_ids)
        key = (func, rep_args)
        d = self._new_dummy()
        eq = (func, arg_ids, d)
        self._const_to_app[d] = (func, arg_ids)
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
        """
        TODO: add docs. arguably the most important method in this file.
        """
        rep_a, rep_b = self._find_repr(a), self._find_repr(b)
        if rep_a == rep_b:
            if isinstance(label, AppliedPredicate):
                if self._insert_cgraph_edge(a, b, label) is not None:
                    self._n_edges_during_union += 1
            return
        # Ensure |ClassList(a)| <= |ClassList(b)|
        if len(self.classlist[rep_a]) > len(self.classlist[rep_b]):
            rep_a, rep_b = rep_b, rep_a
            a, b = b, a
        self._insert_pf_edge(a, b, label)
        self._n_edges_during_union += 1
        level = self._insert_cgraph_edge(a, b, label)
        self._level[frozenset((a, b))] = level
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

    # ------ classical explain()

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

        # reverse the path, self explanators
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
        I.e return the node that is closes to root.
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
        explain without greedy algorithm.
        TODO: add more docs
        """
        self._aux_parent = {}
        pending_proofs = [(c1, c2)]
        while pending_proofs:
            x, y = pending_proofs.pop()
            c = self._nearest_common_ancestor(x, y)
            self._explain_along_path(x, c, output, pending_proofs)
            self._explain_along_path(y, c, output, pending_proofs)

    # ----------- greedy algorithm

    def _insert_cgraph_edge(self, a, b, label, level=None):
        """
        Insert a-b edge to the c-graph.
        """
        if a == b:
            return None
        if level is None:
            level = self._level_counter
            self._level_counter += 1
        self.adjacency[a].append((b, label, level))
        self.adjacency[b].append((a, label, level))
        return level

    def _estimate_size(self, label, memo):
        """
        Compute the size of the proof of the label
        memo is an optimization, stands for memoization.
        """
        # If label is not congruence, size is 1
        if label is None or isinstance(label, AppliedPredicate):
            return 1
        if label in memo:
            return memo[label]
        # if label is congruence i.e f(a,b) = f(c,d)
        # this should require sub-proof, so we walk through each pair of args
        (_, args1, _), (_, args2, _) = label
        size = 0
        for x, y in zip(args1, args2):
            if x != y:
                size += self._tree_path_size(x, y, memo)[0]
        memo[label] = size
        return size

    def _tree_path_size(self, a, b, memo):
        """
        TODO: perhaps we should divide this into multiple methods?
        memo is an optimization, stands for memoization.
        Returns
        =======
        size: proof tree size. compared to DAG size, this also counts duplications as paper requests.
            i.e size is the esimate of the classical explanation.
        merge_level: the level that the nodes a and b became equivalent.


        Reference
        =========
        [2] Section 3.1, 3.2
        """

        # First, find the nearest common ancestor.
        # this is very similar to _nearest_common_ancestor, although the prior
        # uses the secondary find-union, which we cannot do it here
        ancestors = set()
        cursor = a
        while True:
            ancestors.add(cursor)
            if cursor not in self.pf_parent:
                break
            cursor = self.pf_parent[cursor]
        nearest_common_ancestor = b
        while nearest_common_ancestor not in ancestors:
            nearest_common_ancestor = self.pf_parent[nearest_common_ancestor]

        size = 0
        merge_level = -1

        for cursor in (a, b):
            while cursor != nearest_common_ancestor:
                parent = self.pf_parent[cursor]
                size += self._estimate_size(self.pf_label[cursor], memo)
                merge_level = max(merge_level, self._level[frozenset((cursor, parent))])
                cursor = parent
        return size, merge_level

    def _get_canonical_form(self, d):
        """
        if f(a,b) = d, return f(a', b'). The name is a little bit misleading
        because it does not accept the application itself (maybe should change it?)
        Get the application that d did replace with in _flatten, and rewrite it to
        its canonical form. An application is in canonical form if its arguments
        are representatives. I.e this method replaces args with their representatives.

        Reference
        =========
        Algorithm  discussed in [2], section 2.4

        """
        app = self._const_to_app.get(d)
        if app is None:
            return d
        func, arg_ids = app
        return (func, tuple(self._find_repr(v) for v in arg_ids))

    def _compute_extra_edges(self, rep, memo):
        """

        Reference
        =========
        Algortihm  discussed in [2], section 2.4, 4.2
        Compared to the eager algorithm in the paper, this is a lazy version,
        i.e the extra edges are computed when explain() is called.
        """
        canonical_map = defaultdict(list) # canonical version -> original version
        # if a member has app equivalent, add it into the map
        for member in self.classlist[rep]:
            if member not in self._const_to_app:
                continue
            canonical_map[self._get_canonical_form(member)].append(member)
        for members in canonical_map.values():
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    # Inequality to limit the number of operations to not make this method O(n^2)
                    # See [2] Section 4.2
                    if self._n_edges_extra >= 2 * self._n_edges_during_union:
                        return
                    u, v = members[i], members[j]
                    key = frozenset((u, v))
                    # if we have already done this edge, skip to avoid duplicate
                    if key in self._extra_edges_seen:
                        continue
                    self._extra_edges_seen.add(key)
                    func_u, args_u = self._const_to_app[u]
                    func_v, args_v = self._const_to_app[v]
                    # do not get confused with this, any value smaller than 0 would work.
                    level = -1
                    for x, y in zip(args_u, args_v):
                        if x != y:
                            level = max(level, self._tree_path_size(x, y, memo)[1])
                    label = ((func_u, args_u, u), (func_v, args_v, v))
                    self._n_edges_extra += 1
                    self._insert_cgraph_edge(u, v, label, level=level + 1)

    def _shortest_path(self, a, b, memo, max_level):
        """
        Method discussed in [2] no paper clearly discusses this though.
        Compute the shortest weighted path from a to b in the c-graph.
        TODO: add more docs
        """
        distance = {a: 0}
        predecessor = {}
        # (distance, tiebreak, node)
        heap = [(0, 0, a)]
        tiebreak = 1
        settled = set()
        while heap:
            node_distance, _, node = heappop(heap)
            if node == b:
                break
            if node in settled:
                continue
            settled.add(node)
            for neighbor, label, level in self.adjacency[node]:
                if level > max_level or neighbor in settled:
                    continue
                neighbor_distance = node_distance + self._estimate_size(label, memo)
                if neighbor not in distance or neighbor_distance < distance[neighbor]:
                    distance[neighbor] = neighbor_distance
                    predecessor[neighbor] = (node, label)
                    heappush(heap, (neighbor_distance, tiebreak, neighbor))
                    tiebreak += 1
        if a == b:
            return []
        if b not in predecessor:
            return None
        path = []
        node = b
        while node != a:
            node, label = predecessor[node]
            path.append(label)
        return path

    def explain(self, lhs, rhs):
        """
        Arguably second most important method in this file, so
        TODO: add more docs
        """
        a = self._flatten(lhs)
        b = self._flatten(rhs)
        self._process_pending_unions()
        if self._find_repr(a) != self._find_repr(b):
            return None

        memo = {}
        output = set()
        fuel = self.greedy_fuel
        todo = deque([(a, b)])
        seen_pairs = set()
        extra_done = set()
        while todo:
            x, y = todo.popleft()
            if x == y:
                continue
            key = frozenset((x, y))
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            rep = self._find_repr(x)
            if rep not in extra_done:
                extra_done.add(rep)
                self._compute_extra_edges(rep, memo)
            max_level = self._tree_path_size(x, y, memo)[1]
            path = self._shortest_path(x, y, memo, max_level)
            if path is None:
                self._explain_classical(x, y, output)
                continue
            for label in path:
                if isinstance(label, AppliedPredicate):
                    output.add(label)
                elif label is not None:
                    (_, args1, _), (_, args2, _) = label
                    if fuel > 0:
                        fuel -= 1
                        for p, q in zip(args1, args2):
                            if p != q:
                                todo.append((p, q))
                    else:
                        for p, q in zip(args1, args2):
                            if p != q:
                                self._explain_classical(p, q, output)
        return output
