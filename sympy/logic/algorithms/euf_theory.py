"""
EUF Congruence Closure for Equality with Uninterpreted Functions (EUF)
======================================================================

Implements the congruence-closure algorithm for the ground theory
of Equality with Uninterpreted Functions (EUF) as described in:

    "Congruence Closure with Integer Offsets"
    Robert Nieuwenhuis and Albert Oliveras,
    Technical University of Catalonia, 2003,
    https://www.cs.upc.edu/~oliveras/dpllt.pdf

The core algorithm efficiently decides equalities among variables and
(uninterpreted) function applications, using union-find data structures and
a pending queue to enforce function application congruence.

Supported input:
    - SymPy Symbol objects (as variables/constants)
    - EUFFunction objects (for uninterpreted functions of fixed arity)
    - Apply objects (custom, see below), for ground applications of EUFFunction to arguments.
    - Equations (sympy.Eq), e.g. `Eq(f(x), y)`

Key Classes and Structures:
---------------------------

- EUFFunction:  Uninterpreted function symbol (fixed arity).
- Apply:        Application of an EUFFunction to argument tuple.
- EUFCongruenceClosure: Implements the congruence closure engine.
    - Uses union-find to represent equivalence classes of terms.
    - Uses `lookup` and `use` for efficient congruence propagation.
    - Processes equations and propagates function congruence as per the literature.

Design:
-------

This implementation is suitable for integration with SymPy's SAT-based
assumptions and logic systems. It is ground (quantifier-free), works
directly with SymPy Symbol and Eq objects (no need for special encoding),
and supports efficient incremental equality propagation.

References:
-----------

- R. Nieuwenhuis and A. Oliveras, "Congruence Closure with Integer Offsets", 2003.
- Ganzinger, Nieuwenhuis, Oliveras, Tinelli, "DPLL(T): Fast Decision Procedures", 2004.
"""

# from sympy import Symbol, Eq, Dummy, Lambda, Tuple
# from sympy.core import Basic
# from collections import deque, defaultdict
# from sympy.utilities.iterables import numbered_symbols
# from sympy.utilities.iterables import iterable


# class EUFFunction:
#     """
#     Uninterpreted function symbol for use in EUF reasoning.

#     Parameters
#     ----------

#     name : str
#         The name of the uninterpreted function symbol.
#     arity : int
#         The number of arguments the function accepts.
#     """

#     def __init__(self, name: str, arity: int):
#         self.name = name
#         self.arity = arity

#     def __call__(self, *args):
#         assert len(args) == self.arity, (
#             f"{self.name} expects {self.arity} args, got {len(args)}"
#         )
#         return Apply(self, args)

#     def __repr__(self):
#         return self.name

#     def __eq__(self, other):
#         return (
#             isinstance(other, EUFFunction)
#             and self.name == other.name
#             and self.arity == other.arity
#         )

#     def __hash__(self):
#         return hash((self.name, self.arity))


# class Apply(Basic):
#     """
#     Application of an EUFFunction to arguments, suitable for congruence closure.

#     Parameters
#     ----------

#     function : EUFFunction
#         The uninterpreted function symbol.
#     args : tuple
#         The argument terms (Symbol or other Apply objects).
#     """

#     def __new__(cls, function, args):
#         args = tuple(args)
#         obj = Basic.__new__(cls, function, *args)
#         obj._function = function
#         obj._arguments = args
#         return obj

#     @property
#     def function(self):
#         """Return the EUFFunction applied."""
#         return self._function

#     @property
#     def arguments(self):
#         """Return the tuple of argument terms."""
#         return self._arguments

#     def __repr__(self):
#         return f"{self._function}({', '.join(map(str, self._arguments))})"


# class EUFCongruenceClosure:
#     """
#     Implements O(n log n) ground congruence-closure for EUF.

#     Parameters
#     ----------

#     eqs : list
#         A list of SymPy Eq objects or Apply terms encoding the constraints.

#     Attributes
#     ----------

#     _id_of : dict
#         Maps term objects to unique integer ids.
#     _term_of : list
#         Maps integer ids to term objects.
#     rep : list
#         Union-find representative ids for each class.
#     rank : list
#         Union-find rank structure.
#     cls : list
#         Lists of equivalence class members.

#     lookup : dict
#         Maps (function, tuple of argument reps) to result id (for congruence).
#     use : defaultdict
#         For each arg, tracks relevant function applications for congruence propagation.
#     _pending : deque
#         Pairs of class ids pending union.
#     """

#     def __init__(self, eqs):
#         self._id_of = {}          # term object -> id
#         self._term_of = []        # id -> term
#         self.rep = []
#         self.rank = []
#         self.cls = []

#         self.lookup = {}                  # (func, tuple[rep_ids]) -> id
#         self.use = defaultdict(list)      # arg_id -> applications
#         self._pending = deque()

#         def _cid(term):
#             if term not in self._id_of:
#                 tid = len(self._term_of)
#                 self._id_of[term] = tid
#                 self._term_of.append(term)
#                 self.rep.append(tid)
#                 self.rank.append(0)
#                 self.cls.append([tid])
#             return self._id_of[term]

#         for eq in eqs:
#             if isinstance(eq, Eq):
#                 lhs = self._flatten(eq.lhs, _cid)
#                 rhs = self._flatten(eq.rhs, _cid)
#                 self._pending.append((_cid(lhs), _cid(rhs)))
#             else:
#                 self._flatten(eq, _cid)

#         self._process_pending_unions()

#     def _flatten(self, term, _cid):
#         """
#         Recursively flattens terms (Symbol or Apply) for use in congruence closure.
#         Internally registers all terms and argument tuples.
#         """
#         if isinstance(term, (Symbol, int)):
#             _cid(term)
#             return term
#         elif isinstance(term, Apply):
#             flat_args = tuple(self._flatten(arg, _cid) for arg in term.arguments)
#             flat_term = Apply(term.function, flat_args)
#             tid = _cid(flat_term)
#             for arg in flat_args:
#                 aid = _cid(arg)
#                 self.use[aid].append((term.function, flat_args, tid))
#             self.lookup[(term.function, tuple(_cid(arg) for arg in flat_args))] = tid
#             return flat_term
#         else:
#             raise TypeError(f"Unsupported term in congruence closure: {term!r}")

#     def _find(self, i):
#         """
#         Find the representative id for a given id, with path compression.
#         """
#         while self.rep[i] != i:
#             self.rep[i] = self.rep[self.rep[i]]
#             i = self.rep[i]
#         return i

#     def _union(self, i, j):
#         """
#         Merge the equivalence classes of ids i and j.
#         """
#         i, j = self._find(i), self._find(j)
#         if i == j:
#             return
#         if self.rank[i] < self.rank[j]:
#             i, j = j, i
#         self.rep[j] = i
#         if self.rank[i] == self.rank[j]:
#             self.rank[i] += 1
#         self.cls[i].extend(self.cls[j])
#         self.cls[j] = []
#         self._merge_effects(i, j)

#     def _process_pending_unions(self):
#         """
#         Processes all pending unions, applying congruence closure to class merge queue.
#         """
#         while self._pending:
#             i, j = self._pending.popleft()
#             self._union(i, j)

#     def _merge_effects(self, ra, rb):
#         """
#         Propagate congruence after merging equivalence classes ra and rb.
#         Scans all function applications affected by merged arguments, and
#         schedules new unions if two app results (with congruent args) now exist.
#         """
#         for r in (ra, rb):
#             tmp = self.use.pop(r, [])
#             for (func, arg_terms, res_id) in tmp:
#                 rep_args = tuple(self._find(self._id_of[arg]) for arg in arg_terms)
#                 key = (func, rep_args)
#                 if key in self.lookup:
#                     self._pending.append((res_id, self.lookup[key]))
#                 else:
#                     self.lookup[key] = res_id
#                     for aid in rep_args:
#                         self.use[aid].append((func, rep_args, res_id))

#     def add_equality(self, a, b):
#         """
#         Add an equality a == b (incremental), given as Symbol or Apply.

#         Useful for interactive or stepwise usage.
#         """
#         if a not in self._id_of or b not in self._id_of:
#             raise KeyError("Term not registered.")
#         self._pending.append((self._find(self._id_of[a]), self._find(self._id_of[b])))
#         self._process_pending_unions()

#     def are_equal(self, a, b):
#         """
#         Return True iff a and b are currently known to be equal in the congruence closure.

#         If a or b are not present in the current congruence, returns False.
#         """
#         if a not in self._id_of or b not in self._id_of:
#             return False
#         return self._find(self._id_of[a]) == self._find(self._id_of[b])

from sympy import Symbol, Eq, Dummy, Lambda
from sympy.core import Basic
from sympy.core.numbers import Number
from sympy.utilities.iterables import numbered_symbols
from collections import defaultdict, deque

class EUFCongruenceClosure:
    """
    Congruence closure engine for EUF with currying and Dummy constants.

    Implements the algorithm from:
    'Congruence Closure with Integer Offsets' (Nieuwenhuis and Oliveras, 2003).

    Parameters
    ----------
    eqs : iterable of sympy.Eq or SymPy expressions
        Initial equalities or expressions to build the congruence closure.
    """

    def __init__(self, eqs):
        self._flat = {}  # Maps flattened terms to unique Dummy constants
        self._dummies = numbered_symbols('c', Dummy)
        self.rep = {}    # Union-find parent pointers: Dummy or Symbol -> Dummy or Symbol
        self.rank = {}   # Union-find ranks
        self.cls = defaultdict(set)  # Equivalence classes: leader -> set of members
        self.lookup = {} # Map from (func, tuple(args)) to Dummy
        self.use = defaultdict(list)  # For each Dummy/Symbol, list of function applications referencing it
        self.pending = deque()  # Queue of pairs of elements to merge

        # Caches for numbers and other atoms mapped to Dummies
        self._num_dummy_cache = {}
        self._atom_dummy_cache = {}

        # Register all symbols in union-find (as representatives of themselves)
        self._register_symbol_cache = {}

        # Process initial equalities and flatten terms
        for eq in eqs:
            if isinstance(eq, Eq):
                lhs = self._flatten(eq.lhs)
                rhs = self._flatten(eq.rhs)
                self._register_dummy_or_symbol(lhs)
                self._register_dummy_or_symbol(rhs)
                self.pending.append((lhs, rhs))
            else:
                val = self._flatten(eq)
                self._register_dummy_or_symbol(val)

        # Process all pending merges
        self._process_pending()

    def _register_dummy_or_symbol(self, d):
        """Register Dummy or Symbol in union-find data structures if not already."""
        if isinstance(d, Symbol):
            if d not in self._register_symbol_cache:
                self.rep[d] = d
                self.rank[d] = 0
                self.cls[d].add(d)
                self._register_symbol_cache[d] = True
        else:
            # Dummy or anything else expected here
            if d not in self.rep:
                self.rep[d] = d
                self.rank[d] = 0
                self.cls[d].add(d)

    def _flatten(self, expr):
        """
        Recursively flatten the expression to unique Dummies or Symbols.

        - Symbols returned as is, registered in union-find.
        - Numbers replaced by unique Dummies.
        - Lambda expressions are fully curried, body flattened,
          and a unique Dummy assigned to the Lambda itself.
        - Function applications replaced by unique Dummies with congruence info.
        """
        # Already a Dummy (unique constant)
        if isinstance(expr, Dummy):
            return expr

        # Symbol: returned as is and registered
        if isinstance(expr, Symbol):
            self._register_dummy_or_symbol(expr)
            return expr

        # Numbers (int, float, Rational)
        if isinstance(expr, Number):
            if expr not in self._num_dummy_cache:
                d = next(self._dummies)
                self._num_dummy_cache[expr] = d
                self._register_dummy_or_symbol(d)
            return self._num_dummy_cache[expr]

        # Other atoms (e.g., BooleanTrue etc.) convert to Dummy
        if hasattr(expr, "is_Atom") and expr.is_Atom:
            if expr not in self._atom_dummy_cache:
                d = next(self._dummies)
                self._atom_dummy_cache[expr] = d
                self._register_dummy_or_symbol(d)
            return self._atom_dummy_cache[expr]

        # Handle Lambda: curry fully, flatten body, assign unique Dummy to Lambda
        if isinstance(expr, Lambda):
            lam = expr
            if len(lam.variables) != 1:
                lam = lam.curry()
            body_flat = self._flatten(lam.expr)
            lam_curried = Lambda(lam.variables[0], body_flat)

            # Assign unique Dummy for this Lambda expression
            if lam_curried not in self._flat:
                d = next(self._dummies)
                self._flat[lam_curried] = d
                self._register_dummy_or_symbol(d)
            return self._flat[lam_curried]

        # Function application and other expressions
        # Flatten function part
        func = expr.func
        if isinstance(func, Basic) and not isinstance(expr, Lambda):
            func_flat = self._flatten(func)
        else:
            func_flat = func

        args_flat = tuple(self._flatten(arg) for arg in expr.args)

        key = (func_flat, args_flat)

        if key not in self._flat:
            d = next(self._dummies)
            self._flat[key] = d
            self._register_dummy_or_symbol(d)
            self.lookup[key] = d
            for arg in args_flat:
                self.use[arg].append((func_flat, args_flat, d))

        return self._flat[key]

    def _find(self, d):
        """Union-Find find with path compression."""
        if d not in self.rep:
            raise KeyError(f"Element {d} not registered in union-find structure.")
        if self.rep[d] != d:
            self.rep[d] = self._find(self.rep[d])
        return self.rep[d]

    def _union(self, a, b):
        """Union-Find union with rank heuristic and congruence propagation."""
        ra = self._find(a)
        rb = self._find(b)
        if ra == rb:
            return

        # Union by rank
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra

        self.rep[rb] = ra
        self.cls[ra].update(self.cls[rb])
        del self.cls[rb]

        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

        self._propagate(ra, rb)

    def _propagate(self, ra, rb):
        """Propagate congruence closure after union."""
        for leader in (ra, rb):
            uses = self.use.pop(leader, [])
            for (func, argdummies, res) in uses:
                rep_args = tuple(self._find(arg) for arg in argdummies)
                key = (func, rep_args)
                if key in self.lookup:
                    other = self.lookup[key]
                    self.pending.append((res, other))
                else:
                    self.lookup[key] = res
                    for arg in rep_args:
                        self.use[arg].append((func, rep_args, res))

    def _process_pending(self):
        """Process all pending merges until fixpoint."""
        while self.pending:
            a, b = self.pending.popleft()
            self._union(a, b)

    def add_equality(self, a, b):
        """
        Add an equality between expressions a and b incrementally.

        Flatten internally.
        """
        a_flat = self._flatten(a)
        b_flat = self._flatten(b)
        self._register_dummy_or_symbol(a_flat)
        self._register_dummy_or_symbol(b_flat)
        self.pending.append((a_flat, b_flat))
        self._process_pending()

    def are_equal(self, a, b):
        """
        Check if expressions a and b are equal under current congruence closure.

        Returns True or False.
        """
        a_flat = self._flatten(a)
        b_flat = self._flatten(b)
        if a_flat not in self.rep or b_flat not in self.rep:
            return False
        return self._find(a_flat) == self._find(b_flat)
