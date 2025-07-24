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

from sympy import Symbol, Eq
from sympy.core import Basic
from collections import deque, defaultdict


class EUFFunction:
    """
    Uninterpreted function symbol for use in EUF reasoning.

    Parameters
    ----------

    name : str
        The name of the uninterpreted function symbol.
    arity : int
        The number of arguments the function accepts.
    """

    def __init__(self, name: str, arity: int):
        self.name = name
        self.arity = arity

    def __call__(self, *args):
        assert len(args) == self.arity, (
            f"{self.name} expects {self.arity} args, got {len(args)}"
        )
        return Apply(self, args)

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        return (
            isinstance(other, EUFFunction)
            and self.name == other.name
            and self.arity == other.arity
        )

    def __hash__(self):
        return hash((self.name, self.arity))


class Apply(Basic):
    """
    Application of an EUFFunction to arguments, suitable for congruence closure.

    Parameters
    ----------

    function : EUFFunction
        The uninterpreted function symbol.
    args : tuple
        The argument terms (Symbol or other Apply objects).
    """

    def __new__(cls, function, args):
        args = tuple(args)
        obj = Basic.__new__(cls, function, *args)
        obj._function = function
        obj._arguments = args
        return obj

    @property
    def function(self):
        """Return the EUFFunction applied."""
        return self._function

    @property
    def arguments(self):
        """Return the tuple of argument terms."""
        return self._arguments

    def __repr__(self):
        return f"{self._function}({', '.join(map(str, self._arguments))})"


class EUFCongruenceClosure:
    """
    Implements O(n log n) ground congruence-closure for EUF.

    Parameters
    ----------

    eqs : list
        A list of SymPy Eq objects or Apply terms encoding the constraints.

    Attributes
    ----------

    _id_of : dict
        Maps term objects to unique integer ids.
    _term_of : list
        Maps integer ids to term objects.
    rep : list
        Union-find representative ids for each class.
    rank : list
        Union-find rank structure.
    cls : list
        Lists of equivalence class members.

    lookup : dict
        Maps (function, tuple of argument reps) to result id (for congruence).
    use : defaultdict
        For each arg, tracks relevant function applications for congruence propagation.
    _pending : deque
        Pairs of class ids pending union.
    """

    def __init__(self, eqs):
        self._id_of = {}          # term object -> id
        self._term_of = []        # id -> term
        self.rep = []
        self.rank = []
        self.cls = []

        self.lookup = {}                  # (func, tuple[rep_ids]) -> id
        self.use = defaultdict(list)      # arg_id -> applications
        self._pending = deque()

        def _cid(term):
            if term not in self._id_of:
                tid = len(self._term_of)
                self._id_of[term] = tid
                self._term_of.append(term)
                self.rep.append(tid)
                self.rank.append(0)
                self.cls.append([tid])
            return self._id_of[term]

        for eq in eqs:
            if isinstance(eq, Eq):
                lhs = self._flatten(eq.lhs, _cid)
                rhs = self._flatten(eq.rhs, _cid)
                self._pending.append((_cid(lhs), _cid(rhs)))
            else:
                self._flatten(eq, _cid)

        self._process_pending_unions()

    def _flatten(self, term, _cid):
        """
        Recursively flattens terms (Symbol or Apply) for use in congruence closure.
        Internally registers all terms and argument tuples.
        """
        if isinstance(term, (Symbol, int)):
            _cid(term)
            return term
        elif isinstance(term, Apply):
            flat_args = tuple(self._flatten(arg, _cid) for arg in term.arguments)
            flat_term = Apply(term.function, flat_args)
            tid = _cid(flat_term)
            for arg in flat_args:
                aid = _cid(arg)
                self.use[aid].append((term.function, flat_args, tid))
            self.lookup[(term.function, tuple(_cid(arg) for arg in flat_args))] = tid
            return flat_term
        else:
            raise TypeError(f"Unsupported term in congruence closure: {term!r}")

    def _find(self, i):
        """
        Find the representative id for a given id, with path compression.
        """
        while self.rep[i] != i:
            self.rep[i] = self.rep[self.rep[i]]
            i = self.rep[i]
        return i

    def _union(self, i, j):
        """
        Merge the equivalence classes of ids i and j.
        """
        i, j = self._find(i), self._find(j)
        if i == j:
            return
        if self.rank[i] < self.rank[j]:
            i, j = j, i
        self.rep[j] = i
        if self.rank[i] == self.rank[j]:
            self.rank[i] += 1
        self.cls[i].extend(self.cls[j])
        self.cls[j] = []
        self._merge_effects(i, j)

    def _process_pending_unions(self):
        """
        Processes all pending unions, applying congruence closure to class merge queue.
        """
        while self._pending:
            i, j = self._pending.popleft()
            self._union(i, j)

    def _merge_effects(self, ra, rb):
        """
        Propagate congruence after merging equivalence classes ra and rb.
        Scans all function applications affected by merged arguments, and
        schedules new unions if two app results (with congruent args) now exist.
        """
        for r in (ra, rb):
            tmp = self.use.pop(r, [])
            for (func, arg_terms, res_id) in tmp:
                rep_args = tuple(self._find(self._id_of[arg]) for arg in arg_terms)
                key = (func, rep_args)
                if key in self.lookup:
                    self._pending.append((res_id, self.lookup[key]))
                else:
                    self.lookup[key] = res_id
                    for aid in rep_args:
                        self.use[aid].append((func, rep_args, res_id))

    def add_equality(self, a, b):
        """
        Add an equality a == b (incremental), given as Symbol or Apply.

        Useful for interactive or stepwise usage.
        """
        if a not in self._id_of or b not in self._id_of:
            raise KeyError("Term not registered.")
        self._pending.append((self._find(self._id_of[a]), self._find(self._id_of[b])))
        self._process_pending_unions()

    def are_equal(self, a, b):
        """
        Return True iff a and b are currently known to be equal in the congruence closure.

        If a or b are not present in the current congruence, returns False.
        """
        if a not in self._id_of or b not in self._id_of:
            return False
        return self._find(self._id_of[a]) == self._find(self._id_of[b])
