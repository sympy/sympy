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

from sympy import Symbol, Eq, Dummy, Lambda
from sympy.core import Basic
from sympy.core.numbers import Number
from sympy.utilities.iterables import numbered_symbols
from collections import defaultdict, deque


class EUFCongruenceClosure:
    """
    Congruence closure engine for Equality with Uninterpreted Functions (EUF),
    using currying and Dummy symbols for flattening.

    This implementation follows the algorithm described in the paper
    "Congruence Closure with Integer Offsets" (Nieuwenhuis and Oliveras, 2003).

    It maintains equivalence classes of expressions using a union-find
    data structure, supports incremental addition of equalities,
    and enforces congruence closure by propagating equalities through
    function applications.

    Parameters
    ----------

    eqs : iterable of sympy.Eq or SymPy expressions
        An iterable of initial equalities or expressions to initialize the
        congruence closure engine. These are used to build equivalence classes
        and propagate congruences on construction.

    Attributes
    ----------

    _flat : dict
        Maps flattened expressions (tuples of function and argument Dummies/Symbols)
        or curried Lambdas to unique Dummy symbols representing those terms.
    _dummies : generator
        Generator yielding fresh Dummy symbols used to represent flattened terms.
    rep : dict
        Union-find parent pointers mapping terms (Symbols or Dummies) to their current
        class representative.
    cls : collections.defaultdict[set]
        Maps each representative to the set of all terms in its equivalence class.
    lookup : dict
        Maps congruence keys of the form (function, tuple of argument representatives)
        to the Dummy symbol representing the result of the function application.
    use : collections.defaultdict[list]
        Maps each Dummy/Symbol (argument representative) to the list of function applications
        (tuples of function, args, and application Dummy) that use it.
    pending : collections.deque
        Queue of pairs of terms (Dummies/Symbols) pending union, to be processed iteratively.

    Methods
    -------

    add_equality(a, b)
        Incrementally add an equality between two SymPy expressions `a` and `b`,
        updating the congruence classes accordingly.
    are_equal(a, b)
        Check if two SymPy expressions `a` and `b` are equivalent under the current
        congruence closure.
    """

    def __init__(self, eqs):
        self._flat = {}
        self._dummies = numbered_symbols('c', Dummy)
        self.rep = {}                # term -> current representative (as key)
        self.cls = defaultdict(set)  # representative -> set of class members
        self.lookup = {}             # (function, args) -> result term
        self.use = defaultdict(list) # arg -> [(function, args, result term), ...]
        self.pending = deque()

        self._num_dummy_cache = {}
        self._atom_dummy_cache = {}
        self._register_symbol_cache = {}

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

        self._process_pending()

    def _register_dummy_or_symbol(self, d):
        """
        Register a Dummy or Symbol `d` in the union-find data structures,
        initializing its equivalence class if not already present.

        Parameters
        ----------
        d : sympy.Symbol or sympy.Dummy
            The term to register.
        """
        if isinstance(d, Symbol):
            if d not in self._register_symbol_cache:
                self.rep[d] = d
                self.cls[d].add(d)
                self._register_symbol_cache[d] = True
        else:
            if d not in self.rep:
                self.rep[d] = d
                self.cls[d].add(d)

    def _flatten(self, expr):
        """
        Recursively flatten an expression `expr` into a unique Dummy symbol or
        Symbol representing that expression in the congruence closure.

        The flattening performs:
          - Atoms and Symbols remain unchanged or registered directly.
          - Numbers are replaced by unique Dummy symbols for canonical representation.
          - Other atomic expressions (like BooleanTrue) are assigned unique Dummies.
          - Lambda expressions are curried, their bodies flattened, and assigned unique Dummies.
          - Function applications are curried and replaced by unique Dummies, with congruence info maintained.

        Parameters
        ----------
        expr : sympy.core.expr.Expr
            The input SymPy expression to flatten.

        Returns
        -------
        sympy.Symbol or sympy.Dummy
            The canonical representative Dummy or Symbol corresponding to the expression.
        """
        if isinstance(expr, Dummy):
            return expr
        if isinstance(expr, Symbol):
            self._register_dummy_or_symbol(expr)
            return expr
        if isinstance(expr, Number):
            if expr not in self._num_dummy_cache:
                d = next(self._dummies)
                self._num_dummy_cache[expr] = d
                self._register_dummy_or_symbol(d)
            return self._num_dummy_cache[expr]
        if hasattr(expr, "is_Atom") and expr.is_Atom:
            if expr not in self._atom_dummy_cache:
                d = next(self._dummies)
                self._atom_dummy_cache[expr] = d
                self._register_dummy_or_symbol(d)
            return self._atom_dummy_cache[expr]
        if isinstance(expr, Lambda):
            lam = expr
            if len(lam.variables) != 1:
                lam = lam.curry()
            body_flat = self._flatten(lam.expr)
            lam_curried = Lambda(lam.variables[0], body_flat)
            if lam_curried not in self._flat:
                d = next(self._dummies)
                self._flat[lam_curried] = d
                self._register_dummy_or_symbol(d)
            return self._flat[lam_curried]
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
        """
        Find the representative of the equivalence class containing `d`
        using path compression.

        If `d` is not registered, it is initialized as a singleton class.

        Parameters
        ----------
        d : sympy.Symbol or sympy.Dummy
            The term to find the representative for.

        Returns
        -------
        sympy.Symbol or sympy.Dummy
            The canonical representative of the class containing `d`.
        """
        if d not in self.rep:
            self.rep[d] = d
            self.cls[d].add(d)
            if isinstance(d, Symbol):
                self._register_symbol_cache[d] = True
        while d != self.rep[d]:
            self.rep[d] = self.rep[self.rep[d]]
            d = self.rep[d]
        return d

    def _union(self, a, b):
        """
        Merge the equivalence classes of `a` and `b` ensuring the smaller
        class is merged into the larger one for efficiency.

        This updates representatives, merges equivalence classes, merges
        use lists, and propagates congruences by enqueuing new equalities into pending.

        Parameters
        ----------
        a, b : sympy.Symbol or sympy.Dummy
            The terms whose classes are to be merged.
        """
        ra = self._find(a)
        rb = self._find(b)
        if ra == rb:
            return
        if len(self.cls[ra]) > len(self.cls[rb]):
            ra, rb = rb, ra

        # Move all elements and uses from ra into rb
        for c in list(self.cls[ra]):
            self.rep[c] = rb
            self.cls[rb].add(c)
            # Merge the use list for each element c
            uses = self.use.pop(c, [])
            if uses:
                self.use[rb].extend(uses)
        del self.cls[ra]

        # Now process congruence for all applications in use[rb]
        for (func, argdummies, res) in self.use[rb]:
            rep_args = tuple(self._find(arg) for arg in argdummies)
            res_rep = self._find(res)
            key = (func, rep_args)
            if key in self.lookup:
                other = self.lookup[key]
                other_rep = self._find(other)
                if res_rep != other_rep:
                    self.pending.append((res_rep, other_rep))
            self.lookup[key] = res_rep

    def _process_pending(self):
        """
        Process all pending equality merges in a loop until no further
        merges remain.

        This ensures the congruence closure reaches a fixpoint.
        """
        while self.pending:
            a, b = self.pending.popleft()
            self._union(a, b)

    def add_equality(self, a, b):
        """
        Incrementally add equality between two SymPy expressions `a` and `b`,
        updating internal congruence classes.

        Parameters
        ----------
        a, b : sympy.core.expr.Expr
            Expressions declared equal.
        """
        a_flat = self._flatten(a)
        b_flat = self._flatten(b)
        self._register_dummy_or_symbol(a_flat)
        self._register_dummy_or_symbol(b_flat)
        self.pending.append((a_flat, b_flat))
        self._process_pending()

    def are_equal(self, a, b):
        """
        Check whether two expressions `a` and `b` are equal under the current
        congruence closure state.

        Parameters
        ----------
        a, b : sympy.core.expr.Expr
            Expressions to compare.

        Returns
        -------
        bool
            True if `a` and `b` are in the same equivalence class, False otherwise.
        """
        a_flat = self._flatten(a)
        b_flat = self._flatten(b)
        if a_flat not in self.rep or b_flat not in self.rep:
            return False
        return self._find(a_flat) == self._find(b_flat)
