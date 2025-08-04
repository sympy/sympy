import sys
from collections import defaultdict, deque
import sympy as sp
from sympy import Symbol, Eq, Lambda
from sympy.core.numbers import Number
from sympy.core import Basic
from sympy.utilities.iterables import numbered_symbols

# Allow very deep terms
sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))

class EUFCongruenceClosure:
    """
    Congruence closure engine for ground EUF (equality with uninterpreted functions).
    Implements the algorithm from Nieuwenhuis and Oliveras (2003).
    """

    def __init__(self, eqs):
        self._flat = {}               # Maps (func, args) or Lambda to Dummy
        self.lookup = {}              # Congruence lookup for function applications
        self.use = defaultdict(list) # Maps argument to list of applications using it

        self._dummies = numbered_symbols('c', sp.Dummy)
        self.rep = {}                # Union-find parent pointers
        self.cls = defaultdict(set)  # Maps representatives to class members

        self._num_dummy_cache = {}
        self._atom_dummy_cache = {}
        self._register_symbol_cache = {}

        self.pending = deque()

        # Flatten all initial equalities or terms and queue unions
        for eq in eqs:
            if isinstance(eq, Eq):
                lhs = self._flatten(eq.lhs)
                rhs = self._flatten(eq.rhs)
                self._register(lhs)
                self._register(rhs)
                self.pending.append((lhs, rhs))
            else:
                val = self._flatten(eq)
                self._register(val)

        self._process_pending()

    def _register(self, t):
        if t not in self.rep:
            self.rep[t] = t
            self.cls[t].add(t)

    def _new_dummy(self):
        d = next(self._dummies)
        self._register(d)
        return d

    def _flatten(self, expr):
        """
        Recursively flatten expr into canonical Dummy or Symbol.
        Handles Lambdas with currying and caches flattened terms.
        """
        if isinstance(expr, (sp.Dummy, Symbol)):
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

        if isinstance(expr, Lambda):
            lam = expr if len(expr.variables) == 1 else expr.curry()
            body_flat = self._flatten(lam.expr)
            lam_curried = Lambda(lam.variables[0], body_flat)
            if lam_curried not in self._flat:
                d = self._new_dummy()
                self._flat[lam_curried] = d
            return self._flat[lam_curried]

        # Flatten function application
        func = expr.func
        if isinstance(func, Basic) and not isinstance(expr, Lambda):
            func_flat = self._flatten(func)
        else:
            func_flat = func

        args_flat = tuple(self._find(self._flatten(arg)) for arg in expr.args)
        key = (func_flat, args_flat)

        if key not in self._flat:
            d = self._new_dummy()
            self._flat[key] = d
            self.lookup[key] = d
            for arg in args_flat:
                self.use[arg].append((func_flat, args_flat, d))
        return self._flat[key]

    def _find(self, t):
        # Register new elements before find to avoid KeyError
        if t not in self.rep:
            self._register(t)
        while t != self.rep[t]:
            self.rep[t] = self.rep[self.rep[t]]
            t = self.rep[t]
        return t

    def _union(self, a, b):
        ra, rb = self._find(a), self._find(b)
        if ra == rb:
            return
        # Union by size heuristic
        if len(self.cls[ra]) > len(self.cls[rb]):
            ra, rb = rb, ra

        for member in list(self.cls[ra]):
            self.rep[member] = rb
            self.cls[rb].add(member)
            self.use[rb].extend(self.use.pop(member, []))
        del self.cls[ra]

        # Propagate congruence through use list of rb
        for func, arg_ids, res in self.use[rb]:
            rep_args = tuple(self._find(arg) for arg in arg_ids)
            rep_res = self._find(res)
            key = (func, rep_args)
            if key in self.lookup:
                other = self._find(self.lookup[key])
                if other != rep_res:
                    self.pending.append((rep_res, other))
            self.lookup[key] = rep_res
            self._flat[key] = rep_res

    def _process_pending(self):
        while self.pending:
            self._union(*self.pending.popleft())

    def add_equality(self, lhs, rhs):
        self.pending.append((self._flatten(lhs), self._flatten(rhs)))
        self._process_pending()

    def are_equal(self, lhs, rhs):
        lid, rid = self._flatten(lhs), self._flatten(rhs)
        return self._find(lid) == self._find(rid)
