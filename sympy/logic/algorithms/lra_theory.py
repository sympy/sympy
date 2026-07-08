"""Implements "A Fast Linear-Arithmetic Solver for DPLL(T)"

The LRASolver class defined in this file can be used
in conjunction with a SAT solver to check the
satisfiability of formulas involving inequalities.

Here's an example of how that would work:

    Suppose you want to check the satisfiability of
    the following formula:

    >>> from sympy.core.relational import Eq
    >>> from sympy.abc import x, y
    >>> f = ((x > 0) | (x < 0)) & (Eq(x, 0) | Eq(y, 1)) & (~Eq(y, 1) | Eq(1, 2))

    First a preprocessing step should be done on f. During preprocessing,
    f should be checked for any predicates such as `Q.prime` that can't be
    handled. Also unequality like `~Eq(y, 1)` should be split.

    I should mention that the paper says to split both equalities and
    unequality, but this implementation only requires that unequality
    be split.

    >>> f = ((x > 0) | (x < 0)) & (Eq(x, 0) | Eq(y, 1)) & ((y < 1) | (y > 1) | Eq(1, 2))

    Then an LRASolver instance needs to be initialized with this formula.

    >>> from sympy.assumptions.cnf import CNF, EncodedCNF
    >>> from sympy.assumptions.ask import Q
    >>> from sympy.logic.algorithms.lra_theory import LRASolver
    >>> cnf = CNF.from_prop(f)
    >>> enc = EncodedCNF()
    >>> enc.add_from_cnf(cnf)
    >>> lra, conflicts = LRASolver.from_encoded_cnf(enc)

    Any immediate one-lital conflicts clauses will be detected here.
    In this example, `~Eq(1, 2)` is one such conflict clause. We'll
    want to add it to `f` so that the SAT solver is forced to
    assign Eq(1, 2) to False.

    >>> f = f & ~Eq(1, 2)

    Now that the one-literal conflict clauses have been added
    and an lra object has been initialized, we can pass `f`
    to a SAT solver. The SAT solver will give us a satisfying
    assignment such as:

    (1 = 2): False
    (y = 1): True
    (y < 1): True
    (y > 1): True
    (x = 0): True
    (x < 0): True
    (x > 0): True

    Next you would pass this assignment to the LRASolver
    which will be able to determine that this particular
    assignment is satisfiable or not.

    Note that since EncodedCNF is inherently non-deterministic,
    the int each predicate is encoded as is not consistent. As a
    result, the code below likely does not reflect the assignment
    given above.

    >>> lra.assert_lit(-1) #doctest: +SKIP
    >>> lra.assert_lit(2) #doctest: +SKIP
    >>> lra.assert_lit(3) #doctest: +SKIP
    >>> lra.assert_lit(4) #doctest: +SKIP
    >>> lra.assert_lit(5) #doctest: +SKIP
    >>> lra.assert_lit(6) #doctest: +SKIP
    >>> lra.assert_lit(7) #doctest: +SKIP
    >>> is_sat, conflict_or_assignment = lra.check()

    As the particular assignment suggested is not satisfiable,
    the LRASolver will return unsat and a conflict clause when
    given that assignment. The conflict clause will always be
    minimal, but there can be multiple minimal conflict clauses.
    One possible conflict clause could be `~(x < 0) | ~(x > 0)`.

    We would then add whatever conflict clause is given to
    `f` to prevent the SAT solver from coming up with an
    assignment with the same conflicting literals. In this case,
    the conflict clause `~(x < 0) | ~(x > 0)` would prevent
    any assignment where both (x < 0) and (x > 0) were both
    true.

    The SAT solver would then find another assignment
    and we would check that assignment with the LRASolver
    and so on. Eventually either a satisfying assignment
    that the SAT solver and LRASolver agreed on would be found
    or enough conflict clauses would be added so that the
    boolean formula was unsatisfiable.


This implementation is based on [1]_, which includes a
detailed explanation of the algorithm and pseudocode
for the most important functions.

[1]_ also explains how backtracking and theory propagation
could be implemented to speed up the current implementation,
but these are not currently implemented.

TODO:
 - Handle non-rational real numbers
 - Handle positive and negative infinity
 - Implement backtracking and theory propagation

References
==========

.. [1] Dutertre, B., de Moura, L.:
       A Fast Linear-Arithmetic Solver for DPLL(T)
       https://link.springer.com/chapter/10.1007/11817963_11
"""
from __future__ import annotations
from sympy.solvers.solveset import linear_eq_to_matrix
from sympy.matrices.dense import eye
from sympy.assumptions import Predicate
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.ask import Q
from sympy.core import Dummy
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.relational import Eq, Ge, Gt, Le, Lt
from sympy.core.sympify import sympify
from sympy.core.singleton import S
from sympy.core.numbers import Rational, oo
from sympy.matrices.dense import Matrix
from sympy.utilities.iterables import sift
import math


class UnhandledInput(Exception):
    """
    Raised while creating an LRASolver if non-linearity
    or non-rational numbers are present.
    """

# predicates that LRASolver understands and makes use of
ALLOWED_PRED = {Q.eq: Eq, Q.gt: Gt, Q.lt: Lt, Q.le: Le, Q.ge: Ge}

# if true ~Q.gt(x, y) implies Q.le(x, y)
HANDLE_NEGATION = True

class LRASolver():
    """
    Linear Arithmetic Solver for DPLL(T) implemented with an algorithm based on
    the Dual Simplex method. Uses Bland's pivoting rule to avoid cycling.

    References
    ==========

    .. [1] Dutertre, B., de Moura, L.:
           A Fast Linear-Arithmetic Solver for DPLL(T)
           https://link.springer.com/chapter/10.1007/11817963_11
    """

    def __init__(self, A, slack_variables, nonslack_variables,
                 atom_id_to_boundaries, s_subs, testing_mode):
        """
        Use the "from_encoded_cnf" method to create a new LRASolver.
        """
        self.run_checks = testing_mode
        self.s_subs = s_subs  # used only for test_lra_theory.test_random_problems

        if any(not isinstance(a, Rational) for a in A):
            raise UnhandledInput("Non-rational numbers are not handled")
        if not all(isinstance(b.bound, Rational)
               for bs in atom_id_to_boundaries.values() for b in bs):
            raise UnhandledInput("Non-rational numbers are not handled")
        m, n = len(slack_variables), len(slack_variables)+len(nonslack_variables)
        if m != 0:
            assert A.shape == (m, n)
        if self.run_checks:
            assert A[:, n-m:] == -eye(m)

        self.atom_id_to_boundaries = atom_id_to_boundaries
        self.A = A
        self.slack = slack_variables
        self.nonslack = nonslack_variables
        self.all_var = nonslack_variables + slack_variables

        self.slack_set = set(slack_variables)

        self.is_sat = True  # While True, all constraints asserted so far are satisfiable
        self.result = None  # always one of: (True, assignment), (False, conflict clause), None

        self.bound_history = []
        self.last_assign_snapshot = {var: var.assign for var in self.all_var}

    @staticmethod
    def from_encoded_cnf(encoded_cnf, testing_mode=False):
        """
        Creates an LRASolver from an EncodedCNF object
        and a list of conflict clauses for propositions
        that can be simplified to True or False.

        Parameters
        ==========

        encoded_cnf : EncodedCNF

        testing_mode : bool
            Setting testing_mode to True enables some slow assert statements
            and sorting to reduce nonterministic behavior.

        Returns
        =======

        (lra, conflicts)

        lra : LRASolver

        conflicts : list
            Contains a one-literal conflict clause for each proposition
            that can be simplified to True or False.

        Example
        =======

        >>> from sympy.core.relational import Eq
        >>> from sympy.assumptions.cnf import CNF, EncodedCNF
        >>> from sympy.assumptions.ask import Q
        >>> from sympy.logic.algorithms.lra_theory import LRASolver
        >>> from sympy.abc import x, y, z
        >>> phi = (x >= 0) & ((x + y <= 2) | (x + 2 * y - z >= 6))
        >>> phi = phi & (Eq(x + y, 2) | (x + 2 * y - z > 4))
        >>> phi = phi & Q.gt(2, 1)
        >>> cnf = CNF.from_prop(phi)
        >>> enc = EncodedCNF()
        >>> enc.from_cnf(cnf)
        >>> lra, conflicts = LRASolver.from_encoded_cnf(enc, testing_mode=True)
        >>> lra #doctest: +SKIP
        <sympy.logic.algorithms.lra_theory.LRASolver object at 0x7fdcb0e15b70>
        >>> conflicts #doctest: +SKIP
        [[4]]
        """
        # This function has three main jobs:
        # - raise errors if the input formula is not handled
        # - preprocesses the formula into a matrix and single variable constraints
        # - create one-literal conflict clauses from predicates that are always True
        #   or always False such as Q.gt(3, 2)
        #
        # See the preprocessing section of "A Fast Linear-Arithmetic Solver for DPLL(T)"
        # for an explanation of how the formula is converted into a matrix
        # and a set of single variable constraints.

        atom_id_to_boundaries = {}
        A = []

        basic = []
        s_count = 0
        s_subs = {}
        nonbasic = []
        atom_vars = set()

        if testing_mode:
            # sort to reduce nondeterminism
            encoded_cnf_items = sorted(encoded_cnf.encoding.items(),
                                       key=lambda x: str(x))
        else:
            encoded_cnf_items = encoded_cnf.encoding.items()

        empty_var = Dummy()
        var_to_lra_var = {}
        conflicts = []

        for prop, atom_id in encoded_cnf_items:
            if isinstance(prop, Predicate):
                prop = prop(empty_var)
            if not isinstance(prop, AppliedPredicate):
                if prop == True:
                    conflicts.append([atom_id])
                    continue
                if prop == False:
                    conflicts.append([-atom_id])
                    continue

                raise ValueError(f"Unhandled Predicate: {prop}")

            assert prop.function in ALLOWED_PRED
            if prop.lhs == S.NaN or prop.rhs == S.NaN:
                raise ValueError(f"{prop} contains nan")
            if prop.lhs.is_imaginary or prop.rhs.is_imaginary:
                raise UnhandledInput(f"{prop} contains an imaginary component")
            if prop.lhs == oo or prop.rhs == oo:
                raise UnhandledInput(f"{prop} contains infinity")

            expr = prop.lhs - prop.rhs
            pred = ALLOWED_PRED[prop.function](expr, S.Zero)
            if pred == True:
                conflicts.append([atom_id])
                continue
            if pred == False:
                conflicts.append([-atom_id])
                continue
            if not expr.free_symbols:
                raise UnhandledInput(f"{prop} could not be simplified")

            if prop.function in [Q.ge, Q.gt]:
                expr = -expr

            # Example: 2x + 3y, 2 <- _sep_const_terms(2x + 3y + 2)
            vars, const = _sep_const_terms(expr)
            # Examples:
            # x, 2 <- _sep_const_coeff(2x)
            # 2x + 3y, 1 <- _sep_const_coeff(2x + 3y + 2)
            vars, var_coeff = _sep_const_coeff(vars)
            const = const / var_coeff
            # Example: [2x, 3y] <- Add.make_args(2x + 3y)
            terms = Add.make_args(vars)
            for term in terms:
                term, _ = _sep_const_coeff(term)
                assert len(term.free_symbols) > 0
                if term not in var_to_lra_var:
                    var_to_lra_var[term] = LRAVariable(term)
                    nonbasic.append(term)

            if len(terms) > 1:
                if vars not in s_subs:
                    s_count += 1
                    d = Dummy(f"s{s_count}")
                    var_to_lra_var[d] = LRAVariable(d)
                    basic.append(d)
                    s_subs[vars] = d
                    A.append(vars - d)
                var = s_subs[vars]
            else:
                var = terms[0]

            atom_vars.add(var)

            assert var_coeff != 0

            equality = prop.function == Q.eq
            strict = prop.function in [Q.gt, Q.lt]
            if equality:
                b1 = Boundary(var_to_lra_var[var], -const, True, False)  # x <= c
                b2 = Boundary(var_to_lra_var[var], -const, False, False) # x >= c
                atom_id_to_boundaries[atom_id] = [b1, b2]
            else:
                upper = var_coeff > 0
                b = Boundary(var_to_lra_var[var], -const, upper, strict)
                atom_id_to_boundaries[atom_id] = [b]

        fs = [v.free_symbols for v in nonbasic + basic]
        assert all(len(syms) > 0 for syms in fs)
        fs_count = sum(len(syms) for syms in fs)
        if len(fs) > 0 and  len(set.union(*fs)) < fs_count:
            raise UnhandledInput("Nonlinearity is not handled")

        A, _ = linear_eq_to_matrix(A, nonbasic + basic)
        # matrix A is guaranteed to able to be simplified
        # by removing the non-basic (e.g original) non-atom variables from it
        # these removed variables will be replaced by linear equation of existing variables.
        nonatom_vars = {i for i in nonbasic if i not in atom_vars}
        A, basic, nonbasic = _reduce_matrix(A, basic, nonbasic, nonatom_vars)
        nonbasic = [var_to_lra_var[nb] for nb in nonbasic]
        basic = [var_to_lra_var[b] for b in basic]
        for idx, var in enumerate(nonbasic + basic):
            var.col_idx = idx

        solver = LRASolver(A, basic, nonbasic, atom_id_to_boundaries,
                           s_subs, testing_mode)
        return solver, conflicts

    def reset_bounds(self):
        """
        Resets the state of the LRASolver to before
        anything was asserted.
        """
        self.result = None
        for var in self.all_var:
            var.initialize()

    def assert_lit(self, literal):
        """
        Assert a literal representing a constraint
        and update the internal state accordingly.

        Note that due to peculiarities of this implementation
        asserting ~(x > 0) will assert (x <= 0) but asserting
        ~Eq(x, 0) will not do anything.

        Parameters
        ==========

        literal : int
            A mapping of IDs to constraints
            can be found in `self.atom_id_to_boundaries`.

        Returns
        =======

        None or (False, explanation)

        explanation : set of ints
            A conflict clause that "explains" why
            the literals asserted so far are unsatisfiable.
        """
        if abs(literal) not in self.atom_id_to_boundaries:
            return None

        if not HANDLE_NEGATION and literal < 0:
            return None

        boundaries = self.atom_id_to_boundaries[abs(literal)]
        is_literal_negated = literal < 0

        if len(boundaries) > 1 and is_literal_negated:
            # Negated equality is not handled and should only appear in
            # conflict clauses.
            return None

        res = None
        for boundary in boundaries:
            res = self._assert_bound(boundary, literal)
            if res and res[0] is False:
                break

        if self.is_sat and all(b.var not in self.slack_set for b in boundaries):
            self.is_sat = res is None
        else:
            self.is_sat = False

        return res

    def _assert_bound(self, boundary, literal):
        """
        Adjusts the upper or lower bound on variable xi if the new bound is
        more limiting. The assignment of variable xi is adjusted to be
        within the new bound if needed.

        Also calls `self._update` to update the assignment for slack variables
        to keep all equalities satisfied.
        """
        if self.result:
            assert self.result[0] != False
        self.result = None

        xi = boundary.var
        ci, upper = boundary.to_rational(is_negated=literal < 0)

        s = 1 if upper else -1
        target_bound = xi.upper if upper else xi.lower
        opposing_bound = xi.lower if upper else xi.upper
        conflicting_lit = xi.lower_literal if upper else xi.upper_literal

        # If asserting lower bound, convert to equivalent upper bound situation
        # to simplify logic.
        c_norm = ci * s
        target_norm = target_bound * s
        opposing_norm = opposing_bound * s

        # Return `None` if new constraint is weaker than existing constraint.
        if c_norm >= target_norm:
            return None

        # Return conflict if new constraint directly conflicts with opposing constraint.
        if c_norm < opposing_norm:
            assert (opposing_bound.d * s >= 0) is True
            assert (ci.d * s <= 0) is True

            self.result = False, [-conflicting_lit, -literal]
            return self.result

        self.bound_history.append((xi, target_bound, upper))

        xi.set_bound(boundary, literal)

        if xi in self.nonslack and xi.assign * s > c_norm:
            self._update(xi, ci)

        if self.run_checks and all(not math.isinf(v.assign.q)
                                   for v in self.all_var):
            M = self.A
            X = Matrix([v.assign.q for v in self.all_var])
            assert all(abs(val) < 10 ** (-10) for val in M * X)

        return None

    def _update(self, xi, v):
        """
        Updates all slack variables that have equations that contain
        variable xi so that they stay satisfied given xi is equal to v.
        """
        i = xi.col_idx
        assert i is not None
        for j, b in enumerate(self.slack):
            aji = self.A[j, i]
            b.assign = b.assign + (v - xi.assign)*aji
        xi.assign = v

    def check(self):
        """
        Searches for an assignment that satisfies all constraints
        or determines that no such assignment exists and gives
        a minimal conflict clause that "explains" why the
        constraints are unsatisfiable.

        Returns
        =======

        (True, assignment) or (False, explanation)

        assignment : dict of LRAVariables to values
            Assigned values are tuples that represent a rational number
            plus some infinatesimal delta.

        explanation : set of ints
        """
        if self.is_sat:
            self.last_assign_snapshot = {var: var.assign for var in self.all_var}
            return True, self.last_assign_snapshot
        if self.result:
            return self.result

        from sympy.matrices.dense import Matrix
        M = self.A.copy()
        basic = {s: i for i, s in enumerate(self.slack)}  # contains the row index associated with each basic variable
        nonbasic = set(self.nonslack)
        while True:
            if self.run_checks:
                # nonbasic variables must always be within bounds
                assert all(((nb.assign >= nb.lower) == True) and ((nb.assign <= nb.upper) == True) for nb in nonbasic)

                # assignments for x must always satisfy Ax = 0
                # probably have to turn this off when dealing with strict ineq
                if all(not math.isinf(v.assign.q) for v in self.all_var):
                    X = Matrix([v.assign.q for v in self.all_var])
                    assert all(abs(val) < 10**(-10) for val in M*X)

                # check upper and lower match this format:
                # x <= rat + delta iff x < rat
                # x >= rat - delta iff x > rat
                # this wouldn't make sense:
                # x <= rat - delta
                # x >= rat + delta
                assert all(x.upper.d <= 0 for x in self.all_var)
                assert all(x.lower.d >= 0 for x in self.all_var)

            cand = [b for b in basic if b.assign < b.lower or b.assign > b.upper]

            if len(cand) == 0:
                self.last_assign_snapshot = {var: var.assign for var in self.all_var}
                return True, self.last_assign_snapshot

            xi = min(cand, key=lambda v: v.col_idx) # Bland's rule
            i = basic[xi]

            if xi.assign < xi.lower:
                cand = [nb for nb in nonbasic
                        if (M[i, nb.col_idx] > 0 and nb.assign < nb.upper)
                        or (M[i, nb.col_idx] < 0 and nb.assign > nb.lower)]
                if len(cand) == 0:
                    N_plus = [nb for nb in nonbasic if M[i, nb.col_idx] > 0]
                    N_minus = [nb for nb in nonbasic if M[i, nb.col_idx] < 0]

                    conflict = []
                    conflict += [nb.upper_literal for nb in N_plus]
                    conflict += [nb.lower_literal for nb in N_minus]
                    conflict.append(xi.lower_literal)
                    conflict = [-conflicting_lit for conflicting_lit in conflict]
                    return False, conflict
                xj = min(cand, key=str)
                M = self._pivot_and_update(M, basic, nonbasic, xi, xj, xi.lower)

            if xi.assign > xi.upper:
                cand = [nb for nb in nonbasic
                        if (M[i, nb.col_idx] < 0 and nb.assign < nb.upper)
                        or (M[i, nb.col_idx] > 0 and nb.assign > nb.lower)]

                if len(cand) == 0:
                    N_plus = [nb for nb in nonbasic if M[i, nb.col_idx] > 0]
                    N_minus = [nb for nb in nonbasic if M[i, nb.col_idx] < 0]

                    conflict_bounds = []
                    conflict_bounds += [nb.upper_literal for nb in N_minus]
                    conflict_bounds += [nb.lower_literal for nb in N_plus]
                    conflict_bounds.append(xi.upper_literal)

                    conflict = [-conflicting_lit for conflicting_lit in conflict_bounds]
                    return False, conflict
                xj = min(cand, key=lambda v: v.col_idx)
                M = self._pivot_and_update(M, basic, nonbasic, xi, xj, xi.upper)

    def _pivot_and_update(self, M, basic, nonbasic, xi, xj, v):
        """
        Pivots basic variable xi with nonbasic variable xj,
        and sets value of xi to v and adjusts the values of all basic variables
        to keep equations satisfied.
        """
        i, j = basic[xi], xj.col_idx
        assert j is not None
        assert M[i, j] != 0
        theta = (v - xi.assign)*(1/M[i, j])
        xi.assign = v
        xj.assign = xj.assign + theta
        for xk in basic:
            if xk != xi:
                k = basic[xk]
                akj = M[k, j]
                xk.assign = xk.assign + theta*akj
        # pivot
        basic[xj] = basic[xi]
        del basic[xi]
        nonbasic.add(xi)
        nonbasic.remove(xj)
        return self._pivot(M, i, j)

    @staticmethod
    def _pivot(M, i, j):
        """
        Performs a pivot operation about entry i, j of M by performing
        a series of row operations on a copy of M and returning the result.
        The original M is left unmodified.

        Conceptually, M represents a system of equations and pivoting
        can be thought of as rearranging equation i to be in terms of
        variable j and then substituting in the rest of the equations
        to get rid of other occurrences of variable j.

        Example
        =======

        >>> from sympy.matrices.dense import Matrix
        >>> from sympy.logic.algorithms.lra_theory import LRASolver
        >>> from sympy import var
        >>> Matrix(3, 3, var('a:i'))
        Matrix([
        [a, b, c],
        [d, e, f],
        [g, h, i]])

        This matrix is equivalent to:
        0 = a*x + b*y + c*z
        0 = d*x + e*y + f*z
        0 = g*x + h*y + i*z

        >>> LRASolver._pivot(_, 1, 0)
        Matrix([
        [ 0, -a*e/d + b, -a*f/d + c],
        [-1,       -e/d,       -f/d],
        [ 0,  h - e*g/d,  i - f*g/d]])

        We rearrange equation 1 in terms of variable 0 (x)
        and substitute to remove x from the other equations.

        0 = 0 + (-a*e/d + b)*y + (-a*f/d + c)*z
        0 = -x + (-e/d)*y + (-f/d)*z
        0 = 0 + (h - e*g/d)*y + (i - f*g/d)*z
        """
        Mij = M[i, j]
        if Mij == 0:
            raise ZeroDivisionError("Tried to pivot about zero-valued entry.")
        A = M.copy()
        A[i, :] = -A[i, :]/Mij
        for row in range(M.shape[0]):
            if row != i:
                A[row, :] = A[row, :] + A[row, j] * A[i, :]

        return A

    def backtrack(self):
        """
        Revert the most recent bound update to resolve a conflict.

        Pops the last state from the ``bound_history`` stack and restores the
        variable's previous upper or lower bound. It also reverts all variable
        assignments to their previous valid state using a dictionary,
        thus clearing the current conflict and restoring satisfiability.

        Raises
        ======

        ValueError
            If called when the ``bound_history`` stack is empty, indicating
            the solver's internal state is out of sync.
        """
        if not self.bound_history:
            raise ValueError("Cannot backtrack, bound_history stack is empty")

        xi, old_bound, upper = self.bound_history.pop()

        if upper:
            xi.upper = old_bound
        else:
            xi.lower = old_bound

        for var in self.all_var:
            var.assign = self.last_assign_snapshot[var]

        self.is_sat = True
        self.result = None

def _sep_const_coeff(expr):
    """
    Example
    =======

    >>> from sympy.logic.algorithms.lra_theory import _sep_const_coeff
    >>> from sympy.abc import x, y
    >>> _sep_const_coeff(2*x)
    (x, 2)
    >>> _sep_const_coeff(2*x + 3*y)
    (2*x + 3*y, 1)
    """
    if isinstance(expr, Add):
        return expr, sympify(1)
    const, var = sift(Mul.make_args(expr),
                      lambda c: len(sympify(c).free_symbols) == 0,
                      binary=True)
    return Mul(*var), Mul(*const)


def _sep_const_terms(expr):
    """
    Example
    =======

    >>> from sympy.logic.algorithms.lra_theory import _sep_const_terms
    >>> from sympy.abc import x, y
    >>> _sep_const_terms(2*x + 3*y + 2)
    (2*x + 3*y, 2)
    """
    const, var = sift(Add.make_args(expr),
                      lambda t: len(t.free_symbols) == 0,
                      binary=True)
    return Add(*var), Add(*const)


def _reduce_matrix(A, basic, nonbasic, nonatom_vars):
    """
    Remove every non-atom variable from the tableu A. This is discussed in
    Preprocessing part of the paper [1]_ as the "Gaussian Eliminaton".

    The idea is that, all non-atom variables are dependent of atom variables,
    which consistent-wise means that solving for atom variables should directly
    give solutions for non-atom variables.

    Therefore, any information about dependent, or to be more precise, non atom variables
    in the matrix A is not necessary and can be safely discarded without any correctness worries.

    E.g in,

        x >= 0 & x+y >= 1 -> Phi' := (x >= 0 & s1 >= 1), Phi_A := x + y = s1

    Since y is dependent, solving Phi' alone is enough, and _reduce_matrix should reduce Phi_A
    into collapsed matrix since it stores no useful information.

    Example
    =======

    Consider the formula:

        x >= 0 & z <= 1 & (x + y <= 5 | z + y >= 2)

    Here y is the only non-atom variable, so only y is removed, s1 = x+y, s2 = z+y.
    >>> from sympy.abc import x, y, z
    >>> from sympy import symbols
    >>> from sympy.solvers.solveset import linear_eq_to_matrix
    >>> from sympy.logic.algorithms.lra_theory import _reduce_matrix
    >>> s1, s2 = symbols('s1 s2')
    >>> nonbasic, basic = [x, y, z], [s1, s2]
    >>> A, _ = linear_eq_to_matrix([x + y - s1, z + y - s2], nonbasic + basic)
    >>> A
    Matrix([
    [1, 1, 0, -1,  0],
    [0, 1, 1,  0, -1]])
    >>> A, basic, nonbasic = _reduce_matrix(A, basic, nonbasic, {y})
    >>> basic, nonbasic
    ([s1], [x, z, s2])

    Notice that s2 became nonbasic.

    >>> A
    Matrix([[1, -1, 1, -1]])

    It is possible for the matrix A to collapse entirely, which happens when
    all the remaining terms are linearly independent. Or in other terms, The matrix A
    is no longer "stores" information about variables as there are no information to store.
    E.g,

         (x >= 0) & ((x + y <= 2) | (x + 2 * y - z >= 6)) & (Eq(x + y, 2) | (x + 2 * y - z > 4))

    only x is the atom variable so only y and z is removed, s1 = x+y and s2 = x+2*y-z.
    >>> nonbasic, basic = [x, y, z], [s1, s2]
    >>> A, _ = linear_eq_to_matrix([x + y - s1, x + 2 * y - z - s2], nonbasic + basic)
    >>> A
    Matrix([
    [1, 1,  0, -1,  0],
    [1, 2, -1,  0, -1]])
    >>> A, basic, nonbasic = _reduce_matrix(A, basic, nonbasic, {y, z})
    >>> basic, nonbasic
    ([], [x, s1, s2])

    Basic is empty, which in result should mean A has collapsed.
    >>> A.shape
    (0,3)
    """
    if not nonatom_vars:
        return A, basic, nonbasic

    kept_nonbasic = [v for v in nonbasic if v not in nonatom_vars]
    # order starts with the variables we want to eliminate
    # in rref, these variables will become pivots
    order = list(nonatom_vars) + basic + kept_nonbasic
    col_of = {v: i for i, v in enumerate(nonbasic + basic)}
    # reorder the matrix for the rref
    A = A[:, [col_of[v] for v in order]]

    B, pivots = A.rref()

    keep_rows = [r for r, pc in enumerate(pivots) if pc >= len(nonatom_vars)]
    new_basic = [order[pivots[r]] for r in keep_rows]
    basic_set = set(new_basic)
    new_nonbasic = [v for v in kept_nonbasic + basic if v not in basic_set]

    order_pos = {v: i for i, v in enumerate(order)}
    A = -B[keep_rows, [order_pos[v] for v in new_nonbasic + new_basic]]
    return A, new_basic, new_nonbasic


class Boundary:
    """
    Represents an upper or lower bound between a symbol
    and some constant.

    Example
    =======

    >>> from sympy.logic.algorithms.lra_theory import Boundary, LRAVariable
    >>> from sympy.abc import x
    >>> var = LRAVariable(x)
    >>> # x <= 5
    >>> b1 = Boundary(var, 5, upper=True, strict=False)
    >>> b1.get_inequality()
    x <= 5
    >>> # x > 10 (represented as a lower bound with strict=True)
    >>> b2 = Boundary(var, 10, upper=False, strict=True)
    >>> b2.get_inequality()
    x > 10
    """
    def __init__(self, var, const, upper, strict=None):
        self.var = var
        if isinstance(const, tuple):
            s = const[1] != 0
            if strict is not None:
                assert s == strict
            self.bound = const[0]
            self.strict = s
        else:
            self.bound = const
            self.strict = strict
        self.upper = upper
        assert self.strict is not None

    def to_rational(self, is_negated):
        """
        Return the LRARational bound and effective direction (upper=True)
        considering whether the boundary is negated.
        """
        upper = self.upper != is_negated
        delta = 0
        if self.strict != is_negated:
            delta = -1 if upper else 1
        return LRARational(self.bound, delta), upper

    def get_inequality(self):
        if self.upper and self.strict:
            return self.var.var < self.bound
        elif not self.upper and self.strict:
            return self.var.var > self.bound
        elif self.upper:
            return self.var.var <= self.bound
        else:
            return self.var.var >= self.bound

    def __repr__(self):
        return repr("Boundary(" + repr(self.get_inequality()) + ")")

    def __eq__(self, other):
        if not isinstance(other, Boundary):
            return NotImplemented
        return ((self.var, self.bound, self.strict, self.upper)
            == (other.var, other.bound, other.strict, other.upper))

    def __hash__(self):
        return hash((self.var, self.bound, self.strict, self.upper))


class LRARational():
    """
    Represents a rational plus or minus some amount
    of arbitrary small deltas.
    """
    def __init__(self, rational, delta):
        self.value = (rational, delta)

    @property
    def q(self):
        return self.value[0]

    @property
    def d(self):
        return self.value[1]

    def __lt__(self, other):
        return self.value < other.value

    def __le__(self, other):
        return self.value <= other.value

    def __eq__(self, other):
        if not isinstance(other, LRARational):
            return NotImplemented
        return self.value == other.value

    def __add__(self, other):
        return LRARational(self.q + other.q, self.d + other.d)

    def __sub__(self, other):
        return LRARational(self.q - other.q, self.d - other.d)

    def __mul__(self, other):
        assert not isinstance(other, LRARational)
        return LRARational(self.q * other, self.d * other)

    def __getitem__(self, index):
        return self.value[index]

    def __repr__(self):
        return repr(self.value)


class LRAVariable():
    """
    Object to keep track of upper and lower bounds
    on `self.var`.
    """
    def __init__(self, var):
        self.initialize()
        self.var = var
        self.col_idx = None

    def initialize(self):
        self.upper = LRARational(float("inf"), 0)
        self.upper_literal = None
        self.lower = LRARational(-float("inf"), 0)
        self.lower_literal = None
        self.assign = LRARational(0,0)

    def __repr__(self):
        return repr(self.var)

    def set_bound(self, boundary, literal):
        """
        Set the upper or lower bound and record its source.

        Example
        =======

        >>> from sympy.logic.algorithms.lra_theory import LRAVariable, Boundary
        >>> from sympy.abc import x
        >>> v = LRAVariable(x)
        >>> b = Boundary(v, 10, upper=False, strict=False)
        >>> # Asserting a lower bound x >= 10 using literal 5
        >>> v.set_bound(b, 5)
        >>> v.lower
        (10, 0)
        >>> v.lower_literal
        5
        """
        is_negated = literal < 0
        ci, upper = boundary.to_rational(is_negated)
        if upper:
            self.upper = ci
            self.upper_literal = literal
        else:
            self.lower = ci
            self.lower_literal = literal

    def __eq__(self, other):
        if not isinstance(other, LRAVariable):
            return False
        return other.var == self.var

    def __hash__(self):
        return hash(self.var)
