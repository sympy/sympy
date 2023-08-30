from sympy.solvers.solveset import linear_eq_to_matrix
from sympy.matrices.dense import eye
from sympy.assumptions import Predicate
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.ask import Q
from sympy.core import Dummy
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.relational import Eq, Ne
from sympy.core.sympify import sympify
from sympy.core.singleton import S
from sympy import SYMPY_DEBUG
from sympy.core.numbers import Rational, oo
from sympy.matrices.dense import Matrix
from sympy.matrices.common import MatrixKind
from sympy.core.kind import NumberKind


def sep_const_coeff(expr):
    if isinstance(expr, Add):
        return expr, sympify(1)

    if isinstance(expr, Mul):
        coeffs = expr.args
    else:
        coeffs = [expr]

    var, const = [], []
    for c in coeffs:
        c = sympify(c)
        if len(c.free_symbols)==0:
            const.append(c)
        else:
            var.append(c)
    return Mul(*var), Mul(*const)


def list_terms(expr):
    if not isinstance(expr, Add):
        return [expr]

    return expr.args


def sep_const_terms(expr):
    if isinstance(expr, Add):
        terms = expr.args
    else:
        terms = [expr]

    var, const = [], []
    for t in terms:
        if len(t.free_symbols) == 0:
            const.append(t)
        else:
            var.append(t)
    return sum(var), sum(const)


def standardize_binrel(prop):
    assert prop.function in [Q.ge, Q.gt, Q.le, Q.lt, Q.eq]

    expr = prop.lhs - prop.rhs
    if prop.function in [Q.ge, Q.gt]:
        expr = -expr
    var, const = sep_const_terms(expr)

    if prop.function == Q.eq:
        return Q.eq(var, const)
    if prop.function in [Q.gt, Q.lt]:
        return Q.lt(var, const)
    else:
        return Q.le(var, const)

def _eval_binrel(binrel):
    if not (len(binrel.lhs.free_symbols) == 0 and len(binrel.rhs.free_symbols) == 0):
        return binrel
    if binrel.function == Q.lt:
        res = binrel.lhs < binrel.rhs
    if binrel.function == Q.gt:
        res = binrel.lhs > binrel.rhs
    if binrel.function == Q.le:
        res = binrel.lhs <= binrel.rhs
    if binrel.function == Q.ge:
        res = binrel.lhs >= binrel.rhs
    if binrel.function == Q.eq:
        res = Eq(binrel.lhs, binrel.rhs)
    if binrel.function == Q.ne:
        res = Ne(binrel.lhs, binrel.rhs)

    if res == True or res == False:
        return res
    else:
        return None


class Boundry:
    """
    Represents an upper or lower bound or an equality between a symbol
    and some constant.
    """
    def __init__(self, var, const, upper, equality, strict=None):
        if not equality in [True, False]:
            assert equality in [True, False]


        self.var = var
        if isinstance(const, tuple):
            s = const[1] != 0
            if strict:
                assert s == strict
            self.bound = const[0]
            self.strict = s
        else:
            self.bound = const
            self.strict = strict
        self.upper = upper if not equality else None
        self.equality = equality
        self.strict = strict
        assert self.strict is not None

    @staticmethod
    def from_upper(var):
        neg = -1 if var.upper_from_neg else 1
        b = Boundry(var, var.upper[0], True, var.upper_from_eq, var.upper[1] != 0)
        if neg < 0:
            b = b.get_negated()
        return b, neg

    @staticmethod
    def from_lower(var):
        neg = -1 if var.lower_from_neg else 1
        b = Boundry(var, var.lower[0], False, var.lower_from_eq, var.lower[1] != 0)
        if neg < 0:
            b = b.get_negated()
        return b, neg

    def get_negated(self):
        return Boundry(self.var, self.bound, not self.upper, self.equality, not self.strict)

    def get_inequality(self):
        if self.equality:
            return Eq(self.var.var, self.bound)
        elif self.upper and self.strict:
            return self.var.var < self.bound
        elif not self.upper and self.strict:
            return self.var.var > self.bound
        elif self.upper:
            return self.var.var <= self.bound
        else:
            return self.var.var >= self.bound

    def __repr__(self):
        return repr("Boundry(" + repr(self.get_inequality()) + ")")

    def __eq__(self, other):
        other = (other.var, other.bound, other.strict, other.upper, other.equality)
        return (self.var, self.bound, self.strict, self.upper, self.equality) == other

    def __hash__(self):
        return hash((self.var, self.bound, self.strict, self.upper, self.equality))


class LRARational():
    """
    Represents a rational plus or minus some amount
    of arbitrary small deltas.
    """
    def __init__(self, rational, delta):
        self.value = (rational, delta)

    def __lt__(self, other):
        return self.value < other.value

    def __le__(self, other):
        return self.value <= other.value

    def __eq__(self, other):
        return self.value == other.value

    def __add__(self, other):
        return LRARational(self.value[0] + other.value[0], self.value[1] + other.value[1])

    def __sub__(self, other):
        return LRARational(self.value[0] - other.value[0], self.value[1] - other.value[1])

    def __mul__(self, other):
        assert not isinstance(other, LRARational)
        return LRARational(self.value[0] * other, self.value[1] * other)

    def __getitem__(self, index):
        return self.value[index]

    def __repr__(self):
        return repr(self.value)


class LRAVariable():
    def __init__(self, var):
        self.upper = LRARational(float("inf"), 0)
        self.upper_from_eq = False
        self.upper_from_neg = False
        self.lower = LRARational(-float("inf"), 0)
        self.lower_from_eq = False
        self.lower_from_neg = False
        self.assign = LRARational(0,0)
        self.var = var
        self.col_idx = None

    def __repr__(self):
        return repr(self.var)

    def __eq__(self, other):
        return other.var == self.var

    def __hash__(self):
        return hash(self.var)

class UnhandledNumber(Exception):
    pass


ALLOWED_PRED = {Q.eq, Q.gt, Q.lt, Q.le, Q.ge}
HANDLE_NEGATION = True

class LRASolver():
    """
    Linear Arithmatic Solver for DPLL(T) implemented with algorithm based on
    the Dual Simplex method and Bland's pivoting rule.

    TODO:
     - Implement and utilize backtracking
     - Handle non-rational real numbers
     - Handle positive and negative infinity

    References
    ==========

    .. [1] Dutertre, B., de Moura, L.:
           A Fast Linear-Arithmetic Solver for DPLL(T)
           https://link.springer.com/chapter/10.1007/11817963_11
    """

    def __init__(self, A, slack_variables, nonslack_variables, enc_to_boundry, x_subs, s_subs, testing_mode):
        """
        Use the "from_encoded_cnf" method to create a new LRASolver.
        """
        if any(not isinstance(a, Rational) for a in A):
            raise UnhandledNumber
        if any(not isinstance(b.bound, Rational) and (b.bound != oo) and (b.bound != -oo) for b in enc_to_boundry.values()):
            raise UnhandledNumber

        self.run_checks = testing_mode # set to True to turn on assert statements
        m, n = len(slack_variables), len(slack_variables)+len(nonslack_variables)
        if m != 0:
            assert A.shape == (m, n)

        if self.run_checks:
            assert A[:, n-m:] == -eye(m)

        self.enc_to_boundry = enc_to_boundry
        self.boundry_to_enc = {value: key for key, value in enc_to_boundry.items()}

        self.A = A # TODO: Row reduce A
        self.slack = slack_variables
        self.nonslack = nonslack_variables
        self.all_var = nonslack_variables + slack_variables

        self.x_subs = x_subs
        self.s_subs = s_subs

        # if previously we knew this was sat
        # if no assertions about slack variables have been made since
        # and changing bounds on nonslack hasn't lead to unsat
        self.is_sat = True
        self.slack_set = set(slack_variables) # used for checking if var is slack


        # always one of: (True, Assignment), (False, ConflictClause), None
        self.result = None

    # @staticmethod
    # def _remove_uneeded_variables(A, nonbasic, num_unused):
    #     """
    #     WIP speed up. See the following cs stack exchange for more information.
    #     cs.stackexchange.com/questions/161709/a-fast-linear-arithmetic-solver-how-can-gaussian-elimination-be-used-to-simplif/161713#161713
    #     """
    #     M = A.copy()
    #     if num_unused == 0:
    #         return M, nonbasic
    #
    #     rref, pivs = M.rref()
    #
    #     removed = set()
    #     for col in pivs:
    #         if col >= num_unused:
    #             break
    #
    #         # col in terms of other variables
    #         sol = rref[col, :]
    #         if sol*M[col, col] == M[col, :]:
    #             continue
    #
    #         removed.add(col)
    #         for row in range(M.shape[0]):
    #             M[row, :] = M[row, :]- sol*M[row, col]
    #
    #     not_removed = [col for col in range(M.shape[1]) if col not in removed]
    #     M = M[:,not_removed]
    #     nonbasic = [nonbasic[nr] for nr in not_removed if nr < len(nonbasic)]
    #
    #     return M, nonbasic

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

        encoding = {}  # maps int to Boundry
        A = []

        basic = []
        s_count = 0
        s_subs = {}
        nonbasic = []
        x_count = 0
        x_subs = {}

        if testing_mode:
            encoded_cnf_items = sorted(encoded_cnf.encoding.items(), key=lambda x: str(x))
        else:
            encoded_cnf_items = encoded_cnf.encoding.items()


        empty_var = Dummy()
        var_to_lra_var = {}
        conflicts = []

        for prop, enc in encoded_cnf_items:
            if isinstance(prop, Predicate):
                prop = prop(empty_var)
            if not isinstance(prop, AppliedPredicate):
                continue

            assert prop.function in ALLOWED_PRED
            if prop.lhs.kind == MatrixKind(NumberKind) or prop.rhs.kind == MatrixKind(NumberKind):
                raise ValueError(f"{prop} contains matrix variables")
            if prop.lhs == S.NaN or prop.rhs == S.NaN:
                raise ValueError(f"{prop} contains nan")
            if prop.lhs.is_imaginary or prop.rhs.is_imaginary:
                raise UnhandledNumber(f"{prop} contains an imaginary component")
            if prop.lhs == oo or prop.rhs == oo:
                raise UnhandledNumber(f"{prop} contains infinity")

            prop = _eval_binrel(prop)  # simplify to True / False if possible
            if prop == True:
                conflicts.append([enc])
                continue
            elif prop == False:
                conflicts.append([-enc])
                continue
            elif prop is None:
                raise UnhandledNumber(f"{prop} contains no variables and could not be evaluated as True or False")

            expr = prop.lhs - prop.rhs
            if prop.function in [Q.ge, Q.gt]:
                expr = -expr

            var, const = sep_const_terms(expr)
            var, var_coeff = sep_const_coeff(var)
            const = const / var_coeff

            # replace each term in expr with dummy _xi variable
            # e.g.: -x -2*y + z --> [_x3, -_x1, -2*_x2]
            terms = []
            for term in list_terms(var):
                assert not isinstance(term, Add)
                term, term_coeff = sep_const_coeff(term)
                if term not in x_subs:
                    x_count += 1
                    x_subs[term] = Dummy(f"x{x_count}")
                    var_to_lra_var[x_subs[term]] = LRAVariable(x_subs[term])
                    nonbasic.append(x_subs[term])
                terms.append(term_coeff * x_subs[term])
            # If there are multiple variable terms, replace them with a dummy _si variable.
            # If needed (no other expr has this sum of variable terms), create a new Dummy _si varaible.
            if len(terms) > 1:
                var = sum(terms)
                if var not in s_subs:
                    s_count += 1
                    d = Dummy(f"s{s_count}")
                    var_to_lra_var[d] = LRAVariable(d)
                    basic.append(d)
                    s_subs[var] = d
                    A.append(var - d)
                var = s_subs[var]
            else:
                var = terms[0]

            assert var_coeff != 0

            equality = prop.function == Q.eq
            upper = var_coeff > 0 if not equality else None
            strict = prop.function in [Q.gt, Q.lt]
            b = Boundry(var_to_lra_var[var], -const, upper, equality, strict)
            encoding[enc] = b

        A, _ = linear_eq_to_matrix(A, nonbasic + basic)
        nonbasic = [var_to_lra_var[nb] for nb in nonbasic]
        basic = [var_to_lra_var[b] for b in basic]
        for idx, var in enumerate(nonbasic + basic):
            var.col_idx = idx

        return LRASolver(A, basic, nonbasic, encoding, x_subs, s_subs, testing_mode), conflicts

    def reset_bounds(self):
        """
        Resets the state of the LRASolver to before
        anything was asserted.
        """
        self.result = None
        for var in self.all_var:
            var.lower = LRARational(-float("inf"), 0)
            var.lower_from_eq = False
            var.lower_from_neg = False
            var.upper = LRARational(float("inf"), 0)
            var.upper_from_eq= False
            var.lower_from_neg = False
            var.assign = LRARational(0, 0)

    def assert_lit(self, enc_constraint):
        """
        Assert a literal representing a constraint
        and update the interal state accordingly.

        Parameters
        ==========

        enc_constraint : int
            A mapping of encodings to constraints
            can be found in `self.enc_to_boundry`.

        Returns
        =======

        None or (False, explanation)

        explanation : set of ints
            A conflict clause that "explains" why
            the literals asserted so far are unsatisfiable.
        """
        if abs(enc_constraint) not in self.enc_to_boundry:
            raise ValueError("Tried to assert a literal with no encoding")

        if not HANDLE_NEGATION and enc_constraint < 0:
            return None

        boundry = self.enc_to_boundry[abs(enc_constraint)]
        sym, c, negated = boundry.var, boundry.bound, enc_constraint < 0

        if boundry.equality and negated:
            return None # negated equality is not handled and should only appear as conflict clause

        upper = boundry.upper != negated
        if boundry.strict != negated:
            delta = -1 if upper else 1
            c = LRARational(c, delta)
        else:
            c = LRARational(c, 0)

        if boundry.equality:
            res1 = self._assert_lower(sym, c, from_equality=True, from_neg=negated)
            if res1 and res1[0] == False:
                res = res1
            else:
                res2 = self._assert_upper(sym, c, from_equality=True, from_neg=negated)
                res =  res2
        elif upper:
            res = self._assert_upper(sym, c, from_neg=negated)
        else:
            res = self._assert_lower(sym, c, from_neg=negated)

        if self.is_sat and sym not in self.slack_set:
            self.is_sat = res is None
        else:
            self.is_sat = False

        return res

    def _assert_upper(self, xi, ci, from_equality=False, from_neg=False):
        if self.result:
            assert self.result[0] != False
        self.result = None
        if ci >= xi.upper:
            return None
        if ci < xi.lower:
            assert (xi.lower[1] >= 0) is True
            assert (ci[1] <= 0) is True

            lit1, neg1 = Boundry.from_lower(xi)

            lit2 = Boundry(var=xi, const=ci[0], strict=ci[1] != 0, upper=True, equality=from_equality)
            if from_neg:
                lit2 = lit2.get_negated()
            neg2 = -1 if from_neg else 1

            conflict = [-neg1*self.boundry_to_enc[lit1], -neg2*self.boundry_to_enc[lit2]]
            self.result = False, conflict
            return self.result
        xi.upper = ci
        xi.upper_from_eq = from_equality
        xi.upper_from_neg = from_neg
        if xi in self.nonslack and xi.assign > ci:
            self._update(xi, ci)

        if self.run_checks and all(v.assign[0] != float("inf") and v.assign[0] != -float("inf")
                                   for v in self.all_var):
            M = self.A
            X = Matrix([v.assign[0] for v in self.all_var])
            assert all(abs(val) < 10 ** (-10) for val in M * X)

        return None

    def _assert_lower(self, xi, ci, from_equality=False, from_neg=False):
        if self.result:
            assert self.result[0] != False
        self.result = None
        if ci <= xi.lower:
            return None
        if ci > xi.upper:
            assert (xi.upper[1] <= 0) is True
            assert (ci[1] >= 0) is True

            lit1, neg1 = Boundry.from_upper(xi)

            lit2 = Boundry(var=xi, const=ci[0], strict=ci[1] != 0, upper=False, equality=from_equality)
            if from_neg:
                lit2 = lit2.get_negated()
            neg2 = -1 if from_neg else 1

            conflict = [-neg1*self.boundry_to_enc[lit1],-neg2*self.boundry_to_enc[lit2]]
            self.result = False, conflict
            return self.result
        xi.lower = ci
        xi.lower_from_eq = from_equality
        xi.lower_from_neg = from_neg
        if xi in self.nonslack and xi.assign < ci:
            self._update(xi, ci)

        if self.run_checks and all(v.assign[0] != float("inf") and v.assign[0] != -float("inf")
                                   for v in self.all_var):
            M = self.A
            X = Matrix([v.assign[0] for v in self.all_var])
            assert all(abs(val) < 10 ** (-10) for val in M * X)

        return None

    def _update(self, xi, v):
        i = xi.col_idx
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
            return True, {var: var.assign for var in self.all_var}

        if self.result:
            return self.result

        def _debug_internal_state_printer1(iteration, A, bas, variables):
            if not SYMPY_DEBUG:
                return
            import sys
            from sympy.matrices.dense import Matrix
            from sympy import pprint

            bvar = [None]*len(bas)
            for v, idx in bas.items():
                bvar[idx] = v

            r1 = Matrix([variables])
            c1 = Matrix(bvar)
            corner = Matrix([[iteration]])

            tableau = Matrix([[corner, r1], [c1, A]])
            pprint(tableau)
            sys.stderr.write("\n")
            sys.stderr.write(f"{self.assign}\n")
            for v in self.upper:
                sys.stderr.write(str(self.lower[v]) + " <= " + str(v) + " <= " + str(self.upper[v]) + "\n")

        def _debug_internal_state_printer2(xi, xj):
            if not SYMPY_DEBUG:
                return
            import sys
            sys.stderr.write(f"\npivoting {xi} with {xj}\n\n")

        from sympy.matrices.dense import Matrix
        M = self.A.copy()
        basic = {s: i for i, s in enumerate(self.slack)}  # contains the row index associated with each basic variable
        nonbasic = set(self.nonslack)
        iteration = 0
        while True:
            iteration += 1; _debug_internal_state_printer1(iteration, M, basic, self.all_var)

            if self.run_checks:
                # nonbasic variables must always be within bounds
                assert all(((nb.assign >= nb.lower) == True) and ((nb.assign <= nb.upper) == True) for nb in nonbasic)

                # assignments for x must always satisfy Ax = 0
                # probably have to turn this off when dealing with strict ineq
                if all(v.assign[0] != float("inf") and v.assign[0] != -float("inf")
                                   for v in self.all_var):
                    X = Matrix([v.assign[0] for v in self.all_var])
                    assert all(abs(val) < 10**(-10) for val in M*X)

                # check upper and lower match this format:
                # x <= rat + delta iff x < rat
                # x >= rat - delta iff x > rat
                # this wouldn't make sense:
                # x <= rat - delta
                # x >= rat + delta
                assert all(x.upper[1] <= 0 for x in self.all_var)
                assert all(x.lower[1] >= 0 for x in self.all_var)

            cand = [b for b in basic if b.assign < b.lower or b.assign > b.upper]

            if len(cand) == 0:
                return True, {var: var.assign for var in self.all_var}

            xi = sorted(cand, key=lambda v: v.col_idx)[0] # Bland's rule
            i = basic[xi]

            if xi.assign < xi.lower:
                cand = [nb for nb in nonbasic
                        if (M[i, nb.col_idx] > 0 and nb.assign < nb.upper)
                        or (M[i, nb.col_idx] < 0 and nb.assign > nb.lower)]
                if len(cand) == 0:
                    N_plus = [nb for nb in nonbasic if M[i, nb.col_idx] > 0]
                    N_minus = [nb for nb in nonbasic if M[i, nb.col_idx] < 0]

                    conflict = []
                    conflict += [Boundry.from_upper(nb) for nb in N_plus]
                    conflict += [Boundry.from_lower(nb) for nb in N_minus]
                    conflict.append(Boundry.from_lower(xi))
                    conflict = [-neg*self.boundry_to_enc[c] for c, neg in conflict]
                    return False, conflict
                xj = sorted(cand, key=lambda v: str(v))[0]
                _debug_internal_state_printer2(xi, xj)
                M = self._pivot_and_update(M, basic, nonbasic, xi, xj, xi.lower)

            if xi.assign > xi.upper:
                cand = [nb for nb in nonbasic
                        if (M[i, nb.col_idx] < 0 and nb.assign < nb.upper)
                        or (M[i, nb.col_idx] > 0 and nb.assign > nb.lower)]

                if len(cand) == 0:
                    N_plus = [nb for nb in nonbasic if M[i, nb.col_idx] > 0]
                    N_minus = [nb for nb in nonbasic if M[i, nb.col_idx] < 0]

                    conflict = []
                    conflict += [Boundry.from_upper(nb) for nb in N_minus]
                    conflict += [Boundry.from_lower(nb) for nb in N_plus]
                    conflict.append(Boundry.from_upper(xi))

                    conflict = [-neg*self.boundry_to_enc[c] for c, neg in conflict]
                    return False, conflict
                xj = sorted(cand, key=lambda v: v.col_idx)[0]
                _debug_internal_state_printer2(xi, xj)
                M = self._pivot_and_update(M, basic, nonbasic, xi, xj, xi.upper)

    def _pivot_and_update(self, M, basic, nonbasic, xi, xj, v):
        i, j = basic[xi], xj.col_idx
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
        a series of row operations on a copy of M and returing the result.
        The original M is left unmodified.

        Conceptually, M represents a system of equations and pivoting
        can be thought of as rearranging equation i to be in terms of
        variable j and then substituting in the rest of the equations
        to get rid of other occurances of variable j.

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

        We rearage equation 1 in terms of variable 0 (x)
        and substitute to remove x from the other equations.

        0 = 0 + (-a*e/d + b)*y + (-a*f/d + c)*z
        0 = -x + (-e/d)*y + (-f/d)*z
        0 = 0 + (h - e*g/d)*y + (i - f*g/d)*z
        """
        _, _, Mij = M[i, :], M[:, j], M[i, j]
        if Mij == 0:
            raise ZeroDivisionError("Tried to pivot about zero-valued entry.")
        A = M.copy()
        A[i, :] = -A[i, :]/Mij
        for row in range(M.shape[0]):
            if row != i:
                A[row, :] = A[row, :] + A[row, j] * A[i, :]

        return A
