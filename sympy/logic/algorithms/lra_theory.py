from sympy.solvers.solveset import linear_eq_to_matrix
from sympy.matrices.dense import eye
from sympy.logic.boolalg import BooleanFunction
from sympy.assumptions.relation.binrel import AppliedBinaryRelation
from sympy.assumptions.ask import Q
from sympy.core import Dummy
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.relational import Relational, Eq, Ge, Lt, Le
from sympy import SYMPY_DEBUG


def sep_const_terms(expr):
    if isinstance(expr, Add):
        terms = expr.args
    else:
        terms = [expr]

    var, const = [], []
    for t in terms:
        if t.is_constant():
            const.append(t)
        else:
            var.append(t)
    return sum(var), sum(const)
def sep_const_coeff(expr):
    if isinstance(expr, Add):
        return expr, 1

    if isinstance(expr, Mul):
        coeffs = expr.args
    else:
        coeffs = [expr]

    var, const = [], []
    for c in coeffs:
        if c.is_constant():
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
        if t.is_constant():
            const.append(t)
        else:
            var.append(t)
    return sum(var), sum(const)


class Boundry:
    """
    Represents an upper or lower bound or an equality between a symbol and some constant.
    """
    def __init__(self, var, const, upper, equality, strict=False):
        self.var = var
        self.bound = const
        self.upper = upper
        self.equality = equality
        self.strict = strict

class LRASolver():
    """
    Linear Arithmatic Solver for DPLL(T) implemented with algorithm based on
    the Dual Simplex method and Bland's pivoting rule.

    References
    ==========

    .. [1] Dutertre, B., de Moura, L.: A Fast Linear-Arithmetic Solver for DPLL(T)
           https://link.springer.com/chapter/10.1007/11817963_11
    """

    def __init__(self, A, slack_variables, nonslack_variables, boundry_enc=None):
        self.run_checks = False # set to True to turn on assert statements
        m, n = len(slack_variables), len(slack_variables)+len(nonslack_variables)
        if m != 0:
            assert A.shape == (m, n)

        if self.run_checks:
            assert A[:, n-m:] == -eye(m)

        self.boundry_enc = boundry_enc

        self.A = A # TODO: Row reduce A
        self.slack = slack_variables
        self.nonslack = nonslack_variables
        self.all_var = nonslack_variables + slack_variables
        self.col_index = {v: i for i, v in enumerate(self.all_var)}

        self.lower = {}
        self.upper = {}
        self.assign = {}

        self.result = None # always one of: (True, Assignment), (False, Clause), None

        for var in self.all_var:
            self.lower[var] = (-float("inf"), 0)
            self.upper[var] = (float("inf"), 0)
            self.assign[var] = (0, 0)

    @staticmethod
    def from_encoded_cnf(encoded_cnf):
        # TODO: Preprecessing needs to be done to encoded_cnf
        # x - y > 0 should be the same as x > y
        encoding = {}  # maps int to Boundry
        A = []

        basic = []
        s_count = 0
        s_subs = {}
        nonbasic = []
        x_count = 0
        x_subs = {}

        # sorted here to remove nondeterminism
        for prop, enc in sorted(encoded_cnf.encoding.items(), key=lambda x: str(x)):
            if not isinstance(prop, AppliedBinaryRelation):
                continue
            assert prop.function in [Q.le, Q.ge, Q.eq, Q.gt, Q.lt]

            expr = prop.lhs - prop.rhs
            if prop.function in [Q.ge, Q.gt]:
                expr = -expr

            var, const = sep_const_terms(expr)
            var, var_coeff = sep_const_coeff(var)
            const = const / var_coeff

            terms = []
            for term in list_terms(var):
                assert not isinstance(term, Add)
                term, term_coeff = sep_const_coeff(term)
                if term not in x_subs:
                    x_count += 1
                    x_subs[term] = Dummy(f"x{x_count}")
                    nonbasic.append(x_subs[term])
                terms.append(term_coeff * x_subs[term])

            if len(terms) > 1:
                var = sum(terms)
                if var not in s_subs:
                    s_count += 1
                    d = Dummy(f"s{s_count}")
                    basic.append(d)
                    s_subs[var] = d
                    A.append(Eq(var, d))
                var = s_subs[var]
            else:
                var = terms[0]

            assert var_coeff != 0

            equality = prop.function == Q.eq
            upper = var_coeff > 0 if not equality else None
            strict = prop.function in [Q.gt, Q.lt]
            b = Boundry(var, -const, upper, equality, strict)
            encoding[enc] = b
            print(enc, prop, "-->", "var:", b.var, "bound:", b.bound, "upper:", b.upper)

        A, _ = linear_eq_to_matrix(A, nonbasic + basic)
        print(A, x_subs, s_subs)
        return LRASolver(A, basic, nonbasic, boundry_enc=encoding), x_subs, s_subs

    @staticmethod
    def preprocess(BF, variables):
        """
        Note that this function is only for testing.

        It's a bit easier to use than from_encoded_cnf, but it has some bugs which from_encoded_cnf doesn't have.
        Also, from_encoded_cnf is faster.
        """
        from sympy.matrices.dense import Matrix
        equations = []
        count = 0
        sub = {}
        set_variables = set(variables) # this should fix time complexity problems

        def to_standard_form(ineq, variables):
            expr = ineq.args[0] - ineq.args[1]
            A, B = linear_eq_to_matrix(expr, variables)
            lhs = (A * Matrix(variables))[0, 0]
            rhs = B[0, 0]
            return ineq.func(lhs, rhs)

        def do(b):
            if isinstance(b, BooleanFunction):
                return b.func(*[do(arg) for arg in b.args])
            elif isinstance(b, Relational):
                b = to_standard_form(b, variables)
                expr, const = b.args
                if isinstance(expr, Add):
                    if expr not in sub:
                        nonlocal count
                        count += 1
                        sub[expr] = Dummy(f"s{count}")
                        equations.append(Eq(sub[expr], expr))
                    return b.func(sub[expr], const)
                else:
                    return b
            else:
                return b

        res = do(BF)

        nonbasic = variables
        basic = list(sub.values())

        A, _ = linear_eq_to_matrix(equations, variables + basic)
        A = -A # identity matrix should be negative

        return res, LRASolver(A, basic, nonbasic), sub


    def assert_enc_boundry(self, enc_boundry):
        boundry = self.boundry_enc[enc_boundry]
        sym, c = boundry.var, boundry.bound

        if boundry.strict:
            delta = -1 if boundry.upper else 1
            c = (c, delta)
        else:
            c = (c, 0)

        if boundry.equality:
            res1 = self._assert_lower(sym, c)
            if res1[0] == "UNSAT":
                return res1
            res2 = self._assert_upper(sym, c)
            if res1[0] == "UNSAT":
                return res2
            if "OK" in [res1[0], res2[0]]:
                return "OK"
            return res2 # SAT
        elif boundry.upper:
            return self._assert_upper(sym, c)
        else:
            return self._assert_lower(sym, c)





    def assert_con(self, atom):
        if isinstance(atom, AppliedBinaryRelation):
            sym, c = atom.arguments
            assert sym.is_symbol
            assert c.is_constant()
            if atom.function == Q.ge:
                return self._assert_lower(sym, c)
            if atom.function == Q.le:
                return self._assert_upper(sym, c)
        elif isinstance(atom, Relational):
            sym, c = atom.args
            assert sym.is_symbol
            assert c.is_constant()

            if atom.func == Ge:
                return self._assert_lower(sym, c)
            if atom.func == Le:
                return self._assert_upper(sym, c)

        raise ValueError(f"{atom} could not be asserted.")
    def _assert_upper(self, xi, ci):
        self.result = None
        if ci >= self.upper[xi]:
            self.result = "SAT", self.assign
            return self.result
        if ci < self.lower[xi]:
            self.result = "UNSAT", "WIP"#{(xi, self.lower[xi]), (xi <= ci)}
            return self.result
        self.upper[xi] = ci
        if xi in self.nonslack and self.assign[xi] > ci:
            self._update(xi, ci)

        return "OK", None

    def _assert_lower(self, xi, ci):
        self.result = None
        if ci <= self.lower[xi]:
            self.result = "SAT", self.assign
            return self.result
        if ci > self.upper[xi]:
            self.result = "UNSAT", "WIP"#{xi <= self.upper[xi], xi >= ci}
            return self.result
        self.lower[xi] = ci
        if xi in self.nonslack and self.assign[xi] < ci:
            self._update(xi, ci)

        return "OK", None

    def _update(self, xi, v):
        i = self.col_index[xi]
        for j, b in enumerate(self.slack):
            aji = self.A[j, i]
            lhs = self.assign[b][0] + aji*(v[0] - self.assign[xi][0])
            rhs = self.assign[b][1] + aji*(v[1] - self.assign[xi][1])
            self.assign[b] = (lhs, rhs)
        self.assign[xi] = v

    def get_assignment(self, xi):
        pass # not sure if it's possible to get a rational assignment when involving strict inequalities
        # low, high, assign = self.lower[xi], self.high[xi], self.assign
        # if low[0] == float("inf") and high[0] == float("inf"):
        #     return 0
        # if low[0] == float("inf"):
        #     return high[0] - 1
        # if high[0] == float("inf"):
        #     return low[0] + 1



    def backtrack(self):
        pass

    def check(self):
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
                # nonbasic variables always must be within bounds
                assert all(((self.assign[nb] >= self.lower[nb]) == True) and ((self.assign[nb] <= self.upper[nb]) == True) for nb in nonbasic)

                # assignments for x must always satisfy Ax = 0
                # probably have to turn this off when dealing with strict ineq
                X = Matrix([self.assign[v][0] for v in self.col_index])
                assert all(abs(val) < 10**(-10) for val in M*X)

            cand = [b for b in basic
             if self.assign[b] < self.lower[b]
             or self.assign[b] > self.upper[b]]

            if len(cand) == 0:
                return "SAT", self.assign

            xi = sorted(cand, key=lambda v: str(v))[0] # TODO: Do Bland'S rule better
            i = basic[xi]

            if self.assign[xi] < self.lower[xi]:
                cand = [nb for nb in nonbasic
                        if (M[i, self.col_index[nb]] > 0 and self.assign[nb] < self.upper[nb])
                        or (M[i, self.col_index[nb]] < 0 and self.assign[nb] > self.lower[nb])]
                if len(cand) == 0:
                    N_plus = {nb for nb in nonbasic if M[i, self.col_index[nb]] > 0}
                    N_minus = {nb for nb in nonbasic if M[i, self.col_index[nb]] < 0}
                    conflict = set()
                    #conflict |= {nb <= self.upper[nb] for nb in N_plus}
                    #conflict |= {nb >= self.lower[nb] for nb in N_minus}
                    #conflict.add(xi >= self.lower[xi])
                    return "UNSAT", "WIP"#conflict
                xj = sorted(cand, key=lambda v: str(v))[0]
                _debug_internal_state_printer2(xi, xj)
                M = self._pivot_and_update(M, basic, nonbasic, xi, xj, self.lower[xi])

            if self.assign[xi] > self.upper[xi]:
                cand = [nb for nb in nonbasic
                        if (M[i, self.col_index[nb]] < 0 and self.assign[nb] < self.upper[nb])
                        or (M[i, self.col_index[nb]] > 0 and self.assign[nb] > self.lower[nb])]

                if len(cand) == 0:
                    N_plus = {nb for nb in nonbasic if M[i, self.col_index[nb]] > 0}
                    N_minus = {nb for nb in nonbasic if M[i, self.col_index[nb]] < 0}
                    conflict = set()
                    conflict |= {nb <= self.upper[nb] for nb in N_minus}
                    conflict |= {nb >= self.lower[nb] for nb in N_plus}
                    conflict.add(xi <= self.upper[xi])
                    return "UNSAT", conflict
                xj = sorted(cand, key=lambda v: str(v))[0]
                _debug_internal_state_printer2(xi, xj)
                M = self._pivot_and_update(M, basic, nonbasic, xi, xj, self.upper[xi])

    def _pivot_and_update(self, M, basic, nonbasic, xi, xj, v):
        i, j = basic[xi], self.col_index[xj]
        assert M[i, j] != 0
        theta_lhs = (v[0] - self.assign[xi][0])/M[i, j]
        theta_rhs = (v[1] - self.assign[xi][1])/M[i, j]
        self.assign[xi] = v
        self.assign[xj] = self.assign[xj] + (theta_lhs, theta_rhs)
        for xk in basic:
            if xk != xi:
                k = basic[xk]
                akj = M[k, j]
                self.assign[xk] = self.assign[xk] + (akj*theta_lhs, akj*theta_rhs)
        # pivot
        basic[xj] = basic[xi]
        del basic[xi]
        nonbasic.add(xi)
        nonbasic.remove(xj)
        return self._pivot(M, i, j)

    @staticmethod
    def _pivot(M, i, j):
        Mi, Mj, Mij = M[i, :], M[:, j], M[i, j]
        if Mij == 0:
            raise ZeroDivisionError("Tried to pivot about zero-valued entry.")
        A = M.copy()
        A[i, :] = -A[i, :]/Mij
        for row in range(M.shape[0]):
            if row != i:
                A[row, :] = A[row, :] + A[row, j] *A[i, :]


        assert A.rref() == M.rref()
        return A
