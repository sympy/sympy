from sympy.solvers.solveset import linear_eq_to_matrix
from sympy.matrices.dense import eye
from sympy.assumptions.relation.binrel import AppliedBinaryRelation
from sympy.assumptions.ask import Q
from sympy.core import Dummy
from sympy.core.mul import Mul
from sympy.core.add import Add
from sympy.core.relational import Eq
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


class Boundry:
    """
    Represents an upper or lower bound or an equality between a symbol and some constant.
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

    def get_inequality(self):
        if self.equality:
            return Eq(self.var, self.bound)
        elif self.upper and self.strict:
            return self.var < self.bound
        elif not self.upper and self.strict:
            return self.var > self.bound
        elif self.upper:
            return self.var <= self.bound
        else:
            return self.var >= self.bound

    def __eq__(self, other):
        other = (other.var, other.bound, other.strict, other.upper, other.equality)
        return (self.var, self.bound, self.strict, self.upper, self.equality) == other

    def __hash__(self):
        return hash((self.var, self.bound, self.strict, self.upper, self.equality))



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
        self.boundry_rev_enc = {value: key for key, value in boundry_enc.items()}

        self.A = A # TODO: Row reduce A
        self.slack = slack_variables
        self.nonslack = nonslack_variables
        self.all_var = nonslack_variables + slack_variables
        self.col_index = {v: i for i, v in enumerate(self.all_var)}

        self.lower = {}
        self.upper = {}
        self.low_origin = {}
        self.up_origin = {}
        self.assign = {}

        self.result = None # always one of: (True, Assignment), (False, ConflictClause), None

        for var in self.all_var:
            self.lower[var] = (-float("inf"), 0)
            self.low_origin[var] = False
            self.upper[var] = (float("inf"), 0)
            self.up_origin[var] = False
            self.assign[var] = (0, 0)

    @staticmethod
    def from_encoded_cnf(encoded_cnf):
        """
        Creates an LRASolver from an EncodedCNF object.

        Example
        -------

        This example comes from the example in section 3 of Dutertre's and de Moura's paper.

        >>> from sympy.core.relational import Eq
        >>> from sympy.matrices.dense import Matrix
        >>> from sympy.assumptions.cnf import CNF, EncodedCNF
        >>> from sympy.assumptions.ask import Q
        >>> from sympy.logic.algorithms.lra_theory import LRASolver
        >>> from sympy.abc import a, x, y, z
        >>> phi = Q.prime(a) & (x >= 0) & ((x + y <= 2) | (x + 2 * y - z >= 6)) & (Eq(x + y, 2) | (x + 2 * y - z > 4))
        >>> cnf = CNF.from_prop(phi)
        >>> enc = EncodedCNF()
        >>> enc.from_cnf(cnf)
        >>> enc.data #doctest: +SKIP
        [{1, 5}, {3}, {4}, {2, 6}]
        >>> enc.encoding #doctest: +SKIP
        {Q.gt(x + 2*y - z, 4): 1,
         Q.le(x + y, 2): 2,
         Q.prime(a): 3,
         Q.ge(x, 0): 4,
         Q.eq(x + y, 2): 5,
         Q.ge(x + 2*y - z, 6): 6}
        >>> lra, x_subs, s_subs = LRASolver.from_encoded_cnf(enc)

        Each nonslack variable gets replaced with a dummy variable. Here x, y, and z get replaced
        with _x1, _x2, _x3 respectively.

        >>> x_subs
        {x: _x1, y: _x2, z: _x3}

        We convert any constraints with multiple nonslack variables into single variable constraints
        by substituting with slack variables. Here x + y gets substituted with _s1 and x + 2 * y - z gets
        substituted with -_s2.

        >>> s_subs
        {_x1 + _x2: _s1, -_x1 - 2*_x2 + _x3: _s2}

        Each row of the matrix A represents an equallity between a slack variable and some nonslack varaibles.

        >>> lra.A
        Matrix([
        [ 1,  1, 0, -1,  0],
        [-1, -2, 1,  0, -1]])

        To make it very clear what those equalities are, we can multiply A by a vector containing each varaiable.
        The result will be a list of quantities that are equal to zero.

        >>> lra.A * Matrix(lra.all_var)
        Matrix([
        [        -_s1 + _x1 + _x2],
        [-_s2 - _x1 - 2*_x2 + _x3]])

        By substituting terms with multiple constraints with slack variables, each constraint in phi
        is transformed into an upper or lower bound or equality between a single variable and some
        constant. Rather than returning a new encoded cnf object with a new encoding, the new lra
        object has its own encoding stored in `lra.boundry_enc` and `lra.boundry_rev_enc`.

        As boundry objects can't be printed nicely, here's what that looks like if the boundries are
        converted into inequalities.

        >>> {key: value.get_inequality() for key, value in lra.boundry_enc.items()} #doctest: +SKIP
        {5: Eq(_s1, 2), 6: _s2 <= -6, 4: _x1 >= 0, 1: _s2 < -4, 2: _s1 <= 2}

        Notice that there are no encodings for 3. This is because predicates such as
        Q.prime which the LRASolver has no understanding of are ignored.
        """

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
        # TODO: get rid of this to speed things up
        ordered_encoded_cnf = sorted(encoded_cnf.encoding.items(), key=lambda x: str(x))

        # check that preprocessing has been done
        # TODO: get rid of this to speed things up
        #assert all(standardize_binrel(prop) == prop for prop, enc in encoding)

        for prop, enc in ordered_encoded_cnf:
            if not isinstance(prop, AppliedBinaryRelation):
                continue
            assert prop.function in [Q.le, Q.ge, Q.eq, Q.gt, Q.lt]

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
                    nonbasic.append(x_subs[term])
                terms.append(term_coeff * x_subs[term])

            # If there are multiple variable terms, replace them with a dummy _si variable.
            # If needed (no other expr has this sum of variable terms), create a new Dummy _si varaible.
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

        A, _ = linear_eq_to_matrix(A, nonbasic + basic)
        return LRASolver(A, basic, nonbasic, boundry_enc=encoding), x_subs, s_subs

    def assert_enc_boundry(self, enc_boundry):
        """
        Assert an upper or lower bound or equality between
        a variable and a constant. Update the state
        accordingly.

        Parameters
        ==========

        enc_boundry : int
        """
        boundry = self.boundry_enc[enc_boundry]
        sym, c = boundry.var, boundry.bound

        if boundry.strict:
            delta = -1 if boundry.upper else 1
            c = (c, delta)
        else:
            c = (c, 0)

        if boundry.equality:
            res1 = self._assert_lower(sym, c,from_equality=True)
            if res1[0] == "UNSAT":
                return res1
            res2 = self._assert_upper(sym, c,from_equality=True)
            if res1[0] == "UNSAT":
                return res2
            if "OK" in [res1[0], res2[0]]:
                return "OK"
            return res2 # SAT
        elif boundry.upper:
            return self._assert_upper(sym, c)
        else:
            return self._assert_lower(sym, c)

    def _assert_upper(self, xi, ci, from_equality=False):
        if self.result:
            assert self.result[0] != "UNSAT"
        self.result = None
        if ci >= self.upper[xi]:
            return "OK", None
        if ci < self.lower[xi]:
            assert (self.lower[xi][1] >= 0) is True
            assert (ci[1] <= 0) is True


            lit1 = Boundry(var=xi, const=self.lower[xi][0], strict=self.lower[xi][1] != 0, upper=False,
                           equality=self.low_origin[xi])
            lit2 = Boundry(var=xi, const=ci[0], strict=ci[1] != 0, upper=True, equality=from_equality)

            conflict = {-self.boundry_rev_enc[lit1], -self.boundry_rev_enc[lit2]}
            self.result = "UNSAT", conflict
            assert lit1 in self.boundry_rev_enc
            assert lit2 in self.boundry_rev_enc
            return self.result
        self.upper[xi] = ci
        self.up_origin[xi] = from_equality
        if xi in self.nonslack and self.assign[xi] > ci:
            self._update(xi, ci)

        return "OK", None

    def _assert_lower(self, xi, ci, from_equality=False):
        if self.result:
            assert self.result[0] != "UNSAT"
        self.result = None
        if ci <= self.lower[xi]:
            return "OK", None
        if ci > self.upper[xi]:
            assert (self.upper[xi][1] <= 0) is True
            assert (ci[1] >= 0) is True

            lit1 = Boundry(var=xi, const=self.upper[xi][0], strict=self.upper[xi][1] != 0, upper=True,
                           equality=self.up_origin[xi])
            lit2 = Boundry(var=xi, const=ci[0], strict=ci[1] != 0, upper=False, equality=from_equality)

            conflict = {-self.boundry_rev_enc[lit1],-self.boundry_rev_enc[lit2]}
            self.result = "UNSAT", conflict
            return self.result
        self.lower[xi] = ci
        self.low_origin[xi] = from_equality
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
        pass

    def backtrack(self):
        pass

    def check(self):
        """
        Searches for an assignment for all variables that satisfies all their
        upper bounds or determines that such an assignment does not exist.
        """
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
                assert all(((self.assign[nb] >= self.lower[nb]) == True) and ((self.assign[nb] <= self.upper[nb]) == True) for nb in nonbasic)

                # assignments for x must always satisfy Ax = 0
                # probably have to turn this off when dealing with strict ineq
                X = Matrix([self.assign[v][0] for v in self.col_index])
                assert all(abs(val) < 10**(-10) for val in M*X)

                # check upper and lower match this format:
                # x <= rat + delta iff x < rat
                # x >= rat - delta iff x > rat
                # this wouldn't make sense:
                # x <= rat - delta
                # x >= rat + delta
                assert all(self.upper[x][1] <= 0 for x in self.upper)
                assert all(self.lower[x][1] >= 0 for x in self.upper)

            cand = [b for b in basic
             if self.assign[b] < self.lower[b]
             or self.assign[b] > self.upper[b]]

            if len(cand) == 0:
                return "SAT", self.assign

            xi = sorted(cand, key=lambda v: str(v))[0] # TODO: Do Bland's rule better
            i = basic[xi]

            if self.assign[xi] < self.lower[xi]:
                cand = [nb for nb in nonbasic
                        if (M[i, self.col_index[nb]] > 0 and self.assign[nb] < self.upper[nb])
                        or (M[i, self.col_index[nb]] < 0 and self.assign[nb] > self.lower[nb])]
                if len(cand) == 0:
                    N_plus = {nb for nb in nonbasic if M[i, self.col_index[nb]] > 0}
                    N_minus = {nb for nb in nonbasic if M[i, self.col_index[nb]] < 0}
                    upper = [(nb, self.upper[nb]) for nb in N_plus]
                    lower = [(nb, self.lower[nb]) for nb in N_minus]

                    conflict = set()
                    conflict |= {Boundry(nb, up[0], True, self.up_origin[nb], up[1] != 0)
                            for nb, up in upper}
                    conflict |= {Boundry(nb, lo[0], False, self.low_origin[nb], lo[1] != 0)
                                 for nb, lo in lower}
                    conflict.add(Boundry(xi, self.lower[xi][0], False, self.low_origin[xi], self.lower[xi][1] != 0))
                    conflict = set(-self.boundry_rev_enc[c] for c in conflict)
                    return "UNSAT", conflict
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
                    upper = [(nb, self.upper[nb]) for nb in N_minus]
                    lower = [(nb, self.lower[nb]) for nb in N_plus]

                    conflict = set()
                    conflict |= {Boundry(nb, up[0], True, self.up_origin[nb], up[1] != 0)
                                 for nb, up in upper}
                    conflict |= {Boundry(nb, lo[0], False, self.low_origin[nb], lo[1] != 0)
                                 for nb, lo in lower}
                    conflict.add(Boundry(xi, self.upper[xi][0], True, self.up_origin[xi], self.upper[xi][1] != 0))

                    conflict = set(-self.boundry_rev_enc[c] for c in conflict)
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
        self.assign[xj] = (self.assign[xj][0] + theta_lhs, self.assign[xj][1] + theta_rhs)
        for xk in basic:
            if xk != xi:
                k = basic[xk]
                akj = M[k, j]
                self.assign[xk] = (self.assign[xk][0] + akj*theta_lhs, self.assign[xk][1] + akj*theta_rhs)
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
        variable j and then substiututing in the rest of the equations
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
        Mi, Mj, Mij = M[i, :], M[:, j], M[i, j]
        if Mij == 0:
            raise ZeroDivisionError("Tried to pivot about zero-valued entry.")
        A = M.copy()
        A[i, :] = -A[i, :]/Mij
        for row in range(M.shape[0]):
            if row != i:
                A[row, :] = A[row, :] + A[row, j] *A[i, :]

        return A
