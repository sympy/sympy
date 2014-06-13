from __future__ import print_function, division

__all__ = ['Linearizer']

from sympy import Matrix, eye, zeros, Dummy
from sympy.utilities.iterables import flatten
from sympy.physics.vector import dynamicsymbols
from sympy.physics.mechanics.functions import _subs_keep_derivs
from collections import namedtuple

class Linearizer(object):
    """This object holds the general model form used as an intermediate
    for the linearize functions."""

    def __init__(self, f_0, f_1, f_2, f_3, f_4, f_c, f_v, f_a, q, u,
            q_ind=None, q_dep=None, u_ind=None, u_dep=None, r=None, lams=None):
        # Generalized equation form
        self.f_0 = f_0
        self.f_1 = f_1
        self.f_2 = f_2
        self.f_3 = f_3
        self.f_4 = f_4
        self.f_c = f_c
        self.f_v = f_v
        self.f_a = f_a

        # Generalized equation variables
        self.q = Matrix(q)
        self.u = Matrix(u)
        none_handler = lambda x: Matrix(x) if x else Matrix()
        self.q_ind = none_handler(q_ind)
        self.q_dep = none_handler(q_dep)
        self.u_ind = none_handler(u_ind)
        self.u_dep = none_handler(u_dep)
        self.r = none_handler(r)
        self.lams = none_handler(lams)

        # Derivatives of generalized equation variables
        t = dynamicsymbols._t
        self.qd = self.q.diff(t)
        self.ud = self.u.diff(t)
        # If the user doesn't actually use generalized variables, and the
        # qd and u vectors have any intersecting variables, this can cause
        # problems. We'll fix this with some hackery, and Dummy variables
        dup_vars = set(self.qd).intersection(self.u)
        self.qd_dup = Matrix([var if var not in dup_vars else Dummy()
            for var in self.qd])

        # Derive dimesion terms
        l = len(self.f_c)
        m = len(self.f_v)
        n = len(self.q)
        o = len(self.u)
        s = len(self.r)
        k = len(self.lams)
        dims = namedtuple('dims', ['l', 'm', 'n', 'o', 's', 'k'])
        self.dims = dims(l, m, n, o, s, k)

        # Calculates here only need to be run once:
        self._form_permutation_matrices()
        self._form_block_matrices()
        self._form_coefficient_matrices()

    def _form_permutation_matrices(self):
        """Form the permutation matrices Pq and Pu."""

        # Extract dimension variables
        l, m, n, o, s, k = self.dims
        # Compute permutation matrices
        if n != 0:
            self.Pq = permutation_matrix(self.q, Matrix([self.q_ind,
                self.q_dep]))
            if l > 0:
                self.Pqi = self.Pq[:, :-l]
                self.Pqd = self.Pq[:, -l:]
            else:
                self.Pqi = self.Pq
                self.Pqd = Matrix([])
        if o != 0:
            self.Pu = permutation_matrix(self.u, Matrix([self.u_ind,
                self.u_dep]))
            if m > 0:
                self.Pui = self.Pu[:, :-m]
                self.Pud = self.Pu[:, -m:]
            else:
                self.Pui = self.Pu
                self.Pud = Matrix([])
        # Compute combination permutation matrix for computing A and B
        P_col1 = Matrix([self.Pqi, zeros(o + k, n - l)])
        P_col2 = Matrix([zeros(n, o - m), self.Pui, zeros(k, o - m)])
        if P_col1:
            if P_col2:
                self.P_prime = P_col1.row_join(P_col2)
            else:
                self.P_prime = P_col1
        else:
            self.P_prime = P_col2

    def _form_coefficient_matrices(self):
        """Form the coefficient matrices C_0, C_1, and C_2."""

        # Extract dimension variables
        l, m, n, o, s, k = self.dims
        # Build up the coefficient matrices C_0, C_1, and C_2
        # If there are configuration constraints (l > 0), form C_0 as normal.
        # If not, C_0 is I_(nxn). Note that this works even if n=0
        if l > 0:
            f_c_jac_q = self.f_c.jacobian(self.q)
            self.C_0 = (eye(n) - self.Pqd * (f_c_jac_q * self.Pqd).inv() *
                    f_c_jac_q) * self.Pqi
        else:
            self.C_0 = eye(n)
        # If there are motion constraints (m > 0), form C_1 and C_2 as normal.
        # If not, C_1 is 0, and C_2 is I_(oxo). Note that this works even if
        # o = 0.
        if m > 0:
            f_v_jac_u = self.f_v.jacobian(self.u)
            temp = self.Pud * (f_v_jac_u * self.Pud).inv()
            if n != 0:
                f_v_jac_q = self.f_v.jacobian(self.q)
                self.C_1 = -temp * f_v_jac_q
            else:
                self.C_1 = 0
            self.C_2 = (eye(o) - temp * f_v_jac_u) * self.Pui
        else:
            self.C_1 = 0
            self.C_2 = eye(o)

    def _form_block_matrices(self):
        """Form the block matrices for composing M, A, and B."""

        # Extract dimension variables
        l, m, n, o, s, k = self.dims
        # Block Matrix Definitions. These are only defined if under certain
        # conditions. If undefined, an empty matrix is used instead
        if n != 0:
            self.M_qq = self.f_0.jacobian(self.qd)
            self.A_qq = -(self.f_0 + self.f_1).jacobian(self.q)
        else:
            self.M_qq = Matrix([])
            self.A_qq = Matrix([])
        if n != 0 and m != 0:
            self.M_uqc = self.f_a.jacobian(self.qd_dup)
            self.A_uqc = -self.f_a.jacobian(self.q)
        else:
            self.M_uqc = Matrix([])
            self.A_uqc = Matrix([])
        if n != 0 and o - m + k != 0:
            self.M_uqd = self.f_3.jacobian(self.qd_dup)
            self.A_uqd = -(self.f_2 + self.f_3 + self.f_4).jacobian(self.q)
        else:
            self.M_uqd = Matrix([])
            self.A_uqd = Matrix([])
        if o != 0 and m != 0:
            self.M_uuc = self.f_a.jacobian(self.ud)
            self.A_uuc = -self.f_a.jacobian(self.u)
        else:
            self.M_uuc = Matrix([])
            self.A_uuc = Matrix([])
        if o != 0 and o - m + k != 0:
            self.M_uud = self.f_2.jacobian(self.ud)
            self.A_uud = -(self.f_2 + self.f_3).jacobian(self.u)
        else:
            self.M_uud = Matrix([])
            self.A_uud = Matrix([])
        if o != 0 and n != 0:
            self.A_qu = -self.f_1.jacobian(self.u)
        else:
            self.A_qu = Matrix([])
        if k != 0 and o - m + k != 0:
            self.M_uld = self.f_4.jacobian(self.lams)
        else:
            self.M_uld = Matrix([])
        if s != 0 and o - m + k != 0:
            self.B_u = -self.f_3.jacobian(self.r)
        else:
            self.B_u = Matrix([])

    def linearize(self, q_op=None, u_op=None, qd_op=None, ud_op=None,
            r_op=None, lam_op=None, A_and_B=False):
        """Linearize the system about the trim conditions. Note that
        q_op, u_op, qd_op, ud_op must satisfy the equations of motion.
        These may be either symbolic or numeric.

        Parameters
        ==========
        q_op, u_op, qd_op, ud_op, r_op, lam_op : dict
            Dictionaries of the trim conditions. These will be substituted in
            to the linearized system before the linearization is complete.
            Leave blank if you want a completely symbolic form. Note that any
            reduction in symbols (whether substituted for numbers or
            expressions with a common parameter) will result in faster runtime.

        A_and_B : bool
            If set to True, A and B for forming dx = [A]x + [B]r will returned,
            where x = [q_indnd, u_indnd]^T. If set to False, M, A, and B for
            forming [M]x = [A]x + [B]r will be returned. Default is False."""

        # Compose dicts of the trim condition for q, u, and qd, ud
        op_compose = lambda kw: dict() if not kw else kw
        q_op = op_compose(q_op)
        u_op = op_compose(u_op)
        qd_op = op_compose(qd_op)
        ud_op = op_compose(ud_op)
        r_op = op_compose(r_op)
        lam_op = op_compose(lam_op)
        trim = q_op
        trim.update(u_op)
        trim.update(r_op)
        trim.update(lam_op)
        dtrim = qd_op
        dtrim.update(ud_op)

        # Extract dimension variables
        l, m, n, o, s, k = self.dims

        # Rename terms to shorten expressions
        M_qq = self.M_qq
        M_uqc = self.M_uqc
        M_uqd = self.M_uqd
        M_uuc = self.M_uuc
        M_uud = self.M_uud
        M_uld = self.M_uld
        A_qq = self.A_qq
        A_uqc = self.A_uqc
        A_uqd = self.A_uqd
        A_qu = self.A_qu
        A_uuc = self.A_uuc
        A_uud = self.A_uud
        B_u = self.B_u
        C_0 = self.C_0
        C_1 = self.C_1
        C_2 = self.C_2

        # Build up Mass Matrix
        #     |M_qq    0_nxo   0_nxk|
        # M = |M_uqc   M_uuc   0_mxk|
        #     |M_uqd   M_uud   M_uld|
        if o != 0:
            col2 = Matrix([zeros(n, o), M_uuc, M_uud])
        if k != 0:
            col3 = Matrix([zeros(n + m, k), M_uld])
        if n != 0:
            col1 = Matrix([M_qq, M_uqc, M_uqd])
            if o != 0 and k != 0:
                M = col1.row_join(col2).row_join(col3)
            elif o != 0:
                M = col1.row_join(col2)
            else:
                M = col1
        elif k != 0:
            M = col2.row_join(col3)
        else:
            M = col2
        M_eq = _subs_keep_derivs(M.subs(dtrim), trim)
        M_eq.simplify()

        # Build up state coefficient matrix A
        #     |(A_qq + A_qu*C_1)*C_0       A_qu*C_2|
        # A = |(A_uqc + A_uuc*C_1)*C_0    A_uuc*C_2|
        #     |(A_uqd + A_uud*C_1)*C_0    A_uud*C_2|
        # Col 1 is only defined if n != 0
        if n != 0:
            r1c1 = A_qq
            if o != 0:
                r1c1 += A_qu*C_1
            r1c1 = r1c1 * C_0
            if m != 0:
                r2c1 = A_uqc
                if o != 0:
                    r2c1 += A_uuc*C_1
                r2c1 = r2c1 * C_0
            else:
                r2c1 = Matrix([])
            if o - m + k != 0:
                r3c1 = A_uqd
                if o != 0:
                    r3c1 += A_uud*C_1
                r3c1 = r3c1 * C_0
            else:
                r3c1 = Matrix([])
            col1 = Matrix([r1c1, r2c1, r3c1])
        else:
            col1 = Matrix([])
        # Col 2 is only defined if o != 0
        if o != 0:
            if n != 0:
                r1c2 = A_qu * C_2
            else:
                r1c2 = Matrix([])
            if m != 0:
                r2c2 = A_uuc * C_2
            else:
                r2c2 = Matrix([])
            if o - m + k != 0:
                r3c2 = A_uud * C_2
            else:
                r3c2 = Matrix([])
            col2 = Matrix([r1c2, r2c2, r3c2])
        else:
            col2 = Matrix([])
        if col1:
            if col2:
                Amat = col1.row_join(col2)
            else:
                Amat = col1
        else:
            Amat = col2
        Amat_eq = _subs_keep_derivs(Amat.subs(dtrim), trim)
        Amat_eq.simplify()

        # Build up the B matrix if there are forcing variables
        #     |0_(n + m)xs|
        # B = |B_u        |
        if s != 0 and o - m + k != 0:
            Bmat = zeros(n + m, s).col_join(B_u)
            Bmat_eq = _subs_keep_derivs(Bmat.subs(dtrim), trim)
        else:
            Bmat_eq = Matrix([])

        # kwarg A_and_B indicates to return  A, B for forming the equation
        # dx = [A]x + [B]r, where x = [q_indnd, u_indnd]^T,
        if A_and_B:
            Minv = M_eq.inv()
            A_cont = self.P_prime.T * Minv * Amat_eq
            A_cont.simplify()
            if Bmat_eq:
                B_cont = self.P_prime.T * Minv * Bmat_eq
                B_cont.simplify()
            else:
                B_cont = Bmat_eq
            return A_cont, B_cont
        # Otherwise return M, A, B for forming the equation
        # [M]dx = [A]x + [B]r, where x = [q_indnd, u_indnd]^T
        else:
            return M_eq, Amat_eq, Bmat_eq


def permutation_matrix(orig_vec, per_vec):
    """Compute the permutation matrix to change order of
    orig_vec into order of per_vec.

    Parameters
    ==========
    orig_vec : Vector or List
        (n x 1) of original symbol ordering

    per_vec : Vector or List
        (n x 1) of desired symbol ordering"""
    if not isinstance(orig_vec, (list, tuple)):
        orig_list = flatten(orig_vec)
    per_list = [orig_list.index(i) for i in per_vec]
    p_matrix = zeros(len(orig_list))
    for i, j in enumerate(per_list):
        p_matrix[i, j] = 1
    return p_matrix
