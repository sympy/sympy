from __future__ import print_function, division

__all__ = ['Linearizer']

from sympy import symbols, Matrix, eye, zeros, solve, cos, sqrt, Poly, flatten
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point,\
        dot, cross, mprint, inertia

class Linearizer(object):
    """ This object holds the general model form used as an intermediate
    for the linearize functions"""

    def __init__(self, f_0, f_1, f_2, f_3, f_c, f_v, f_a, q, u, q_i=[], q_d=[], u_i=[], u_d=[], r=[]):
        # Generalized equation form
        self.f_0 = f_0
        self.f_1 = f_1
        self.f_2 = f_2
        self.f_3 = f_3
        self.f_c = f_c
        self.f_v = f_v
        self.f_a = f_a

        # Generalized equation variables
        self.q = Matrix(q)
        self.u = Matrix(u)
        self.q_i = Matrix(q_i)
        self.q_d = Matrix(q_d)
        self.u_i = Matrix(u_i)
        self.u_d = Matrix(u_d)
        self.r = Matrix(r)

        # Derivatives of generalized equation variables
        t = symbols('t')
        self.qd = self.q.diff(t)
        self.ud = self.u.diff(t)
        

    def check_setpoint(self, q_trim, u_trim, r_trim):
        pass

    def linearize(self, eq_ud, eq_u, eq_qd, qd_sym, eq_q, qd_sym_inv):
        """ Perform the actual linearization"""

        # Rename terms to shorten expressions
        q = self.q
        u = self.u
        qd = self.qd
        ud = self.ud
        r = self.r

        q_i = self.q_i
        q_d = self.q_d
        u_i = self.u_i
        u_d = self.u_d

        f_0 = self.f_0
        f_1 = self.f_1
        f_2 = self.f_2
        f_3 = self.f_3
        f_a = self.f_a
        f_v = self.f_v
        f_c = self.f_c

        # Dimension terms
        n = len(q)
        o = len(u)
        s = len(r)
        l = len(f_c)
        m = len(f_a)

        # Compute permutation matrices
        Pq = permutation_matrix(q, q_i.col_join(q_d))
        Pqi = Pq[:, :-l]
        Pqd = Pq[:, -l:]
        Pu = permutation_matrix(u, u_i.col_join(u_d))
        Pui = Pu[:, :-m]
        Pud = Pu[:, -m:]

        # Block Matrix Definitions (eq 57)
        M_qq = f_0.jacobian(qd)
        M_uqc = f_a.jacobian(qd)
        M_uuc = f_a.jacobian(ud)
        M_uqd = f_2.jacobian(qd)
        M_uud = f_2.jacobian(ud)
        A_qq = -(f_0 + f_1).jacobian(q)
        A_qu = -f_1.jacobian(u)
        A_uqc = -f_a.jacobian(q)
        A_uuc = -f_a.jacobian(u)
        A_uqd = -(f_2 + f_3).jacobian(q)
        A_uud = -f_3.jacobian(u)

        # Build up Mass Matrix 
        #     |M_qq    0_nxo|
        # M = |M_uqc   M_uuc|
        #     |M_uqd   M_uud|
        row1 = M_qq.row_join(zeros(n, o))
        row2 = M_uqc.row_join(M_uuc)
        row3 = M_uqd.row_join(M_uud)
        M = row1.col_join(row2).col_join(row3)
        M_eq = M.subs(eq_ud).subs(eq_u).subs(eq_qd).subs(qd_sym).subs(eq_q).subs(qd_sym_inv)
        M_eq.simplify()

        # Jacobian of constraint matrices
        f_c_jac_q = f_c.jacobian(q)
        f_v_jac_q = f_v.jacobian(q)
        f_v_jac_u = f_v.jacobian(u)

        # Matrices defined in equations (59) and equation (60)
        C_0 = (eye(n) - Pqd * (f_c_jac_q * Pqd).inv() * f_c_jac_q) * Pqi
        C_1 = -Pud * (f_v_jac_u * Pud).inv() * f_v_jac_q
        C_2 = (eye(o) - Pud * (f_v_jac_u * Pud).inv() * f_v_jac_u) * Pui

        # Build up state coefficient matrix A
        #     |(A_qq + A_qu*C_1)*C_0      A_qu*C_2|
        # A = |(A_uqc + A_uuc*C_1)*C_0    A_uuc*C_2|
        #     |(A_uqd + A_uud*C_1)*C_0    A_uud*C_2|
        row1 = ((A_qq + A_qu * C_1) * C_0).row_join(A_qu * C_2)
        row2 = ((A_uqc + A_uuc * C_1) * C_0).row_join(A_uuc * C_2)
        row3 = ((A_uqd + A_uud * C_1) * C_0).row_join(A_uud * C_2)
        Amat = row1.col_join(row2).col_join(row3)
        Amat_eq = Amat.subs(eq_ud).subs(eq_u).subs(eq_qd).subs(qd_sym).subs(eq_q).subs(qd_sym_inv)
        Amat_eq.simplify()

        # Build up the B matrix if there are forcing variables
        #     |0_(n + m)xs|
        # B = |B_u        |
        if s > 0:
            B_u = -f_3.jacobian(r)
            Bmat = zeros(n + m, s).col_join(B_u)
            Bmat_eq = Bmat.subs(eq_ud).subs(eq_u).subs(eq_qd).subs(qd_sym).subs(eq_q).subs(qd_sym_inv)
        else:
            Bmat_eq = Matrix([])

        return M_eq, Amat_eq, Bmat_eq


def permutation_matrix(orig_vec, per_vec):
    """ Compute the permutation matrix to change order of
    orig_vec into order of per_vec 
    
    Parameters
    ==========
    orig_vec : Vector or List
        (n x 1) of original symbol ordering

    per_vec : Vector or List
        (n x 1) of desired symbol ordering
    """
    if not isinstance(orig_vec, (list, tuple)):
        orig_list = flatten(orig_vec)
    per_list = [orig_list.index(i) for i in per_vec]
    p_matrix = zeros(len(orig_list))
    for i, j in enumerate(per_list):
        p_matrix[i, j] = 1
    return p_matrix










