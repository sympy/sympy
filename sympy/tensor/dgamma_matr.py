from sympy import Symbol, S, I
from sympy.combinatorics import Permutation
from sympy.tensor.tensor import (TensorIndexType, tensor_indices,
  TensorSymmetry, get_symmetric_group_sgs, TensorType, tensor_mul, Tensor,
  TensAdd, TensorHead, tensorlist_contract_metric)


class GammaMatrices(object):
    """
    Gamma matrices in dimensional regularization
    with G5 naively anticommuting with all gamma matrices (NDR)

    NDR is inconsistent but widely used.
    Notice in particular that the cyclic property of the
    trace in presence of G5 does not hold.

    The original dimensional regularization scheme by 't Hooft and Veltman
    has `G5` anticommuting only with `G(m)` for `m` in 0,1,2,3,
    and commuting for `m > 3`. It is consistent, and it is used to
    compute the chiral anomaly; however it breaks gauge invariance
    involving G5 even when there are no anomalies, which complicates
    a lot maintaining the Ward (or Slavnov-Taylor) identities.

    TODO: introduce the option of using the 't Hooft and Veltman scheme.
    """

    def __init__(self, typ, gctr=4, g5c=I):
        """
        typ  TensorIndexType (Lorentz or Eulidean)

        gct3  `tr(G(m)*G(n)) = gctr*g(m,n)`
        in D dimensions, when D is an integer, the gamma
        matrices are represented with matrices or rank `2**(D//2)`,
        so `gctr = 2**(D//2)`;
        in dimensional regularization one uses usually
        `gctr = 4*g(m,n)`, so we choose `gctr=4` as default.

        g5c  `tr(G(m0)*G(m1)*G(m2)*G(m3)) = gctr*g5c*epsilon(m0,m1,m2,m3)`
        with Lorentz signature `g5c = I`, which is the default;
        in Euclidean space `g5c = 1`

        """
        self.g = typ.metric
        if not typ.dim:
            raise ValueError('Dimension not assigned')
        self.D = typ.dim
        self.typ = typ
        sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        S1 = TensorType([typ], sym1)
        self.G = S1('G', None)
        sym0 = TensorSymmetry(([], [Permutation(1)]))
        S0 = TensorType([], sym0)
        self.Gamma5 = S0('G5', None)
        self.G5 = Tensor(S.One, [self.Gamma5], [], [])
        self.g5c = g5c
        self.epsilon = typ.epsilon
        self.eps_dim = typ.eps_dim
        self.gctr = gctr

    def G5_to_right(self, t):
        """
        move G5 to the right
        """
        if not t.is_Tensor:
            return t
        if t.is_TensAdd:
            a = [self.G5_to_right(x) for x in t.args]
            return TensAdd(*a)
        components = t.components
        if not t.components:
            return t
        ncomps = len(components)
        G = self.G
        Gamma5 = self.Gamma5
        G5 = self.G5
        # [g5,g,g,g5,g,g]
        numG = 0
        vposG5 = []
        for i in range(ncomps):
            if components[i] == Gamma5:
                vposG5.append(i)
                for j in range(i + 1, ncomps):
                    if components[j] == G:
                        numG += 1
        ct = t.coeff if 0 in vposG5 else S.One
        a = t.split()
        a1 = []
        for i in range(ncomps):
            if not i in vposG5:
                a1.append(a[i])
        for i in range(ncomps - 1, -1, -1):
            if components[i] == G:
                break
        if len(vposG5) % 2:
            a1 = a1[:i + 1] + [G5] + a1[i + 1:ncomps]
        else:
            a1 = a1[:i + 1] + a1[i + 1:ncomps]
        t1 = ((-1)**numG*ct)*tensor_mul(*a1)
        return t1

    def match1_gamma(self, t, n):
        #t = t.sorted_components()
        components = t.components
        ncomps = len(components)
        G = self.G
        for i in range(ncomps):
            if not components[i] == G:
                continue
            for j in range(i + 1, ncomps):
                if components[j] == G and abs(j - i) == n + 1 and \
                        (0, 0, i, j) in t.dum:
                    return i, j


    def rule1_gamma(self, t, n, r=None):
        """
        simplify products of gamma matrices `G(m)*G(m_1)...*G(m_n)*G(-m)`

        t   Gamma matrix monomial
        n   apply rule for `G(m)*G(m_1)...*G(m_n)*G(-m)`

        Examples
        ========

        >>> from sympy import Symbol, S
        >>> from sympy.tensor.tensor import (TensorIndexType, tensor_indices,\
            TensorSymmetry, get_symmetric_group_sgs, TensorType, Tensor)
        >>> from sympy.tensor.dgamma_matr import GammaMatrices
        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        >>> S1 = TensorType([Lorentz], sym1)
        >>> GM = GammaMatrices(Lorentz)
        >>> G = GM.G
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> GM.gamma_trace(G(m0)*G(m1)*G(-m0)*G(m3))
        (-4*D + 8)*metric(m1, m3)
        >>> t = G(m1)*G(m0)*G(-m1)*G(-m0)
        >>> GM.rule1_gamma(t, 1)
        (-D + 2)*G(L_0)*G(-L_0)
        """
        if not r:
            r = self.match1_gamma(t, n)
            if not r:
                return t
        i, j = r
        a = t.split()
        za = zip(a, range(len(a)))
        tc = t.coeff if i == 0 else S.One
        a1 = [x for x, y in za if y not in r]
        D = self.D
        if n == 0:
            # G(m)*G(-m) = d
            t1 = tensor_mul(*a1)
            t2 = Tensor(D*t1.coeff*tc, t1.components, t1.free, t1.dum)
            return t2
        if n == 1:
            # G(m)*G(n)*G(-m) = (2 - d)*G(n)
            t1 = tensor_mul(*a1)
            t2 = Tensor((2-D)*t1.coeff*tc, t1.components, t1.free, t1.dum)
            return t2
        if n == 2:
            # G(m)*G(n1)*G(n2)*G(-m) = (d-4)*G(n1)*G(n2) + 4*delta(n1, n2)
            t1 = ((D-4)*tc)*tensor_mul(*a1)
            a2 = a[:i] + a[j + 1:]
            ind1 = a[r[0] + 1].free[0][0]
            ind2 = a[r[1] - 1].free[0][0]
            if a2:
                a2 = tensorlist_contract_metric(a2, self.g(ind1, ind2))
                t2 = (4*tc)*tensor_mul(*a2)
            else:
                t2 = (4*tc)*self.g(ind1, ind2)
            return t1 + t2
        if n == 3:
            # G(m)*G(n1)*G(n2)*G(n3)*G(-m) = (4-d)*G(n1)*G(n2)*G(n3) - \
            #  2*G(n3)*G(n2)*G(n1)
            t1 = ((4-D)*tc)*tensor_mul(*a1)
            a2a = a[i + 1:j]
            a2a.reverse()
            a2 = a[:i] + a2a + a[j + 1:]
            t2 = (-2*tc)*tensor_mul(*a2[:])
            return t1 + t2
        if n > 3:
            # G(m0)*G(m_1)*...*G(m_k)*G(-m0) =
            # 2*G(m_k)*G(m_1)*...*G(m_(k-1)) -
            # G(m0)*G(m_1)*...*G(m_(k-1))G(-m0)*G(m_k))
            a1 = a[:i] + [a[j-1]] + a[i+1:j-1] + a[j + 1:]
            t1 = (2*tc)*tensor_mul(*a1)
            a2 = a[:j-1] + [a[j], a[j-1]] + a[j + 1:]
            t2 = -tensor_mul(*a2)
            return t1 + t2


    def do_rule1_gamma(self, t, nmax=4, doall=False):
        """
        simplify products of gamma matrices `G(m)*G(m_1)...*G(m_n)*G(-m)`

        `t` Gamma matrix expression

        nmax  maximum number for which `rule1\_gamma(t, n)` is applied

        doall  if true apply `rule1\_gamma` till the expression does not change

        Examples
        ========

        >>> from sympy import Symbol, S
        >>> from sympy.tensor.tensor import (TensorIndexType, tensor_indices,\
            TensorSymmetry, get_symmetric_group_sgs, TensorType, Tensor)
        >>> from sympy.tensor.dgamma_matr import GammaMatrices
        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        >>> S1 = TensorType([Lorentz], sym1)
        >>> GM = GammaMatrices(Lorentz)
        >>> G = GM.G
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> t = G(m1)*G(m0)*G(-m1)*G(-m0)
        >>> GM.do_rule1_gamma(t, doall=True)
        D*(-D + 2)
        """
        if not t.is_Tensor:
            return t
        if t.is_TensMul:
            for n in range(nmax + 1):
                #print 'DB10 t=%s n=%d' %(t, n)
                r = self.match1_gamma(t, n)
                if not r:
                    continue
                t1 = self.rule1_gamma(t, n, r)
                t1 = t1.contract_metric(self.g, contract_all=True)
                if doall:
                    if t1 != t:
                        t1 = self.do_rule1_gamma(t1, nmax, doall)
                return t1
            return t
        if t.is_TensAdd:
            a = [self.do_rule1_gamma(tx, nmax, doall) for tx in t.args]
            t1 = TensAdd(*a)
            return t1.contract_metric(self.g, contract_all=True)

    def match2_gamma(self, t, n):
        components = t.components
        ncomps = len(components)
        G = self.G
        # G(a)*G(b)*...*p(-a)*p(-b)
        for i in range(ncomps - 1 - n):
            if not components[i] == G:
                continue
            if not components[i + n + 1] == G:
                continue
            p_pos1 = p_pos2 = 0
            for dx in t.dum:
                if dx[2] == i:
                    p_pos1 = dx[3]
                if dx[2] == i + n + 1:
                    p_pos2 = dx[3]
            if p_pos1 > 0 and p_pos2 > 0:
                if components[p_pos1] == components[p_pos2] \
                    and components[p_pos1].commuting == 0:
                    return i, p_pos1, p_pos2
        return None

    def rule2_gamma(self, t, n):
        """
        simplify `G(m0)*p(-m0)*G(i_1)...*G(i_n)*G(m1)*p(-m1)`

        `t`    Gamma matrix monomial
        `n` as in  `G(m0)*p(-m0)*G(i_1)...*G(i_n)*G(m1)*p(-m1)`

        Examples
        ========

        >>> from sympy import Symbol, S
        >>> from sympy.tensor.tensor import (TensorIndexType, tensor_indices,\
            TensorSymmetry, get_symmetric_group_sgs, TensorType, Tensor)
        >>> from sympy.tensor.dgamma_matr import GammaMatrices
        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        >>> S1 = TensorType([Lorentz], sym1)
        >>> GM = GammaMatrices(Lorentz)
        >>> G = GM.G
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> p = S1('p')
        >>> ps = G(m0)*p(-m0)
        >>> t = G(m0)*ps*ps*G(-m0)
        >>> GM.rule2_gamma(t, 0)
        G(L_0)*G(-L_0)*p(L_1)*p(-L_1)
        """
        t = t.canon_bp()
        r = self.match2_gamma(t, n)
        if not r:
            return t
        i0, i1, i2 = r
        if n == 0:
            a = t.split()
            ct = t.coeff if i0 == 0 else S.One
            p = a[i1]
            ind_p = p.free[0][0]
            # G(m0)*G(m1)*p(-m0)*p(-m1) = p(m0)*p(-m0)
            a1 = a[:i0] + a[i0 + 2:i1] + a[i1 + 1:i2] + a[i2 + 1:]
            a2 = [p.substitute_indices((ind_p, -ind_p)), p]
            a3 = a1 + a2
            t = tensor_mul(*a3)*ct
            return t
        elif n == 1:
            # G(m0)*G(ind1)*G(m1)*p(-m0)*p(-m1) =
            # 2*G(m0)*p(-m0)*p(ind1) - G(ind1)*p(m0)*p(-m0)
            a = t.split()
            ind1 = a[i0 + 1].free[0][0]
            ct = t.coeff if i0 == 0 else S.One
            p = a[i1]
            ind_p = p.free[0][0]
            a1 = a[:i0] + [a[i0 + 1]] + a[i0 + 3:i1] + a[i1 + 1:i2] + a[i2 + 1:]
            a2 = [p.substitute_indices((ind_p, -ind_p)), p]
            a3 = a1 + a2
            t1 = tensor_mul(*a3)*(-ct)
            args = [t1]
            a1 = a[:i0 + 1] + a[i0 + 3:i2] + a[i2 + 1:]
            p = a[i2]
            ind_p = p.free[0][0]
            a2 = [p.substitute_indices((ind_p, ind1))]
            a3 = a1 + a2
            t2 = tensor_mul(*a3)*2
            t3 = t1 + t2
            return t3
        else:
            raise NotImplementedError

    def do_rule2_gamma(self, t, nmax=1):
        """
        simplify `G(m0)*p(-m0)*G(m1)*p(-m1)` to `p(m)*p(-m)`

        `t`    Gamma matrix expression

        Examples
        ========

        >>> from sympy import Symbol, S
        >>> from sympy.tensor.tensor import (TensorIndexType, tensor_indices,\
            TensorSymmetry, get_symmetric_group_sgs, TensorType, Tensor)
        >>> from sympy.tensor.dgamma_matr import GammaMatrices
        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        >>> S1 = TensorType([Lorentz], sym1)
        >>> GM = GammaMatrices(Lorentz)
        >>> G = GM.G
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> p = S1('p')
        >>> ps = G(m0)*p(-m0)
        >>> M = Symbol('M')
        >>> t = G(m0)*(ps + M)*(ps - M)*G(-m0)
        >>> GM.do_rule2_gamma(t)
        (-M**2)*G(L_0)*G(-L_0) + G(L_0)*G(-L_0)*p(L_1)*p(-L_1)
        """
        if not t.is_Tensor:
            return t
        if t.is_TensMul:
            for n in range(nmax + 1):
                t = self.rule2_gamma(t, n)
        if t.is_TensAdd:
            a = [self.do_rule2_gamma(tx) for tx in t.args]
            t = TensAdd(*a)
        return t


    def gamma_trace(self, t):
        """
        Trace of gamma matrices

        Examples
        ========

        >>> from sympy import Symbol, S
        >>> from sympy.tensor.tensor import (TensorIndexType, tensor_indices,\
                TensorSymmetry, get_symmetric_group_sgs, TensorType, Tensor)
        >>> from sympy.tensor.dgamma_matr import GammaMatrices
        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        >>> S1 = TensorType([Lorentz], sym1)
        >>> GM = GammaMatrices( Lorentz)
        >>> gamma_trace = GM.gamma_trace
        >>> G = GM.G
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> gamma_trace(G(m0)*G(m1)*G(-m0)*G(m3))
        (-4*D + 8)*metric(m1, m3)

        >>> p, q = S1('p,q')
        >>> t = G(m0)*G(m1)*G(m2)*G(m3)*p(-m1)*q(-m3)
        >>> gamma_trace(t)
        -4*metric(m0, m2)*p(L_0)*q(-L_0) + 4*p(m0)*q(m2) + 4*p(m2)*q(m0)
        """
        if not t.is_Tensor:
            return self.gctr*t
        t = self.G5_to_right(t)
        if t.is_TensAdd:
            a = [self.gamma_trace(x) for x in t.args]
            return TensAdd(*a)
        components = t.components
        ncomps = len(components)
        G = self.G
        G5 = self.G5
        Gamma5 = self.Gamma5
        g = self.g
        if all(x != G for x in components):
            if any(x == Gamma5 for x in components):
                return Tensor(S.Zero, [], [], [])
            return self.gctr*t
        withG5 = False
        t = t.canon_bp()
        components = t.components
        for i in range(ncomps):
            if components[i] == G:
                break
        for j in range(i + 1, ncomps):
            if not components[j] == G:
                break
        else:
            j = ncomps
        numG = j - i
        if any(x == Gamma5 for x in components):
            withG5 = True
        if withG5:
            if numG < 4 or numG % 2 != self.eps_dim % 2:
                return Tensor(S.Zero, [], [], [])
            if numG == 4:
                a = t.split()
                indices = [x.free[0][0] for x in a[i:j]]
                ct = t.coeff if i == 0 else S.One
                t1 = (self.gctr*ct*self.g5c)*self.epsilon(*indices)
                a2 = a[:i] + a[j + 1:]
                t2 = tensor_mul(*a2)
                res = t1*t2
                res = res.canon_bp()
                return res
            if numG > 4:
                prev = Tensor(S.Zero, [], [], [])
                t2 = t
                while True:
                    t2 = self.do_rule1_gamma(t2)
                    t2 = t2.contract_metric(g, contract_all=True)
                    t2 = self.do_rule2_gamma(t2)
                    if t2 == prev:
                        break
                    prev = t2
                if t == t2:
                    a = t.split()
                    args = []
                    for k1 in range(i, j):
                        ind1 = a[k1].free[0][0]
                        for k2 in range(k1 + 1, j):
                            ind2 = a[k2].free[0][0]
                            sign = 1 if (k1 - k2) % 2 == 1 else -1
                            a1 = a[:k1] + a[k1 + 1:k2] + a[k2 + 1:]
                            a1 = tensorlist_contract_metric(a1, g(ind1, ind2))
                            ct = t.coeff if k1 == 0 else S.One
                            t1 = (sign*ct)*tensor_mul(*a1)
                            args.append(t1)
                    t2 = TensAdd(*args)
                    return self.gamma_trace(t2)
                else:
                    return self.gamma_trace(t2)
        else:
            # without G5
            if numG % 2 == 1:
                return Tensor(S.Zero, [], [], [])
            if numG > 4:
                prev = Tensor(S.Zero, [], [], [])
                t2 = t
                while True:
                    t2 = self.do_rule1_gamma(t2)
                    t2 = t2.contract_metric(g, contract_all=True)
                    t2 = self.do_rule2_gamma(t2)
                    if t2 == prev:
                        break
                    prev = t2
                if t == t2:
                    a = t.split()
                    # G are from i to j1 included; anticommute G(j1) through
                    ind1 = a[i].free[0][0]
                    ind2 = a[i + 1].free[0][0]
                    aa = a[:i] + a[i + 2:]
                    aa = tensorlist_contract_metric(aa, g(ind1, ind2))
                    ct = t.coeff if i == 0 else S.One
                    t1 = ct*tensor_mul(*aa)
                    args = [t1]
                    sign = 1
                    for k in range(i + 2, j):
                        sign = -sign
                        ind2 = a[k].free[0][0]
                        aa = a[:i] + a[i + 1:k] + a[k + 1:]
                        aa = tensorlist_contract_metric(aa, g(ind1, ind2))
                        t2 = sign*tensor_mul(*aa)
                        args.append(t2)
                    t3 = TensAdd(*args)
                    return self.gamma_trace(t3)
                return self.gamma_trace(t2)

            a = t.split()
            t1 = self.gamma_trace1(*a[i:j])
            a2 = a[:i] + a[j:]
            t2 = tensor_mul(*a2)
            ct = t.coeff if i == 0 else S.One
            t3 = ct*t1*t2
            if not t3:
                return t3
            t3 = t3.contract_metric(g, contract_all=True)
            return t3


    def gamma_trace1(self, *a):
        if not a:
            return self.gctr
        n = len(a)
        if n%2 == 1:
            return Tensor(S.Zero, [], [], [])
        if n == 2:
            ind0 = a[0].free[0][0]
            ind1 = a[1].free[0][0]
            return self.gctr*self.g(ind0, ind1)
        if n == 4:
            ind0 = a[0].free[0][0]
            ind1 = a[1].free[0][0]
            ind2 = a[2].free[0][0]
            ind3 = a[3].free[0][0]
            return self.gctr*(self.g(ind0, ind1)*self.g(ind2, ind3) - \
                self.g(ind0, ind2)*self.g(ind1, ind3) + self.g(ind0, ind3)*self.g(ind1, ind2))
        else:
            raise NotImplementedError
