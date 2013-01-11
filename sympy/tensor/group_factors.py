from sympy import I, S
from sympy.tensor.tensor import (TensorIndexType, tensor_indices, \
TensorSymmetry, get_symmetric_group_sgs, TensorType, TensMul, tensorlist_contract_metric, tensor_mul, TensAdd, TensExpr, tensorsymmetry)



def match_CijCij(t, th):
    """
    match for two tensors components ``th`` with two or more indices contracted
    """
    components = t._components
    rdum = [tuple(sorted([cpos1, cpos2])) for ipos1, ipos2, cpos1, cpos2 \
            in t._dum if components[cpos1] == th and components[cpos2] == th]
    rdum.sort()
    vcount = [rdum.count(x) for x in rdum]
    repeated = {}
    for i in range(len(rdum)):
        if vcount[i] > 1 and rdum[i] not in repeated:
            repeated[rdum[i]] = vcount[i]
    if not repeated:
        return None
    repeated = repeated.items()
    res = max(repeated, key=lambda x: x[1])
    return res

def match_C_triangle(t, th):
    """
    match for three tensor components ``th`` with triangle of contractions
    """
    components = t._components
    rdum = [sorted([cpos1, cpos2]) for ipos1, ipos2, cpos1, cpos2 in \
            t._dum if components[cpos1] == th and components[cpos2] == th]
    rdum.sort()
    for cpos1, cpos2 in rdum:
        for i in range(cpos1 + 1, len(rdum)):
            if [cpos1, i] in rdum:
                if [cpos2, i] in rdum:
                    return cpos1, cpos2, i
                elif [i, cpos2] in rdum:
                    return cpos1, i, cpos2

def match_C_square(t, th):
    """
    match for four tensor components ``th`` with square of contractions

    It is assumed that there are not smaller cycles
    """
    components = t._components
    rdum = [sorted([cpos1, cpos2]) for ipos1, ipos2, cpos1, cpos2 in \
            t._dum if components[cpos1] == th and components[cpos2] == th]
    rdum.sort()
    for i in range(len(rdum)):
        cpos1, cpos2 = rdum[i]
        for j in range(i + 1, len(rdum)):
            cpos3, cpos4 = rdum[j]
            if cpos3 in rdum[i] or cpos4 in rdum[i]:
                continue
            if [cpos1, cpos3] in rdum:
                if [cpos2, cpos4] in rdum or [cpos4, cpos2] in rdum:
                    return cpos1, cpos2, cpos4, cpos3
            if [cpos1, cpos4] in rdum:
                if [cpos2, cpos3] in rdum or [cpos3, cpos2] in rdum:
                    return cpos1, cpos2, cpos3, cpos4


class SuNGroupFactors(object):
    """
    SU(N) generators in any representation satisfy the
    Lie algebra ``[T(i), T(j)] = I*C(i,j, -k)*T(k)``

    ``tr(rep, T(i)*T(j)) = a_rep*g(i, j)`` where
    ``a_rep`` depends on the representation.

    ``g(i, j)`` is the metric in the adjoint representation space.

    For ``N > 2`` there is no metric in the defining representation.

    Let ``T(i, -a, b)`` be the generator in the defining representation ``N``
    normalized with ``a_rep = 1``

    Define ``theta(i_1, ..., i_n) = tr(T(i_1)*...*T^(i_n))

    Multiplying the Lie algebra relation by ``T(-k)`` and tracing one gets
    ``C_(i j k) = -I*theta(i, j, k) + I*theta(j, i, k)``

    For ``SU(N)`` one has the relation
    ``T(i,-d,a)*T(-i,-b,c) =``
    ``delta(-b,a)*delta(-d,c) - 1/N*delta(-d,a)*delta(-b,c)``

    Using this relation one obtains easily relations among contracted
    ``theta`` tensors.
    """
    def __init__(self, N):
        self.N = N
        self.DSuN = TensorIndexType('DSuN', metric=None, dim=N, dummy_fmt='D')
        self.ASuN = TensorIndexType('ASuN', metric=0, dim = N**2-1, dummy_fmt='A')
        self.delta = self.DSuN.delta
        self.g = self.ASuN.metric
        symT = tensorsymmetry(*[[1]]*3)
        ST = TensorType([self.ASuN, self.DSuN, self.DSuN], symT)
        self.T = ST('T')
        symA = tensorsymmetry([3])
        SC = TensorType([self.ASuN]*3, symA)
        self.C = SC('C')
        self._theta = [None]*3 + [TensorType([self.ASuN]*i, \
            tensorsymmetry(*[[1]]*i))('theta') for i in range(3, 6)]
        i0, i1, i2 = tensor_indices('i0,i1,i2', self.ASuN)
        self.tC = self.C(i0,i1,i2)
        self.tCs = -I*self.theta(i0, i1, i2) + I*self.theta(i1, i0, i2)

    def theta(self, *indices):
        n = len(indices)
        if n == 1:
            return S.Zero
        if n == 2:
            return self.g(*indices)
        try:
            return self._theta[n](*indices)
        except IndexError:
            for i in range(len(self._theta), n + 1):
                self._theta.append(TensorType([self.ASuN]*i, tensorsymmetry(*[[1]]*i))('theta'))
            return self._theta[n](*indices)

    def match_thetaii(self, t):
        """
        match a ``theta`` with two indices contracted

        Returns the slot position of the contracted indices and the position
        of the ``theta``.
        """
        components = t._components
        dum = t._dum
        for ipos1, ipos2, cpos1, cpos2 in dum:
            if cpos1 == cpos2 and \
              components[cpos1].name == 'theta' and \
              components[cpos1].types == [self.ASuN]:
                return ipos1, ipos2, cpos1

    def match_thetaithetai(self, t):
        """
        match for two component tensors ``theta`` with
        an index contracted with each other

        Return the positions of the contracted indices and the positions
        of the two tensors.
        """
        components = t._components
        dum = t._dum
        for ipos1, ipos2, cpos1, cpos2 in dum:
            if cpos1 != cpos2 and \
               components[cpos1].name == 'theta' and \
               components[cpos2].name == 'theta' and \
               components[cpos1].types == components[cpos1].types == [self.ASuN]:
               return ipos1, ipos2, cpos1, cpos2


    def rule_C_triangle(self, t):
        """
        evaluate a triangle of C's
        """
        if isinstance(t, TensAdd):
            args = t.args
            args = [self.rule_C_triangle(x) for x in args]
            return TensAdd(*args)
        r = match_C_triangle(t, self.C)
        if not r:
            return t
        i0, i1, i2 = sorted(r)
        a = t.split()
        a1 = a[:i0] + a[i0 + 1: i1] + a[1 + 1:i2] + a[i2 + 1:]
        t1 = tensor_mul(*a1)
        t2 = tensor_mul(a[i0], a[i1], a[i2])
        t2 = self.rule_C2theta_all(t2)
        t = t1*t2
        prev = S.Zero
        while t != prev:
            prev = t
            t = self.rule_thetaithetai(t)
            t = t.contract_metric(self.g, True)
            t = self.rule_thetaii(t)
        return t

    def rule_C_square(self, t):
        """
        evaluate a square of C's
        """
        if isinstance(t, TensAdd):
            args = t.args
            args = [self.rule_C_square(x) for x in args]
            return TensAdd(*args)
        r = match_C_square(t, self.C)
        if not r:
            return t
        i0, i1, i2, i3 = sorted(r)
        a = t.split()
        a1 = a[:i0] + a[i0 + 1: i1] + a[1 + 1:i2] + a[i2 + 1:i3] + a[i3 + 1:]
        t1 = tensor_mul(*a1)
        t2 = tensor_mul(a[i0], a[i1], a[i2], a[i3])
        t2 = self.rule_C2theta_all(t2)
        t = t1*t2
        prev = S.Zero
        while t != prev:
            prev = t
            t = self.rule_thetaithetai(t)
            t = t.contract_metric(self.g, True)
            t = self.rule_thetaii(t)
        return t

    def rule_thetaii(self, t):
        """
        theta(k,i_0,..,i_m,-k,j_0,..,j_n) =
        theta(i_0,..,i_m)*theta(j_0,..,j_n) - 1/N*theta(i_0,..,i_m,j_0,..,j_n)

        """
        if isinstance(t, TensAdd):
            args = t.args
            args = [self.rule_thetaii(x) for x in args]
            return TensAdd(*args)
        r = self.match_thetaii(t)
        if not r:
            return t
        ipos1, ipos2, cpos1 = r
        a = t.split()
        ct = t._coeff if cpos1 == 0 else S.One
        a1 = a[:cpos1] + a[cpos1 + 1:]
        t1 = a[cpos1]
        indices = t1.get_indices()
        if ipos1 < ipos2:
            indices = indices[ipos1:] + indices[:ipos1]
            i = ipos2 - ipos1
        else:
            indices = indices[ipos2:] + indices[:ipos2]
            i = ipos1 - ipos2
        indices1 = indices[1:i]
        indices2 = indices[i + 1:]
        theta = self.theta
        N = self.N
        if indices1 == []:
            t2 = (N - 1/N)*theta(*indices2)
        elif indices2 == []:
            t2 = (N - 1/N)*theta(*indices1)
        else:
            t2 = theta(*indices1)*theta(*indices2) - 1/N*theta(*(indices1 + indices2))
        t3 = tensor_mul(*a1)
        return t2*t3*ct

    def rule_thetaithetai(self, t):
        """
        theta(k,i_0,..,i_m)*theta(k,j_0,...,j_n) =
        theta(i_0,...,i_m,j_0,...,j_n) - 1/N*theta(i_0,..,i_m)*theta(j_0,...,j_n)
        """
        if isinstance(t, TensAdd):
            args = t.args
            args = [self.rule_thetaithetai(x) for x in args]
            return TensAdd(*args)

        r = self.match_thetaithetai(t)
        if not r:
            return t
        ipos1, ipos2, cpos1, cpos2 = r
        a = t.split()
        if cpos1 < cpos2:
            a1 = a[:cpos1] + a[cpos1 + 1:cpos2] + a[cpos2 + 1:]
        else:
            a1 = a[:cpos2] + a[cpos2 + 1:cpos1] + a[cpos1 + 1:]
        if cpos1 == 0:
            ct = t._coeff
        elif cpos2 == 0:
            ct = t._coeff
        else:
            ct = S.One
        t1 = a[cpos1]
        t2 = a[cpos2]
        indices1 = t1.get_indices()
        indices2 = t2.get_indices()
        indices1 = indices1[ipos1 + 1:] + indices1[:ipos1]
        indices2 = indices2[ipos2 + 1:] + indices2[:ipos2]
        theta = self.theta
        N = self.N
        t3 = theta(*(indices1 + indices2))
        t4 = theta(*indices1)*theta(*indices2)
        t5 = theta(*(indices1 + indices2)) - 1/N*theta(*indices1)*theta(*indices2)
        res = tensor_mul(*a1)*t5*ct
        return res

    def rule_C2theta_all(self, t):
        """
        convert all ``C`` to ``theta`` and apply simplifying rules.

        Examples
        ========

        >>> from sympy.tensor.group_factors import SuNGroupFactors
        >>> from sympy import Symbol
        >>> from sympy.tensor.tensor import tensor_indices
        >>> N = Symbol('N')
        >>> SN = SuNGroupFactors(N)
        >>> C = SN.C
        >>> i0,i1,i2,i3,i4,i5,i6 = tensor_indices('i0:7', SN.ASuN)
        >>> t = C(i0,i1,i2)*C(-i0,-i2,i3)
        >>> t = SN.rule_C2theta_all(t)
        >>> SN.rule_C2theta_all(t)
        (-2*N)*metric(i1, i3)
        """

        tC = self.tC
        tCs = self.tCs
        t = t.substitute_tensor(tC, tCs, True)
        prev = S.Zero
        while t != prev:
            prev = t
            t = t.substitute_tensor(tC, tCs, True)
            t = self.rule_thetaithetai(t)
            t = t.contract_metric(self.g, True)
            t = self.rule_thetaii(t)
        return t
