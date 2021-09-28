"""Modules in number fields. """

from sympy import sympify
from sympy.core.power import binpow
from sympy.polys.domains import GF, QQ, ZZ
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.numberfields.utilities import is_int


class Module:
    """
    Abstract base class for Modules in number fields.
    Do not instantiate directly.
    """

    def __init__(self, W):
        self.W = W
        self.m = W.rows
        self.n = W.cols

    @property
    def matrix(self):
        return self.W

    @property
    def mult_tab(self):
        raise NotImplementedError

    def __call__(self, v):
        """
        Cast a list v of coefficients as a ModuleElement over this basis.

        Alternatively: v may be an integer 0 <= v < self.n, which we interpret
        as an elementary basis vector with 1 in the vth entry.
        """
        if isinstance(v, int):
            z = [0] * self.n
            z[v] = 1
            v = z
        return ModuleElement(self, v)

    def basis_elements(self):
        return [self(j) for j in range(self.n)]

    def standard_rep(self, elt):
        """
        Given a ModuleElement belonging to this module, recast it as an
        algebraic number, given by a StandardRep.
        """
        raise NotImplementedError

    def with_unity(self):
        """
        Say whether the first generator of this module equals unity.
        """
        return self(0) == 1

    def submodule(self, B, mult_closed=True):
        n = B.cols
        gens = [self(B.col(j).flat()) for j in range(n)]
        M = {}
        # FIXME:
        #  Should `is_square` actually be a condition here?
        #  Or are we supposed to be using sth like A2.3.5 for non-square cases?
        if B.is_square and mult_closed:
            for u in range(n):
                M[u] = {}
                for v in range(u, n):
                    h = gens[u] * gens[v]
                    # TODO: user should be able to control domain where we solve?

                    #x = B.solve(Matrix([h.coeffs]).T)
                    x = B.solve(h.column())

                    M[u][v] = x.flat()
        return Submodule(B, self, M)


class ModuleWithDenominator(Module):

    def __init__(self, W, d):
        super().__init__(W)
        self.d = d


class Submodule(Module):

    def __init__(self, W, container, mult_tab):
        super().__init__(W)
        self.container = container
        self._mult_tab = mult_tab

    @property
    def mult_tab(self):
        return self._mult_tab

    def discard_before(self, r):
        """
        Produce a new module by discarding all generators before a given index r.
        """
        W = self.W[:, r:]
        s = self.n - r
        mt = self.mult_tab
        M = {}
        for u in range(s):
            M[u] = {}
            for v in range(u, s):
                M[u][v] = mt[r + u][r + v][r:]
        return Submodule(W, self.container, M)

    def recast_in_container(self, elt):
        """
        Recast a ModuleElement belonging to this submodule as a ModuleElement
        belonging to this submodule's containing module.
        """
        return ModuleElement.from_column(self.container, self.W * elt.column())

    def standard_rep(self, elt):
        container_elt = self.recast_in_container(elt)
        return self.container.standard_rep(container_elt)

    def represent(self, elt, modulus=None):
        """
        Given a ModuleElement that formally belongs either to this submodule
        or its container, represent it as a linear combination over this
        submodule's basis.

        Parameters
        ----------
        elt: ModuleElement instance, whose `module` property should
          equal this module, or its `container` property.
        modulus: (optional) If a rational prime p is provided, we interpret our
          matrix of coefficients as being over GF(p), and solve over this field.
          Otherwise, since we need to solve over a field, we work over QQ, but
          require that the final coefficients be in ZZ (and cast them there).

        Returns
        -------
        List of coefficients representing the element over this module's basis.

        Raises
        ------
        ValueError if the linear system cannot be solved in the required domain.
        """
        from sympy.matrices.common import NonInvertibleMatrixError
        if elt.module == self:
            return elt.coeffs
        # TODO: Generalize. Can recursively support an element belonging to
        #  any container among all our ancestors.
        if elt.module != self.container:
            raise ValueError('Mismatched Submodule and ModuleElement')
        b = elt.column()
        dom = QQ
        if modulus is not None:
            if modulus < 2:  # skip primality test for sake of speed
                raise ValueError('If defined, modulus must be positive prime.')
            dom = GF(modulus)
        dW = self.W.to_domain(dom)
        db = b.to_domain(dom)
        try:
            dx = dW.lu_solve(db)
            x = dx.convert_to(ZZ)
        except (NonInvertibleMatrixError, CoercionFailed):
            raise ValueError('Cannot represent module element in submodule.')
        return self(x.flat())


class ModuleElement:
    """
    Element of a Module.
    """

    def __init__(self, module, coeffs):
        self.module = module
        self.coeffs = coeffs
        self.n = self.module.n

    def column(self):
        from sympy.matrices import Matrix
        return Matrix([self.coeffs]).T

    @classmethod
    def from_column(cls, module, col):
        return cls(module, col.flat())

    @classmethod
    def one(cls, module):
        return cls.from_int(module, 1)

    @classmethod
    def from_int(cls, module, a):
        if module.with_unity():
            return cls(module, [a] + [0] * (module.n - 1))
        return NotImplemented

    def clone(self):
        return ModuleElement(self.module, self.coeffs[:])

    def is_compat(self, other):
        return isinstance(other, ModuleElement) and other.module == self.module

    def standard_rep(self):
        return self.module.standard_rep(self)

    def __eq__(self, a):
        if self.is_compat(a):
            return all(self.coeffs[i] == a.coeffs[i] for i in range(self.n))
        if is_int(a):
            a = sympify(a)
            return self.standard_rep() == a
        return NotImplemented

    def __add__(self, other):
        if self.is_compat(other):
            return ModuleElement(self.module, [a + b for a, b in zip(self.coeffs, other.coeffs)])
        if is_int(other):
            other = sympify(other)
            return self + ModuleElement.from_int(self.module, other)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return self * (-1)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, a):
        if self.is_compat(a):
            M = self.module.mult_tab
            C, A = self.coeffs, a.coeffs
            n = self.module.n
            P = [0] * n
            for u in range(n):
                for v in range(u, n):
                    R = M[u][v]
                    b = C[u] * A[v]
                    if v > u:
                        b += C[v] * A[u]
                    if b != 0:
                        for k in range(n):
                            P[k] += b*R[k]
            return ModuleElement(self.module, P)
        if is_int(a):
            a = sympify(a)
            return ModuleElement(self.module, [a * c for c in self.coeffs])
        return NotImplemented

    __rmul__ = __mul__

    def __pow__(self, a):
        if is_int(a):
            if a < 0:
                return NotImplemented
            elif a == 0:
                return ModuleElement.one(self.module)
            elif a == 1:
                return self.clone()
            else:
                return binpow(self, a)
        return NotImplemented


class EndomorphismRing:

    def __init__(self, domain):
        self.domain = domain

    def inner_endomorphism(self, multiplier):
        return InnerEndomorphism(self, multiplier)

    def represent(self, element, modulus=None):
        raise NotImplementedError


class InnerEndomorphism:

    def __init__(self, domain, multiplier):
        self.domain = domain
        self.multiplier = multiplier

    def matrix(self, modulus=None):
        raise NotImplementedError


class HnfEndomorphismRing(EndomorphismRing):
    """
    End(R) where R is given as an HNF.
    """

    def inner_endomorphism(self, multiplier):
        return HnfInnerEndomorphism(self, multiplier)

    def represent(self, element, modulus=None):
        if isinstance(element, InnerEndomorphism) and element.domain.domain == self.domain:
            M = element.matrix(modulus=modulus)
            # The purpose of this method is to produce a list of coefficients,
            # (which can be included as a single column vector in a matrix
            # representing a linear map). So must transform. The list should
            # reproduce the columns, one after another.
            return M.T.flat()
        return NotImplemented


class HnfInnerEndomorphism(InnerEndomorphism):

    @property
    def hnf(self):
        return self.domain.domain

    def matrix(self, modulus=None):
        from sympy.matrices import Matrix
        alpha = self.multiplier.standard_rep()
        basis = self.hnf.standard_reps()
        M = Matrix()
        for b in basis:
            # Since self.domain is an HNF, it represents an ideal, and thus is
            # closed under multiplication from outside. Therefore its `represent`
            # method should always succeed.
            c = self.hnf.represent(alpha * b, modulus=modulus)
            M = M.row_join(Matrix([c]).T)
        return M


class ModuleHomomorphism:

    def __init__(self, domain, codomain, mapping, modulus=None):
        self.domain = domain
        self.codomain = codomain
        self.mapping = mapping
        self.modulus = modulus

    def matrix(self):
        from sympy.matrices import Matrix
        basis = self.domain.basis_elements()
        M = Matrix()
        for elt in basis:
            c = self.codomain.represent(self.mapping(elt), modulus=self.modulus)
            M = M.row_join(Matrix([c]).T)
        return M

    def kernel(self):
        from sympy.matrices.common import NonInvertibleMatrixError
        M = self.matrix()
        dom = GF(self.modulus) if self.modulus else QQ
        dM = M.to_domain(dom)
        try:
            dK = dM.nullspace()
            K = dK.convert_to(ZZ).to_Matrix().T
        except (NonInvertibleMatrixError, CoercionFailed):
            raise ValueError('Cannot compute kernel in desired domain.')
        return self.domain.submodule(K)


class ModuleEndomorphism(ModuleHomomorphism):

    def __init__(self, domain, mapping, modulus=None):
        super().__init__(domain, domain, mapping, modulus=modulus)
