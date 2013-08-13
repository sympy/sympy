# -*- coding: utf-8 -*-
from .cartan_type import CartanType
from sympy.core import Basic

class RootSystem(Basic):
    """

    Every simple Lie algebra has a unique root system.

    To find the root system, we first consider the Cartan subalgebra of g,
    which is the maximal abelian subalgebra, and consider the adjoint
    action of g on this subalgebra.  There is a root system associated
    with this action.  Now, a root system over a vector space V is a set
    of finite vectors Φ(called roots), which satisfy:
    1.  The roots span V
    2.  The only scalar multiples of x in Φ are x and -x
    3.  For every x in Φ, the set Φ is closed under reflection
        through the hyperplane perpendicular to x.
    4.  If x and y are roots in Φ, then the projection of y onto
        the line through x is a half-integral multiple of x.

    Now, there is a subset of Φ, which we will call Δ, such that:
    1.  Δ is a basis of V
    2.  Each root x in Φ can be written x = Σ k_y y for y in Δ

    The elements of Δ are called the simple roots.
    Therefore, we see that the simple roots span the root space of a given
    simple Lie algebra.

    References: https://en.wikipedia.org/wiki/Root_system
                Lie Algebras and Representation Theory - Humphreys

    """

    def __new__(cls, cartantype):
        """
        Creates a new RootSystem object.  This method assigns an attribute
        called cartan_type to each instance of a RootSystem object.  When
        an instance of RootSystem is called, it needs an argument, which
        should be an instance of a simple Lie algebra.  We then take the
        CartanType of this argument and set it as the cartan_type attribute
        of the RootSystem instance.
        """
        obj = Basic.__new__(cls, cartantype)
        obj.cartan_type = CartanType(cartantype)
        return obj

    def simple_roots(self):
        """
        This method generates and returns the simple roots of the Lie
        algebra.  The rank of the Lie algebra determines the number of
        simple roots that it has.  This method obtains the rank of the
        Lie algebra, and then uses the simple_root method from the Lie
        algebra classes to generate all the simple roots.

        Example
        ====
        >>> from sympy.liealgebras.root_system import RootSystem
        >>> c = RootSystem("A3")
        >>> roots = c.simple_roots()
        >>> roots
        {1: [1, -1, 0, 0], 2: [0, 1, -1, 0], 3: [0, 0, 1, -1]}
        """
        n = self.cartan_type.rank()
        roots = {}
        for i in range(1, n+1):
            root = self.cartan_type.simple_root(i)
            roots[i] = root
        return roots


    def root_space(self):
        """
        The root space is the vector space spanned by the
        simple roots, i.e. it is a vector space with a
        distinguished basis, the simple roots.
        """
        n = self.cartan_type.rank()
        rs = " + ".join("alpha["+str(i) +"]" for i in range(1, n+1))
        return rs

    def cartan_matrix(self):
        """
        Return the Cartan matrix of Lie algebra associated
        with this root system.
        """
        return self.cartan_type.cartan_matrix()

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram of the Lie algebra
        associated with this root system.
        """
        return self.cartan_type.dynkin_diagram()
