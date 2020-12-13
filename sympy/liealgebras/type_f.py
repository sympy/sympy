from .cartan_base import Standard_Cartan
from sympy.core.backend import Rational, Matrix


class TypeF(Standard_Cartan):

    def __new__(cls, n):
        if n != 4:
            raise ValueError("n should be 4")
        return Standard_Cartan.__new__(cls, "F", 4)

    def dimension(self):
        """Dimension of the vector space V underlying the Lie algebra

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("F4")
        >>> c.dimension()
        4
        """

        return 4


    def basic_root(self, i, j):
        """Generate roots with 1 in ith position and -1 in jth position

        """

        n = self.rank
        root = [0]*n
        root[i] = 1
        root[j] = -1
        return Matrix([root])

    def simple_root(self, i):
        """The ith simple root of F_4

        Every lie algebra has a unique root system.
        Given a root system Q, there is a subset of the
        roots such that an element of Q is called a
        simple root if it cannot be written as the sum
        of two elements in Q.  If we let D denote the
        set of simple roots, then it is clear that every
        element of Q can be written as a linear combination
        of elements of D with all coefficients non-negative.

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("F4")
        >>> c.simple_root(3)
        Matrix([[0, 0, 0, 1]])

        """

        if i < 3:
            return self.basic_root(i, i+1)
        if i == 3:
            root = [0]*4
            root[-1] = 1
            return Matrix([root])
        if i == 4:
            root = [Rational(-1, 2)]*4
            root[0]*=-1
            return Matrix([root])

    def roots(self):
        """
        Returns the total number of roots for F_4
        """
        return 48

    def basis(self):
        """
        Returns the number of independent generators of F_4
        """
        return 52

    def dynkin_diagram(self):
        diag = "0---0=>=0---0\n"
        diag += "   ".join(str(i) for i in range(1, 5))
        return diag
