from __future__ import annotations

from sympy.core.backend import Matrix
from .cartan_type import Standard_Cartan


class TypeB(Standard_Cartan):

    def __new__(cls, n: int):
        if n < 2:
            raise ValueError("n cannot be less than 2")
        return Standard_Cartan.__new__(cls, "B", n)

    def dimension(self) -> int:
        """Dimension of the vector space V underlying the Lie algebra

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("B3")
        >>> c.dimension()
        3
        """

        return self.rank()

    def basic_root(self, i: int, j: int) -> Matrix:
        """
        This is a method just to generate roots
        with a 1 iin the ith position and a -1
        in the jth position.

        """
        root = [0]*self.rank()
        root[i] = 1
        root[j] = -1
        return Matrix([root])

    def simple_root(self, i: int) -> Matrix:
        """
        Every lie algebra has a unique root system.
        Given a root system Q, there is a subset of the
        roots such that an element of Q is called a
        simple root if it cannot be written as the sum
        of two elements in Q.  If we let D denote the
        set of simple roots, then it is clear that every
        element of Q can be written as a linear combination
        of elements of D with all coefficients non-negative.

        In B_n the first n-1 simple roots are the same as the
        roots in A_(n-1) (a 1 in the ith position, a -1 in
        the (i+1)th position, and zeroes elsewhere).  The n-th
        simple root is the root with a 1 in the nth position
        and zeroes elsewhere.

        This method returns the ith simple root for the B series.

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("B3")
        >>> c.simple_root(2)
        Matrix([[0, 1, -1]])

        """
        n = self.rank()
        if i < 0 or i > n:
            raise ValueError(f"Lie algebra B{n} has {n} simple roots")

        if i < n:
            return Matrix(self.basic_root(i-1, i))
        else:
            root = [0]*n
            root[n-1] = 1
            return Matrix([root])

    def roots(self) -> int:
        """
        Returns the total number of roots for B_n"
        """

        n = self.rank()
        return 2*(n**2)

    def basis(self) -> int:
        """
        Returns the number of independent generators of B_n
        """

        n = self.rank()
        return (n**2 - n) // 2

    def lie_algebra(self) -> str:
        """
        Returns the Lie algebra associated with B_n
        """

        n = self.rank()
        return "so(" + str(2*n) + ")"

    def dynkin_diagram(self) -> str:
        n = self.rank()
        diag = "---".join("0" for _ in range(1, n)) + "=>=0\n"
        diag += "   ".join(str(i) for i in range(1, n+1))
        return diag
