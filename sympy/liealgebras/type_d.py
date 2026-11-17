from __future__ import annotations

from sympy.core.backend import Matrix
from .cartan_type import Standard_Cartan


class TypeD(Standard_Cartan):

    def __new__(cls, n: int):
        if n < 3:
            raise ValueError("n cannot be less than 3")
        return Standard_Cartan.__new__(cls, "D", n)


    def dimension(self) -> int:
        """Dimension of the vector space V underlying the Lie algebra

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("D4")
        >>> c.dimension()
        4
        """

        return self.n

    def basic_root(self, i: int, j: int) -> Matrix:
        """
        This is a method just to generate roots
        with a 1 iin the ith position and a -1
        in the jth position.

        """

        n = self.n
        root = [0]*n
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

        In D_n, the first n-1 simple roots are the same as
        the roots in A_(n-1) (a 1 in the ith position, a -1
        in the (i+1)th position, and zeroes elsewhere).
        The nth simple root is the root in which there 1s in
        the nth and (n-1)th positions, and zeroes elsewhere.

        This method returns the ith simple root for the D series.

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("D4")
        >>> c.simple_root(2)
        Matrix([[0, 1, -1, 0]])

        """
        n = self.rank()
        if i < 0 or i > n:
            raise ValueError(f"Lie algebra D{n} has {n} simple roots")

        if i < n:
            return self.basic_root(i-1, i)
        else:
            root = [0]*n
            root[n-2] = 1
            root[n-1] = 1
            return Matrix([root])

    def roots(self) -> int:
        """
        Returns the total number of roots for D_n"
        """

        n = self.n
        return 2*n*(n-1)

    def basis(self) -> int:
        """
        Returns the number of independent generators of D_n
        """
        n = self.n
        return n*(n-1) // 2

    def lie_algebra(self) -> str:
        """
        Returns the Lie algebra associated with D_n"
        """

        n = self.n
        return "so(" + str(2*n) + ")"

    def dynkin_diagram(self) -> str:
        n = self.n
        diag = " "*4*(n-3) + str(n-1) + "\n"
        diag += " "*4*(n-3) + "0\n"
        diag += " "*4*(n-3) +"|\n"
        diag += " "*4*(n-3) + "|\n"
        diag += "---".join("0" for i in range(1,n)) + "\n"
        diag += "   ".join(str(i) for i in range(1, n-1)) + "   "+str(n)
        return diag
