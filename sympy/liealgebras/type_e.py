from __future__ import annotations

from typing import Literal

from .cartan_type import Standard_Cartan
from sympy.core.backend import Rational, Matrix

class TypeE(Standard_Cartan):

    def __new__(cls, n: Literal[6, 7, 8]):
        if n < 6 or n > 8:
            raise ValueError("Invalid value of n")
        return Standard_Cartan.__new__(cls, "E", n)

    def dimension(self) -> int:
        """Dimension of the vector space V underlying the Lie algebra

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("E6")
        >>> c.dimension()
        8
        """

        return 8

    def basic_root(self, i: int, j: int) -> Matrix:
        """
        This is a method just to generate roots
        with a -1 in the ith position and a 1
        in the jth position.

        """

        root = [0]*8
        root[i] = -1
        root[j] = 1
        return Matrix([root])

    def simple_root(self, i: int) -> Matrix:
        """
        Every Lie algebra has a unique root system.
        Given a root system Q, there is a subset of the
        roots such that an element of Q is called a
        simple root if it cannot be written as the sum
        of two elements in Q.  If we let D denote the
        set of simple roots, then it is clear that every
        element of Q can be written as a linear combination
        of elements of D with all coefficients non-negative.

        This method returns the i-th simple root for E_n. Simple
        roots are chosen in the basis that allows for easy subgrouping
        of E6 and E7 from E8 simple roots.

        Sources
        =======
        - https://en.wikipedia.org/wiki/E8_(mathematics)#Simple_roots
        - https://en.wikipedia.org/wiki/E7_(mathematics)#Root_system
        - https://en.wikipedia.org/wiki/E6_(mathematics)#Roots_of_E6

        Sources
        =======
        - https://en.wikipedia.org/wiki/E8_(mathematics)#Simple_roots
        - https://en.wikipedia.org/wiki/E7_(mathematics)#Root_system
        - https://en.wikipedia.org/wiki/E6_(mathematics)#Roots_of_E6

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("E6")
        >>> c.simple_root(2)
        Matrix([[-1, 1, 0, 0, 0, 0, 0, 0]])
        """
        n = self.rank()
        if i < 0 or i > n:
            raise ValueError(f"Lie algebra E{n} has {n} simple roots")

        e8 = [[Rational(1, 2), Rational(-1, 2), Rational(-1, 2),
            Rational(-1, 2), Rational(-1, 2), Rational(-1, 2),
            Rational(-1, 2), Rational(1, 2)],
            [-1, 1, 0, 0, 0, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0, 0],
            [0, 0, -1, 1, 0, 0, 0, 0],
            [0, 0, 0, -1, 1, 0, 0, 0],
            [0, 0, 0, 0, -1, 1, 0, 0],
            [0, 0, 0, 0, 0, -1, 1, 0],
            [1, 1, 0, 0, 0, 0, 0, 0]]

        roots = e8[:self.rank() - 1] + [e8[-1]]
        return Matrix([roots[i-1]])

    def roots(self) -> int:
        """
        Returns the total number of roots of E_n
        """

        n = self.rank()
        if n == 6:
            return 72
        if n == 7:
            return 126
        return 240

    def basis(self) -> int:
        """
        Returns the number of independent generators of E_n
        """

        n = self.rank()
        if n == 6:
            return 78
        if n == 7:
            return 133
        return 248

    def dynkin_diagram(self) -> str:
        n = self.rank()
        diag = " "*8 + str(2) + "\n"
        diag += " "*8 + "0\n"
        diag += " "*8 + "|\n"
        diag += " "*8 + "|\n"
        diag += "---".join("0" for _ in range(1, n)) + "\n"
        diag += "1   " + "   ".join(str(i) for i in range(3, n+1))
        return diag
