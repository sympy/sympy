from .cartan_base import Standard_Cartan
from sympy.core.backend import eye, ones, sqrt, Matrix


class TypeE(Standard_Cartan):

    def __new__(cls, n):
        if n < 6 or n > 8:
            raise ValueError("Invalid value of n")
        return Standard_Cartan.__new__(cls, "E", n)

    def dimension(self):
        """Dimension of the vector space V underlying the Lie algebra

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("E6")
        >>> c.dimension()
        6
        """

        return self.rank

    def basic_root(self, i, j):
        """
        This is a method just to generate roots
        with a -1 in the ith position and a 1
        in the jth position.

        """

        root = [0]*8
        root[i] = -1
        root[j] = 1
        return Matrix([root])

    def simple_root(self, i):
        """
        Every lie algebra has a unique root system.
        Given a root system Q, there is a subset of the
        roots such that an element of Q is called a
        simple root if it cannot be written as the sum
        of two elements in Q.  If we let D denote the
        set of simple roots, then it is clear that every
        element of Q can be written as a linear combination
        of elements of D with all coefficients non-negative.

        This method returns the ith simple root for E_n. Simple
        roots are chosen in the basis that allows for easy subgrouping
        of E6 and E7 from E8 simple roots.

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
        Matrix([[0, 1, -1, 0, 0, 0]])

        """
        n = self.rank
        e = eye(8)

        if i == n:
            root =  e.row(n-3) - e.row(n-2)
        elif i == n - 1:
            root = - ones(1,8) / 2
        elif i == n - 2:
            root = e.row(n-3) + e.row(n-2)
        else:
            root = e.row(i-1) - e.row(i)
        root = root.tolist()[0]
        if n != 8:
            # in order to simplify the sum, we just use a static
            # calculation and check whether this simple root has a
            # nonzero index in the last index.
            root[n-1] = sqrt(9-n) / 2 * int(root[n-1] != 0)

        return Matrix([root[:n]])


    def roots(self):
        """
        Returns the total number of roots of E_n
        """

        n = self.rank
        if n == 6:
            return 72
        if n == 7:
            return 126
        if n == 8:
            return 240

    def basis(self):
        """
        Returns the number of independent generators of E_n
        """

        n = self.rank
        if n == 6:
            return 78
        if n == 7:
            return 133
        if n == 8:
            return 248

    def dynkin_diagram(self):
        n = self.rank
        diag = " "*8 + str(n) + "\n"
        diag += " "*8 + "0\n"
        diag += " "*8 + "|\n"
        diag += " "*8 + "|\n"
        diag += "---".join("0" for i in range(1, n)) + "\n"
        diag += "1   " + "   ".join(str(i) for i in range(2,n))
        return diag
