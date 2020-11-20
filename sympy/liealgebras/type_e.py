from .cartan_base import Standard_Cartan
from sympy.core.backend import eye, Rational, ones, sqrt


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

        return self.n

    def basic_root(self, i, j):
        """
        This is a method just to generate roots
        with a -1 in the ith position and a 1
        in the jth position.

        """

        root = [0]*8
        root[i] = -1
        root[j] = 1
        return root

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

        This method returns the ith simple root for E_n.

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("E6")
        >>> c.simple_root(2)
        [0, 1, -1, 0, 0, 0]

        """
        n = self.n
        def _dim8_weights(i,n):
            """
            Returns the roots in a basis
            that makes it easy to subset
            e7 and e6.

            """
            e = eye(8)
            if i == n:
                return (e.row(n-3) - e.row(n-2)).T
            if i == n - 1:
                return  - ones(8,1) / 2
            if i == n - 2:
                return (e.row(n-3) + e.row(n-2)).T
            else:
                return (e.row(i-1) - e.row(i)).T

        _root = _dim8_weights(i, n).T
        if i != 8:
            _root[n-1] = sqrt(sum(_root[j]**2 for j in range(n-1, 8)))            
        return _root[:n]



    def positive_roots(self):
        """
        This method generates all the positive roots of
        A_n.  This is half of all of the roots of E_n;
        by multiplying all the positive roots by -1 we
        get the negative roots.

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("A3")
        >>> c.positive_roots()
        {1: [1, -1, 0, 0], 2: [1, 0, -1, 0], 3: [1, 0, 0, -1], 4: [0, 1, -1, 0],
                5: [0, 1, 0, -1], 6: [0, 0, 1, -1]}

        """
        n = self.n
        if n == 6:
            posroots = {}
            k = 0
            for i in range(n-1):
                for j in range(i+1, n-1):
                    k += 1
                    root = self.basic_root(i, j)
                    posroots[k] = root
                    k += 1
                    root = self.basic_root(i, j)
                    root[i] = 1
                    posroots[k] = root

            root = [Rational(1, 2), Rational(1, 2), Rational(1, 2), Rational(1, 2), Rational(1, 2),
                    Rational(-1, 2), Rational(-1, 2), Rational(1, 2)]
            for a in range(0, 2):
                for b in range(0, 2):
                    for c in range(0, 2):
                        for d in range(0, 2):
                            for e in range(0, 2):
                                if (a + b + c + d + e)%2 == 0:
                                    k += 1
                                    if a == 1:
                                        root[0] = Rational(-1, 2)
                                    if b == 1:
                                        root[1] = Rational(-1, 2)
                                    if c == 1:
                                        root[2] = Rational(-1, 2)
                                    if d == 1:
                                        root[3] = Rational(-1, 2)
                                    if e == 1:
                                        root[4] = Rational(-1, 2)
                                    posroots[k] = root

            return posroots
        if n == 7:
            posroots = {}
            k = 0
            for i in range(n-1):
                for j in range(i+1, n-1):
                    k += 1
                    root = self.basic_root(i, j)
                    posroots[k] = root
                    k += 1
                    root = self.basic_root(i, j)
                    root[i] = 1
                    posroots[k] = root

            k += 1
            posroots[k] = [0, 0, 0, 0, 0, 1, 1, 0]
            root = [Rational(1, 2), Rational(1, 2), Rational(1, 2), Rational(1, 2), Rational(1, 2),
                    Rational(-1, 2), Rational(-1, 2), Rational(1, 2)]
            for a in range(0, 2):
                for b in range(0, 2):
                    for c in range(0, 2):
                        for d in range(0, 2):
                            for e in range(0, 2):
                                for f in range(0, 2):
                                    if (a + b + c + d + e + f)%2 == 0:
                                        k += 1
                                        if a == 1:
                                            root[0] = Rational(-1, 2)
                                        if b == 1:
                                            root[1] = Rational(-1, 2)
                                        if c == 1:
                                            root[2] = Rational(-1, 2)
                                        if d == 1:
                                            root[3] = Rational(-1, 2)
                                        if e == 1:
                                            root[4] = Rational(-1, 2)
                                        if f == 1:
                                            root[5] = Rational(1, 2)
                                        posroots[k] = root

            return posroots
        if n == 8:
            posroots = {}
            k = 0
            for i in range(n):
                for j in range(i+1, n):
                    k += 1
                    root = self.basic_root(i, j)
                    posroots[k] = root
                    k += 1
                    root = self.basic_root(i, j)
                    root[i] = 1
                    posroots[k] = root

            root = [Rational(1, 2), Rational(1, 2), Rational(1, 2), Rational(1, 2), Rational(1, 2),
                    Rational(-1, 2), Rational(-1, 2), Rational(1, 2)]
            for a in range(0, 2):
                for b in range(0, 2):
                    for c in range(0, 2):
                        for d in range(0, 2):
                            for e in range(0, 2):
                                for f in range(0, 2):
                                    for g in range(0, 2):
                                        if (a + b + c + d + e + f + g)%2 == 0:
                                            k += 1
                                            if a == 1:
                                                root[0] = Rational(-1, 2)
                                            if b == 1:
                                                root[1] = Rational(-1, 2)
                                            if c == 1:
                                                root[2] = Rational(-1, 2)
                                            if d == 1:
                                                root[3] = Rational(-1, 2)
                                            if e == 1:
                                                root[4] = Rational(-1, 2)
                                            if f == 1:
                                                root[5] = Rational(1, 2)
                                            if g == 1:
                                                root[6] = Rational(1, 2)
                                            posroots[k] = root

            return posroots



    def roots(self):
        """
        Returns the total number of roots of E_n
        """

        n = self.n
        if n == 6:
            return 72
        if n == 7:
            return 126
        if n == 8:
            return 240


    def cartan_matrix(self):
        """
        Returns the Cartan matrix for G_2
        The Cartan matrix matrix for a Lie algebra is
        generated by assigning an ordering to the simple
        roots, (alpha[1], ...., alpha[l]).  Then the ijth
        entry of the Cartan matrix is (<alpha[i],alpha[j]>).

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType('A4')
        >>> c.cartan_matrix()
        Matrix([
        [ 2, -1,  0,  0],
        [-1,  2, -1,  0],
        [ 0, -1,  2, -1],
        [ 0,  0, -1,  2]])

        """

        return super().cartan_matrix()


    def basis(self):
        """
        Returns the number of independent generators of E_n
        """

        n = self.n
        if n == 6:
            return 78
        if n == 7:
            return 133
        if n == 8:
            return 248

    def dynkin_diagram(self):
        n = self.n
        diag = " "*8 + str(2) + "\n"
        diag += " "*8 + "0\n"
        diag += " "*8 + "|\n"
        diag += " "*8 + "|\n"
        diag += "---".join("0" for i in range(1, n)) + "\n"
        diag += "1   " + "   ".join(str(i) for i in range(3, n+1))
        return diag
