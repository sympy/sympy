from .cartan_base import Standard_Cartan
from sympy.core.backend import Matrix,ImmutableMatrix, Rational
from sympy.utilities.iterables import multiset_permutations

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
        Matrix([[-1, 1, 0, 0, 0, 0, 0, 0]])

        """
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


        roots = e8[:self.rank - 1] + [e8[-1]]
        return ImmutableMatrix([roots[i-1]])



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

    def rootsystem(self, **kwargs):
        # doc string inherited
        if self._is_default_basis:
            return _e_root_system(self.rank, self.basic_root, self._orbit_sorter_lambda())
        return super().rootsystem(**kwargs)


def _e_root_system(n, basic_root, sort_key):
    """Returns the TypeE root system based based on enumeration"""
    def mirror_add(l, r):
        l.extend([ImmutableMatrix(r),ImmutableMatrix(-r)])

    proots = []
    # enumerated integer roots
    for i in range(n-1):
        top = n-1 if n != 8 else n
        for j in range(i+1, top):
            mirror_add(proots, basic_root(i, j))
            r = basic_root(i, j)
            r[i] = 1
            mirror_add(proots, r)


    # all possible plus/minus combos of rationals for n-1 entries
    # some are not possible as they have non integer values for dynkin coeffecients (aka omega basis)
    m=n-1
    for combo in set(tuple(x[:m]) for x in multiset_permutations([Rational(1, 2)]*m + [Rational(-1, 2)]*m)):
        tail = [Rational(-1, 2), Rational(-1, 2), Rational(1, 2)][n-6:]
        frac_root = ImmutableMatrix([list(combo) + tail])

        _, t = sort_key(frac_root)
        if all(x.is_integer for x in t):
            proots.append(ImmutableMatrix(frac_root))
            proots.append(ImmutableMatrix(-frac_root))

    if n == 7:
        mirror_add(proots, ImmutableMatrix([[0, 0, 0, 0, 0, 0, -1, 1]]))

    proots += [ImmutableMatrix(basic_root(0,0) * 0)] * n

    proots = sorted(proots, key=sort_key)

    return proots
