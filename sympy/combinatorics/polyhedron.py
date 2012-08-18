from sympy.core import Basic, Tuple, FiniteSet
from sympy.core.sympify import sympify
from sympy.combinatorics import Permutation
from sympy.utilities.misc import default_sort_key

from random import choice

class Polyhedron(Basic):
    """
    Represents the Polyhedral symmetry group.

    It is one of the symmetry groups of the Platonic solids.
    There are three polyhedral groups: the tetrahedral group
    of order 12, the octahedral group of order 24, and the
    icosahedral group of order 60.

    All doctests have been given in the docstring of the
    constructor of the object.

    Reference:
    [1] http://mathworld.wolfram.com/PolyhedralGroup.html
    """
    _edges = None

    def __new__(cls, corners, faces=[], pgroups=[]):
        """
        The constructor of the Polyhedron group object.

        It takes up to three parameters: the corners, faces, and
        allowed transformations.

        The corners/vertices are entered as a list of arbitrary
        expressions that are used to identify each vertex.

        The faces are entered as a list of tuples of indices; a tuple
        of indices identifies the vertices which define the face. They
        should be entered in a cw or ccw order. If no faces are given
        then no edges will be computed.

        The allowed transformations are entered as allowable permutations
        of the vertices for the polyhedron. Instance of Permutations
        (as with faces) should refer to the supplied vertices by index.

        Examples
        ========
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.polyhedron import Polyhedron
        >>> from sympy.abc import w, x, y, z

        Here we construct the Polyhedron object for a tetrahedron.

        >>> corners = [w, x, y, z]
        >>> faces = [(0,1,2), (0,2,3), (0,3,1), (1,2,3)]

        Next, allowed transformations of the polyhedron must be given. This
        is given as permutations of vertices.

        Although the vertices of a tetrahedron can be numbered in 24 (4!)
        different ways, there are only 12 different orientations for a
        physical tetrahedron. The following permutations, applied once or
        twice, will generate all 12 of the orientations. (The identity
        permutation, Permutation(range(4)), is not included since it does
        not change the orientation of the vertices.)

        >>> pgroups = [Permutation([[0,1,2], [3]]),\
                      Permutation([[0,1,3], [2]]),\
                      Permutation([[0,2,3], [1]]),\
                      Permutation([[1,2,3], [0]]),\
                      Permutation([[0,1], [2,3]]),\
                      Permutation([[0,2], [1,3]]),\
                      Permutation([[0,3], [1,2]])]

        The Polyhedron is now constructed and demonstrated:

        >>> tetra = Polyhedron(corners, faces, pgroups)
        >>> tetra.size
        4
        >>> tetra.edges
        {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)}
        >>> tetra.corners
        (w, x, y, z)

        It can be rotated with an arbitrary permutation of vertices, e.g.
        the following permutation is not in the pgroups:

        >>> tetra.rotate(Permutation([0, 1, 3, 2]))
        >>> tetra.corners
        (w, x, z, y)

        An allowed permutation of the vertices can be constructed by
        repeatedly applying permutations from the pgroups to the vertices.
        Here is a demonstration that applying p and p**2 for every p in
        pgroups generates all the orientations of a tetrahedron and no others:

        >>> all = ( (w, x, y, z) ,\
                    (x, y, w, z) ,\
                    (y, w, x, z) ,\
                    (w, z, x, y) ,\
                    (z, w, y, x) ,\
                    (w, y, z, x) ,\
                    (y, z, w, x) ,\
                    (x, z, y, w) ,\
                    (z, y, x, w) ,\
                    (y, x, z, w) ,\
                    (x, w, z, y) ,\
                    (z, x, w, y) )

        >>> got = []
        >>> for p in (pgroups + [p**2 for p in pgroups]):
        ...     h = Polyhedron(corners)
        ...     h.rotate(p)
        ...     got.append(h.corners)
        ...
        >>> set(got) == set(all)
        True

        The make_perm method will randomly pick permutations given in
        pgroups, multiply them together, and return the permutation that
        can be applied to the polyhedron to give the orientation produced
        by those individual permutations.

        Here, 3 permutations are used:

        >>> tetra.make_perm(3) # doctest: +SKIP
        Permutation([0, 3, 1, 2])

        To select the permutations that should be used, supply a list
        of indices to the permutations in pgroups:

        >>> use = [0, 0, 2]
        >>> tetra.make_perm(3, use)
        Permutation([1, 0, 3, 2])
        >>> saved = _

        Apply them one at a time:

        >>> h = Polyhedron(corners)
        >>> for i in use:
        ...     h.rotate(pgroups[i])
        ...
        >>> h.vertices
        (x, w, z, y)
        >>> sequentially = _

        Apply the saved "composite" permutation:

        >>> h = Polyhedron(corners)
        >>> h.rotate(saved)
        >>> h.corners
        (x, w, z, y)
        >>> _ in all and _ == sequentially
        True

        For convenience, the vertices and faces are defined for the following
        standard solids (but the allowed transformations are not provided).

        - tetrahedron
        - square
        - octahedron
        - dodecahedron
        - icosahedron

        >>> from sympy.combinatorics.polyhedron import square
        >>> Polyhedron(*square).edges
        {(0, 1), (0, 3), (0, 5), '...', (4, 7), (5, 6), (6, 7)}

        If you want to use letters or other names for the corners you
        can still use the pre-calculated faces:

        >>> _, faces = square
        >>> corners = list('abcdefgh')
        >>> Polyhedron(corners, faces).corners
        (a, b, c, d, e, f, g, h)

        References
        ==========
        [1] www.ocf.berkeley.edu/~wwu/articles/platonicsolids.pdf

        """
        args = [Tuple(*a) for a in (corners, faces, pgroups)]
        obj = Basic.__new__(cls, *args)
        obj._corners = tuple(args[0]) # in order given
        obj._faces = FiniteSet([tuple(a) for a in args[1]])
        if args[2] and args[2][0].size != len(obj.corners):
            raise ValueError("Permutation size unequal to number of corners.")
        if len(set([a.size for a in args[2]])) > 1:
            raise ValueError("All permutations must be of the same size.")
        obj._pgroups = args[2]
        return obj

    @property
    def corners(self):
        """
        Get the corners of the Polyhedron.

        The method ``vertices`` is an alias for ``corners``.

        Examples
        ========
        >>> from sympy.combinatorics import Polyhedron
        >>> from sympy.abc import a, b, c, d
        >>> p = Polyhedron(list('abcd'))
        >>> p.corners == p.vertices == (a, b, c, d)
        True

        """
        return self._corners
    vertices = corners

    @property
    def size(self):
        """
        Get the number of corners of the Polyhedron.
        """
        return len(self._corners)

    @property
    def faces(self):
        """
        Get the faces of the Polyhedron.
        """
        return self._faces

    @property
    def pgroups(self):
        """
        Get the permutations of the Polyhedron.
        """
        return self._pgroups

    @property
    def edges(self):
        """
        Given the faces of the polyhedra we can get the edges.

        Examples
        ========
        >>> from sympy.combinatorics import Polyhedron
        >>> from sympy.abc import a, b, c
        >>> corners = (a, b, c)
        >>> faces = [(0, 1, 2)]
        >>> Polyhedron(corners, faces).edges
        {(0, 1), (0, 2), (1, 2)}

        """
        if self._edges is None:
            output = set()
            n = len(self.corners)
            for face in self.faces:
                for i in xrange(len(face)):
                    edge = tuple(sorted([face[i], face[i - 1]]))
                    if edge not in output:
                        output.add(edge)
            self._edges = FiniteSet(*output)
        return self._edges

    def make_perm(self, n, seed=None):
        """
        Multiply ``n`` randomly selected permutations from
        pgroups together, starting with the identity
        permutation.

        ``seed`` is used to set the seed for the random selection
        of permutations from pgroups. If this is a list of integers,
        the corresponding permutations from pgroups will be selected
        in the order give.

        Examples
        ========
        >>> from sympy.combinatorics import Permutation, Polyhedron
        >>> pgroups = [Permutation([1, 0, 3, 2]), Permutation([1, 3, 0, 2])]
        >>> h = Polyhedron(list('abcd'), pgroups=pgroups)
        >>> h.make_perm(1, [0])
        Permutation([1, 0, 3, 2])
        >>> h.make_perm(3, [0, 1, 0])
        Permutation([2, 0, 3, 1])

        """
        from sympy.utilities.randtest import _randrange
        randrange = _randrange(seed)

        # start with the identity permutation
        result = Permutation(range(len(self.corners)))
        m = len(self.pgroups)
        for i in xrange(n):
            result *= self.pgroups[randrange(m)]
        return result

    def rotate(self, perm):
        """
        Apply permutation to corners of a polyhedron *in place*.

        This is an operation that is analogous to rotation about
        an axis by a fixed increment.

        Examples
        ========
        >>> from sympy.combinatorics import Polyhedron, Permutation
        >>> h = Polyhedron(list('abcde'))
        >>> h.rotate(Permutation([3, 0, 1, 2, 4]))
        >>> h.corners
        (d, a, b, c, e)

        """
        if perm.size != self.size:
            raise ValueError("The size of the permutation and polyhedron must match.")
        temp = []
        for i in xrange(len(self.corners)):
            temp.append(self.corners[perm.array_form[i]])
        self._corners = tuple(temp)

tetrahedron = ([0, 1, 2, 3], [[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 2, 3]])
square = ([0, 1, 2, 3, 4, 5, 6, 7], [[0, 1, 2, 3], [0, 1, 4, 5], [1, 2, 5, 6],
[2, 3, 6, 7], [3, 0, 7, 4], [4, 5, 6, 7]])
octahedron = ([0, 1, 2, 3, 4, 5], [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
[1, 2, 5], [2, 3, 5], [3, 4, 5], [4, 1, 5]])
dodecahedron = ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
18, 19], [[0, 1, 2, 3, 4], [0, 1, 5, 6, 10], [1, 2, 6, 7, 11], [2, 3, 7,
8, 12], [3, 4, 8, 9, 13], [4, 0, 9, 5, 14], [5, 10, 11, 15, 16], [6,
11, 12, 16, 17], [7, 12, 13, 17, 18], [8, 13, 14, 18, 19], [9, 14, 10,
19, 15], [15, 16, 17, 18, 19]])
icosahedron = ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], [[0, 1, 2], [0, 2, 3],
[0, 3, 4], [0, 4, 5], [0, 5, 1], [[1, 2, 6], [2, 3, 6], [3, 4, 6], [4, 5, 6],
[5, 1, 6]], [[1, 6, 7], [1, 7, 8], [1, 8, 9], [1, 9, 10], [1, 10, 6]],
[[1, 2, 7], [2, 3, 7], [3, 4, 7], [4, 5, 7], [5, 1, 7]], [[2, 6, 7],
[2, 7, 8], [2, 8, 9], [2, 9, 10], [2, 10, 6]], [[1, 2, 8], [2, 3, 8],
[3, 4, 8], [4, 5, 8], [5, 1, 8]], [[3, 6, 7], [3, 7, 8], [3, 8, 9],
[3, 9, 10], [3, 10, 6]], [[1, 2, 9], [2, 3, 9], [3, 4, 9], [4, 5, 9],
[5, 1, 9]], [[4, 6, 7], [4, 7, 8], [4, 8, 9], [4, 9, 10], [4, 10, 6]],
[[1, 2, 10], [2, 3, 10], [3, 4, 10], [4, 5, 10], [5, 1, 10]], [[5, 6,
7], [5, 7, 8], [5, 8, 9], [5, 9, 10], [5, 10, 6]], [6, 7, 11], [7, 8,
11], [8, 9, 11], [9, 10, 11], [10, 6, 11]])
