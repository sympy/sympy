from sympy.core import Basic, Tuple, FiniteSet
from sympy.core.sympify import sympify
from sympy.combinatorics import Permutation
from sympy.utilities.misc import default_sort_key
from sympy.utilities.iterables import (rotate_left, has_variety,
    is_sequence, minlex)
from sympy.utilities.randtest import _randrange

lmul = Permutation.lmul

class Polyhedron(Basic):
    """
    Represents the polyhedral symmetry group (PSG).

    The PSG is one of the symmetry groups of the Platonic solids.
    There are three polyhedral groups: the tetrahedral group
    of order 12, the octahedral group of order 24, and the
    icosahedral group of order 60.

    All doctests have been given in the docstring of the
    constructor of the object.

    References
    ==========

    http://mathworld.wolfram.com/PolyhedralGroup.html
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
        should be entered in a cw or ccw order; they will be standardized
        by reversal and rotation to be give the lowest lexical ordering.
        If no faces are given then no edges will be computed.

            >>> from sympy.combinatorics.polyhedron import Polyhedron
            >>> Polyhedron(list('abc'), [(1, 2, 0)]).faces
            {(0, 1, 2)}
            >>> Polyhedron(list('abc'), [(1, 0, 2)]).faces
            {(0, 1, 2)}

        The allowed transformations are entered as allowable permutations
        of the vertices for the polyhedron. Instance of Permutations
        (as with faces) should refer to the supplied vertices by index.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
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

        Notes
        =====

        It is not necessary to enter any permutations, nor is necessary to
        enter a complete set of transforations. In fact, two permutations
        (corresponding to a rotation on an axis through a vertex and face
        and another for the rotation through a different vertex or from
        one edge to the opposite edge) are sufficient to generate all
        orientations. For simplicity of presentation, consider a square --
        not a cube -- with vertices 1, 2, 3, and 4:

        1-----2  We could think of axes of rotation being:
        |     |  1) through the face
        |     |  2) from midpoint 1-2 to 3-4 or 1-3 to 2-4
        3-----4  3) lines 1-4 or 2-3


        To determine how to write the permutations, imagine 4 cameras, one at
        each corner:

        A       B          A       B
         1-----2            1-----3             vertex index:
         |     |            |     |                 1   0
         |     |            |     |                 2   1
         3-----4            2-----4                 3   2
        C       D          C       D                4   3

        original           after rotation
                           along 1-4

        A diagonal and a face axis will be chosen for the "permutation group"
        from which any orientation can be constructed.

        >>> pgroup = []

        Imagine rotating clockwise when viewing 1-4 from camera A. The new
        orientation is (in camera-order): 1, 3, 2, 4 so the permutation is
        given using the *indices* of the vertices as:

        >>> pgroup.append(Permutation((0, 2, 1, 3)))

        Now imagine rotating clockwise when looking down an axis entering the
        center of the square as viewed. The new camera-order would be
        3, 1, 4, 2 so the permutation is (using indices):

        >>> pgroup.append(Permutation((2, 0, 3, 1)))

        The square can now be constructed:
            ** use real-world labels for the vertices, entering them in
               camera order
            ** for the faces we use zero-based indices of the vertices
               in *edge-order* as the face is traversed; neither the
               direction nor the starting point matter -- the faces are
               only used to define edges (if so desired).

        >>> square = Polyhedron((1, 2, 3, 4), [(0, 1, 3, 2)], pgroup)

        To rotate the square a single permutation we can do:

        >>> sq = square.copy()
        >>> sq.rotate(square.pgroups[0]); sq.corners
        (1, 3, 2, 4)

        To use more than one permutation (or to use one permutation more
        than once) it is more convenient to use the make_perm method:

        >>> p011 = square.make_perm([0,1,1]) # diag flip and 2 rotations
        >>> sq = square.copy()
        >>> sq.rotate(p011); sq.corners
        (4, 2, 3, 1)


        Predefined Polyhedra
        ====================

        For convenience, the vertices and faces are defined for the following
        standard solids along with a permutation group for transformations.
        When the polyhedron is oriented as indicated below, the vertices in
        a given horizontal plane are numbered in ccw direction, starting from
        the vertex that will give the lowest indices in a given face. (In the
        net of the vertices, indices preceded by "-" indicate replication of
        the lhs index in the net.)

        tetrahedron, tetrahedron_faces
        ------------------------------

            4 vertices (vertex up) net:

                 0 0-0
                1 2 3-1

            4 faces:

            (0,1,2) (0,2,3) (0,3,1) (1,2,3)

        cube, cube_faces
        ----------------

            8 vertices (face up) net:

                0 1 2 3-0
                4 5 6 7-4

            6 faces:

            (0,1,2,3)
            (0,1,5,4) (1,2,6,5) (2,3,7,6) (0,3,7,4)
            (4,5,6,7)

        octahedron, octahedron_faces
        ----------------------------

            6 vertices (vertex up) net:

                 0 0 0-0
                1 2 3 4-1
                 5 5 5-5

            8 faces:

            (0,1,2) (0,2,3) (0,3,4) (0,1,4)
            (1,2,5) (2,3,5) (3,4,5) (1,4,5)

        dodecahedron, dodecahedron_faces
        --------------------------------

            20 vertices (vertex up) net:

                  0  1  2  3  4 -0
                  5  6  7  8  9 -5
                14 10 11 12 13-14
                15 16 17 18 19-15

            12 faces:

            (0,1,2,3,4)
            (0,1,6,10,5) (1,2,7,11,6) (2,3,8,12,7) (3,4,9,13,8) (0,4,9,14,5)
            (5,10,16,15,14) (6,10,16,17,11) (7,11,17,18,12) (8,12,18,19,13) (9,13,19,15,14)
            (15,16,17,18,19)

        icosahedron, icosahedron_faces
        ------------------------------

            12 vertices (face up) net:

                 0  0  0  0 -0
                1  2  3  4  5 -1
                 6  7  8  9  10 -6
                  11 11 11 11 -11

            20 faces:

            (0,1,2) (0,2,3) (0,3,4) (0,4,5) (0,1,5)
            (1,2,6) (2,3,7) (3,4,8) (4,5,9) (1,5,10)
            (2,6,7) (3,7,8) (4,8,9) (5,9,10) (1,6,10)
            (6,7,11,) (7,8,11) (8,9,11) (9,10,11) (6,10,11)

        >>> from sympy.combinatorics.polyhedron import cube
        >>> cube.edges
        {(0, 1), (0, 3), (0, 4), '...', (4, 7), (5, 6), (6, 7)}

        If you want to use letters or other names for the corners you
        can still use the pre-calculated faces:

        >>> corners = list('abcdefgh')
        >>> Polyhedron(corners, cube.faces).corners
        (a, b, c, d, e, f, g, h)

        References
        ==========

        [1] www.ocf.berkeley.edu/~wwu/articles/platonicsolids.pdf

        """
        faces = [minlex(f, directed=False) for f in faces]
        corners, faces, pgroups = args = \
            [Tuple(*a) for a in (corners, faces, pgroups)]
        obj = Basic.__new__(cls, *args)
        obj._corners = tuple(corners) # in order given
        obj._faces = FiniteSet(faces)
        if pgroups and pgroups[0].size != len(corners):
            raise ValueError("Permutation size unequal to number of corners.")
        if has_variety(a.size for a in pgroups) > 1:
            raise ValueError("All permutations must be of the same size.")
        obj._pgroups = pgroups
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
                for i in range(len(face)):
                    edge = tuple(sorted([face[i], face[i - 1]]))
                    if edge not in output:
                        output.add(edge)
            self._edges = FiniteSet(*output)
        return self._edges

    def make_perm(self, n, seed=None):
        """
        Multiply ``n`` randomly selected permutations from
        pgroups together, starting with the identity
        permutation. If ``n`` is a list of integers, those
        integers will be used to select the permutations.

        ``seed`` is used to set the seed for the random selection
        of permutations from pgroups. If this is a list of integers,
        the corresponding permutations from pgroups will be selected
        in the order give. This is mainly used for testing purposes.

        Examples
        ========

        >>> from sympy.combinatorics import Permutation, Polyhedron
        >>> pgroups = [Permutation([1, 0, 3, 2]), Permutation([1, 3, 0, 2])]
        >>> h = Polyhedron(list('abcd'), pgroups=pgroups)
        >>> h.make_perm(1, [0])
        Permutation([1, 0, 3, 2])
        >>> h.make_perm(3, [0, 1, 0])
        Permutation([2, 0, 3, 1])
        >>> h.make_perm([0, 1, 0])
        Permutation([2, 0, 3, 1])

        """
        if is_sequence(n):
            if is_sequence(seed):
                raise ValueError('If n is a sequence, seed should be None')
            n, seed = len(n), n
        randrange = _randrange(seed)

        # start with the identity permutation
        result = Permutation(range(len(self.corners)))
        m = len(self.pgroups)
        for i in range(n):
            result = lmul(result, self.pgroups[randrange(m)])
        return result

    def rotate(self, perm):
        """
        Apply permutation to corners of a polyhedron *in place*.

        This is an operation that is analogous to rotation about
        an axis by a fixed increment.

        Examples
        ========

        >>> from sympy.combinatorics import Polyhedron, Permutation
        >>> shadow = h = Polyhedron(list('abcde'))
        >>> p = Permutation([3, 0, 1, 2, 4])
        >>> h.rotate(p)
        >>> h.corners
        (d, a, b, c, e)
        >>> _ == shadow.corners
        True
        >>> copy = h.copy()
        >>> h.rotate(p)
        >>> h.corners == copy.corners
        False
        """
        if perm.size != self.size:
            raise ValueError("The size of the permutation and polyhedron must match.")
        temp = []
        for i in range(len(self.corners)):
            temp.append(self.corners[perm.array_form[i]])
        self._corners = tuple(temp)

def _pgroup_calcs():
    """
    Although only 2 permutations are needed for a polyhedron in order to
    generate all the possible orientations, it is customary to give a
    list of permutations (P0, P1, ...) such that powers of them alone are
    able to generate the orientations, e.g. P0, P0**2, P0**3, P1, P1**2,
    etc..., instead of mixed permutations (P0*P1**2*P0). The following
    work was used to calculate the permutation group of the polyhedra.
    """
    def _pgroups_of_double(polyh, ordered_faces, pgroup):
        from sympy.utilities import unflatten, flatten
        from sympy.ntheory.residue_ntheory import int_tested
        n = len(ordered_faces[0])
        # the vertices of the double which sits inside a give polyhedron
        # can be found by tracking the faces of the outer polyhedron.
        # A map between face and the vertex of the double is made so that
        # after rotation the position of the vertices can be located
        fmap = dict(zip(ordered_faces,
                         range(len(ordered_faces))))
        flat_faces = flatten(ordered_faces)
        new_pgroup = []
        for i, p in enumerate(pgroup):
            h = polyh.copy()
            h.rotate(p)
            c = h.corners
            # reorder corners in the order they should appear when
            # enumerating the faces
            reorder = unflatten([c[j] for j in flat_faces], n)
            # make them canonical
            reorder = [tuple(int_tested(minlex(f, directed=False)))
                for f in reorder]
            # map face to vertex: the resulting list of vertices are the
            # permutation that we seek for the double
            new_pgroup.append(Permutation([fmap[f] for f in reorder]))
        return new_pgroup

    tetrahedron_faces = [(0, 1, 2), (0, 2, 3), (0, 3, 1), (1, 2, 3)]

    _t_pgroups = [
        Permutation([[0,1,2], [3]]),\
        Permutation([[0,1,3], [2]]),\
        Permutation([[0,2,3], [1]]),\
        Permutation([[1,2,3], [0]]),\
        Permutation([[0,1], [2,3]]),\
        Permutation([[0,2], [1,3]]),\
        Permutation([[0,3], [1,2]])]

    tetrahedron = Polyhedron(
        range(4),
        tetrahedron_faces,
        _t_pgroups)

    cube_faces = [
    (0, 1, 2, 3),
    (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (0, 3, 7, 4),
    (4, 5, 6, 7)]

    _c_pgroups = [Permutation(p) for p in
        [[1,2,3,0,5,6,7,4],
        [4,0,3,7,5,1,2,6],
        [4,5,1,0,7,6,2,3],

        [1,0,4,5,2,3,7,6],
        [6,2,1,5,7,3,0,4],
        [6,7,3,2,5,4,0,1],
        [3,7,4,0,2,6,5,1],
        [6,5,4,7,2,1,0,3],
        [4,7,6,5,0,3,2,1],

        [0,3,7,4,1,2,6,5],
        [5,1,0,4,6,2,3,7],
        [5,6,2,1,4,7,3,0],
        [7,4,0,3,6,5,1,2]]]

    cube = Polyhedron(
        range(8),
        cube_faces,
        _c_pgroups)

    octahedron_faces = [
        (0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 1, 4),
        (1, 2, 5), (2, 3, 5), (3, 4, 5), (1, 4, 5)]

    octahedron = Polyhedron(
        range(6),
        octahedron_faces,
        _pgroups_of_double(cube, cube_faces, _c_pgroups))

    dodecahedron_faces = [
         (0,1,2,3,4),
                         (0,1,6,10,5), (1,2,7,11,6), (2,3,8,12,7),
                         (3,4,9,13,8), (0,4,9,14,5),
                         (5,10,16,15,14), (6,10,16,17,11), (7,11,17,18,12),
                         (8,12,18,19,13), (9,13,19,15,14),
         (15,16,17,18,19)]

    def _string_to_perm(s):
        rv = [Permutation(range(20))]
        p = None
        for si in s:
            if si not in '01':
                count = int(si) - 1
            else:
                count = 1
                if si == '0':
                    p = _f0
                elif si == '1':
                    p = _f1
            rv.extend([p]*count)
        return Permutation.lmul(*rv)

    # top face cw
    _f0 = Permutation([
        1, 2, 3, 4, 0, 6, 7, 8, 9, 5, 11,
        12, 13, 14, 10, 16, 17, 18, 19, 15])
    # front face cw
    _f1 = Permutation([
        5, 0, 4, 9, 14, 10, 1, 3, 13, 15,
        6, 2, 8, 19, 16, 17, 11, 7, 12, 18])
    # the strings below, like 0104 are shorthand for F0*F1*F0**4 and are
    # the remaining 4 face rotations, 15 edge permutations, and the
    # 10 vertex rotations.
    _dodeca_pgroups = [_f0, _f1] + [_string_to_perm(s) for s in '''
    0104 140 014 0410
    010 1403 03104 04103 102
    120 1304 01303 021302 03130
    0412041 041204103 04120410 041204104 041204102
    10 01 1402 0140 04102 0412 1204 1302 0130 03120'''.strip().split()]

    dodecahedron = Polyhedron(
        range(20),
        dodecahedron_faces,
        _dodeca_pgroups)

    icosahedron_faces = [
        [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 1, 5],
        [1, 6, 7], [1, 2, 7], [2,7,  8], [2,3,8  ], [3, 8, 9 ],
        [3, 4, 9], [4,9,10 ], [4, 5,10], [5, 6, 10], [1, 5, 6 ],
        [6, 7, 11], [7, 8, 11], [8, 9, 11], [9, 10, 11], [6, 10, 11]]

    icosahedron = Polyhedron(
        range(12),
        icosahedron_faces,
        _pgroups_of_double(dodecahedron, dodecahedron_faces, _dodeca_pgroups))

    return (tetrahedron, cube, octahedron, dodecahedron, icosahedron,
        tetrahedron_faces, cube_faces, octahedron_faces,
        dodecahedron_faces, icosahedron_faces)

(tetrahedron, cube, octahedron, dodecahedron, icosahedron,
tetrahedron_faces, cube_faces, octahedron_faces,
dodecahedron_faces, icosahedron_faces) = _pgroup_calcs()
