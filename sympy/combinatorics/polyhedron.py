from sympy.core import Basic
from sympy.combinatorics import Permutation

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
    _perm_size = None
    _faces = []
    _pgroups = []
    _corners = []

    @property
    def corners(self):
        return self._corners

    @property
    def size(self):
        return self._perm_size

    @property
    def faces(self):
        return self._faces

    @property
    def edges(self):
        """
        Given the faces of the polyhedra we can get the edges.
        """
        if self._edges is None:
            output = []
            for face in self.faces:
                for i in xrange(len(face)):
                    edge = tuple(sorted([face[i], face[i - 1]]))
                    if edge not in output:
                        output.append(edge)
            self._edges = output
        return self._edges

    @property
    def pgroups(self):
        """
        Get the permutations of the Polyhedron.
        """
        return self._pgroups

    def make_perm(self, n):
        """
        Multiply randomly selected permutations from
        pgroups together, starting with the identity
        permutation.  Randomly choose and multiply n times
        """
        # For polyhedra, the last permutation is the identity permutation
        result = self.pgroups[-1]
        for i in xrange(n):
            result *= choice(self.pgroups)
        return result

    def rotate(self, perm):
        """
        Apply permutation to corners of a polyhedron -- an
        operation analogous to rotation about an axis by a
        fixed increment
        """
        if perm.size != self.size:
            raise ValueError("The size of the permutation and polyhedron must match.")
        temp = []
        for i in xrange(len(self.corners)):
            temp.append(self.corners[perm.array_form[i]])
        self._corners = temp

    def __new__(cls, *args):
        """
        The constructor of the Polyhedron group object.
        It takes three parameters, a representation of
        the corners, the faces and a representation of
        the edges.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.polyhedron import Polyhedron
        >>> from sympy.abc import w, x, y, z
        >>> pgroup = [Permutation([[0,1,2], [3]]),\
                      Permutation([[0,1,3], [2]]),\
                      Permutation([[0,2,3], [1]]),\
                      Permutation([[1,2,3], [0]]),\
                      Permutation([[0,1], [2,3]]),\
                      Permutation([[0,2], [1,3]]),\
                      Permutation([[0,3], [1,2]]),\
                      Permutation([[0, 1, 2, 3]])]
        >>> corners = [w, x, y, z]
        >>> faces = [(w,x,y),(w,y,z),(w,z,x),(x,y,z)]
        >>> tetra = Polyhedron(corners, faces, pgroup)
        >>> tetra.size
        4
        >>> tetra.edges
        [(w, y), (w, x), (x, y), (w, z), (y, z), (x, z)]
        >>> tetra.corners
        [w, x, y, z]
        >>> tetra.rotate(Permutation([3,2,1,0]))
        >>> tetra.corners
        [z, y, x, w]
        """
        ret_obj = Basic.__new__(cls, *args)
        ret_obj._corners = args[0]
        ret_obj._faces = args[1]
        ret_obj._perm_size = args[2][0].size
        ret_obj._pgroups = [args[2][0]] + filter(lambda x: x.size == \
                                                ret_obj._perm_size,
                                                args[2][1:])
        if len(ret_obj._pgroups) != len(args[2]):
            raise ValueError("All permutations must be of the same size.")
        return ret_obj
