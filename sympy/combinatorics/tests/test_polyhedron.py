from sympy import Symbol
from sympy.combinatorics.polyhedron import Polyhedron
from sympy.combinatorics.permutations import Permutation

C1 = Symbol('C1')
C2 = Symbol('C2')
C3 = Symbol('C3')
C4 = Symbol('C4')
C5 = Symbol('C5')
C6 = Symbol('C6')
C7 = Symbol('C7')
C8 = Symbol('C8')

def test_polyhedron():
    pgroup = [Permutation([[0,7,2,5],[6,1,4,3]]),\
              Permutation([[0,7,1,6],[5,2,4,3]]),\
              Permutation([[3,6,0,5],[4,1,7,2]]),\
              Permutation([[7,4,5],[1,3,0],[2],[6]]),\
              Permutation([[1,3,2],[7,6,5],[4],[0]]),\
              Permutation([[4,7,6],[2,0,3],[1],[5]]),\
              Permutation([[1,2,0],[4,5,6],[3],[7]]),\
              Permutation([[4,2],[0,6],[3,7],[1,5]]),\
              Permutation([[3,5],[7,1],[2,6],[0,4]]),\
              Permutation([[2,5],[1,6],[0,4],[3,7]]),\
              Permutation([[4,3],[7,0],[5,1],[6,2]]),\
              Permutation([[4,1],[0,5],[6,2],[7,3]]),\
              Permutation([[7,2],[3,6],[0,4],[1,5]]),\
              Permutation([0,1,2,3,4,5,6,7])]

    faces = ((C1,C8,C3,C6),(C1,C8,C2,C7),(C2,C5,C3,C8),
             (C2,C5,C4,C7),(C4,C7,C1,C6),(C3,C5,C4,C6))


    corners = ['A','B','C','D','E','F','G','H']
    cube = Polyhedron(corners, faces, pgroup)

    assert cube.size == 8
    assert cube.edges == [(C1, C6), (C1, C8), (C3, C8), (C3, C6), \
                          (C1, C7), (C2, C8), (C2, C7), (C2, C5), \
                          (C3, C5), (C4, C5), (C4, C7), (C4, C6)]

    for i in xrange(3):    #  add 180 degree face rotations
        cube.rotate(cube.pgroups[i]**2)

    assert cube.corners == ['A', 'B', 'F', 'C', 'E', 'D', 'G', 'H']

    for i in range(3,7):  # add 240 degree axial corner rotations
        cube.rotate(cube.pgroups[i]**2)

    assert cube.corners == ['A', 'F', 'C', 'B', 'E', 'G', 'H', 'D']

