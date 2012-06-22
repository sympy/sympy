"""
This module contains the functionality to arrange the nodes of a
diagram on an abstract grid, and then to produce a graphical
representation of the grid.

The currently supported back-ends are Xy-pic [Xypic]

[Xypic] http://www.tug.org/applications/Xy-pic/
"""

from sympy.core import Basic, FiniteSet
from sympy.categories import CompositeMorphism, IdentityMorphism

class _GrowableGrid:
    """
    Holds a growable grid of objects.

    It is possible to append or prepend a row or a column to the grid
    using the corresponding methods.  Prepending rows or columns has
    the effect of changing the coordinates of the already existing
    elements.

    This class currently represents a naive implementation of the
    functionality with little attempt at optimisation.
    """
    def __init__(self, width, height):
        self._width = width
        self._height = height

        self._array = [[None for j in xrange(width)] for i in xrange(height)]

    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height


    def __getitem__(self, (i, j)):
        """
        Returns the element located at in the i-th line and j-th
        column.
        """
        return self._array[i][j]

    def __setitem__(self, (i, j), newvalue):
        """
        Sets the element located at in the i-th line and j-th
        column.
        """
        self._array[i][j] = newvalue

    def append_row(self):
        """
        Appends an empty row to the grid.
        """
        self._height += 1
        self._array.append([None for j in xrange(self._width)])

    def append_column(self):
        """
        Appends an empty column to the grid.
        """
        self._width += 1
        for i in xrange(self._height):
            self._array[i].append(None)

    def prepend_row(self):
        """
        Prepends the grid with an empty row.
        """
        self._height += 1
        self._array.insert(0, [None for j in xrange(self._width)])

    def prepend_column(self):
        """
        Prepends the grid with an empty column.
        """
        self._width += 1
        for i in xrange(self._height):
            self._array[i].insert(0, None)

class DiagramGrid(Basic):
    """
    Constructs and holds the fitting of the diagram into a grid.

    The constructor of this class takes a :class:`Diagram` as an
    argument and fits it into a grid.  It is then possible to find the
    with and height of the grid using the properties ``width`` and
    ``height``, as well as accessing the objects located at certain
    coordinates and finding out which moprhisms connect the object
    with other objects in the grid.  TODO: Explain how to do that.
    """
    @staticmethod
    def _simplify_morphisms(morphisms):
        """
        Given a dictionary mapping morphisms to their properties,
        returns a new dictionary in which there are no morphisms which
        do not have properties, and which are compositions of other
        morphisms included in the dictionary.  Identities are dropped
        as well.
        """
        newmorphisms = {}
        for morphism, props in morphisms.items():
            if isinstance(morphism, CompositeMorphism) and not props:
                continue
            elif isinstance(morphism, IdentityMorphism):
                continue
            else:
                newmorphisms[morphism] = props
        return newmorphisms

    @staticmethod
    def _merge_premises_conclusions(premises, conclusions):
        """
        Given two dictionaries of morphisms and their properties,
        produces a single dictionary which includes elements from both
        dictionaries.  If a morphism has some properties in premises
        and also in conclusions, the properties in conclusions take
        priority.
        """
        merged = premises
        for morphism, props in conclusions.items():
            merged[morphism] = props
        return merged

    @staticmethod
    def _juxtapose_edges(edge1, edge2):
        """
        If ``edge1`` and ``edge2`` have a common extremity, returns an
        edge which would form a triangle with ``edge1`` and ``edge2``.

        If ``edge1`` and ``edge2`` don't have a common extremity,
        returns ``None``.

        If ``edge1`` and ``edge`` are the same edge, returns ``None``.
        """
        if (edge1 == edge2) or ((edge1[0] == edge2[1]) and (edge1[1] == edge2[0])):
            return None

        for i in [0, 1]:
            for j in [0, 1]:
                if edge1[i] == edge2[j]:
                    # Some extremities match, return the other two.
                    return (edge1[i ^ 1], edge2[j ^ 1])

        # No extremities match, return None.
        return None

    @staticmethod
    def _add_edge_append(dictionary, edge, elem):
        """
        If ``edge`` is not ``dictionary``, adds ``edge`` to the
        dictionary and sets its value to ``[elem]``.  Otherwise
        appends ``elem`` to the value of existing entry.

        Note that edges are undirected, thus `(A, B) = (B, A)`.
        """
        (A, B) = edge
        if (A, B) in dictionary:
            dictionary[(A, B)].append(elem)
        elif (B, A) in dictionary:
             dictionary[(B, A)].append(elem)
        else:
            dictionary[(A, B)] = [elem]

    @staticmethod
    def _build_skeleton(morphisms):
        """
        Creates a dictionary which maps edges to corresponding
        morphisms.  Thus for a morphism `f:A\rightarrow B`, the edge
        `(A, B)` will be associated with `f`.  This function also adds
        to the list those edges which are formed by juxtaposition of
        two edges already in the list.  These new edges are not
        associated with any morphism and are only added to assure that
        the diagram can be decomposed into triangles.
        """
        edges = {}
        # Create edges for morphisms.
        for morphism in morphisms:
            DiagramGrid._add_edge_append(
                edges, (morphism.domain,morphism.codomain), morphism)

        # Create new edges by juxtaposing existing edges.
        edges1 = dict(edges)
        for w in edges1:
            for v in edges1:
                wv = DiagramGrid._juxtapose_edges(w, v)
                if wv:
                    (A, B) = wv
                    if ((A, B) not in edges) and ((B, A) not in edges):
                        edges[(A, B)] = []

        return edges

    @staticmethod
    def _list_triangles(edges):
        """
        Builds the set of triangles formed by the supplied edges.  The
        triangles are arbitrary and need not be commutative.  A
        triangle is a set contains all three sides.
        """
        triangles = []

        for w in edges:
            for v in edges:
                wv = DiagramGrid._juxtapose_edges(w, v)
                if wv:
                    (A, B) = wv
                    if (A, B) in edges:
                        triangle = FiniteSet(w, v, (A, B))
                        triangles.append(triangle)
                    elif (B, A) in edges:
                        triangle = FiniteSet(w, v, (B, A))
                        triangles.append(triangle)

        return FiniteSet(triangles)

    @staticmethod
    def _drop_redundant_triangles(triangles, skeleton):
        """
        Returns a list which contains only those triangles who have
        morphisms associated with at least two edges.
        """
        return [tri for tri in triangles if \
                len([e for e in tri if skeleton[e]]) >= 2]

    @staticmethod
    def _morphism_length(morphism):
        """
        Returns the length of a morphism.  The length of a morphism is
        the number of components it consists of.  A non-composite
        morphism is of length 1.
        """
        if isinstance(morphism, CompositeMorphism):
            return len(morphism.components)
        else:
            return 1

    @staticmethod
    def _compute_triangle_min_sizes(triangles, edges):
        """
        Returns a dictionary mapping triangles to their minimal sizes.
        The minimal size of a triangle is the sum of maximal lengths
        of morphisms associated to the sides of the triangle.  The
        length of a morphism is the number of components it consists
        of.  A non-composite morphism is of length 1.
        """
        triangle_sizes = {}
        for triangle in triangles:
            size = 0
            for e in triangle:
                morphisms = edges[e]
                if morphisms:
                    size += max([DiagramGrid._morphism_length(m) \
                                 for m in morphisms])
            triangle_sizes[triangle] = size
        return triangle_sizes

    @staticmethod
    def _triangle_objects(triangle):
        """
        Given a triangle, returns the objects included in it.
        """
        return FiniteSet(sum(triangle, () ))

    @staticmethod
    def _other_vertex(triangle, edge):
        """
        Given a triangle and an edge of it, returns the vertex which
        opposes the edge.
        """
        return (DiagramGrid._triangle_objects(triangle) - FiniteSet(edge)).args[0]

    @staticmethod
    def _find_triangle_welding(triangle, fringe, grid):
        """
        Finds if possible an edge in the fringe to which the supplied
        triangle could be attached and returns the index of the
        corresponding edge in the fringe.

        This function relies on the fact that objects are unique in
        the diagram.
        """
        for (a, b) in fringe:
            if ((grid[a], grid[b]) in triangle) or \
               ((grid[b], grid[a]) in triangle):
                return (a, b)
        return None

    @staticmethod
    def _fix_degerate(edge, pt1, pt2, fringe, grid):
        """
        Checks if either of the points ``pt1`` or ``pt2`` forms a
        perpendicular edge on ``edge``, which is in the fringe.  If
        this is indeed so, replaces the two edges with the
        corresponding diagonal and returns ``True``.  Otherwise
        returns ``False``.
        """
        (a, b) = edge

        if (a, pt1) in fringe:
            fringe.remove((a, b))
            fringe.remove((a, pt1))
            fringe.append((pt1, b))
            return True
        elif (pt1, a) in fringe:
            fringe.remove((a, b))
            fringe.remove((pt1, a))
            fringe.append((pt1, b))
            return True
        elif (b, pt2) in fringe:
            fringe.remove((a, b))
            fringe.remove((b, pt2))
            fringe.append((a, pt2))
            return True
        elif (pt2, b) in fringe:
            fringe.remove((a, b))
            fringe.remove((pt2, b))
            fringe.append((a, pt2))
            return True

        return False

    @staticmethod
    def _empty_point(pt, grid):
        """
        Checks if the cell at coordinates ``pt`` is either empty or
        out of the bounds of the grid.
        """
        if (pt[0] < 0) or (pt[1] < 0) or \
           (pt[0] >= grid.height) or (pt[1] >= grid.width):
            return True
        return grid[pt] is None

    @staticmethod
    def _weld_triangle(triangles, fringe, grid, skeleton):
        """
        Welds a triangle to the fringe and returns ``True``, if
        possible.  Otherwise returns ``False``.
        """
        for tri in triangles:
            welding_edge = DiagramGrid._find_triangle_welding(
                tri, fringe, grid)
            if not welding_edge:
                continue

            a, b = welding_edge
            target_cell = None

            obj = DiagramGrid._other_vertex(tri, (grid[a], grid[b]))

            # We now have a triangle and an wedge where it can be
            # welded to the fringe.  Decide where to place the
            # other vertex of the triangle and check for
            # degenerate situations en route.

            if (abs(a[0] - b[0]) == 1) and (abs(a[1] - b[1]) == 1):
                # A diagonal edge.
                target_cell = (a[0], b[1])
                if grid[target_cell]:
                    # That cell is already occupied.
                    target_cell = (b[0], a[1])

                    if grid[target_cell]:
                        # Degenerate situation, this edge is not
                        # on the actual fringe.  Correct the
                        # fringe and go on.
                        fringe.remove((a, b))
                        break
            elif a[0] == b[0]:
                # A horizontal edge.  Suppose a triangle can be
                # built in the upward direction.

                up_left = a[0] - 1, a[1]
                up_right = a[0] - 1, b[1]

                if DiagramGrid._fix_degerate(
                    (a, b), up_left, up_right, fringe, grid):
                    # The fringe has just been corrected.  Restart.
                    break

                if DiagramGrid._empty_point(up_left, grid):
                    target_cell = up_left
                elif DiagramGrid._empty_point(up_right, grid):
                    target_cell = up_right
                else:
                    # No room above this edge.  Check below.
                    down_left = a[0] + 1, a[1]
                    down_right = a[0] + 1, b[1]

                    if DiagramGrid._fix_degerate(
                        (a, b), up_left, up_right, fringe, grid):
                        # The fringe has just been corrected.
                        # Restart.
                        break

                    if DiagramGrid._empty_point(up_left, grid):
                        target_cell = down_left
                    elif DiagramGrid._empty_point(up_right, grid):
                        target_cell = down_right
                    else:
                        # This edge is not in the fringe, remove it
                        # and restart.
                        fringe.remove((a, b))
                        break

            elif a[1] == b[1]:
                # A vertical edge.  Suppose a triangle can be built to
                # the left.
                left_up = a[0], a[1] - 1
                left_down = b[0], a[1] - 1

                if DiagramGrid._fix_degerate(
                    (a, b), left_up, left_down, fringe, grid):
                    # The fringe has just been corrected.  Restart.
                    break

                if DiagramGrid._empty_point(left_up, grid):
                    target_cell = left_up
                elif DiagramGrid._empty_point(left_down, grid):
                    target_cell = left_down
                else:
                    # No room to the left.  See what's to the right.
                    right_up = a[0], a[1] + 1
                    right_down = b[0], a[1] + 1

                    if DiagramGrid._fix_degerate(
                        (a, b), right_up, right_down, fringe, grid):
                        # The fringe has just been corrected.  Restart.
                        break

                    if DiagramGrid._empty_point(right_up, grid):
                        target_cell = right_up
                    elif DiagramGrid._empty_point(right_down, grid):
                        target_cell = right_down
                    else:
                        # This edge is not in the fringe, remove it
                        # and restart.
                        fringe.remove((a, b))
                        break

    def __new__(cls, diagram):
        premises = DiagramGrid._simplify_morphisms(diagram.premises)
        conclusions = DiagramGrid._simplify_morphisms(diagram.conclusions)
        merged_morphisms = DiagramGrid._merge_premises_conclusions(
            premises, conclusions)

        skeleton = DiagramGrid._build_skeleton(merged_morphisms)

        triangles = DiagramGrid._list_triangles(skeleton)
        triangles = DiagramGrid._drop_redundant_triangles(triangles, skeleton)
        triangle_sizes = DiagramGrid._compute_triangle_min_sizes(
            triangles, skeleton)

        triangles = sorted(triangles, key=lambda triangle: triangle_sizes[triangle])

        grid = _GrowableGrid(2, 1)

        # Place the first edge on the grid.
        root_edge = [e for e in triangles[0] if skeleton[e]][0]
        grid[0, 0], grid[0, 1] = root_edge
        fringe = [((0,0), (0, 1))]

        while DiagramGrid._weld_triangle(triangles, fringe, grid, skeleton):
            pass
