"""
This module contains the functionality to arrange the nodes of a
diagram on an abstract grid, and then to produce a graphical
representation of the grid.

The currently supported back-ends are Xy-pic [Xypic]

[Xypic] http://www.tug.org/applications/Xy-pic/
"""

from sympy.core import Basic, FiniteSet
from sympy.categories import CompositeMorphism, IdentityMorphism
from sympy.utilities import default_sort_key

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

class DiagramGrid:
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
    def _put_object(coords, obj, grid, fringe):
        """
        Places an object at the coordinate ``cords`` in ``grid``,
        growing the grid and updating ``fringe``, if necessary.
        Returns (0, 0) if no row or column has been prepended, (1, 0)
        if a row was prepended, (0, 1) if a column was prepended and
        (1, 1) if both a column and a row were prepended.
        """
        (i, j) = coords
        offset = (0, 0)
        if i == -1:
            grid.prepend_row()
            i = 0
            offset = (1, offset[1])
            for k in xrange(len(fringe)):
                ((i1, j1), (i2, j2)) = fringe[k]
                fringe[k] = ((i1 + 1, j1), (i2 + 1, j2))
        elif i == grid.height:
            grid.append_row()

        if j == -1:
            j = 0
            offset = (offset[0], 1)
            grid.prepend_column()
            for k in xrange(len(fringe)):
                ((i1, j1), (i2, j2)) = fringe[k]
                fringe[k] = ((i1, j1 + 1), (i2, j2 + 1))
        elif j == grid.width:
            grid.append_column()

        grid[i, j] = obj
        return offset

    @staticmethod
    def _choose_target_cell(pt1, pt2, edge, obj, skeleton, grid):
        """
        Given two points, ``pt1`` and ``pt2``, and the welding edge
        ``edge``, chooses one of the two points to place the opposing
        vertex ``obj`` of the triangle.  If neither of this points
        fits, returns ``None``.
        """
        pt1_empty = DiagramGrid._empty_point(pt1, grid)
        pt2_empty = DiagramGrid._empty_point(pt2, grid)

        if pt1_empty and pt2_empty:
            # Both cells are empty.  Of them two, choose that cell
            # which will assure that a visible edge of the triangle
            # will be drawn perpendicularly to the current welding
            # edge.

            A = grid[edge[0]]
            B = grid[edge[1]]

            if skeleton.get((A, obj)) or skeleton.get((obj, A)):
                return pt1
            else:
                return pt2
        if pt1_empty:
            return pt1
        elif pt2_empty:
            return pt2
        else:
            return None

    # The possible return values of ``_weld_triangle``.
    _WELDING_FAILURE = 1
    _RESTART = 2

    @staticmethod
    def _weld_triangle(triangles, fringe, grid, skeleton):
        """
        If possible, welds a triangle to the fringe and returns the
        welded triangle.  Otherwise returns ``_WELDING_FAILURE``.  If
        this method encounters a degenerate situation and corrects the
        fringe such that a restart of the search is required, it
        returns ``_RESTART``.
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
                        return DiagramGrid._RESTART
            elif a[0] == b[0]:
                # A horizontal edge.  Suppose a triangle can be
                # built in the downward direction.

                down_left = a[0] + 1, a[1]
                down_right = a[0] + 1, b[1]

                target_cell = DiagramGrid._choose_target_cell(
                    down_left, down_right, (a, b), obj, skeleton, grid)

                if not target_cell:
                    # No room below this edge.  Check above.
                    up_left = a[0] - 1, a[1]
                    up_right = a[0] - 1, b[1]

                    target_cell = DiagramGrid._choose_target_cell(
                        up_left, up_right, (a, b), obj, skeleton, grid)

                    if not target_cell:
                        # This edge is not in the fringe, remove it
                        # and restart.
                        fringe.remove((a, b))
                        return DiagramGrid._RESTART

            elif a[1] == b[1]:
                # A vertical edge.  Suppose a triangle can be built to
                # the right.
                right_up = a[0], a[1] + 1
                right_down = b[0], a[1] + 1

                target_cell = DiagramGrid._choose_target_cell(
                    right_up, right_down, (a, b), obj, skeleton, grid)

                if not target_cell:
                    # No room to the left.  See what's to the right.
                    left_up = a[0], a[1] - 1
                    left_down = b[0], a[1] - 1

                    target_cell = DiagramGrid._choose_target_cell(
                        left_up, left_down, (a, b), obj, skeleton, grid)

                    if not target_cell:
                        # This edge is not in the fringe, remove it
                        # and restart.
                        fringe.remove((a, b))
                        return DiagramGrid._RESTART

            # We now know where to place the other vertex of the
            # triangle.
            offset = DiagramGrid._put_object(target_cell, obj, grid, fringe)

            # Take care of the displacement of coordinates if a row or
            # a column was prepended.
            target_cell = (target_cell[0] + offset[0],
                           target_cell[1] + offset[1])
            a = (a[0] + offset[0], a[1] + offset[1])
            b = (b[0] + offset[0], b[1] + offset[1])

            fringe.extend([(a, target_cell), (b, target_cell)])

            triangles.remove(tri)

            return tri

        return DiagramGrid._WELDING_FAILURE

    @staticmethod
    def _triangle_key(tri, triangle_sizes):
        """
        Returns a key for the supplied triangle.  It should be the
        same independently of the hash randomisation.
        """
        objects = sorted([x.name for x in DiagramGrid._triangle_objects(tri)])
        return (triangle_sizes[tri], objects)

    @staticmethod
    def _pick_root_edge(tri, skeleton):
        """
        For a given triangle always picks the same root edge.  The
        root edge is the edge that will be placed first on the grid.
        """
        candidates = [e for e in tri if skeleton[e]]
        sorted_candidates = sorted(candidates, key=default_sort_key)
        return sorted_candidates[0]

    @staticmethod
    def _drop_irrelevant_triangles(triangles, placed_objects):
        """
        Returns only those triangles whose set of objects is not
        completely included in ``placed_objects``.
        """
        return [tri for tri in triangles if not placed_objects.subset(
            DiagramGrid._triangle_objects(tri))]

    @staticmethod
    def _grow_pseudopod(triangles, fringe, grid, skeleton, placed_objects):
        """
        Starting from an object in the existing structure on the grid,
        adds an edge to which a triangle from ``triangles`` could be
        welded.  If this method has found a way to do so, it returns
        the object it has just added.

        This method should be applied when ``_weld_triangle`` cannot
        find weldings any more.
        """
        for i in xrange(grid.height):
            for j in xrange(grid.width):
                obj = grid[i, j]
                if not obj:
                    continue

                # Here we need to choose a triangle which has only
                # ``obj`` in common with the existing structure.  The
                # situations when this is not possible should be
                # handled elsewhere.

                def good_triangle(tri):
                    objs = DiagramGrid._triangle_objects(tri)
                    return objs.contains(obj) and \
                           placed_objects & (objs - FiniteSet(obj)) == FiniteSet()

                tris = [tri for tri in triangles if good_triangle(tri)]
                if not tris:
                    # This object is not interesting.
                    continue

                # We have found a triangle which could be attached to
                # the existing structure by a vertex.

                candidates = sorted([e for e in tri if skeleton[e]],
                                    key=default_sort_key)
                edges = [e for e in candidates if obj in e]
                if not edges:
                    # No meaningful edge could be drawn from this
                    # object, so we are not interested.
                    continue

                # Get the object at the other end of the edge.
                edge = edges[0]
                other_obj = edge[0]
                if other_obj == obj:
                    other_obj = edge[1]

                # Now check for free directions.  When checking for
                # free directions, prefer the horizontal and vertical
                # directions.
                neighbours = [(i-1, j), (i, j+1), (i+1, j), (i, j-1)] + \
                             [(i-1,j-1), (i-1, j+1), (i+1,j-1), (i+1, j+1)]

                for pt in neighbours:
                    if DiagramGrid._empty_point(pt, grid):
                        # We have a found a place to grow the
                        # pseudopod into.
                        offset = DiagramGrid._put_object(
                            pt, other_obj, grid, fringe)

                        i += offset[0]
                        j += offset[1]
                        pt = (pt[0] + offset[0], pt[1] + offset[1])
                        fringe.append(((i, j), pt))

                        return other_obj

    def __init__(self, diagram):
        premises = DiagramGrid._simplify_morphisms(diagram.premises)
        conclusions = DiagramGrid._simplify_morphisms(diagram.conclusions)
        merged_morphisms = DiagramGrid._merge_premises_conclusions(
            premises, conclusions)

        skeleton = DiagramGrid._build_skeleton(merged_morphisms)

        all_objects = diagram.objects
        grid = _GrowableGrid(2, 1)

        if len(skeleton) == 1:
            # This diagram contains only one morphism.  Draw it
            # horizontally.
            objects = sorted(all_objects, key=default_sort_key)
            grid[0, 0] = objects[0]
            grid[0, 1] = objects[1]

            self._grid = grid

            return

        triangles = DiagramGrid._list_triangles(skeleton)
        triangles = DiagramGrid._drop_redundant_triangles(triangles, skeleton)
        triangle_sizes = DiagramGrid._compute_triangle_min_sizes(
            triangles, skeleton)

        triangles = sorted(triangles, key=lambda tri:
                           DiagramGrid._triangle_key(tri, triangle_sizes))

        # Place the first edge on the grid.
        root_edge = DiagramGrid._pick_root_edge(triangles[0], skeleton)
        grid[0, 0], grid[0, 1] = root_edge
        fringe = [((0,0), (0, 1))]

        # Record which objects we now have on the grid.
        placed_objects = FiniteSet(root_edge)

        while placed_objects != all_objects:
            triangle = DiagramGrid._weld_triangle(
                triangles, fringe, grid, skeleton)

            if triangle == DiagramGrid._RESTART:
                # ``_weld_triangle`` wants to have the search for
                # welding restarted.
                continue

            if triangle == DiagramGrid._WELDING_FAILURE:
                # No more weldings found.  Try to attach triangles by
                # vertices.
                new_obj = DiagramGrid._grow_pseudopod(
                    triangles, fringe, grid, skeleton, placed_objects)

                placed_objects = placed_objects | FiniteSet(new_obj)

                # Now, hopefully, a new welding will be found.
            else:
                placed_objects = placed_objects | \
                                 DiagramGrid._triangle_objects(triangle)

            triangles = DiagramGrid._drop_irrelevant_triangles(
                triangles, placed_objects)

        self._grid = grid

    @property
    def width(self):
        """
        Returns the number of columns in this diagram layout.

        Examples
        ========
        TODO: Add examples.
        """
        return self._grid.width

    @property
    def height(self):
        """
        Returns the number of rows in this diagram layout.

        Examples
        ========
        TODO: Add examples.
        """
        return self._grid.height

    def __getitem__(self, (i, j)):
        """
        Returns the object placed in the row ``i`` and column ``j``.
        The indices are 0-based.
        """
        return self._grid[i, j]
