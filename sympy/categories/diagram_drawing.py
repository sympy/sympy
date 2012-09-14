r"""
This module contains the functionality to arrange the nodes of a
diagram on an abstract grid, and then to produce a graphical
representation of the grid.

The currently supported back-ends are Xy-pic [Xypic].

Layout Algorithm
================

This section provides an overview of the algorithms implemented in
:class:`DiagramGrid` to lay out diagrams.

The first step of the algorithm is the removal composite and identity
morphisms which do not have properties in the supplied diagram.  The
premises and conclusions of the diagram are then merged.

The generic layout algorithm begins with the construction of the
"skeleton" of the diagram.  The skeleton is an undirected graph which
has the objects of the diagram as vertices and has an (undirected)
edge between each pair of objects between which there exist morphisms.
The direction of the morphisms does not matter at this stage.  The
skeleton also includes an edge between each pair of vertices `A` and
`C` such that there exists an object `B` which is connected via
a morphism to `A`, and via a morphism to `C`.

The skeleton constructed in this way has the property that every
object is a vertex of a triangle formed by three edges of the
skeleton.  This property lies at the base of the generic layout
algorithm.

After the skeleton has been constructed, the algorithm lists all
triangles which can be formed.  Note that some triangles will not have
all edges corresponding to morphisms which will actually be drawn.
Triangles which have only one edge or less which will actually be
drawn are immediately discarded.

The list of triangles is sorted according to the number of edges which
correspond to morphisms, then the triangle with the least number of such
edges is selected.  One of such edges is picked and the corresponding
objects are placed horizontally, on a grid.  This edge is recorded to
be in the fringe.  The algorithm then finds a "welding" of a triangle
to the fringe.  A welding is an edge in the fringe where a triangle
could be attached.  If the algorithm succeeds in finding such a
welding, it adds to the grid that vertex of the triangle which was not
yet included in any edge in the fringe and records the two new edges in
the fringe.  This process continues iteratively until all objects of
the diagram has been placed or until no more weldings can be found.

An edge is only removed from the fringe when a welding to this edge
has been found, and there is no room around this edge to place
another vertex.

When no more weldings can be found, but there are still triangles
left, the algorithm searches for a possibility of attaching one of the
remaining triangles to the existing structure by a vertex.  If such a
possibility is found, the corresponding edge of the found triangle is
placed in the found space and the iterative process of welding
triangles restarts.

When logical groups are supplied, each of these groups is laid out
independently.  Then a diagram is constructed in which groups are
objects and any two logical groups between which there exist morphisms
are connected via a morphism.  This diagram is laid out.  Finally,
the grid which includes all objects of the initial diagram is
constructed by replacing the cells which contain logical groups with
the corresponding laid out grids, and by correspondingly expanding the
rows and columns.

The sequential layout algorithm begins by constructing the
underlying undirected graph defined by the morphisms obtained after
simplifying premises and conclusions and merging them (see above).
The vertex with the minimal degree is then picked up and depth-first
search is started from it.  All objects which are located at distance
`n` from the root in the depth-first search tree, are positioned in
the `n`-th column of the resulting grid.  The sequential layout will
therefore attempt to lay the objects out along a line.

References
==========

[Xypic] http://www.tug.org/applications/Xy-pic/
"""

from sympy.core import Basic, FiniteSet, Dict
from sympy.categories import (CompositeMorphism, IdentityMorphism,
                              NamedMorphism, Diagram)
from sympy.utilities import default_sort_key
from itertools import chain
from sympy.core.compatibility import iterable

class _GrowableGrid(object):
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

class DiagramGrid(object):
    r"""
    Constructs and holds the fitting of the diagram into a grid.

    The mission of this class is to analyse the structure of the
    supplied diagram and to place its objects on a grid such that,
    when the objects and the morphisms are actually drawn, the diagram
    would be "readable", in the sense that there will not be many
    intersections of moprhisms.  This class does not perform any
    actual drawing.  It does strive nevertheless to offer sufficient
    metadata to draw a diagram.

    Consider the following simple diagram.

    >>> from sympy.categories import Object, NamedMorphism
    >>> from sympy.categories import Diagram, DiagramGrid
    >>> from sympy import pprint
    >>> A = Object("A")
    >>> B = Object("B")
    >>> C = Object("C")
    >>> f = NamedMorphism(A, B, "f")
    >>> g = NamedMorphism(B, C, "g")
    >>> diagram = Diagram([f, g])

    The simplest way to have a diagram laid out is the following:

    >>> grid = DiagramGrid(diagram)
    >>> (grid.width, grid.height)
    (2, 2)
    >>> pprint(grid)
    A  B
    <BLANKLINE>
       C

    Sometimes one sees the diagram as consisting of logical groups.
    One can advise ``DiagramGrid`` as to such groups by employing the
    ``groups`` keyword argument.

    Consider the following diagram:

    >>> D = Object("D")
    >>> f = NamedMorphism(A, B, "f")
    >>> g = NamedMorphism(B, C, "g")
    >>> h = NamedMorphism(D, A, "h")
    >>> k = NamedMorphism(D, B, "k")
    >>> diagram = Diagram([f, g, h, k])

    Lay it out with generic layout:

    >>> grid = DiagramGrid(diagram)
    >>> pprint(grid)
    A  B  D
    <BLANKLINE>
       C

    Now, we can group the objects `A` and `D` to have them near one
    another:

    >>> grid = DiagramGrid(diagram, groups=[[A, D], B, C])
    >>> pprint(grid)
    B     C
    <BLANKLINE>
    A  D

    Note how the positioning of the other objects changes.

    Further indications can be supplied to the constructor of
    :class:`DiagramGrid` using keyword arguments.  The currently
    supported hints are explained in the following paragraphs.

    :class:`DiagramGrid` does not automatically guess which layout
    would suit the supplied diagram better.  Consider, for example,
    the following linear diagram:

    >>> E = Object("E")
    >>> f = NamedMorphism(A, B, "f")
    >>> g = NamedMorphism(B, C, "g")
    >>> h = NamedMorphism(C, D, "h")
    >>> i = NamedMorphism(D, E, "i")
    >>> diagram = Diagram([f, g, h, i])

    When laid out with the generic layout, it does not get to look
    linear:

    >>> grid = DiagramGrid(diagram)
    >>> pprint(grid)
    A  B
    <BLANKLINE>
       C  D
    <BLANKLINE>
          E

    To get it laid out in a line, use ``layout="sequential"``:

    >>> grid = DiagramGrid(diagram, layout="sequential")
    >>> pprint(grid)
    A  B  C  D  E

    One may sometimes need to transpose the resulting layout.  While
    this can always be done by hand, :class:`DiagramGrid` provides a
    hint for that purpose:

    >>> grid = DiagramGrid(diagram, layout="sequential", transpose=True)
    >>> pprint(grid)
    A
    <BLANKLINE>
    B
    <BLANKLINE>
    C
    <BLANKLINE>
    D
    <BLANKLINE>
    E

    Separate hints can also be provided for each group.  For an
    example, refer to ``tests/test_drawing.py``, and see the different
    ways in which the five lemma [FiveLemma] can be laid out.

    See Also
    ========
    Diagram

    References
    ==========
    [FiveLemma] http://en.wikipedia.org/wiki/Five_lemma
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
        return dict(chain(premises.items(), conclusions.items()))

    @staticmethod
    def _juxtapose_edges(edge1, edge2):
        """
        If ``edge1`` and ``edge2`` have precisely one common endpoint,
        returns an edge which would form a triangle with ``edge1`` and
        ``edge2``.

        If ``edge1`` and ``edge2`` don't have a common endpoint,
        returns ``None``.

        If ``edge1`` and ``edge`` are the same edge, returns ``None``.
        """
        intersection = edge1 & edge2
        if len(intersection) != 1:
            # The edges either have no common points or are equal.
            return None

        # The edges have a common endpoint.  Extract the different
        # endpoints and set up the new edge.
        return (edge1 - intersection) | (edge2 - intersection)

    @staticmethod
    def _add_edge_append(dictionary, edge, elem):
        """
        If ``edge`` is not in ``dictionary``, adds ``edge`` to the
        dictionary and sets its value to ``[elem]``.  Otherwise
        appends ``elem`` to the value of existing entry.

        Note that edges are undirected, thus `(A, B) = (B, A)`.
        """
        if edge in dictionary:
            dictionary[edge].append(elem)
        else:
            dictionary[edge] = [elem]

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
                edges, frozenset([morphism.domain,morphism.codomain]), morphism)

        # Create new edges by juxtaposing existing edges.
        edges1 = dict(edges)
        for w in edges1:
            for v in edges1:
                wv = DiagramGrid._juxtapose_edges(w, v)
                if wv and wv not in edges:
                    edges[wv] = []

        return edges

    @staticmethod
    def _list_triangles(edges):
        """
        Builds the set of triangles formed by the supplied edges.  The
        triangles are arbitrary and need not be commutative.  A
        triangle is a set that contains all three of its sides.
        """
        triangles = set()

        for w in edges:
            for v in edges:
                wv = DiagramGrid._juxtapose_edges(w, v)
                if wv and wv in edges:
                    triangles.add(frozenset([w, v, wv]))

        return triangles

    @staticmethod
    def _drop_redundant_triangles(triangles, skeleton):
        """
        Returns a list which contains only those triangles who have
        morphisms associated with at least two edges.
        """
        return [tri for tri in triangles
                if len([e for e in tri if skeleton[e]]) >= 2]

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
        r"""
        Returns a dictionary mapping triangles to their minimal sizes.
        The minimal size of a triangle is the sum of maximal lengths
        of morphisms associated to the sides of the triangle.  The
        length of a morphism is the number of components it consists
        of.  A non-composite morphism is of length 1.

        Sorting triangles by this metric attempts to address two
        aspects of layout.  For triangles with only simple morphisms
        in the edge, this assures that triangles with all three edges
        visible will get typeset after triangles with less visible
        edges, which sometimes minimises the necessity in diagonal
        arrows.  For triangles with composite morphisms in the edges,
        this assures that objects connected with shorter morphisms
        will be laid out first, resulting the visual proximity of
        those objects which are connected by shorter morphisms.
        """
        triangle_sizes = {}
        for triangle in triangles:
            size = 0
            for e in triangle:
                morphisms = edges[e]
                if morphisms:
                    size += max(DiagramGrid._morphism_length(m)
                                for m in morphisms)
            triangle_sizes[triangle] = size
        return triangle_sizes

    @staticmethod
    def _triangle_objects(triangle):
        """
        Given a triangle, returns the objects included in it.
        """
        # A triangle is a frozenset of three two-element frozensets
        # (the edges).  This chains the three edges together and
        # creates a frozenset from the iterator, thus producing a
        # frozenset of objects of the triangle.
        return frozenset(chain(*tuple(triangle)))

    @staticmethod
    def _other_vertex(triangle, edge):
        """
        Given a triangle and an edge of it, returns the vertex which
        opposes the edge.
        """
        # This gets the set of objects of the triangle and then
        # subtracts the set of objects employed in ``edge`` to get the
        # vertex opposite to ``edge``.
        return list(DiagramGrid._triangle_objects(triangle) - set(edge))[0]

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
            offset = (1, 0)
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
            # Both cells are empty.  Of these two, choose that cell
            # which will assure that a visible edge of the triangle
            # will be drawn perpendicularly to the current welding
            # edge.

            A = grid[edge[0]]
            B = grid[edge[1]]

            if skeleton.get(frozenset([A, obj])):
                return pt1
            else:
                return pt2
        if pt1_empty:
            return pt1
        elif pt2_empty:
            return pt2
        else:
            return None

    @staticmethod
    def _find_triangle_to_weld(triangles, fringe, grid):
        """
        Finds, if possible, a triangle and an edge in the fringe to
        which the triangle could be attached.  Returns the tuple
        containing the triangle and the index of the corresponding
        edge in the fringe.

        This function relies on the fact that objects are unique in
        the diagram.
        """
        for triangle in triangles:
            for (a, b) in fringe:
                if frozenset([grid[a], grid[b]]) in triangle:
                    return (triangle, (a, b))
        return None

    @staticmethod
    def _weld_triangle(tri, welding_edge, fringe, grid, skeleton):
        """
        If possible, welds the triangle ``tri`` to ``fringe`` and
        returns ``False``.  If this method encounters a degenerate
        situation in the fringe and corrects it such that a restart of
        the search is required, it returns ``True`` (which means that
        a restart in finding triangle weldings is required).

        A degenerate situation is a situation when an edge listed in
        the fringe does not belong to the visual boundary of the
        diagram.
        """
        a, b = welding_edge
        target_cell = None

        obj = DiagramGrid._other_vertex(tri, (grid[a], grid[b]))

        # We now have a triangle and an edge where it can be welded to
        # the fringe.  Decide where to place the other vertex of the
        # triangle and check for degenerate situations en route.

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
                    return True
        elif a[0] == b[0]:
            # A horizontal edge.  We first attempt to build the
            # triangle in the downward direction.

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
                    return True
        elif a[1] == b[1]:
            # A vertical edge.  We will attempt to place the other
            # vertex of the triangle to the right of this edge.
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
                    return True

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

        # No restart is required.
        return False

    @staticmethod
    def _triangle_key(tri, triangle_sizes):
        """
        Returns a key for the supplied triangle.  It should be the
        same independently of the hash randomisation.
        """
        objects = sorted(DiagramGrid._triangle_objects(tri), key=default_sort_key)
        return (triangle_sizes[tri], default_sort_key(objects))

    @staticmethod
    def _pick_root_edge(tri, skeleton):
        """
        For a given triangle always picks the same root edge.  The
        root edge is the edge that will be placed first on the grid.
        """
        candidates = [sorted(e, key=default_sort_key)
                      for e in tri if skeleton[e]]
        sorted_candidates = sorted(candidates, key=default_sort_key)
        # Don't forget to assure the proper ordering of the vertices
        # in this edge.
        return tuple(sorted(sorted_candidates[0], key=default_sort_key))

    @staticmethod
    def _drop_irrelevant_triangles(triangles, placed_objects):
        """
        Returns only those triangles whose set of objects is not
        completely included in ``placed_objects``.
        """
        return [tri for tri in triangles if not placed_objects.issuperset(
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
                    return obj in objs and \
                           placed_objects & (objs - set([obj])) == set()

                tris = [tri for tri in triangles if good_triangle(tri)]
                if not tris:
                    # This object is not interesting.
                    continue

                # Pick the "simplest" of the triangles which could be
                # attached.  Remember that the list of triangles is
                # sorted according to their "simplicity" (see
                # _compute_triangle_min_sizes for the metric).
                #
                # Note that ``tris`` are sequentially built from
                # ``triangles``, so we don't have to worry about hash
                # randomisation.
                tri = tris[0]

                # We have found a triangle which could be attached to
                # the existing structure by a vertex.

                candidates = sorted([e for e in tri if skeleton[e]],
                                    key=lambda e: FiniteSet(e).sort_key())
                edges = [e for e in candidates if obj in e]

                # Note that a meaningful edge (i.e., and edge that is
                # associated with a morphism) containing ``obj``
                # always exists.  That's because all triangles are
                # guaranteed to have at least two meaningful edges.
                # See _drop_redundant_triangles.

                # Get the object at the other end of the edge.
                edge = edges[0]
                other_obj = tuple(edge - frozenset([obj]))[0]

                # Now check for free directions.  When checking for
                # free directions, prefer the horizontal and vertical
                # directions.
                neighbours = [(i-1, j), (i, j+1), (i+1, j), (i, j-1),
                              (i-1,j-1), (i-1, j+1), (i+1,j-1), (i+1, j+1)]

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

        # This diagram is actually cooler that I can handle.  Fail cowardly.
        return None

    @staticmethod
    def _handle_groups(diagram, groups, merged_morphisms, hints):
        """
        Given the slightly preprocessed morphisms of the diagram,
        produces a grid laid out according to ``groups``.

        If a group has hints, it is laid out with those hints only,
        without any influence from ``hints``.  Otherwise, it is laid
        out with ``hints``.
        """
        def lay_out_group(group, local_hints):
            """
            If ``group`` is a set of objects, uses a ``DiagramGrid``
            to lay it out and returns the grid.  Otherwise returns the
            object (i.e., ``group``).  If ``local_hints`` is not
            empty, it is supplied to ``DiagramGrid`` as the dictionary
            of hints.  Otherwise, the ``hints`` argument of
            ``_handle_groups`` is used.
            """
            if isinstance(group, FiniteSet):
                # Set up the corresponding object-to-group
                # mappings.
                for obj in group:
                    obj_groups[obj] = group

                # Lay out the current group.
                if local_hints:
                    groups_grids[group] = DiagramGrid(
                        diagram.subdiagram_from_objects(group), **local_hints)
                else:
                    groups_grids[group] = DiagramGrid(
                        diagram.subdiagram_from_objects(group), **hints)
            else:
                obj_groups[group] = group

        def group_to_finiteset(group):
            """
            Converts ``group`` to a :class:``FiniteSet`` if it is an
            iterable.
            """
            if iterable(group):
                return FiniteSet(group)
            else:
                return group

        obj_groups = {}
        groups_grids = {}

        # We would like to support various containers to represent
        # groups.  To achieve that, before laying each group out, it
        # should be converted to a FiniteSet, because that is what the
        # following code expects.

        if isinstance(groups, dict) or isinstance(groups, Dict):
            finiteset_groups = {}
            for group, local_hints in groups.items():
                finiteset_group = group_to_finiteset(group)
                finiteset_groups[finiteset_group] = local_hints
                lay_out_group(group, local_hints)
            groups = finiteset_groups
        else:
            finiteset_groups = []
            for group in groups:
                finiteset_group = group_to_finiteset(group)
                finiteset_groups.append(finiteset_group)
                lay_out_group(finiteset_group, None)
            groups = finiteset_groups

        new_morphisms = []
        for morphism in merged_morphisms:
            dom = obj_groups[morphism.domain]
            cod = obj_groups[morphism.codomain]
            # Note that we are not really interested in morphisms
            # which do not employ two different groups, because
            # these do not influence the layout.
            if dom != cod:
                # These are essentially unnamed morphisms; they are
                # not going to mess in the final layout.  By giving
                # them the same names, we avoid unnecessary
                # duplicates.
                new_morphisms.append(NamedMorphism(dom, cod, "dummy"))

        # Lay out the new diagram.  Since these are dummy morphisms,
        # properties and conclusions are irrelevant.
        top_grid = DiagramGrid(Diagram(new_morphisms))

        # We now have to substitute the groups with the corresponding
        # grids, laid out at the beginning of this function.  Compute
        # the size of each row and column in the grid, so that all
        # nested grids fit.

        def group_size(group):
            """
            For the supplied group (or object, eventually), returns
            the size of the cell that will hold this group (object).
            """
            if group in groups_grids:
                grid = groups_grids[group]
                return (grid.height, grid.width)
            else:
                return (1, 1)

        row_heights = [max(group_size(top_grid[i, j])[0]
                           for j in xrange(top_grid.width))
                       for i in xrange(top_grid.height)]

        column_widths = [max(group_size(top_grid[i, j])[1]
                             for i in xrange(top_grid.height))
                         for j in xrange(top_grid.width)]

        grid = _GrowableGrid(sum(column_widths), sum(row_heights))

        real_row = 0
        real_column = 0
        for logical_row in xrange(top_grid.height):
            for logical_column in xrange(top_grid.width):
                obj = top_grid[logical_row, logical_column]

                if obj in groups_grids:
                    # This is a group.  Copy the corresponding grid in
                    # place.
                    local_grid = groups_grids[obj]
                    for i in xrange(local_grid.height):
                        for j in xrange(local_grid.width):
                            grid[real_row + i, real_column + j] = local_grid[i, j]
                else:
                    # This is an object.  Just put it there.
                    grid[real_row, real_column] = obj

                real_column += column_widths[logical_column]
            real_column = 0
            real_row += row_heights[logical_row]

        return grid

    @staticmethod
    def _generic_layout(diagram, merged_morphisms):
        """
        Produces the generic layout for the supplied diagram.
        """
        all_objects = set(diagram.objects)
        if len(all_objects) == 1:
            # There only one object in the diagram, just put in on 1x1
            # grid.
            grid = _GrowableGrid(1, 1)
            grid[0, 0] = tuple(all_objects)[0]
            return grid

        skeleton = DiagramGrid._build_skeleton(merged_morphisms)

        grid = _GrowableGrid(2, 1)

        if len(skeleton) == 1:
            # This diagram contains only one morphism.  Draw it
            # horizontally.
            objects = sorted(all_objects, key=default_sort_key)
            grid[0, 0] = objects[0]
            grid[0, 1] = objects[1]

            return grid

        triangles = DiagramGrid._list_triangles(skeleton)
        triangles = DiagramGrid._drop_redundant_triangles(triangles, skeleton)
        triangle_sizes = DiagramGrid._compute_triangle_min_sizes(
            triangles, skeleton)

        triangles = sorted(triangles, key=lambda tri:
                           DiagramGrid._triangle_key(tri, triangle_sizes))

        # Place the first edge on the grid.
        root_edge = DiagramGrid._pick_root_edge(triangles[0], skeleton)
        grid[0, 0], grid[0, 1] = root_edge
        fringe = [((0, 0), (0, 1))]

        # Record which objects we now have on the grid.
        placed_objects = set(root_edge)

        while placed_objects != all_objects:
            welding = DiagramGrid._find_triangle_to_weld(triangles, fringe, grid)

            if welding:
                (triangle, welding_edge) = welding

                restart_required = DiagramGrid._weld_triangle(
                    triangle, welding_edge, fringe, grid, skeleton)
                if restart_required:
                    continue

                placed_objects.update(
                    DiagramGrid._triangle_objects(triangle))
            else:
                # No more weldings found.  Try to attach triangles by
                # vertices.
                new_obj = DiagramGrid._grow_pseudopod(
                    triangles, fringe, grid, skeleton, placed_objects)

                if not new_obj:
                    # No more triangles can be attached, not even by
                    # the edge.  We will set up a new diagram out of
                    # what has been left, laid it out independently,
                    # and then attach it to this one.

                    remaining_objects = all_objects - placed_objects

                    remaining_diagram = diagram.subdiagram_from_objects(
                        FiniteSet(remaining_objects))
                    remaining_grid = DiagramGrid(remaining_diagram)

                    # Now, let's glue ``remaining_grid`` to ``grid``.
                    final_width = grid.width + remaining_grid.width
                    final_height = max(grid.height, remaining_grid.height)
                    final_grid = _GrowableGrid(final_width, final_height)

                    for i in xrange(grid.width):
                        for j in xrange(grid.height):
                            final_grid[i, j] = grid[i, j]

                    start_j = grid.width
                    for i in xrange(remaining_grid.height):
                        for j in xrange(remaining_grid.width):
                            final_grid[i, start_j + j] = remaining_grid[i, j]

                    return final_grid

                placed_objects.add(new_obj)

            triangles = DiagramGrid._drop_irrelevant_triangles(
                triangles, placed_objects)

        return grid

    @staticmethod
    def _get_undirected_graph(objects, merged_morphisms):
        """
        Given the objects and the relevant morphisms of a diagram,
        returns the adjacency lists of the underlying undirected
        graph.
        """
        adjlists = {}
        for obj in objects:
            adjlists[obj] = []

        for morphism in merged_morphisms:
            adjlists[morphism.domain].append(morphism.codomain)
            adjlists[morphism.codomain].append(morphism.domain)

        # Assure that the objects in the adjacency list are always in
        # the same order.
        for obj in adjlists.keys():
            adjlists[obj].sort(key=default_sort_key)

        return adjlists

    @staticmethod
    def _sequential_layout(diagram, merged_morphisms):
        r"""
        Lays out the diagram in "sequential" layout.  This method
        will attempt to produce a result as close to a line as
        possible.  For linear diagrams, the result will actually be a
        line.
        """
        objects = diagram.objects
        sorted_objects = sorted(objects, key=default_sort_key)

        # Set up the adjacency lists of the underlying undirected
        # graph of ``merged_morphisms``.
        adjlists = DiagramGrid._get_undirected_graph(objects, merged_morphisms)

        # Find an object with the minimal degree.  This is going to be
        # the root.
        root = sorted_objects[0]
        mindegree = len(adjlists[root])
        for obj in sorted_objects:
            current_degree = len(adjlists[obj])
            if current_degree < mindegree:
                root = obj
                mindegree = current_degree

        grid = _GrowableGrid(1, 1)
        grid[0, 0] = root

        placed_objects = set([root])
        def place_objects(pt, placed_objects):
            """
            Does depth-first search in the underlying graph of the
            diagram and places the objects en route.
            """
            # We will start placing new objects from here.
            new_pt = (pt[0], pt[1] + 1)

            for adjacent_obj in adjlists[grid[pt]]:
                if adjacent_obj in placed_objects:
                    # This object has already been placed.
                    continue

                DiagramGrid._put_object(new_pt, adjacent_obj, grid, [])
                placed_objects.add(adjacent_obj)
                placed_objects.update(place_objects(new_pt, placed_objects))

                new_pt = (new_pt[0] + 1, new_pt[1])

            return placed_objects

        place_objects((0, 0), placed_objects)

        return grid

    @staticmethod
    def _drop_inessential_morphisms(merged_morphisms):
        r"""
        Removes those morphisms which should appear in the diagram,
        but which have no relevance to object layout.

        Currently this removes "loop" morphisms: the non-identity
        morphisms with the same domains and codomains.
        """
        morphisms = [m for m in merged_morphisms if m.domain != m.codomain]
        return morphisms

    @staticmethod
    def _get_connected_components(objects, merged_morphisms):
        """
        Given a container of morphisms, returns a list of connected
        components formed by these morphisms.  A connected component
        is represented by a diagram consisting of the corresponding
        morphisms.
        """
        component_index = {}
        for o in objects:
            component_index[o] = None

        # Get the underlying undirected graph of the diagram.
        adjlist = DiagramGrid._get_undirected_graph(objects, merged_morphisms)

        def traverse_component(object, current_index):
            """
            Does a depth-first search traversal of the component
            containing ``object``.
            """
            component_index[object] = current_index
            for o in adjlist[object]:
                if component_index[o] is None:
                    traverse_component(o, current_index)

        # Traverse all components.
        current_index = 0
        for o in adjlist:
            if component_index[o] is None:
                traverse_component(o, current_index)
                current_index += 1

        # List the objects of the components.
        component_objects = [[] for i in xrange(current_index)]
        for o, idx in component_index.items():
            component_objects[idx].append(o)

        # Finally, list the morphisms belonging to each component.
        #
        # Note: If some objects are isolated, they will not get any
        # morphisms at this stage, and since the layout algorithm
        # relies, we are essentially going to lose this object.
        # Therefore, check if there are isolated objects and, for each
        # of them, provide the trivial identity morphism.  It will get
        # discarded later, but the object will be there.

        component_morphisms = []
        for component in component_objects:
            current_morphisms = {}
            for m in merged_morphisms:
                if (m.domain in component) and (m.codomain in component):
                    current_morphisms[m] = merged_morphisms[m]

            if len(component) == 1:
                # Let's add an identity morphism, for the sake of
                # surely having morphisms in this component.
                current_morphisms[IdentityMorphism(component[0])] = FiniteSet()

            component_morphisms.append(Diagram(current_morphisms))

        return component_morphisms

    def __init__(self, diagram, groups=None, **hints):
        premises = DiagramGrid._simplify_morphisms(diagram.premises)
        conclusions = DiagramGrid._simplify_morphisms(diagram.conclusions)
        all_merged_morphisms = DiagramGrid._merge_premises_conclusions(
            premises, conclusions)
        merged_morphisms = DiagramGrid._drop_inessential_morphisms(
            all_merged_morphisms)

        # Store the merged morphisms for later use.
        self._morphisms = all_merged_morphisms

        components = DiagramGrid._get_connected_components(
            diagram.objects, all_merged_morphisms)

        if groups and (groups != diagram.objects):
            # Lay out the diagram according to the groups.
            self._grid = DiagramGrid._handle_groups(
                diagram, groups, merged_morphisms, hints)
        elif len(components) > 1:
            # Note that we check for connectedness _before_ checking
            # the layout hints because the layout strategies don't
            # know how to deal with disconnected diagrams.

            # The diagram is disconnected.  Lay out the components
            # independently.
            grids = []

            # Sort the components to eventually get the grids arranged
            # in a fixed, hash-independent order.
            components = sorted(components, key=default_sort_key)

            for component in components:
                grid = DiagramGrid(component, **hints)
                grids.append(grid)

            # Throw the grids together, in a line.
            total_width = sum(g.width for g in grids)
            total_height = max(g.height for g in grids)

            grid = _GrowableGrid(total_width, total_height)
            start_j = 0
            for g in grids:
                for i in xrange(g.height):
                    for j in xrange(g.width):
                        grid[i, start_j + j] = g[i, j]

                start_j += g.width

            self._grid = grid
        elif "layout" in hints:
            if hints["layout"] == "sequential":
                self._grid = DiagramGrid._sequential_layout(
                    diagram, merged_morphisms)
        else:
            self._grid = DiagramGrid._generic_layout(diagram, merged_morphisms)

        if hints.get("transpose"):
            # Transpose the resulting grid.
            grid = _GrowableGrid(self._grid.height, self._grid.width)
            for i in xrange(self._grid.height):
                for j in xrange(self._grid.width):
                    grid[j, i] = self._grid[i, j]
            self._grid = grid

    @property
    def width(self):
        """
        Returns the number of columns in this diagram layout.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> from sympy.categories import Diagram, DiagramGrid
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> diagram = Diagram([f, g])
        >>> grid = DiagramGrid(diagram)
        >>> grid.width
        2

        """
        return self._grid.width

    @property
    def height(self):
        """
        Returns the number of rows in this diagram layout.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> from sympy.categories import Diagram, DiagramGrid
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> diagram = Diagram([f, g])
        >>> grid = DiagramGrid(diagram)
        >>> grid.height
        2

        """
        return self._grid.height

    def __getitem__(self, (i, j)):
        """
        Returns the object placed in the row ``i`` and column ``j``.
        The indices are 0-based.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> from sympy.categories import Diagram, DiagramGrid
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> diagram = Diagram([f, g])
        >>> grid = DiagramGrid(diagram)
        >>> (grid[0, 0], grid[0, 1])
        (Object("A"), Object("B"))
        >>> (grid[1, 0], grid[1, 1])
        (None, Object("C"))

        """
        return self._grid[i, j]

    @property
    def morphisms(self):
        """
        Returns those morphisms (and their properties) which are
        sufficiently meaningful to be drawn.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> from sympy.categories import Diagram, DiagramGrid
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> diagram = Diagram([f, g])
        >>> grid = DiagramGrid(diagram)
        >>> grid.morphisms
        {NamedMorphism(Object("A"), Object("B"), "f"): EmptySet(),
        NamedMorphism(Object("B"), Object("C"), "g"): EmptySet()}

        """
        return self._morphisms

    def __str__(self):
        """
        Produces a string representation of this class.

        This method returns a string representation of the underlying
        list of lists of objects.

        Examples
        ========

        >>> from sympy.categories import Object, NamedMorphism
        >>> from sympy.categories import Diagram, DiagramGrid
        >>> A = Object("A")
        >>> B = Object("B")
        >>> C = Object("C")
        >>> f = NamedMorphism(A, B, "f")
        >>> g = NamedMorphism(B, C, "g")
        >>> diagram = Diagram([f, g])
        >>> grid = DiagramGrid(diagram)
        >>> print grid
        [[Object("A"), Object("B")],
        [None, Object("C")]]

        """
        return repr(self._grid._array)
