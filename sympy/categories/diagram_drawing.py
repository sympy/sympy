"""
This module contains the functionality to arrange the nodes of a
diagram on an abstract grid, and then to produce a graphical
representation of the grid.

The currently supported back-ends are Xy-pic [Xypic]

[Xypic] http://www.tug.org/applications/Xy-pic/
"""

from sympy.core import Basic
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

    def __new__(cls, diagram):
        premises = DiagramGrid._simplify_morphisms(diagram.premises)
        conclusions = DiagramGrid._simplify_morphisms(diagram.conclusions)
        merged_morphisms = DiagramGrid._merge_premises_conclusions(
            premises, conclusions)

        skeleton = DiagramGrid._build_skeleton(merged_morphisms)
