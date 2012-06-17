"""
This module contains the functionality to arrange the nodes of a
diagram on an abstract grid, and then to produce a graphical
representation of the grid.

The currently supported back-ends are Xy-pic [Xypic]

[Xypic] http://www.tug.org/applications/Xy-pic/
"""

from sympy.core import Basic

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
    def __new__(cls, diagram):
        pass
