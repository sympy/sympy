from sympy.printing import sstr
from sympy.geometry import Ellipse, Polygon


class CrossSection(object):
    """
    A class to define the cross-section of any structural element,
    especially beam.
    The CrossSection class supports any arbitrary geometric shape to be defined
    as the cross-section of the structure.
    It gives functionality to provide different properties of the cross-section.

    Attributes
    ==========

    geometry: an object of the geometry module

    Examples
    ========

    >>> from sympy import Polygon
    >>> from sympy.physics.continuum_mechanics.cross_section import CrossSection
    >>> p = CrossSection(Polygon((0, 0), (9, 0), (9, 13), (0, 13)))
    >>> p.shape
    Polygon(Point2D(0, 0), Point2D(9, 0), Point2D(9, 13), Point2D(0, 13))
    >>> p.centroid()
    Point2D(9/2, 13/2)
    >>> p.area()
    117
    >>> p.second_moment()
    (6591/4, 3159/4, 0)
    >>> p.polar_modulus()
    4875/2
    >>> p.section_modulus()
    (507/2, 351/2)

    >>> from sympy import Circle
    >>> c = CrossSection(Circle((15, 15), 15))
    >>> c.centroid()
    Point2D(15, 15)
    >>> c.area()
    225*pi
    >>> c.second_moment()
    (50625*pi/4, 50625*pi/4, 0)
    >>> c.polar_modulus()
    50625*pi/2
    >>> c.section_modulus()
    (3375*pi/4, 3375*pi/4)
    """
    def __init__(self, shape):
        self.shape = shape

    def __str__(self):
        str_1 = 'CrossSection({})'.format(sstr(self._shape))
        return str_1

    __repr__ = __str__

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, shape):
        self._shape = shape

    def area(self):
        """Returns the area of the cross-section geometry"""
        return self.shape.area

    def centroid(self):
        """Returns a point with the centroid coordinates"""
        if isinstance(self.shape, Ellipse):
            return self.shape.center
        return self.shape.centroid

    def second_moment(self):
        """Returns a tuple with 3 elements Ixx, Iyy, Ixy
           about the centroid of the cross-section geometry
        """
        return self.shape.second_moment_of_area()

    def polar_modulus(self):
        """Returns the polar modulus of the cross-section geometry"""
        return self.second_moment()[0] + self.second_moment()[1]

    def section_modulus(self):
        """Returns a tuple with the section modulus of the
           cross-section geometry about the centroid.
        """
        x_min, y_min, x_max, y_max = self.shape.bounds

        y_lower = self.centroid()[1] - y_min
        y_upper = y_max - self.centroid()[1]
        x_left = self.centroid()[0] - x_min
        x_right = x_max - self.centroid()[0]

        y = y_upper if y_upper > y_lower else y_lower
        x = x_right if x_right > x_left else x_left

        S_x, S_y = self.second_moment()[0]/y, self.second_moment()[1]/x
        return S_x, S_y


def tsection(*dimensions):
    """
    Returns a Polygon object of t shape with the required dimensions

    Arguments
    =========

    flange_dimensions: tuple
                       (flange_width, flange_geight)
    web_dimensions: tuple
                    (web_width, web_height)

    Examples
    ========

    >>> from sympy.physics.continuum_mechanics.cross_section import tsection, CrossSection
    >>> t = tsection((5, 3), (3, 8))
    >>> t
    Polygon(Point2D(1, 0), Point2D(4, 0), Point2D(4, 8), Point2D(5, 8),
    Point2D(5, 11), Point2D(0, 11), Point2D(0, 8), Point2D(1, 8))
    >>> s = CrossSection(t)
    >>> s.centroid()
    Point2D(5/2, 159/26)
    >>> s.area()
    39
    >>> s.second_moment()
    (21761/52, 197/4, 0)
    >>> s.polar_modulus()
    12161/26
    >>> s.section_modulus()
    (21761/318, 197/10)
    """
    flange, web = dimensions
    flange_width, flange_height = flange
    web_width, web_height = web
    total_height = flange_height + web_height
    tshape = Polygon((flange_width/2 - web_width/2, 0), (flange_width/2 + web_width/2, 0),
                     (flange_width/2 + web_width/2, web_height), (flange_width, web_height),
                     (flange_width, web_height + flange_height), (0, web_height + flange_height),
                     (0, web_height), (flange_width/2 - web_width/2, web_height))

    return tshape


def isection(*dimensions):
    """
    Returns a Polygon object representing an I-shape with the required
    dimensions

    Arguments
    ==========

    upper_flange_dimensions: tuple
                             (upper_flange_width, upper_flange_height)
    web dimensions: tuple
                    (web_width, web_height))
    lower_flange_dimensions: tuple
                             (lower_flange_width, lower_flange_height)

    Examples
    ========

    >>> from sympy.physics.continuum_mechanics.cross_section import isection, CrossSection
    >>> t = isection((10, 3), (3, 8), (20, 4))
    >>> t
    Polygon(Point2D(0, 0), Point2D(20, 0), Point2D(20, 4), Point2D(23/2, 4),
    Point2D(23/2, 12), Point2D(15, 12), Point2D(15, 15), Point2D(5, 15),
    Point2D(5, 12), Point2D(17/2, 12), Point2D(17/2, 4), Point2D(0, 4))

    >>> s = CrossSection(t)
    >>> s.centroid()
    Point2D(10, 757/134)
    >>> s.area()
    134
    >>> s.second_moment()
    (1328281/402, 8804/3, 0)
    >>> s.polar_modulus()
    2508017/402
    >>> s.section_modulus()
    (1328281/3759, 4402/15)
    """
    upper_flange, web, lower_flange = dimensions
    upper_flange_width, upper_flange_height = upper_flange
    lower_flange_width, lower_flange_height = lower_flange
    web_width, web_height = web

    d = lower_flange_width - upper_flange_width
    total_height = lower_flange_height + web_height + upper_flange_height

    ishape = Polygon((0,0), (lower_flange_width, 0), (lower_flange_width, lower_flange_height),
                     (lower_flange_width/2 + web_width/2, lower_flange_height),
                     (lower_flange_width/2 + web_width/2, lower_flange_height + web_height),
                     (lower_flange_width/2 + upper_flange_width/2, lower_flange_height + web_height),
                     (lower_flange_width/2 + upper_flange_width/2, total_height), (d/2, total_height),
                     (d/2, lower_flange_height + web_height),
                     (lower_flange_width/2 - web_width/2, lower_flange_height + web_height),
                     (lower_flange_width/2 - web_width/2, lower_flange_height), (0, lower_flange_height))

    # in case the upper flange enters into second quadrant
    # this is done inorder to get the correct values of centroid
    if upper_flange_width > lower_flange_width:
        ishape = ishape.translate((upper_flange_width - lower_flange_width)/2, 0)

    return ishape
