from sympy.geometry.polygon import Polygon, Point, Circle
from sympy.core import sympify
from sympy.printing import sstr
from sympy.core import Symbol, S


class CrossSection(object):
    def __init__(self, shape, variable, *dimensions):
        self.shape = shape
        self.dimensions = dimensions
        self.variable = variable

        shape = shape.lower()
        
        if shape == "circular":
            self._circular(self.dimensions)
        elif shape == "triangular":
            self._triangular(self.dimensions)
        elif shape == "rectangular":
            self._rectangular(self.dimensions)
        elif shape == "tsection":
            self._tsection(self.dimensions)
        elif shape == "isection":
            self._isection(self.dimensions)
        else:
            raise ValueError("No such cross-sectional shape is supported")

    def __str__(self):
        str_1 = 'CrossSection({}, {})'.format(sstr(self._shape), sstr(self._dimensions))
        return str_1

    __repr__ = __str__

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, shape):
        self._shape = shape
    
    @property
    def dimensions(self):
        return self._dimensions

    @dimensions.setter
    def dimensions(self, dim):
        self._dimensions = dim
    
    @property
    def variable(self):
        return self._variable
    
    @variable.setter
    def variable(self, var):
        if isinstance(var, Symbol):
            self._variable = var
        else:
            raise TypeError("""The variable should be a Symbol object.""")

    @property
    def area(self):
       return self._area

    @property
    def centroid(self):
       return self._centroid
    
    @property
    def second_moment(self):
       return self._second_moment
            
    @property
    def section_modulus(self):
       return self._section_modulus
               
    @property
    def polar_modulus(self):
       return self._polar_modulus


    def _circular(self, radius):
        radius = radius[0]
        circle = Circle(Point(radius, radius), radius)

        self._area = circle.area
        self._centroid = circle.center
        self._second_moment = circle.second_moment_of_area()
        self._section_modulus = self._second_moment[0]/(self.variable - self.centroid.y)
        self._polar_modulus = self._second_moment[0] + self._second_moment[1]


    def _triangular(self, dimensions):
        width, height = dimensions
        triangle = Polygon((0,0), (width, 0), (width/2, height))

        self._centroid = traingle.centroid
        self._area = traingle.area
        self._second_moment = traingle.second_moment_of_area()
        self._centroid = traingle.centroid
        self._section_modulus = self._second_moment[0]/(self.variable - self._centroid.y)
        self._polar_modulus = self._second_moment[0] + self._second_moment[1]


    def _rectangular(self, dimensions):
        width, height = dimensions
        rectangle = Polygon((0, 0), (width, 0), (width, height), (0, height))

        self._area = rectangle.area
        self._centroid = rectangle.centroid
        self._second_moment = rectangle.second_moment_of_area()
        self._section_modulus = self._second_moment[0]/self.variable - self.centroid.y
        self._polar_modulus = self._second_moment[0] + self._second_moment[1]


    def _tsection(self, dimensions):
        flange, web = dimensions
        flange_width, flange_height = flange
        web_width, web_height = web
        total_height = flange_height + web_height
        tshape = Polygon((flange_width/2 - web_width/2, 0), (flange_width/2 + web_width/2, 0),
                         (flange_width/2 + web_width/2, web_height), (flange_width, web_height),
                         (flange_width, web_height + flange_height), (0, web_height + flange_height),
                         (0, web_height), (flange_width/2 - web_width/2, web_height))

        self._area = tshape.area
        self._centroid = tshape.centroid
        self._second_moment = tshape.second_moment_of_area()
        self._centroid = tshape.centroid
        self._section_modulus = self._second_moment[0]/(self.variable - self._centroid.y)
        self._polar_modulus = self._second_moment[0] + self._second_moment[1]


    def _isection(self, dimensions):
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

        self._area = ishape.area
        self._centroid = ishape.centroid
        self._second_moment = ishape.second_moment_of_area()
        self._section_modulus = self._second_moment[0]/(self.variable - self._centroid.y)
        self._polar_modulus = self._second_moment[0] + self._second_moment[1]
