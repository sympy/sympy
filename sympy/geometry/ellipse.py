from sympy.core.basic import S, C, Basic, sympify
from sympy.simplify import simplify, trigsimp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.geometry.exceptions import GeometryError
from entity import GeometryEntity
from point import Point
from line import LinearEntity, Line

class Ellipse(GeometryEntity):
    """
    An ellipse in space. Constructed from a center and two radii, the
    first being the horizontal radius (along the x-axis) and the second
    being the vertical radius (along the y-axis).

    Notes:
    ======
        - Rotation is currently not supported since an ellipse is defined
          on horizontal/vertical radii

    Example:
    ========
        >>> from sympy.geometry import Ellipse, Point
        >>> e = Ellipse(Point(0, 0), 5, 1)
        >>> e.hradius, e.vradius
        (5, 1)

    Plotting example
    ----------------
    In [1]: c1 = Circle(Point(0,0), 1)

    In [2]: Plot(c1)
    Out[2]: [0]: cos(t), sin(t), 'mode=parametric'

    In [3]: p = Plot()

    In [4]: p[0] = c1

    In [5]: radius = Segment(c1.center, c1.random_point())

    In [6]: p[1] = radius

    In [7]: p
    Out[7]:
    [0]: cos(t), sin(t), 'mode=parametric'
    [1]: t*cos(1.546086215036205357975518382),
    t*sin(1.546086215036205357975518382), 'mode=parametric'
    """

    def __new__(cls, center, hradius, vradius, **kwargs):
        hradius = sympify(hradius)
        vradius = sympify(vradius)
        if not isinstance(center, Point):
            raise TypeError("center must be be a Point")

        if hradius == vradius:
            return Circle(center, hradius, **kwargs)
        return GeometryEntity.__new__(cls, center, hradius, vradius, **kwargs)

    @property
    def center(self):
        """The center of the ellipse."""
        return self.__getitem__(0)

    @property
    def hradius(self):
        """The horizontal radius of the ellipse."""
        return self.__getitem__(1)

    @property
    def vradius(self):
        """The vertical radius of the ellipse."""
        return self.__getitem__(2)

    @property
    def area(self):
        """The area of the ellipse."""
        return simplify(S.Pi * self.hradius * self.vradius)

    @property
    def circumference(self):
        """The circumference of the ellipse."""
        # TODO It's fairly complicated, but we could use Ramanujan's
        #      approximation.
        raise NotImplementedError

    @property
    def foci(self):
        """The foci of the ellipse, if the radii are numerical."""
        c = self.center
        if self.hradius == self.vradius:
            return c

        hr, vr = self.hradius, self.vradius
        if hr.atoms(C.Symbol) or vr.atoms(C.Symbol):
            raise ValueError("foci can only be determined on non-symbolic radii")

        v = sqrt(abs(vr**2 - hr**2))
        if hr < vr:
            return (c+Point(0, -v), c+Point(0, v))
        else:
            return (c+Point(-v, 0), c+Point(v, 0))

    def tangent_line(self, p):
        """
        If p is on the ellipse, returns the tangent line through point p.
        Otherwise, returns the tangent line(s) from p to the ellipse, or
        None if no tangent line is possible (e.g., p inside ellipse).

        Example:

        In [1]: e = Ellipse(Point(0,0), 3, 2)

        In [2]: t = e.tangent_line(e.random_point())

        In [3]: p = Plot()

        In [4]: p[0] = e

        In [5]: p[1] = t

        The above will plot an ellipse together with a tangent line.
        """
        if p in self:
            rise = (self.vradius ** 2)*(self.center[0] - p[0])
            run = (self.hradius ** 2)*(p[1] - self.center[1])
            p2 = Point(simplify(p[0] + run),
                       simplify(p[1] + rise))
            return Line(p, p2)
        else:
            # TODO If p is not on the ellipse, attempt to create the
            #      tangent(s) from point p to the ellipse..?
            raise NotImplementedError("Cannot find tangent lines when p is not on the ellipse")

    def is_tangent(self, o):
        """Returns True if o is tangent to the ellipse, False otherwise."""
        inter = None
        if isinstance(o, Ellipse):
            inter = self.intersection(o)
            return (inter is not None and isinstance(inter[0], Point) and len(inter) == 1)
        elif isinstance(o, LinearEntity):
            inter = self._do_line_intersection(o)
            if (inter is not None) and len(inter) == 1:
                return (inter[0] in o)
            else:
                return False
        elif isinstance(o, Polygon):
            c = 0
            for seg in o.sides:
                inter = self._do_line_intersection(seg)
                c += len([True for point in inter if point in seg])
            return (c == 1)
        else:
            raise NotImplementedError("Unknown argument type")

    def arbitrary_point(self, parameter_name='t'):
        """Returns a symbolic point that is on the ellipse."""
        t = C.Symbol(parameter_name, real=True)
        return Point(
                self.center[0] + self.hradius*C.cos(t),
                self.center[1] + self.vradius*C.sin(t))

    def plot_interval(self, parameter_name='t'):
        """Returns a typical plot interval used by the plotting module."""
        t = C.Symbol(parameter_name, real=True)
        return [t, -S.Pi, S.Pi]

    def random_point(self):
        """Returns a random point on the ellipse."""
        from random import random
        from sys import maxint
        t = C.Symbol('t', real=True)
        p = self.arbitrary_point('t')
        # get a random value in [-pi, pi)
        subs_val = float(S.Pi)*(2*random() - 1)
        return Point(p[0].subs(t, subs_val), p[1].subs(t, subs_val))

    def equation(self, x='x', y='y'):
        """
        Returns the equation of the ellipse.

        Optional parameters x and y can be used to specify symbols, or the
        names of the symbols used in the equation.
        """
        if isinstance(x, basestring):   x = C.Symbol(x, real=True)
        if isinstance(y, basestring):   y = C.Symbol(y, real=True)
        t1 = ((x - self.center[0]) / self.hradius)**2
        t2 = ((y - self.center[1]) / self.vradius)**2
        return t1 + t2 - 1

    def _do_line_intersection(self, o):
        """
        Find the intersection of a LinearEntity and the ellipse. Makes no
        regards to what the LinearEntity is because it assumes a Line. To
        ensure correct intersection results one must invoke intersection()
        to remove bad results.
        """
        def dot(p1, p2):
            sum = 0
            for ind in xrange(0, len(p1)):
                sum += p1[ind] * p2[ind]
            return simplify(sum)

        hr_sq = self.hradius ** 2
        vr_sq = self.vradius ** 2
        lp = o.points

        ldir = lp[1] - lp[0]
        diff = lp[0] - self.center
        mdir = (ldir[0] / hr_sq, ldir[1] / vr_sq)
        mdiff = (diff[0] / hr_sq, diff[1] / vr_sq)

        a = dot(ldir, mdir)
        b = dot(ldir, mdiff)
        c = dot(diff, mdiff) - 1
        det = simplify(b*b - a*c);

        result = []
        if det == 0:
            t = -b / a
            result.append( lp[0] + (lp[1] - lp[0]) * t )
        else:
            is_good = True
            try:
                is_good = (det > 0)
            except NotImplementedError: #symbolic, allow
                is_good = True

            if is_good:
                root = sqrt(det)
                t_a = (-b - root) / a
                t_b = (-b + root) / a
                result.append( lp[0] + (lp[1] - lp[0]) * t_a )
                result.append( lp[0] + (lp[1] - lp[0]) * t_b )
        return result

    def intersection(self, o):
        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return []
        elif isinstance(o, LinearEntity):
            # LinearEntity may be a ray/segment, so check the points
            # of intersection for coincidence first
            result = self._do_line_intersection(o)
            if result is not None:
                for ind in xrange(len(result)-1, -1, -1):
                    if result[ind] not in o:
                        del result[ind]
            return result
        elif isinstance(o, Ellipse):
            if o == self:
                return self
            else:
                # TODO This is a bit more complicated
                pass

        raise NotImplementedError()

    def __eq__(self, o):
        return ((self.center == o.center) \
                    and (self.hradius == o.hradius) \
                    and (self.vradius == o.vradius))

    def __contains__(self, o):
        if isinstance(o, Point):
            x = C.Symbol('x', real=True, dummy=True)
            y = C.Symbol('y', real=True, dummy=True)

            res = self.equation(x, y).subs({x: o[0], y: o[1]})
            res = trigsimp(simplify(res))
            return res == 0
        elif isinstance(o, Ellipse):
            return (self == o)
        return False


class Circle(Ellipse):
    """
    A circle in space. Consturcted simply from a center and a radius, or
    from three non-collinear points.

    Example:
    ========
        >>> from sympy.geometry import Point, Circle
        >>> c1 = Circle(Point(0, 0), 5)
        >>> c1.hradius, c1.vradius, c1.radius
        (5, 5, 5)
        >>> c2 = Circle(Point(0, 0), Point(1, 1), Point(1, 0))
        >>> c2.hradius, c2.vradius, c2.radius, c2.center
        (2**(1/2)/2, 2**(1/2)/2, 2**(1/2)/2, Point(1/2, 1/2))

    """
    def __new__(cls, *args, **kwargs):
        c, r = None, None
        if len(args) == 3 and isinstance(args[0], Point):
            from polygon import Triangle
            t = Triangle(args[0], args[1], args[2])
            if t.area == 0:
                raise GeometryError("Cannot construct a circle from three collinear points")
            c = t.circumcenter
            r = t.circumradius
        elif len(args) == 2:
            # Assume (center, radius) pair
            c = args[0]
            r = sympify(args[1])

        if not (c is None or r is None):
            return GeometryEntity.__new__(cls, c, r, **kwargs)

        raise GeometryError("Circle.__new__ received unknown arguments")

    @property
    def hradius(self):
        """The horizontal radius of the ellipse."""
        return self.__getitem__(1)

    @property
    def vradius(self):
        """The vertical radius of the ellipse."""
        return self.__getitem__(1)

    @property
    def radius(self):
        """The radius of the circle."""
        return self.__getitem__(1)

    @property
    def circumference(self):
        """The circumference of the circle."""
        return 2 * S.Pi * self.radius

    def equation(self, x='x', y='y'):
        """
        Returns the equation of the circle.

        Optional parameters x and y can be used to specify symbols, or the
        names of the symbols used in the equation.
        """
        if isinstance(x, basestring):   x = C.Symbol(x, real=True)
        if isinstance(y, basestring):   y = C.Symbol(y, real=True)
        t1 = (x - self.center[0])**2
        t2 = (y - self.center[1])**2
        return t1 + t2 - self.hradius**2

    def intersection(self, o):
        if isinstance(o, Circle):
            dx,dy = o.center - self.center
            d = sqrt( simplify(dy**2 + dx**2) )
            a = simplify((self.radius**2 - o.radius**2 + d**2) / (2*d))

            x2 = self.center[0] + (dx * a/d)
            y2 = self.center[1] + (dy * a/d)

            h = sqrt( simplify(self.radius**2 - a**2) )
            rx = -dy * (h/d)
            ry =  dx * (h/d)

            xi_1 = simplify(x2 + rx)
            xi_2 = simplify(x2 - rx)
            yi_1 = simplify(y2 + ry)
            yi_2 = simplify(y2 - ry)

            ret = [Point(xi_1, yi_1)]
            if xi_1 != xi_2 or yi_1 != yi_2:
                ret.append(Point(xi_2, yi_2))
            return ret
        elif isinstance(o, Ellipse):
            a, b, r = o.hradius, o.vradius, self.radius
            x = a*sqrt(simplify((r**2 - b**2)/(a**2 - b**2)))
            y = b*sqrt(simplify((a**2 - r**2)/(a**2 - b**2)))
            return list(set([Point(x,y), Point(x,-y), Point(-x,y), Point(-x,-y)]))

        return Ellipse.intersection(self, o)
