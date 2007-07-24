from sympy import *
from sympy.modules.simplify import simplify#, trigsimp
from entity import GeometryEntity
from point import Point
from line import LinearEntity, Line

class Ellipse(GeometryEntity):
    """An ellipse in space."""

    def __init__(self, center, h_radius, v_radius, **kwargs):
        GeometryEntity.__init__(self, [center, h_radius, v_radius], **kwargs)
        if not isinstance(center, Point):
            raise TypeError("center must be be a Point")
        self._c = center
        self._hr = Basic.sympify(h_radius)
        self._vr = Basic.sympify(v_radius)

    @property
    def horizontal_radius(self):
        """Returns the horizontal radius of the ellipse."""
        return self._hr

    @property
    def vertical_radius(self):
        """Returns the vertical radius of the ellipse."""
        return self._vr

    @property
    def center(self):
        """Returns the center of the ellipse."""
        return self._c

    @property
    def area(self):
        """Returns the area of the ellipse."""
        return simplify(pi * self._hr * self._vr)

    @property
    def circumference(self):
        """Returns the circumference of the ellipse."""
        # TODO It's fairly complicated, but we could use Ramanujan's approximation.
        raise NotImplementedError

    @property
    def foci(self):
        """Returns the foci of the ellipse, if the radii are numerical."""
        c = self._c
        if self._hr == self._vr:
            return c

        if self._hr.atoms(type=Symbol) or self._vr.atoms(type=Symbol):
            raise Exception("foci can only be determined on numerical radii")

        elif self._hr < self._vr:
            v = sqrt(self._vr**2 - self._hr**2)
            return (c+Point(0, -v), c+Point(0, v))
        else:
            v = sqrt(self._hr**2 - self._vr**2)
            return (c+Point(-v, 0), c+Point(v, 0))

    def tangent_line(self, p):
        """
        If p is on the ellipse, returns the tangent line through point p.
        Otherwise, returns the tangent line(s) from p to the ellipse, or
        None if no tangent line is possible (e.g., p inside ellipse).
        """
        if p in self:
            rise = (self._hr ** Rational(2))*(self._c[0] - p[0])
            run = (self._vr ** Rational(2))*(p[1] - self._c[1])
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
                inter = self._do_line_intersection(o)
                if (inter is not None and inter[0] in o):
                    c += len(inter)
            return (c == 1)
        else:
            raise NotImplementedError("Unknown argument type")

    def arbitrary_point(self, parameter_name='t'):
        """Returns a symbolic point that is on the ellipse."""
        t = Symbol(parameter_name)
        return Point(self._c[0] + self._hr*cos(t), self._c[1] + self._vr*sin(t))

    def random_point(self):
        """Returns a random point on the ellipse."""
        from random import randint
        from sys import maxint
        t = Symbol('t')
        p = self.arbitrary_point('t')
        subs_val = randint(-maxint-1, maxint)
        return Point(p[0].subs(t, subs_val), p[1].subs(t, subs_val))

    def equation(self, xaxis_name='x', yaxis_name='y'):
        """Returns the equation of the ellipse."""
        x = Symbol(xaxis_name)
        y = Symbol(yaxis_name)
        t1 = ((x - self._c[0]) / self._hr)**Rational(2)
        t2 = ((y - self._c[1]) / self._vr)**Rational(2)
        return t1 + t2 - Rational(1)

    def _do_line_intersection(self, o):
        """
        Find the intersection of a LinearEntity and the ellipse. Makes no regards
        to what the LinearEntity is because it assumes a Line. To ensure correct
        intersection results one must invoke intersection() to remove bad results.
        """
        def dot(p1, p2):
            sum = Rational(0)
            for ind in xrange(0, len(p1)):
                sum += p1[ind] * p2[ind]
            return simplify(sum)

        hr_sq = self._hr ** Rational(2)
        vr_sq = self._vr ** Rational(2)
        lp = o.points

        ldir = lp[1] - lp[0]
        diff = lp[0] - self._c
        mdir = (ldir[0] / hr_sq, ldir[1] / vr_sq)
        mdiff = (diff[0] / hr_sq, diff[1] / vr_sq)

        a = dot(ldir, mdir)
        b = dot(ldir, mdiff)
        c = dot(diff, mdiff) - Rational(1)
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

    def _intersection(self, o):
        """
        Returns the intersection of the ellipse and another entity, or None if
        there is no intersection.
        """
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
        return ((self._c == o._c) and
                (self._hr == o._hr) and (self._vr == o._vr))

    def __contains__(self, o):
        if isinstance(o, Point):
            x = Symbol('x')
            y = Symbol('y')
            res = self.equation('x', 'y').subs_dict({x: o[0], y: o[1]})
            #res = trigsimp(simplify(res)) 
            res = simplify(res)
            return bool(res == 0)
        elif isinstance(o, Ellipse):
            return (self == o)
        else:
            return False

    def __str__(self):
        return "Ellipse(%s, %s, %s)" % (repr(self._c), repr(self._hr), repr(self._vr))


class Circle(Ellipse):
    """A circle in space."""

    def __init__(self, *args, **kwargs):
        if len(args) == 2:
            # Assume (center, radius) pair
            Ellipse.__init__(self, args[0], args[1], args[1], **kwargs)
        elif len(args) == 3 and isinstance(args[0], Point):
            from polygon import Triangle
            t = Triangle(args[0], args[1], args[2])
            if t.area == 0:
                raise Exception("Given points are not concyclic")
            c = t.circumcenter
            r = t.circumradius
            Ellipse.__init__(self, c, r, r, **kwargs)
        else:
            raise Exception("Unknown set of arguments")

    @property
    def radius(self):
        return self._hr

    @property
    def circumference(self):
        return Rational(2) * pi * self.radius

    def equation(self, xaxis_name='x', yaxis_name='y'):
        """Returns the equation of the circle."""
        x = Symbol(xaxis_name)
        y = Symbol(yaxis_name)
        t1 = (x - self._c[0])**Rational(2)
        t2 = (y - self._c[1])**Rational(2)
        return t1 + t2 - self._hr**Rational(2)

    def _intersection(self, o):
        if isinstance(o, Circle):
            dx,dy = o._c - self._c
            d = sqrt( simplify(dy**2 + dx**2) )
            a = simplify((self._hr**2 - o._hr**2 + d**2) / (Rational(2) * d))

            x2 = self._c[0] + (dx * a/d)
            y2 = self._c[1] + (dy * a/d)

            h = sqrt( simplify(self._hr**2 - a**2) )
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
            a,b,r = o.horizontal_radius,o.vertical_radius,self._hr
            x = a*sqrt(simplify((r**2 - b**2)/(a**2 - b**2)))
            y = b*sqrt(simplify((a**2 - r**2)/(a**2 - b**2)))
            return list(set([Point(x,y), Point(x,-y), Point(-x,y), Point(-x,-y)]))

        return Ellipse._intersection(self, o)

    def __str__(self):
        return "Circle(%s, %s)" % (repr(self._c), repr(self._hr))