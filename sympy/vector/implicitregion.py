from sympy import S
from sympy.core import Basic, Tuple, diff, expand, Eq
from sympy.core.symbol import _symbol
from sympy.solvers import solveset, nonlinsolve
from sympy.polys import total_degree
from sympy.geometry import Point


class ImplicitRegion(Basic):
    """
    Represents an implicit region in space.

    Example
    =======

    >>> from sympy import Eq
    >>> from sympy.abc import x, y, z, t
    >>> from sympy.vector import ImplicitRegion

    >>> ImplicitRegion((x, y), x**2 + y**2 - 4)
    ImplicitRegion((x, y), x**2 + y**2 - 4)
    >>> ImplicitRegion((x, y), Eq(y*x, 1))
    ImplicitRegion((x, y), x*y - 1)

    >>> parabola = ImplicitRegion((x, y), y**2 - 4*x)
    >>> parabola.degree
    2
    >>> parabola.equation
    -4*x + y**2
    >>> parabola.rational_parametrization(t)
    (4/t**2, 4/t)

    >>> r = ImplicitRegion((x, y, z), Eq(z, x**2 + y**2))
    >>> r.variables
    (x, y, z)
    >>> r.singular_points()
    FiniteSet((0, 0, 0))
    >>> r.regular_point()
    (-10, -10, 200)

    Parameters
    ==========

    variables: tuple to map variables in implicit equation to base scalars.

    equation: An expression or Eq denoting the implicit equation of the region.

    """
    def __new__(cls, variables, equation):
        if not isinstance(variables, Tuple):
            variables = Tuple(*variables)

        if isinstance(equation, Eq):
            equation = equation.lhs - equation.rhs

        return super().__new__(cls, variables, equation)

    @property
    def variables(self):
        return self.args[0]

    @property
    def equation(self):
        return self.args[1]

    @property
    def degree(self):
        return total_degree(self.equation)

    def regular_point(self):
        """
        Returns a point on the implicit region.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy.vector import ImplicitRegion
        >>> circle = ImplicitRegion((x, y), (x + 2)**2 + (y - 3)**2 - 16)
        >>> circle.regular_point()
        (-6, 3)

        """

        # TODO: implement the algorithm to find regular point on a Conic.
        # Now it iterates over a range and tries to find a point on the region.
        equation = self.equation

        if len(self.variables) == 1:
            return (list(solveset(equation, self.variables[0], domain=S.Reals))[0],)
        elif len(self.variables) == 2:
            x, y = self.variables

            for x_reg in range(-100, 100):
                if len(solveset(equation.subs(x, x_reg), self.variables[1], domain=S.Reals)) > 0:
                    return (x_reg, list(solveset(equation.subs(x, x_reg), domain=S.Reals))[0])
        elif len(self.variables) == 3:
            x, y, z = self.variables

            for x_reg in range(-10, 10):
                for y_reg in range(-10, 10):
                    if len(solveset(equation.subs({x: x_reg, y: y_reg}), self.variables[2], domain=S.Reals)) > 0:
                        return (x_reg, y_reg, list(solveset(equation.subs({x: x_reg, y: y_reg})))[0])

        if len(self.singular_points()) != 0:
            return list[self.singular_points()][0]

        raise NotImplementedError()

    def singular_points(self):
        """
        Returns a set of singular points of the region.

        The singular points are those points on the region
        where all partial derivatives vanish.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy.vector import ImplicitRegion
        >>> I = ImplicitRegion((x, y), (y-1)**2 -x**3 + 2*x**2 -x)
        >>> I.singular_points()
        FiniteSet((1, 1))

        """
        eq_list = [self.equation]
        for var in self.variables:
            eq_list += [diff(self.equation, var)]

        return nonlinsolve(eq_list, list(self.variables))

    def multiplicity(self, point):
        """
        Returns the multiplicity of a singular point on the region.

        A singular point (x,y) of region is said to be of multiplicity m
        if all the partial derivatives off to order m - 1 vanish there.

        Examples
        ========

        >>> from sympy.abc import x, y, z
        >>> from sympy.vector import ImplicitRegion
        >>> I = ImplicitRegion((x, y, z), x**2 + y**3 - z**4)
        >>> I.singular_points()
        FiniteSet((0, 0, 0))
        >>> I.multiplicity((0, 0, 0))
        2

        """
        if isinstance(point, Point):
            point = point.args

        modified_eq = self.equation

        for i, var in enumerate(self.variables):
            modified_eq = modified_eq.subs(var, var + point[i])
        modified_eq = expand(modified_eq)

        if len(modified_eq.args) != 0:
            terms = modified_eq.args
            m = min([total_degree(term) for term in terms])
        else:
            terms = modified_eq
            m = total_degree(terms)

        return m

    def rational_parametrization(self, parameter1='t', parameter2='s', reg_point=None):
        """
        Returns the rational parametrization of implict region.

        Examples
        ========

        >>> from sympy import Eq
        >>> from sympy.abc import x, y, z, s, t
        >>> from sympy.vector import ImplicitRegion

        >>> parabola = ImplicitRegion((x, y), y**2 - 4*x)
        >>> parabola.rational_parametrization()
        (4/t**2, 4/t)

        >>> circle = ImplicitRegion((x, y), Eq(x**2 + y**2, 4))
        >>> circle.rational_parametrization()
        (-2 + 4/(t**2 + 1), 4*t/(t**2 + 1))

        >>> I = ImplicitRegion((x, y), x**3 + x**2 - y**2)
        >>> I.rational_parametrization()
        (t**2 - 1, t*(t**2 - 1))

        >>> cubic_curve = ImplicitRegion((x, y), x**3 + x**2 - y**2)
        >>> cubic_curve.rational_parametrization(t)
        (t**2 - 1, t*(t**2 - 1))

        >>> sphere = ImplicitRegion((x, y, z), x**2 + y**2 + z**2 - 4)
        >>> sphere.rational_parametrization(parameter1=t, parameter2=s)
        (-2 + 4/(s**2 + t**2 + 1), 4*s/(s**2 + t**2 + 1), 4*t/(s**2 + t**2 + 1))

        For some conics, regular_points() is unable to find a point on curve.
        To calulcate the parametric representation in such cases, user need
        to determine a point on the region and pass it using reg_point.

        >>> c = ImplicitRegion((x, y), (x  - 1/2)**2 + (y)**2 - (1/4)**2)
        >>> c.rational_parametrization(reg_point=(3/4, 0))
        (0.75 - 0.5/(t**2 + 1), -0.5*t/(t**2 + 1))

        Refrences
        =========

        - Christoph M. Hoffmann, Conversion Methods between Parametric and
        Implicit Curves and Surfaces.

        """
        equation = self.equation
        degree = self.degree

        if degree == 1:
            if len(self.variables) == 1:
                return (equation,)
            elif len(self.variables) == 2:
                x, y = self.variables
                y_par = list(solveset(equation, y))[0]
                return x, y_par
            else:
                raise NotImplementedError()

        point = ()

        # Finding the (n - 1) fold point of the monoid of degree
        if degree == 2:
            # For degree 2 curves, either a regular point or a singular point can be used.
            if reg_point is not None:
                # Using point provided by the user as regular point
                if isinstance(point, Point):
                    point = point.args
                point = reg_point
            else:
                if len(self.singular_points()) != 0:
                    point = list(self.singular_points())[0]
                else:
                    point = self.regular_point()

        if len(self.singular_points()) != 0:
            singular_points = self.singular_points()
            for spoint in singular_points:
                if self.multiplicity(spoint) == degree - 1:
                    point = spoint
                    break

        if len(point) == 0:
            # The region in not a monoid
            raise NotImplementedError()

        modified_eq = equation

        # Shifting the region such that fold point moves to origin
        for i, var in enumerate(self.variables):
            modified_eq = modified_eq.subs(var, var + point[i])
        modified_eq = expand(modified_eq)

        hn = hn_1 = 0
        for term in modified_eq.args:
            if total_degree(term) == degree:
                hn += term
            else:
                hn_1 += term

        hn_1 = -1*hn_1

        if len(self.variables) == 2:

            if parameter1 == 's':
                # To avoid name conflict between parameters
                s = _symbol('s_', real=True)
            else:
                s = _symbol('s', real=True)
            t = _symbol(parameter1, real=True)

            hn = hn.subs({self.variables[0]: s, self.variables[1]: t})
            hn_1 = hn_1.subs({self.variables[0]: s, self.variables[1]: t})

            x_par = (s*(hn_1/hn)).subs(s, 1) + point[0]
            y_par = (t*(hn_1/hn)).subs(s, 1) + point[1]

            return x_par, y_par

        elif len(self.variables) == 3:

            if parameter1 == 'r' or parameter2 == 'r':
                # To avoid name conflict between parameters
                r = _symbol('r_', real=True)
            else:
                r = _symbol('r', real=True)
            s = _symbol(parameter2, real=True)
            t = _symbol(parameter1, real=True)

            hn = hn.subs({self.variables[0]: r, self.variables[1]: s, self.variables[2]: t})
            hn_1 = hn_1.subs({self.variables[0]: r, self.variables[1]: s, self.variables[2]: t})

            x_par = (r*(hn_1/hn)).subs(r, 1) + point[0]
            y_par = (s*(hn_1/hn)).subs(r, 1) + point[1]
            z_par = (t*(hn_1/hn)).subs(r, 1) + point[2]

            return x_par, y_par, z_par

        raise NotImplementedError()
