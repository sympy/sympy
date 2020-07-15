from sympy import S
from sympy.core import Basic, Tuple, Eq, diff, expand
from sympy.core.symbol import _symbol
from sympy.solvers import solveset, nonlinsolve
from sympy.polys.polytools import total_degree


class ImplicitRegion(Basic):
    """
    Represents an implicit region in space. 

    Example
    =======

    >>> from sympy.abc import x, y
    >>> from sympy.vector import ImplicitRegion, parametric_region_list, vector_integrate

    >>> circle1 = ImplicitRegion((x, y), Eq(x**2 + y**2 -1))
    >>> parametric_region_list(circle)
    [ParametricRegion((sin(2*theta), -cos(2*theta)), (theta, 0, pi))]

    >>> circle2 = ImplicitRegion((x, y), (x - 3)**2 + (y + 2)**2 - 16)
    >>> parametric_region_list(circle2)
    [ParametricRegion((2*(-sqrt(7)*tan(theta) + 3)*cos(theta)**2, 3*sin(2*theta) + sqrt(7)*cos(2*theta) - 2), (theta, 0, pi))]
    >>> vector_integrate(1, circle2)
    8*pi

    >>> ellipse = ImplicitRegion(x**2/4 + y**2/16 - 1)
    >>> parametric_region_list(ellipse)
    [ParametricRegion((8*tan(theta)/(tan(theta)**2 + 4), -4 + 8*tan(theta)**2/(tan(theta)**2 + 4)), (theta, 0, pi))]
    >>> vector_integrate(2, ellipse)
    ### It gets stuck

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
        Returns a smooth point on the implicit region.
        
        >>> circle = ImplicitRegion((x, y), (x + 2)**2 + (y - 3)**2 - 16)
        >>> circle.regular_point
        (-6, 3)
         
        """
        equation = self.equation

        if len(self.variables) == 1:
            return list(solveset(equation, self.variables[0], domain=S.Reals))[0]
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
        
        raise NotImplementedError()
        
    def singular_points(self):
        """
        Returns a set of singular points of the region.
        
        The singular points are those points on the region
        where all partial derivatives vanish.
        
        Examples
        ========
        
        >>> I = ImplicitRegion((x, y), (y-1)**2 -x**3 + 2*x**2 -x)
        >>> I.singular_points()
        FiniteSet((1, 1))
        
        """
        eq_list = [self.equation]
        for var in self.variables:
            eq_list += [diff(self.equation, var)]
        
        print(eq_list)    
        return nonlinsolve(eq_list, list(self.variables))

    def multiplicity(self, point):
        """
        Returns the multiplicity of a singular point on the region.
        
        A singular point (x,y) of region is said to be of multiplicity m 
        if all the partial derivatives off to order m - 1 vanish there.
        
        Examples
        ========
        
        >>> I = ImplicitRegion((x, y), x**2 + y**3 - z**4)
        >>> I.singular_points()
        FiniteSet((0, 0))
        >>> I.multiplicity((0, 0))
        2
        
        """
        if point not in self.singular_points():
            raise ValueError()

        modified_eq = self.equation

        for i, var in enumerate(self.variables):
            modified_eq = modified_eq.subs(var, var + point[i])
        modified_eq = expand(modified_eq)
        
        terms = modified_eq.args
        m = min([total_degree(term) for term in terms])
        return m
        
    def rational_parametrization(self, singular_point=None, parameter='theta'):
        """
        Returns a rational parametrization of implict region.
        
        Examples
        ========
        
        >>> parabola = ImplicitRegion((x, y), y**2 - 4*x)
        >>> parabola.rational_parametrization()
        (4/t**2, 4/t)
        
        >>> circle = ImplicitRegion((x, y), Eq(x**2 + y**2, 4))
        >>> circle.rational_parametrization()
        (2 + 4/(t**2 + 1), 4*t/(t**2 + 1))
        
        >>> I = ImplicitRegion((x, y), x**3 + x**2 - y**2)
        >>> I.rational_parametrization()
        (t**2 - 1, t*(t**2 - 1))
        
        >>> cubic_curve = ImplicitRegion((x, y), x**3 + x**2 - y**2)
        >>> cubic_curve.rational_parametrization()
        (t**2 - 1, t*(t**2 - 1))
        
        >>> sphere = ImplicitRegion((x, y, z), x**2 + y**2 + z**2 - 4) 
        >>> sphere.rational_parametrization()
        (-2 + 4/(s**2 + t**2 + 1), 4*s/(s**2 + t**2 + 1), 4*t/(s**2 + t**2 + 1))
         
        """
        equation = self.equation
        degree = self.degree
        
        if degree == 1:
            if len(self.variables) == 1:
                return equation
            elif len(self.variables) == 2:
                x, y = self.variables
                y_par = list(solveset(equation, x))[0]
                return x, y_par 
            else:
                raise NotImplementedError()
        
        point = ()
        if len(self.singular_points()) != 0:
            singular_points = self.singular_points()
            for spoint in singular_points:
                if self.multiplicity(spoint) == degree - 1:
                    point = spoint                                
        elif degree == 2 and len(point) == 0:
                point = self.regular_point()
                    
        modified_eq = equation
        
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

            s = _symbol('s', real=True)
            t = _symbol('t', real=True)
            
            hn = hn.subs({self.variables[0]: s, self.variables[1]: t})
            hn_1 = hn_1.subs({self.variables[0]: s, self.variables[1]: t})
            
            x_par = (s*(hn_1/hn)).subs(s, 1) + point[0]
            y_par = (t*(hn_1/hn)).subs(s, 1) + point[1]
            
            return x_par, y_par
        
        elif len(self.variables) == 3:
            r = _symbol('r', real=True)
            s = _symbol('s', real=True)
            t = _symbol('t', real=True)
            
            hn = hn.subs({self.variables[0]: r, self.variables[1]: s, self.variables[2]: t})
            hn_1 = hn_1.subs({self.variables[0]: r, self.variables[1]: s, self.variables[2]: t})
            
            x_par = (r*(hn_1/hn)).subs(r, 1) + point[0]
            y_par = (s*(hn_1/hn)).subs(r, 1) + point[1]
            z_par = (t*(hn_1/hn)).subs(r, 1) + point[2]
            print(hn, hn_1, "eah")
            return x_par, y_par, z_par
            
        raise NotImplementedError()
