from sympy.core.basic import Basic
from sympy.vector import CoordSys3D, Vector
from sympy.vector.operators import _get_coord_sys_from_expr
from sympy.integrals import Integral, integrate
from sympy.core.function import diff

class ParametricIntegral(Basic):
    """
    Represents integral of a scalar or vector field
    over a Parametric Region

    Examples
    ========

    >>> from sympy.vector import CoordSys3D, ParametricRegion
    >>> from sympy.abc import t

    >>> C = CoordSys3D('C)
    >>> curve = ParametricRegion(t, (3*t - 2, t + 1), {t: (1, 2)})

    >>> ParametricIntegral(C.x, curve)
    0
    >>> semisphere = ParametricRegion((phi, theta), (2*sin(phi)*cos(theta), 2*sin(phi)*sin(theta), 2*cos(phi))
                            {theta: (0, 2*pi), phi: (0, pi/2)})
    >>> ParametricIntegral(C.z, semisphere)
    8*pi

    >>> ParametricIntegral(C.j - C.k, ParametricRegion((r, theta), (r*cos(theta), r*sin(theta))))
    ParametricIntegral(C.j - C.k, ParametricRegion((r, theta), (r*cos(theta), r*sin(theta))))

    """
    
    def __new__(cls, field, parametricregion):
        obj = super().__new__(cls, field, parametricregion)
        return obj

    def eval(self):
        parametricregion = self.parametricregion
        field = self.field
        
        coord_sys = _get_coord_sys_from_expr(field)
        
        if len(coord_sys) > 1:
            raise ValueError
        
        coord_sys = next(iter(coord_sys))

        if parametricregion.dimension == 1:
            parameter = parametricregion.parameters[0]
            
            r_diff_tuple = [diff(comp, parameter) for comp in parametricregion.definition]
            
            base_vectors = coord_sys.base_vectors()
            base_scalars = coord_sys.base_scalars()
            r_diff = Vector.zero
            parametricfield = field
            
            for i in range(len(parametricregion.definition)):
                r_diff += base_vectors[i]*r_diff_tuple[i]        
            
            for i in range(len(parametricregion.definition)):
                parametricfield = parametricfield.subs(base_scalars[i], parametricregion.definition[i])
            
            lower, upper = parametricregion.limits[parameter][0], parametricregion.limits[parameter][1]
            
            if isinstance(parametricfield, Vector):
                result = integrate(r_diff.dot(parametricfield), (parameter, lower, upper))
            else:
                result = integrate(r_diff.magnitude()*parametricfield, (parameter, lower, upper))
                
            return result            
            
    @property
    def field(self):
        return self.args[0]
    
    @property
    def parametricregion(self):
        return self.args[1]  
