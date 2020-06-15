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
    >>> semisphere = ParametricRegion((theta, phi), (2*sin(phi)*cos(theta), 2*sin(phi)*sin(theta), 2*cos(phi))
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
        base_vectors = coord_sys.base_vectors()
        base_scalars = coord_sys.base_scalars()

        parametricfield = field

        r = Vector.zero
        for i in range(len(parametricregion.definition)):
            r += base_vectors[i]*parametricregion.definition[i]
    
        for i in range(len(parametricregion.definition)):
            parametricfield = parametricfield.subs(base_scalars[i], parametricregion.definition[i])

        if parametricregion.dimension == 1:
            parameter = parametricregion.parameters[0]

            r_diff = diff(r, parameter)        

            lower, upper = parametricregion.limits[parameter][0], parametricregion.limits[parameter][1]
            
            if isinstance(parametricfield, Vector):
                result = integrate(r_diff.dot(parametricfield), (parameter, lower, upper))
            else:
                result = integrate(r_diff.magnitude()*parametricfield, (parameter, lower, upper))
                
        elif parametricregion.dimension == 2:
            u, v = parametricregion.parameters[0], parametricregion.parameters[1]
            
            r_u = diff(r, u)
            r_v = diff(r, v)
            normal_vector = r_u.cross(r_v)
            
            lower_u, upper_u = parametricregion.limits[u][0], parametricregion.limits[u][1]
            lower_v, upper_v = parametricregion.limits[v][0], parametricregion.limits[v][1]
            
            if isinstance(parametricfield, Vector):
                result = integrate(parametricfield.dot(normal_vector), (u, lower_u, upper_u),\
                            (v, lower_v, upper_v))
            else:
                result = integrate(parametricfield*normal_vector.magnitude(), (u, lower_u, upper_u), (v, lower_v, upper_v))
        
        return result   
            
    @property
    def field(self):
        return self.args[0]
    
    @property
    def parametricregion(self):
        return self.args[1]  
