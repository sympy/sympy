from sympy.core import (Basic, Expr, Symbol, symbols, sympify, diff, Pow, Mul,
                        Add, S)
from sympy.integrals.integrals import integrate, Integral
from sympy.vector.vector import (_all_coordinate_systems, _all_base_scalars,
                                 _vect_add, _vect_mul, _coord_sys_scalar_list,
                                 _has_base_scalar)
from sympy.vector.vector import (dot, cross, express, grad, div, curl,
                                 laplacian)
from sympy.vector.vector import (BaseScalar, Vector, BaseVector, VectMul,
                                 ZeroVector)


class ParamRegion(object):
    """
    A class to represent parametric region in space
    """
    def __init__(self, params, coord_sys, definition, bounds):
        """
        Create a ParamRegion object to represent a parametrically defined
        region in space.
        params : A tuple of length <=2. Both elements of the tuple are symbols
        that act as parameters.
        coord_sys : an instance of subclass of the CoordSys class.
        definition : a tuple of length 3. Each element corresponds to the
        paramentric definition of the BaseScalars for the coord_sys
        bounds : a tuple of 2 tuples. Each inner tuple has 2 elements which
        act as bounds for the Symbols in params, respectively.
        """
        # sanity check
        if not len(params) <= 2:
            raise ValueError("params should be a tuple of length 2")
        if(not isinstance(params[0], Symbol) or
           not isinstance(params[1], Symbol)):
            raise ValueError("all elements of params should be SymPy Symbols")
        if not len(definition) == 3:
            raise ValueError("definition is a tuple of length 3")

        # Let's call the first parameter 'u' and the second 'v'
        if len(params) == 1:
            self._u = params[0]
        elif len(params) == 2:
            self._u = params[0]
            self._v = params[1]

        self.coord_sys = coord_sys
        self.definition = definition

    @property
    def u(self):
        return self._u

    @property
    def v(self):
        return self._v


class VectIntegral(object):
    """
    Base class for representing different type of integrals of vector fields.
    Not to be initialized directly by the user.
    """
    def __init__(self, vect, param_region):
        self.vect = vect
        self.param_region = param_region

    # More methods that are common to all vector integral classes will go here.


class LineVectIntegral(VectIntegral):
    """
    A container for holding line intergrals of vector fields.
    Holds the vector field and the ParamRegion objects.

    Not to be initialized directly. An instance of this class is returned
    by the integrate method on vectors.
    """
    def __init__(self, vect, param_region, coord_sys):
        # The integral must be of the type F.dl
        super(LineVectIntegral, self).__init__(vect, param_region)

    def eval(self):
        """
        Evaluate the integral symbolically, if possible. Else, return the
        LineVectIntegral object back.
        """
        vect = self.vect.express(self.param_region.coord_sys)
        # calculate the differential length element for coord_sys
        dl = self.param_region.coord_sys.h_list
        base_scalars = self.param_region.coord_sys.base_scalars
        # Calculate differntials for each of the elements of dl using
        # the definitions in the param_region object
        # Also note that because this is a F.dl type integral, therefore,
        # there is only one parameter; self.param_region.
        for i in range(len(dl)):
            dl[i] = diff(self.param_region.definition[i], self.param_region.u)

        # Now since vect and dl are in the same coorfinates, we have:
        comp = vect.components
        expr = S.Zero
        for i in range(len(dl)):
            expr = expr + (comp[i] * dl[i])
        # Now, expr has been completely express in coord_sys. We just need to
        # change everything to the parametric variable, 'u'
        for i, scalar in enumerate(base_scalars):
            subs_dict = {scalar: self.param_region.definition[i]}
            expr = expr.subs(subs_dict)

        # Now, we can integrate the expression just as any regular integral
        bounds = self.param_region.bounds
        res = integrate(expr,
                        (self.param_region.u, bounds[0][0], bounds[0][1]))
        if not isinstance(res, Integral):
            return res
        return self


class SurfaceVectIntegral(VectIntegral):
    """
    A container for holding surface intergrals of vector fields.
    Holds the vector field and the ParamRegion objects.

    Not to be initialized directly. An instance of this class is returned
    by the integrate method on vectors.
    """
    def __init__(self, vect, param_region):
        # The integral must be of the type F.dl
        super(LineVectIntegral, self).__init__(vect, param_region)

    def eval(self):
        vect = self.vect.express(self.param_region.coord_sys)
        # First, we need to calculate dS (differntial surface element vector)
        # if r(u, v) is the position vector of any point on the surface
        # normal = (d/du)r x (d/dv)r
        r = ZeroVector
        base_vectors = self.param_region.coord_sys.base_vectors
        for i in range(len(base_vectors)):
            r = r + self.param_region.definitions[i] * base_vectors[i]

        r = r.expand()
        r_u = r.diff(self.param_region.u)
        r_v = r.diff(self.param_region.v)

        # TODO : Implement normalize
        n = cross(r_u, r_v)
        # http://en.wikipedia.org/wiki/Surface_integral
        vect_dot_n = dot(vect, n, self.param_region.coord_sys.param_region)

        # Now we need to calculate double_integral( vect_dot_n du dv )
        bounds = self.param_region.bounds

        # Check the bounds to ascertain the order
        # Case 1 : Both bounds are constants - can perform either order
        # Case 2 : Bounds for u are in terms of v : integrate wrt to u first
        # Case 3 : Bounds for v are in terms of u : integrate wrt to v first
        case = _bounds_case(bounds, self.param_region.u, self.param_region.v)
        if case == 1 or case == 2:
            res = integrate(vect_dot_n,
                            (self.param_region.u, bounds[0][0], bounds[0][1]),
                            (self.param_region.v, bounds[1][0], bounds[1][1]))
            if not isinstance(res, Integral):
                return res
            elif case == 2:
                return self

        if case == 1 or case == 3:
            res = integrate(vect_dot_n,
                            (self.param_region.v, bounds[0][0], bounds[0][1]),
                            (self.param_region.u, bounds[1][0], bounds[1][1]))
            if not isinstance(res, Integral):
                return res

        # Integral couldn't be evaluated
        return self

def _bounds_case(bounds, u, v):
    """
    Check whether the given bounds depend on variables of integration and
    return the correct case.
    """
    # Check bounds of u for v
    lower = bounds[0][0]
    upper = bounds[0][1]

    atoms_lower = lower.atoms()
    atoms_upper = upper.atoms()
    if atoms_lower.issuperset(set([v])) or atoms_upper.issuperset(set([v])):
        return 2

    # Check bounds of v for u
    lower = bounds[1][0]
    upper = bounds[1][1]

    atoms_lower = lower.atoms()
    atoms_upper = upper.atoms()
    if atoms_lower.issuperset(set([u])) or atoms_upper.issuperset(set([u])):
        return 3

    # No bounds occur in the other. Bounds are constants.
    return 1

def vect_integrate(vect, param_region):
    """
    Takes an object with is_Vector == True and a ParamRegion object and
    retuns an instance of the appropriate subclass of VectIntegral class.
    vect : an is_Vector == True object
    param_region : a ParamRegion object
    coord_sys : an instance of subclass of CoordSys class
    """
    # First check if type of the integral
    if hasattr(param_region, '_u') and not hasattr(param_region, '_v'):
        # Integral is of type F.dl
        obj = LineVectIntegral(vect, param_region)
    elif hasattr(param_region, '_u') and hasattr(param_region, '_v'):
        # Integral is a of type F.dS (surface integral)
        obj = SurfaceVectIntegral(vect, param_region)
    # We can also consider cases with scalar*vector type integral
    # TODO: Implement the above comment
    return obj
