from sympy.core import (Basic, Expr, Dummy, Function, Symbol, symbols,
                        sympify, diff, Pow, Mul, Add, S, AtomicExpr)

class ParamRegion(object):
    """
    A class to represent parametric region in space
    """
    def __init__(self, params, coord_sys, definition):
        """
        Create a ParamRegion object to represent a parametrically defined
        region in space.
        params : A tuple of length 2. Both elements of the tuple are symbols
        that act as parameters.
        coord_sys : an instance of subclass of the CoordSys class.
        definition : a tuple of length 3. Each element corresponds to the
        paramentric definition of the BaseScalars for the coord_sys
        """
        # sanity check
        if not len(params) == 2:
            raise ValueError("params should be a tuple of length 2")
        if(not isinstance(params[0], Symbol) or
           not isinstance(params[1], Symbol)):
            raise ValueError("all elements of params should be SymPy Symbols")
        if not len(definition) == 3:
            raise ValueError("definition is a tuple of length 3")
