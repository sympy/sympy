from sympy.assumptions import Q
from sympy.assumptions.handlers import CommonHandler
from sympy.logic.boolalg import conjuncts

class AskAssociativeHandler(CommonHandler):
    """
    Handler for key 'associative'
    """
    @staticmethod
    def Map(expr, assumptions):
        assumps = conjuncts(assumptions)
        if Q.associative(expr) in assumps:
            return True
        if ~Q.associative(expr) in assumps:
            return False
        if expr.is_associative is not None:
            return expr.is_associative

class AskCommutativeMapHandler(CommonHandler):
    """
    Handler for key 'commutative_map'.

    Explanation
    ===========

    `Q.commutative` deals with whether the expression is commutative
    with respect to multiplication; e.g. ``x*y == y*x``.  
    `Q.commutative_map` deals with whether the map commutates its arguments;
    e.g. ``f(x, y) == f(y, x)``.

    See Also
    ========

    assumptions.handlers.common.AskCommutativeHandler

    """
    @staticmethod
    def Map(expr, assumptions):
        # Unlike other handler, which implies that the object is commutative
        # with respect to multiplication, this handler implies that this map
        # commutes the arguments.
        assumps = conjuncts(assumptions)
        if Q.commutative_map(expr) in assumps:
            return True
        if ~Q.commutative_map(expr) in assumps:
            return False
        return expr.is_commutative_map
