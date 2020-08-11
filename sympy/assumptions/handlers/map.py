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
        if expr.associative is not None:
            return expr.associative
