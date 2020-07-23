from sympy.assumptions.handlers import CommonHandler

class AskAssociativeHandler(CommonHandler):
    """
    Handler for key 'associative'
    """
    @staticmethod
    def BinaryOperator(expr, assumptions):
        result = getattr(expr, 'is_associative', None)
        if result is not None:
            return result
