from sympy.assumptions.handlers import CommonHandler

class AskAssociativeHandler(CommonHandler):
    """
    Handler for key 'associative'
    """
    @staticmethod
    def BinaryOperator(expr, assumptions):
        result = getattr(expr, 'associative', None)
        if result is not None:
            return result

    @staticmethod
    def LeftDivision(expr, assumptions):
        return False

    @staticmethod
    def RightDivision(expr, assumptions):
        return False

class AskLeftDivisibleHandler(CommonHandler):
    r"""
    Handler for key 'left_divisible'

    Represents that a left division operator $\backslash$ exists for binary
    operator *expr*.
    """
    @staticmethod
    def BinaryOperator(expr, assumptions):
        result = getattr(expr, 'left_divisible', None)
        if result is not None:
            return result

class AskRightDivisibleHandler(CommonHandler):
    """
    Handler for key 'right_divisible'

    Represents that a right division operator $/$ exists for binary
    operator *expr*.
    """
    @staticmethod
    def BinaryOperator(expr, assumptions):
        result = getattr(expr, 'right_divisible', None)
        if result is not None:
            return result
