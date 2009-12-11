from sympy.logic.boolalg import conjuncts
from sympy.assumptions import Q, ask

class AskHandler(object):
    """Base class that all Ask Handlers must inherit"""
    pass

class CommonHandler(AskHandler):
    """Defines some useful methods common to most Handlers """

    @staticmethod
    def NaN(expr, assumptions):
        return False

class AskCommutativeHandler(CommonHandler):
    """
    Handler for key 'commutative'
    """

    @staticmethod
    def Symbol(expr, assumptions):
        """Objects are expected to be commutative unless otherwise stated"""
        if assumptions is True: return True
        for assump in conjuncts(assumptions):
            if assump.expr == expr and assump.key == 'commutative':
                return assump.value
        return True

    @staticmethod
    def Basic(expr, assumptions):
        for arg in expr.args:
            if not ask(arg, Q.commutative, assumptions):
                return False
        return True

    @staticmethod
    def Number(expr, assumptions):
        return True

    @staticmethod
    def NaN(expr, assumptions):
        return True
