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
        assumps = conjuncts(assumptions)
        if Q.commutative(expr) in assumps:
            return True
        elif ~Q.commutative(expr) in assumps:
            return False
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

class TautologicalHandler(AskHandler):
    """Wrapper allowing to query the truth value of a boolean expression."""

    @staticmethod
    def bool(expr, assumptions):
        return expr

    @staticmethod
    def Assume(expr, assumptions):
        return ask(expr.expr, expr.key, assumptions)

    @staticmethod
    def Not(expr, assumptions):
        value = ask(expr.args[0], Q.is_true, assumptions=assumptions)
        if value in (True, False):
            return not value
        else:
            return None


    @staticmethod
    def Or(expr, assumptions):
        result = False
        for arg in expr.args:
            p = ask(arg, Q.is_true, assumptions=assumptions)
            if p == True:
                return True
            if p == None:
                result = None
        return result

    @staticmethod
    def And(expr, assumptions):
        result = True
        for arg in expr.args:
            p = ask(arg, Q.is_true, assumptions=assumptions)
            if p == False:
                return False
            if p == None:
                result = None
        return result

    @staticmethod
    def Implies(expr, assumptions):
        p, q = expr.args
        return ask(~p | q, Q.is_true, assumptions=assumptions)

    @staticmethod
    def Equivalent(expr, assumptions):
        p, q = expr.args
        pt = ask(p, Q.is_true, assumptions=assumptions)
        if pt == None:
            return None
        qt = ask(q, Q.is_true, assumptions=assumptions)
        if qt == None:
            return None
        return pt == qt
