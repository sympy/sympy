def test_custom_AskHandler():
    from sympy.assumptions import register_handler, ask, Q
    from sympy.assumptions.handlers import AskHandler
    from sympy.logic.boolalg import conjuncts
    from sympy import Symbol

    class MersenneHandler(AskHandler):
        @staticmethod
        def Integer(expr, assumptions):
            from sympy import log
            if ask(Q.integer(log(expr + 1, 2))):
                return True
        @staticmethod
        def Symbol(expr, assumptions):
            if expr in conjuncts(assumptions):
                return True
    register_handler('mersenne', MersenneHandler)

    n = Symbol('n', integer=True)
    assert ask(Q.mersenne(n), Q.mersenne(n))
