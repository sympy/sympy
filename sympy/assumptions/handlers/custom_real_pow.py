from sympy.assumptions.handlers import AskHandler

class CustomRealPowHandler(AskHandler):
    def _eval_ask(self, expr, assumptions):
        base, exp = expr.args

        # Special case: 0 ** -1 is not real
        if base.is_zero and exp.is_negative:
            return False

        # Let SymPy handle everything else
        return None
