
class Printer(object):
    """
    """
    def __init__(self):
        self._depth = -1
        self._str = str
        self.emptyPrinter = str

    def doprint(self, expr):
        """Returns the pretty representation for expr (as a string)"""
        return self._str(self._print(expr))

    def _print(self, expr):
        self._depth += 1

        # See if the class of expr is known, or if one of its super
        # classes is known, and use that pretty function
        res = None
        for cls in expr.__class__.__mro__:
            if hasattr(self, '_print_'+cls.__name__):
                res = getattr(self, '_print_'+cls.__name__)(expr)
                break

        # Unknown object, just use its string representation
        if res is None:
            res = self.emptyPrinter(expr)

        self._depth -= 1
        return res