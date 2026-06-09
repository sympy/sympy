from .jscode import AbstractJavaFamilyPrinter, known_functions as _known_functions

_known_constants = {
    "Infinity": "Double.POSITIVE_INFINITY"
}

class JavaPrinter(AbstractJavaFamilyPrinter):

    printmethod = '_java'
    language = 'Java'

    def __init__(self, settings=None):
        super().__init__(settings)
        # Known functions and constants handler
        self.known_functions = dict(_known_functions, **(settings or {}).get(
            'user_functions', {}))
        self.known_constants = dict(_known_constants, **(settings or {}).get(
            'user_constants', {}))

    def _print_isnan(self, e):
        arg, = e.args
        return f"Double.isNaN({self._print(arg)})"

    def _print_Raise(self, e):
        arg, = e.args
        return f"throw new {self._print(arg)}"

    def _print_RuntimeError_(self, e):
        arg, = e.args
        return f"Exception({self._print(arg)})"

    def _print_Print(self, e):
        return "System.out." + super()._print_Print(e)


def javacode(expr, assign_to=None, **settings):
    return JavaPrinter(settings).doprint(expr, assign_to)
