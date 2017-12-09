from sympy.external import import_module

LaTeXSyntaxError = import_module('sympy.parsing._latex',
                                 __import__kwargs={'fromlist': ['_latex']}).LaTeXSyntaxError


def parse_latex(latex_str):
    _latex = import_module('sympy.parsing._latex',
                           __import__kwargs={'fromlist': ['_latex']})

    if _latex is not None:
        return _latex.parse_latex(latex_str)
