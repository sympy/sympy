from sympy.printing.codeprinter import CodePrinter

_kw_py2and3 = {
    'and', 'as', 'assert', 'break', 'class', 'continue', 'def', 'del', 'elif',
    'else', 'except', 'finally', 'for', 'from', 'global', 'if', 'import', 'in',
    'is', 'lambda', 'not', 'or', 'pass', 'raise', 'return', 'try', 'while',
    'with', 'yield', 'None'  # 'None' is actually not in Python 2's keyword.kwlist
}
_kw_only_py2 = {'exec', 'print'}
_kw_only_py3 = {'False', 'nonlocal', 'True'}

_known_functions = {
    'Abs': 'abs'
}

class PythonCodePrinter(CodePrinter):
    printmethod = "_pythoncode"
    language = "Python"
    standard = "python3"
    reserved_words = _kw_py2and3.union(_kw_only_py3)
    _kf = _known_functions
    _operators = {'and': 'and', 'or': 'or', 'not': 'not'}
    _indentation = '    '
    _default_settings = dict(CodePrinter._default_settings, precision=17)

    def __init__(self, settings=None):
        super(PythonCodePrinter, self).__init__(settings)
        self.known_functions = dict(self._kf, **(settings or {}).get(
            'user_functions', {}))

    def _format_code(self, lines):
        return lines

    def _get_comment(self, text):
        return "  # {0}".format(text)

    def _print_Mod(self, expr):
        return self.parenthesize('{0} % {1}'.format(*map(self._print, expr.args)))

    def _print_Piecewise(self, expr):
        lines = []
        for i, (e, c) in enumerate(expr.args):
            code0 = None
            if i == 0:
                lines.append("if %s:" % self._print(c))
            elif i == len(expr.args) - 1 and c == True:
                lines.append('else:')
            else:
                lines.append('elif %s' % self._print(c))
            lines.append(_indentation + (code0 or self._print(e)))
            if i == len(expr.args) - 1 and c != True:
                lines.append('else:')
                lines.append(_indentation + "raise NotImplementedError('Unhandled condition')")
