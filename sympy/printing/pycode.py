from sympy.printing.codeprinter import CodePrinter

_kw_py2and3 = {
    'and', 'as', 'assert', 'break', 'class', 'continue', 'def', 'del', 'elif',
    'else', 'except', 'finally', 'for', 'from', 'global', 'if', 'import', 'in',
    'is', 'lambda', 'not', 'or', 'pass', 'raise', 'return', 'try', 'while',
    'with', 'yield'
}
_kw_only_py2 = {'exec', 'print'}
_kw_only_py3 = {'False', 'None', 'nonlocal', 'True'}

_known_functions = {
    'Abs': 'abs'
}

class PythonCodePrinter(CodePrinter):
    printmethod = "_pythoncode"
    language = "Python"
    standard = "python3"
    reserved_words = _kw_py2and3.union(_kw_only_py2)
    _default_settings = {
        'user_functions': {},
        'human': True,
        'error_on_reserved': False,
        'reserved_word_suffix': '_',
    }
    _kf = _known_functions

    def __init__(self, settings={}):
        super(PythonCodePrinter, self).__init__(settings)
        self.known_functions = dict(self._kf, **settings.get(
            'user_functions', {}))

    def _format_code(self, lines):
        return lines
