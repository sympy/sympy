from sympy import Basic, Function, Symbol
from sympy.printing.codeprinter import CodePrinter
from sympy.core.function import _coeff_isneg
from sympy.printing.precedence import precedence
from sympy.core.compatibility import string_types, range
from sympy.core import S
from sympy.codegen.ast import Assignment

known_functions = {
    'Abs': 'abs',
    'sin': 'sin',
    'cos': 'cos',
    'tan': 'tan',
    'acos': 'acos',
    'asin': 'asin',
    'atan': 'atan',
    'atan2': 'atan2',
    'ceiling': 'ceil',
    'floor': 'floor',
    'sign': 'sign',
    'exp': 'exp',
    'log': 'log',
}

class GLSLPrinter(CodePrinter):
    """
    Rudimentary GLSL printing tools.  Capable of printing with mul/add/sub functions instead of
    operators.  Prints pow(b,n) instead of b**n.   Goals: support for emulated double/quad/octal
    precision.


    Available settings:
    'use_operators': Boolean (should the printer use operators for +,-,*, or functions?)
    """
    _not_supported = set()
    printmethod = "_ccode"
    language = "C"

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 32,
        'user_functions': {},
        'human': True,
        'contract': True,
        'error_on_reserved': False,
        'reserved_word_suffix': '_',
        'use_operators': True,
        'element_wrap': '%s',
        'add_wrap': 'add(%s, %s)',
        'sub_wrap': 'sub(%s, %s)',
        'mul_wrap': 'mul(%s, %s)',
        'pow_wrap': 'pow(%s, %s)',
        'abs_wrap': 'abs(%s)'
    }

    def __init__(self, settings={}):
        CodePrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        self.known_functions.update(userfuncs)

###########################

    def _rate_index_position(self, p):
        return p*5

    def _get_statement(self, codestring):
        return "%s;" % codestring

    def _get_comment(self, text):
        return "// {0}".format(text)

    def _declare_number_const(self, name, value):
        return "float {0} = {1};".format(name, value)

    def _format_code(self, lines):
        return self.indent_code(lines)

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        tab = "   "
        inc_token = ('{', '(', '{\n', '(\n')
        dec_token = ('}', ')')

        code = [line.lstrip(' \t') for line in code]

        increase = [int(any(map(line.endswith, inc_token))) for line in code]
        decrease = [int(any(map(line.startswith, dec_token))) for line in code]

        pretty = []
        level = 0
        for n, line in enumerate(code):
            if line == '' or line == '\n':
                pretty.append(line)
                continue
            level -= decrease[n]
            pretty.append("%s%s" % (tab*level, line))
            level += increase[n]
        return pretty

    def _traverse_matrix_indices(self, mat):
        rows, cols = mat.shape
        return ((i, j) for i in range(rows) for j in range(cols))

    def _print_MatrixElement(self, expr):
        return "{0}[{1}]".format(expr.parent, expr.j +
                expr.i*expr.parent.shape[1])

    def _get_loop_opening_ending(self, indices):
        open_lines = []
        close_lines = []
        loopstart = "for (int %(varble)s=%(start)s; %(varble)s<%(end)s; %(varble)s++){"
        for i in indices:
            # GLSL arrays start at 0 and end at dimension-1
            open_lines.append(loopstart % {
                'varble': self._print(i.label),
                'start': self._print(i.lower),
                'end': self._print(i.upper + 1)})
            close_lines.append("}")
        return open_lines, close_lines


    # def _print_Exp1(self, expr):
    #     return '2.7182818284590452353602874713527'
    #
    # def _print_Pi(self, expr):
    #     return '3.1415926535897932384626433832795'

    def _print_Piecewise(self, expr):
        if expr.args[-1].cond != True:
            # We need the last conditional to be a True, otherwise the resulting
            # function may not return a result.
            raise ValueError("All Piecewise expressions must contain an "
                             "(expr, True) statement to be used as a default "
                             "condition. Without one, the generated "
                             "expression may not evaluate to anything under "
                             "some condition.")
        lines = []
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s) {" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
                    lines.append("else {")
                else:
                    lines.append("else if (%s) {" % self._print(c))
                code0 = self._print(e)
                lines.append(code0)
                lines.append("}")
            return "\n".join(lines)
        else:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that inline operators will
            # not work for statements that span multiple lines (Matrix or
            # Indexed expressions).
            ecpairs = ["((%s) ? (\n%s\n)\n" % (self._print(c), self._print(e))
                    for e, c in expr.args[:-1]]
            last_line = ": (\n%s\n)" % self._print(expr.args[-1].expr)
            return ": ".join(ecpairs) + last_line + " ".join([")"*len(ecpairs)])

    def _print_Idx(self, expr):
        return self._print(expr.label)

    def _print_Indexed(self, expr):
        # calculate index for 1d array
        dims = expr.shape
        elem = S.Zero
        offset = S.One
        for i in reversed(range(expr.rank)):
            elem += expr.indices[i]*offset
            offset *= dims[i]
        return "%s[%s]" % (self._print(expr.base.label), self._print(elem))

###############################


    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp == -1:
            return '1.0/%s' % (self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5:
            return 'sqrt(%s)' % self._print(expr.base)
        else:
            try:
                e = self._print(float(expr.exp))
            except TypeError:
                e = self._print(expr.exp)
            return self._settings['pow_wrap'] % (self._print(expr.base),e)
    #
    # def _print_Pow(self, expr):
    #     return self._settings['pow_wrap'] % (self._print(expr.base),
    #                              self._print(float(expr.exp)))

    def _print_Abs(self, expr):
        return self._settings['abs_wrap'] % self._print(expr.args[0])

    def _print_int(self, expr):
        return str(float(expr))

    def _print_Add(self, expr, order=None):
        if(self._settings['use_operators']):
            return CodePrinter._print_Add(self,expr,order)

        terms = list(expr.args)

        def partition(p,l):
            return reduce(lambda x, y: (x[0]+[y], x[1]) if p(y) else (x[0], x[1]+[y]), l,  ([], []))
        def add(a,b):
            return self._settings['add_wrap'] % (a,b)
        neg, pos = partition(lambda arg: _coeff_isneg(arg), terms)
        s = pos = reduce(lambda a,b: add(a,b), map(lambda t: self._print(t),pos))
        if(len(neg) > 0):
            # sum the absolute values of the negative terms
            neg = reduce(lambda a,b: add(a,b), map(lambda n: self._print(-n),neg))
            # then subtract them from the positive terms
            s = self._settings['sub_wrap'] % (pos,neg)
        return s

    def _print_Mul(self, expr, order=None):
        if(self._settings['use_operators']):
            return CodePrinter._print_Mul(self,expr)

        terms = list(expr.args)
        def partition(p,l):
            return reduce(lambda x, y: (x[0]+[y], x[1]) if p(y) else (x[0], x[1]+[y]), l,  ([], []))
        def mul(a,b):
            return self._settings['mul_wrap'] % (a,b)
        s = reduce(lambda a,b: mul(a,b), map(lambda t: self._print(t),terms))
        return s

    # def makeGL(self,expr,assign_to=None,):
    #         return self.doprint(expr,assign_to)

def glslcode(expr,assign_to=None,**settings):
    return GLSLPrinter(settings).doprint(expr,assign_to)

def print_glsl(expr, **settings):
    """Prints the GLSL representation of the given expression.

       See GLSLPrinter init function for settings.
    """
    print(glslcode(expr, **settings))
