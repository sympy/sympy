# sympy/galgebra/printing.py

"""
printing.py implements GA_Printer and GA_LatexPrinter classes for
multivectors by derivation from the SymPy Printer and LatexPrinter
classes.  Where we wish to change the behavior of printing for
existing SymPy types the required function have been rewritten
and included in GA_Printer and GA_LatexPrinter.  For example
we have rewritten:

    _print_Derivative
    _print_Function
    _print_Matrix
    _print_Pow
    _print_Symbol

for GA_LatexPrinter and

    _print_Derivative
    _print_Function

for GA_Printer.

There is also and enhanced_print class that allows multivectors,
multivector functions, and multivector derivatives to print out
in different colors on ansi terminals.
"""

from __future__ import print_function

import os
import sys

from sympy.core.compatibility import StringIO, range
from sympy.core.operations import AssocOp
from sympy.core.power import Pow

from sympy import S, Basic, Symbol, Matrix

from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.conventions import split_super_sub

SYS_CMD = {'linux2': {'rm': 'rm', 'evince': 'evince', 'null': ' > /dev/null', '&': '&'},
           'win32': {'rm': 'del', 'evince': '', 'null': ' > NUL', '&': ''},
           'darwin': {'rm': 'rm', 'evince': 'open', 'null': ' > /dev/null', '&': '&'}}

ColorCode = {
    'black':        '0;30',         'bright gray':   '0;37',
    'blue':         '0;34',         'white':         '1;37',
    'green':        '0;32',         'bright blue':   '1;34',
    'cyan':         '0;36',         'bright green':  '1;32',
    'red':          '0;31',         'bright cyan':   '1;36',
    'purple':       '0;35',         'bright red':    '1;31',
    'yellow':       '0;33',         'bright purple': '1;35',
    'dark gray':    '1;30',         'bright yellow': '1;33',
    'normal':       '0'
}

InvColorCode = dict(zip(ColorCode.values(), ColorCode.keys()))

accepted_latex_functions = ['arcsin', 'arccos', 'arctan', 'sin', 'cos', 'tan',
                            'theta', 'beta', 'alpha', 'gamma', 'sinh', 'cosh', 'tanh', 'sqrt',
                            'ln', 'log', 'sec', 'csc', 'cot', 'coth', 're', 'im', 'frac', 'root',
                            'arg', 'zeta']


def find_executable(executable, path=None):
    """
    Try to find 'executable' in the directories listed in 'path' (a
    string listing directories separated by 'os.pathsep'; defaults to
    os.environ['PATH']).  Returns the complete filename or None if not
    found
    """
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if os.name == 'os2':
        (base, ext) = os.path.splitext(executable)
        # executable files on OS/2 can have an arbitrary extension, but
        # .exe is automatically appended if no dot is present in the name
        if not ext:
            executable = executable + ".exe"
    elif sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (base, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    else:
        return None


class enhance_print:
    """
    A class for color coding the string printing going to a terminal.
    """

    normal = ''
    base = ''
    fct = ''
    deriv = ''
    bold = ''

    def __init__(self, base=None, fct=None, deriv=None, on=True, keys=False):
        if on:
            if 'win' in sys.platform:

                if base is None:
                    enhance_print.base = ColorCode['blue']
                else:
                    enhance_print.base = ColorCode[base]
                if fct is None:
                    enhance_print.fct = ColorCode['red']
                else:
                    enhance_print.fct = ColorCode[fct]
                if deriv is None:
                    enhance_print.deriv = ColorCode['cyan']
                else:
                    enhance_print.deriv = ColorCode[deriv]
                enhance_print.normal = '\033[0m'

            else:

                if base is None:
                    enhance_print.base = ColorCode['dark gray']
                else:
                    enhance_print.base = ColorCode[base]
                if fct is None:
                    enhance_print.fct = ColorCode['red']
                else:
                    enhance_print.fct = ColorCode[fct]
                if deriv is None:
                    enhance_print.deriv = ColorCode['cyan']
                else:
                    enhance_print.deriv = ColorCode[deriv]
                enhance_print.normal = '\033[0m'

            if keys:
                print('Enhanced Printing is on:')
                print('Base/Blade color is ' + InvColorCode[enhance_print.base])
                print('Function color is ' + InvColorCode[enhance_print.fct])
                print('Derivative color is ' + InvColorCode[enhance_print.deriv] + '\n')

            enhance_print.base = '\033[' + enhance_print.base + 'm'
            enhance_print.fct = '\033[' + enhance_print.fct + 'm'
            enhance_print.deriv = '\033[' + enhance_print.deriv + 'm'

    @staticmethod
    def enhance_base(s):
        return enhance_print.base + s + enhance_print.normal

    @staticmethod
    def enhance_fct(s):
        return enhance_print.fct + s + enhance_print.normal

    @staticmethod
    def enhance_deriv(s):
        return enhance_print.deriv + s + enhance_print.normal

    @staticmethod
    def strip_base(s):
        new_s = s.replace(enhance_print.base, '')
        new_s = new_s.replace(enhance_print.normal, '')
        return new_s


class GA_Printer(StrPrinter):
    """
    An enhanced string printer that is galgebra-aware.
    """

    function_names = ('acos', 'acosh', 'acot', 'acoth', 'arg', 'asin', 'asinh',
                      'atan', 'atan2', 'atanh', 'ceiling', 'conjugate', 'cos',
                      'cosh', 'cot', 'coth', 'exp', 'floor', 'im', 'log', 're',
                      'root', 'sin', 'sinh', 'sqrt', 'sign', 'tan', 'tanh')

    def _print_Function(self, expr):
        name = expr.func.__name__

        if expr.args:
            if name in GA_Printer.function_names:
                return expr.func.__name__ + "(%s)" % self.stringify(expr.args, ", ")

        return enhance_print.enhance_fct("%s" % (name, ))

    def _print_Derivative(self, expr):
        diff_args = list(map(self._print, expr.args))
        return enhance_print.enhance_deriv('D{%s}' % (diff_args[1], )) + '%s' % (diff_args[0], )

    def _print_MV(self, expr):
        if expr.obj.is_zero:
            return '0'
        else:
            if expr.print_blades:
                expr.base_to_blade()
            ostr = expr.get_normal_order_str()
            return ostr

    def _print_Vector(self, expr):
        if expr.obj.is_zero:
            return '0'
        else:
            ostr = GA_Printer().doprint(expr.obj)
            ostr = ostr.replace(' ', '')
            return ostr

    @staticmethod
    def _on():
        GA_Printer.Basic__str__ = Basic.__str__
        Basic.__str__ = lambda self: GA_Printer().doprint(self)
        return

    @staticmethod
    def _off():
        Basic.__str__ = GA_Printer.Basic__str__
        return

    def __enter__(self):
        GA_Printer._on()
        return self

    def __exit__(self, type, value, traceback):
        GA_Printer._off()


class GA_LatexPrinter(LatexPrinter):
    r"""
    An enhanced LaTeX printer that is galgebra-aware.

    The latex printer is turned on with the function (in ga.py) -

        Format(Fmode=True,Dmode=True,ipy=False)

    where Fmode is the function printing mode that surpresses printing arguments,
    Dmode is the derivative printing mode that does not use fractions, and
    ipy=True is the Ipython notebook mode that does not redirect the print output.

    The latex output is post processed and displayed with the function (in GAPrint.py) -

        xdvi(filename='tmplatex.tex',debug=False)

    where filename is the name of the tex file one would keep for future
    inclusion in documents and debug=True would display the tex file
    immediately.

    There are three options for printing multivectors in latex.  They are
    acessed with the multivector member function -

        A.Fmt(self,fmt=1,title=None)

    where fmt=1, 2, or 3 determines whether the entire multivector A is
    printed entirely on one line, or one grade is printed per line, or
    one base is printed per line.  If title is not None then the latex
    string generated is of the form -

        title+' = '+str(A)

    where it is assumed that title is a latex math mode string. If title
    contains '%' it is treated as a pure latex math mode string.  If it
    does not contain '%' then the following character mappings are applied -

        - 'grad' replaced by '\\bm{\\nabla} '
        - '*' replaced by ''
        - '^' replaced by '\\W '
        - '|' replaced by '\\cdot '
        - '>' replaced by '\\lfloor '
        - '<' replaced by '\\rfloor '

    In the case of a print statement of the form -

        print(title,A)

    everthing in the title processing still applies except that the multivector
    formatting is one multivector per line.

    For print statements of the form -

        print(title)

    where no program variables are printed if title contains '#' then title
    is printed as regular latex line.  If title does not contain '#' then
    title is printed in equation mode. '%' has the same effect in title as
    in the Fmt() member function.
    """

    preamble = \
"""
\\pagestyle{empty}
\\usepackage[latin1]{inputenc}
\\usepackage{amsmath}
\\usepackage{bm}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{amsbsy}
\\usepackage{tensor}
\\usepackage{listings}
\\usepackage{color}
\\definecolor{gray}{rgb}{0.95,0.95,0.95}
\\setlength{\\parindent}{0pt}
\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}
\\newcommand{\\lp}{\\left (}
\\newcommand{\\rp}{\\right )}
\\newcommand{\\half}{\\frac{1}{2}}
\\newcommand{\\llt}{\\left <}
\\newcommand{\\rgt}{\\right >}
\\newcommand{\\abs}[1]{\\left |{#1}\\right | }
\\newcommand{\\pdiff}[2]{\\bfrac{\\partial {#1}}{\\partial {#2}}}
\\newcommand{\\lbrc}{\\left \\{}
\\newcommand{\\rbrc}{\\right \\}}
\\newcommand{\\W}{\\wedge}
\\newcommand{\\prm}[1]{{#1}'}
\\newcommand{\\ddt}[1]{\\bfrac{d{#1}}{dt}}
\\newcommand{\\R}{\\dagger}
\\newcommand{\\deriv}[3]{\\bfrac{d^{#3}#1}{d{#2}^{#3}}}
\\newcommand{\\grade}[1]{\\left < {#1} \\right >}
\\newcommand{\\f}[2]{{#1}\\lp{#2}\\rp}
\\newcommand{\\eval}[2]{\\left . {#1} \\right |_{#2}}
\\usepackage{float}
\\floatstyle{plain} % optionally change the style of the new float
\\newfloat{Code}{H}{myc}
\\lstloadlanguages{Python}

\\begin{document}
"""
    postscript = '\\end{document}\n'

    Dmode = False  # True - Print derivative contracted
    Fmode = False  # True - Print function contracted
    latex_flg = False

    @staticmethod
    def redirect(ipy=False):
        GA_LatexPrinter.ipy = ipy
        GA_LatexPrinter.latex_flg = True
        GA_LatexPrinter.Basic__str__ = Basic.__str__
        GA_LatexPrinter.Matrix__str__ = Matrix.__str__
        Basic.__str__ = lambda self: GA_LatexPrinter().doprint(self)
        Matrix.__str__ = lambda self: GA_LatexPrinter().doprint(self)
        if not ipy:
            GA_LatexPrinter.stdout = sys.stdout
            sys.stdout = StringIO()
        return

    @staticmethod
    def restore():
        GA_LatexPrinter.latex_flg = False
        if not GA_LatexPrinter.ipy:
            sys.stdout = GA_LatexPrinter.stdout
        Basic.__str__ = GA_LatexPrinter.Basic__str__
        Matrix.__str__ = GA_LatexPrinter.Matrix__str__
        return

    def _print_Pow(self, expr):
        base = self._print(expr.base)
        if ('_' in base or '^' in base) and 'cdot' not in base:
            mode = True
        else:
            mode = False

        # Treat x**Rational(1,n) as special case
        if expr.exp.is_Rational and abs(expr.exp.p) == 1 and expr.exp.q != 1:
            expq = expr.exp.q

            if expq == 2:
                tex = r"\sqrt{%s}" % base
            elif self._settings['itex']:
                tex = r"\root{%d}{%s}" % (expq, base)
            else:
                tex = r"\sqrt[%d]{%s}" % (expq, base)

            if expr.exp.is_negative:
                return r"\frac{1}{%s}" % tex
            else:
                return tex
        elif self._settings['fold_frac_powers'] \
            and expr.exp.is_Rational \
                and expr.exp.q != 1:
            base, p, q = self._print(expr.base), expr.exp.p, expr.exp.q
            if mode:
                return r"{\lp %s \rp}^{%s/%s}" % (base, p, q)
            else:
                return r"%s^{%s/%s}" % (base, p, q)

        elif expr.exp.is_Rational and expr.exp.is_negative and expr.base.is_Function:
            # Things like 1/x
            return r"\frac{%s}{%s}" % \
                (1, self._print(Pow(expr.base, -expr.exp)))
        else:
            if expr.base.is_Function:
                return self._print(expr.base, self._print(expr.exp))
            else:
                if expr.is_commutative and expr.exp == -1:
                    """
                    solves issue 4129
                    As Mul always simplify 1/x to x**-1
                    The objective is achieved with this hack
                    first we get the latex for -1 * expr,
                    which is a Mul expression
                    """
                    tex = self._print(S.NegativeOne * expr).strip()
                    # the result comes with a minus and a space, so we remove
                    if tex[:1] == "-":
                        return tex[1:].strip()
                if self._needs_brackets(expr.base):
                    tex = r"\left(%s\right)^{%s}"
                else:
                    if mode:
                        tex = r"{\lp %s \rp}^{%s}"
                    else:
                        tex = r"%s^{%s}"

                return tex % (self._print(expr.base),
                              self._print(expr.exp))

    def _print_Symbol(self, expr):

        def str_symbol(name_str):
            (name, supers, subs) = split_super_sub(name_str)

            # translate name, supers and subs to tex keywords
            greek = set(['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                         'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu',
                         'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon',
                         'phi', 'chi', 'psi', 'omega'])

            greek_translated = {'lamda': 'lambda', 'Lamda': 'Lambda'}

            other = set(['aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth',
                        'hbar', 'hslash', 'mho'])

            def translate(s):
                tmp = s.lower()
                if tmp in greek or tmp in other:
                    return "\\" + s
                if s in greek_translated:
                    return "\\" + greek_translated[s]
                else:
                    return s

            name = translate(name)

            if supers != []:
                supers = list(map(translate, supers))

            if subs != []:
                subs = list(map(translate, subs))

            # glue all items together:
            if len(supers) > 0:
                name += "^{%s}" % " ".join(supers)
            if len(subs) > 0:
                name += "_{%s}" % " ".join(subs)

            return name

        if expr in self._settings['symbol_names']:
            return self._settings['symbol_names'][expr]

        name_str = expr.name

        if '.' in name_str and name_str[0] == '(' and name_str[-1] == ')':
            name_str = name_str[1:-1]
            name_lst = name_str.split('.')
            name_str = r'\lp ' + str_symbol(name_lst[0]) + r'\cdot ' + str_symbol(name_lst[1]) + r'\rp '
            return name_str

        return str_symbol(expr.name)

    def _print_Function(self, expr, exp=None):
        func = expr.func.__name__

        name = func

        if hasattr(self, '_print_' + func):
            return getattr(self, '_print_' + func)(expr, exp)
        else:
            args = [str(self._print(arg)) for arg in expr.args]
            # How inverse trig functions should be displayed, formats are:
            # abbreviated: asin, full: arcsin, power: sin^-1
            inv_trig_style = self._settings['inv_trig_style']
            # If we are dealing with a power-style inverse trig function
            inv_trig_power_case = False
            # If it is applicable to fold the argument brackets
            can_fold_brackets = self._settings['fold_func_brackets'] and \
                len(args) == 1 and \
                not self._needs_function_brackets(expr.args[0])

            inv_trig_table = ["asin", "acos", "atan", "acot"]

            # If the function is an inverse trig function, handle the style
            if func in inv_trig_table:
                if inv_trig_style == "abbreviated":
                    func = func
                elif inv_trig_style == "full":
                    func = "arc" + func[1:]
                elif inv_trig_style == "power":
                    func = func[1:]
                    inv_trig_power_case = True

                    # Can never fold brackets if we're raised to a power
                    if exp is not None:
                        can_fold_brackets = False

            if inv_trig_power_case:
                if func in accepted_latex_functions:
                    name = r"\%s^{-1}" % func
                else:
                    name = r"\operatorname{%s}^{-1}" % func
            elif exp is not None:
                if func in accepted_latex_functions:
                    name = r"\%s^{%s}" % (func, exp)
                else:
                    name = latex(Symbol(func))
                    if '_' in func or '^' in func:
                        name = r'{\lp ' + name + r'\rp }^{' + exp + '}'
                    else:
                        name += '^{' + exp + '}'
            else:
                if func in accepted_latex_functions:
                    name = r"\%s" % func
                else:
                    name = latex(Symbol(func))
                    if exp is not None:
                        if '_' in name or '^' in name:
                            name = r'\lp ' + name + r'\rp^{' + exp + '}'
                        else:
                            name += '^{' + exp + '}'

            if can_fold_brackets:
                if func in accepted_latex_functions:
                    # Wrap argument safely to avoid parse-time conflicts
                    # with the function name itself
                    name += r" {%s}"
                else:
                    if not GA_LatexPrinter.Fmode:
                        name += r"%s"
            else:
                if func in accepted_latex_functions or not GA_LatexPrinter.Fmode:
                    name += r"{\left (%s \right )}"

            if inv_trig_power_case and exp is not None:
                name += r"^{%s}" % exp

            if func in accepted_latex_functions or not GA_LatexPrinter.Fmode:
                if len(args) == 1:
                    return name % args[0]
                else:
                    return name % ",".join(args)
            else:
                return name

    def _print_Derivative(self, expr):
        dim = len(expr.variables)

        imax = 1

        if dim == 1:
            if GA_LatexPrinter.Dmode:
                tex = r"\partial_{%s}" % self._print(expr.variables[0])
            else:
                tex = r"\frac{\partial}{\partial %s}" % self._print(expr.variables[0])
        else:
            multiplicity, i, tex = [], 1, ""
            current = expr.variables[0]
            for symbol in expr.variables[1:]:
                if symbol == current:
                    i = i + 1
                else:
                    multiplicity.append((current, i))
                    current, i = symbol, 1

            else:
                imax = max(imax, i)
                multiplicity.append((current, i))

            if GA_LatexPrinter.Dmode and imax == 1:
                tmp = ''
                dim = 0
                for x, i in multiplicity:
                    dim += i
                    tmp += r"%s" % (self._print(x), )
                tex = r"\partial^{%s}_{" % (dim, ) + tmp + '}'
            else:
                for x, i in multiplicity:
                    if i == 1:
                        tex += r"\partial %s" % self._print(x)
                    else:
                        tex += r"\partial^{%s} %s" % (i, self._print(x))
                tex = r"\frac{\partial^{%s}}{%s} " % (dim, tex)

        if isinstance(expr.expr, AssocOp):
            return r"%s\left(%s\right)" % (tex, self._print(expr.expr))
        else:
            return r"%s %s" % (tex, self._print(expr.expr))

    def _print_MV(self, expr):
        if expr.obj.is_zero:
            return '0 \n'
        else:
            if expr.print_blades:
                expr.base_to_blade()
            ostr = expr.get_latex_normal_order_str()
            return ostr

    def _print_Matrix(self, expr):
        lines = []
        for line in range(expr. rows):  # horrible, should be 'rows'
            lines.append(" & ".join([self._print(i) for i in expr[line, :]]))

        ncols = 0
        for i in expr[line, :]:
            ncols += 1

        out_str = '\\left [ \\begin{array}{' + ncols * 'c' + '} '
        for line in lines[:-1]:
            out_str += line + ' \\\\ '
        out_str += lines[-1] + ' \\end{array}\\right ] '
        return out_str


def latex(expr, **settings):
    "Return the LaTeX representation of the given expression."
    return GA_LatexPrinter(settings).doprint(expr)


def print_latex(expr, **settings):
    "Print the LaTeX representation of the given expression."
    print(latex(expr, **settings))


def xdvi(filename=None, debug=False, paper=(14, 11)):
    """
    Post processes LaTeX output (see comments below), adds preamble and
    postscript, generates tex file, inputs file to latex, displays resulting
    pdf file.
    """
    if GA_LatexPrinter.ipy:
        GA_LatexPrinter.restore(GA_LatexPrinter.ipy)
        return

    sys_cmd = SYS_CMD[sys.platform]

    latex_str = sys.stdout.getvalue()
    GA_LatexPrinter.restore()
    latex_lst = latex_str.split('\n')
    latex_str = ''

    lhs = ''
    code_flg = False

    for latex_line in latex_lst:
        if len(latex_line) > 0 and '##' == latex_line[:2]:
            if code_flg:
                code_flg = False
                latex_line = latex_line[2:]
            else:
                code_flg = True
                latex_line = latex_line[2:]
        elif code_flg:
                    pass
        elif len(latex_line) > 0 and '#' in latex_line:  # Non equation mode output (comment)
            latex_line = latex_line.replace('#', '')
            if '%' in latex_line:  # Equation mode with no variables to print (comment)
                latex_line = latex_line.replace('%', '')
                latex_line = '\\begin{equation*} ' + latex_line + ' \\end{equation*}\n'

        else:
            latex_line = latex_line.replace('.', r' \cdot ')  # For components of metric tensor
            if '=' in latex_line:  # determing lhs of equation/align
                eq_index = latex_line.rindex('=') + 1
                lhs = latex_line[:eq_index]
                latex_line = latex_line.replace(lhs, '')
                if '%' in lhs:  # Do not preprocess lhs of equation/align
                    lhs = lhs.replace('%', '')
                else:  # preprocess lhs of equation/align
                    lhs = lhs.replace('|', r'\cdot ')
                    lhs = lhs.replace('^', r'\W ')
                    lhs = lhs.replace('*', ' ')
                    lhs = lhs.replace('grad', r'\bm{\nabla} ')
                    lhs = lhs.replace('<<', r'\rfloor ')
                    lhs = lhs.replace('<', r'\rfloor ')
                    lhs = lhs.replace('>>', r'\lfloor ')
                    lhs = lhs.replace('>', r'\lfloor ')
                latex_line = lhs + latex_line

            if r'\begin{align*}' in latex_line:  # insert lhs of align environment
                latex_line = latex_line = latex_line.replace(lhs, '')
                latex_line = latex_line.replace(r'\begin{align*}', r'\begin{align*} ' + lhs)
                lhs = ''
            else:  # normal sympy equation
                latex_line = latex_line.strip()
                if len(latex_line) > 0:
                    latex_line = '\\begin{equation*} ' + latex_line + ' \\end{equation*}'

        latex_str += latex_line + '\n'

    latex_str = latex_str.replace('\n\n', '\n')

    if debug:
        print(latex_str)

    if paper == 'letter':
        paper_size = \
"""
\\documentclass[10pt,fleqn]{report}
"""
    else:
        paper_size = \
"""
\\documentclass[10pt,fleqn]{report}
\\usepackage[vcentering]{geometry}
"""
        paper_size += '\\geometry{papersize={' + str(paper[0]) + \
                      'in,' + str(paper[1]) + 'in},total={' + str(paper[0] - 1) + \
                      'in,' + str(paper[1] - 1) + 'in}}\n'

    latex_str = paper_size + GA_LatexPrinter.preamble + latex_str + GA_LatexPrinter.postscript

    if filename is None:
        pyfilename = sys.argv[0]
        rootfilename = pyfilename.replace('.py', '')
        filename = rootfilename + '.tex'

    print('latex file =', filename)

    latex_file = open(filename, 'w')
    latex_file.write(latex_str)
    latex_file.close()

    latex_str = None

    pdflatex = find_executable('pdflatex')

    print('pdflatex path =', pdflatex)

    if pdflatex is not None:
        latex_str = 'pdflatex'
    else:
        return

    if latex_str is not None:
        if debug:  # Display latex excution output for debugging purposes
            os.system(latex_str + ' ' + filename[:-4])
        else:  # Works for Linux don't know about Windows
            os.system(latex_str + ' ' + filename[:-4] + sys_cmd['null'])

        print_cmd = sys_cmd['evince'] + ' ' + filename[:-4] + '.pdf ' + sys_cmd['&']
        print(print_cmd)

        os.system(print_cmd)
        raw_input('!!!!Return to continue!!!!\n')

        if debug:
            os.system(sys_cmd['rm'] + ' ' + filename[:-4] + '.aux ' + filename[:-4] + '.log')
        else:
            os.system(sys_cmd['rm'] + ' ' + filename[:-4] + '.aux ' + filename[:-4] + '.log ' + filename[:-4] + '.tex')
    return


prog_str = ''
off_mode = False


def Get_Program(off=False):
    global prog_str, off_mode
    off_mode = off
    if off_mode:
        return
    prog_file = open(sys.argv[0], 'r')
    prog_str = prog_file.read()
    prog_file.close()
    return


def Print_Function():
    global prog_str, off_mode
    if off_mode:
        return
    fct_name = str(sys._getframe(1).f_code.co_name)
    ifct = prog_str.find('def ' + fct_name)
    iend = prog_str.find('def ', ifct + 4)
    tmp_str = prog_str[ifct:iend - 1]
    fct_name = fct_name.replace('_', ' ')
    if GA_LatexPrinter.latex_flg:
        print('##\\begin{lstlisting}[language=Python,showspaces=false,' +
              'showstringspaces=false,backgroundcolor=\color{gray},frame=single]')
        print(tmp_str)
        print('##\\end{lstlisting}')
        print('#Code Output:')
    else:
        print('\n' + 80 * '*')
        print(tmp_str)
        print('Code output:\n')
    return
