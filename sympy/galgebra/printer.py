#printer.py

import os
import sys
import StringIO
import re
from sympy import Matrix, Basic, S, C, Symbol, Function, Derivative
from itertools import islice
from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter, accepted_latex_functions
from sympy.core.function import _coeff_isneg
from inspect import getouterframes, currentframe
import ga
import mv
import lt
import metric

ip_cmds = \
"""
$\\DeclareMathOperator{\Tr}{Tr}
\\DeclareMathOperator{\Adj}{Adj}
\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}
\\newcommand{\\lp}{\\left (}
\\newcommand{\\rp}{\\right )}
\\newcommand{\\paren}[1]{\\lp {#1} \\rp}
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
\\newcommand{\\eval}[2]{\\left . {#1} \\right |_{#2}}$
"""

SYS_CMD = {'linux2': {'rm': 'rm', 'evince': 'evince', 'null': ' > /dev/null', '&': '&'},
           'win32': {'rm': 'del', 'evince': '', 'null': ' > NUL', '&': ''},
           'darwin': {'rm': 'rm', 'evince': 'open', 'null': ' > /dev/null', '&': '&'}}


def find_executable(executable, path=None):
    """Try to find 'executable' in the directories listed in 'path' (a
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


def ostr(obj, dict_mode=False, indent=True):
    """
    Recursively convert iterated object (list/tuple/dict/set) to string.
    """
    def ostr_rec(obj, dict_mode):
        global ostr_s
        if isinstance(obj, Matrix):
            ostr_s += str(obj)
        elif isinstance(obj, tuple):
            if len(obj) == 0:
                ostr_s += '(),'
            else:
                ostr_s += '('
                for obj_i in obj:
                    ostr_rec(obj_i, dict_mode)
                ostr_s = ostr_s[:-1] + '),'
        elif isinstance(obj, list):
            if len(obj) == 0:
                ostr_s += '[],'
            else:
                ostr_s += '['
                for obj_i in obj:
                    ostr_rec(obj_i, dict_mode)
                ostr_s = ostr_s[:-1] + '],'
        elif isinstance(obj, dict):
            if dict_mode:
                ostr_s += '\n'
                for key in obj.keys():
                    ostr_rec(key, dict_mode)
                    if ostr_s[-1] == ',':
                        ostr_s = ostr_s[:-1]
                    ostr_s += ' -> '
                    ostr_rec(obj[key], dict_mode)
                    if ostr_s[-1] == ',':
                        ostr_s = ostr_s[:-1]
                    ostr_s += '\n'
            else:
                ostr_s += '{'
                for key in obj.keys():
                    ostr_rec(key, dict_mode)
                    if ostr_s[-1] == ',':
                        ostr_s = ostr_s[:-1]
                    ostr_s += ':'
                    ostr_rec(obj[key], dict_mode)
                ostr_s = ostr_s[:-1] + '} '
        elif isinstance(obj, set):
            tmp_obj = list(obj)
            ostr_s += '{'
            for obj_i in tmp_obj:
                ostr_rec(obj_i, dict_mode)
            ostr_s = ostr_s[:-1] + '},'
        else:
            ostr_s += str(obj) + ','
        return
    global ostr_s
    ostr_s = ''
    if isinstance(obj, Matrix):
        ostr_s += '\n' + str(obj)
        return ostr_s
    elif isinstance(obj, (tuple, list, dict, set)):
        ostr_rec(obj, dict_mode)
        ostr_s = ostr_s[:-1]
    else:
        ostr_s = str(obj)
    return ostr_s


def find_functions(expr):
    f_lst = []
    for f in list(expr.atoms(Function)):
        if str(f) not in GaPrinter.function_names:
            f_lst.append(f)
    f_lst += list(expr.atoms(Derivative))
    return f_lst

    #return list(expr.atoms(Function)) + list(expr.atoms(Derivative))


def coef_simplify(expr):
    fcts = find_functions(expr)
    return expr.collect(fcts)


def oprint(*args, **kwargs):
    """
    Debug printing for iterated (list/tuple/dict/set) objects. args is
    of form (title1,object1,title2,object2,...) and prints:

        title1 = object1
        title2 = object2
        ...

    If you only wish to print a title set object = None.
    """

    if 'dict_mode' in kwargs:
        dict_mode = kwargs['dict_mode']
    else:
        dict_mode = False

    if isinstance(args[0], str) or args[0] is None:
        titles = list(islice(args, None, None, 2))
        objs = tuple(islice(args, 1, None, 2))
        if len(args) > 2:
            if objs[0] is None:
                n = 0
            else:
                n = len(titles[0])
            for (title, obj) in zip(titles[1:], objs[1:]):
                if obj is None:
                    if not (dict_mode and isinstance(obj, dict)):
                        n = max(n, len(title))
        else:
            n = len(titles[0])

        for (title, obj) in zip(titles, objs):
            if obj is None:
                print title
            else:
                npad = n - len(title)
                if isinstance(obj, dict):
                    print title + ':' + ostr(obj, dict_mode)
                else:
                    print title + npad * ' ' + ' = ' + ostr(obj, dict_mode)
    else:
        for arg in args:
            print ostr(arg, dict_mode)
    return


class Eprint:

    ColorCode = {'black': '0;30', 'bright gray': '0;37',
                    'blue': '0;34', 'white': '1;37',
                    'green': '0;32', 'bright blue': '1;34',
                    'cyan': '0;36', 'bright green': '1;32',
                    'red': '0;31', 'bright cyan': '1;36',
                    'purple': '0;35', 'bright red': '1;31',
                    'yellow': '0;33', 'bright purple': '1;35',
                    'dark gray': '1;30', 'bright yellow': '1;33',
                    'normal': '0'}

    InvColorCode = dict(zip(ColorCode.values(), ColorCode.keys()))

    normal = ''
    base = ''
    fct = ''
    deriv = ''
    bold = ''

    defaults = {('win', 'base'): 'blue', ('unix', 'base'): 'dark gray',
                ('win', 'fct'): 'red', ('unix', 'fct'): 'red',
                ('win', 'deriv'): 'cyan', ('unix', 'deriv'): 'cyan'}

    def __init__(self, base=None, fct=None, deriv=None, on=True, debug=False):
        if on:
            OS = 'unix'
            if 'win' in sys.platform:
                OS = 'win'

            if base is None:
                Eprint.base = Eprint.ColorCode[Eprint.defaults[(OS, 'base')]]
            else:
                Eprint.base = Eprint.ColorCode[base]
            if fct is None:
                Eprint.fct = Eprint.ColorCode[Eprint.defaults[(OS, 'fct')]]
            else:
                Eprint.fct = Eprint.ColorCode[fct]
            if deriv is None:
                Eprint.deriv = Eprint.ColorCode[Eprint.defaults[(OS, 'deriv')]]
            else:
                Eprint.deriv = Eprint.ColorCode[deriv]
            Eprint.normal = '\033[0m'

            if debug:
                print 'Enhanced Printing is on:'
                print 'Base/Blade color is ' + Eprint.InvColorCode[Eprint.base]
                print 'Function color is ' + Eprint.InvColorCode[Eprint.fct]
                print 'Derivative color is ' + Eprint.InvColorCode[Eprint.deriv] + '\n'

            Eprint.base = '\033[' + Eprint.base + 'm'
            Eprint.fct = '\033[' + Eprint.fct + 'm'
            Eprint.deriv = '\033[' + Eprint.deriv + 'm'

    @staticmethod
    def Base(s):
        return Eprint.base + s + Eprint.normal

    @staticmethod
    def Fct(s):
        return Eprint.fct + s + Eprint.normal

    @staticmethod
    def Deriv(s):
        return Eprint.deriv + s + Eprint.normal

    @staticmethod
    def Strip(s):
        new_s = s.replace(Eprint.base, '')
        new_s = new_s.replace(Eprint.normal, '')
        return new_s


class GaPrinter(StrPrinter):

    function_names = ('acos', 'acosh', 'acot', 'acoth', 'arg', 'asin', 'asinh',
                      'atan', 'atan2', 'atanh', 'ceiling', 'conjugate', 'cos',
                      'cosh', 'cot', 'coth', 'exp', 'floor', 'im', 'log', 're',
                      'root', 'sin', 'sinh', 'sqrt', 'sign', 'tan', 'tanh', 'Abs')

    str_flg = True

    def _print_Function(self, expr):
        name = expr.func.__name__

        if expr.func.nargs is not None:
            if name in GaPrinter.function_names:
                return(expr.func.__name__ + "(%s)" % self.stringify(expr.args, ", "))

        return Eprint.Fct("%s" % (name,))

    def _print_Derivative(self, expr):
        diff_args = map(self._print, expr.args)
        xi = []
        ni = []
        for x in diff_args[1:]:
            if x in xi:
                i = xi.index(x)
                ni[i] += 1
            else:
                xi.append(x)
                ni.append(1)

        s = 'D'
        for (x, n) in zip(xi, ni):
            s += '{' + str(x) + '}'
            if n > 1:
                s += '^' + str(n)
        s += str(diff_args[0])
        return Eprint.Deriv(s)

    def _print_Matrix(self, expr):
        out_str = ostr(list(expr))
        return out_str

    def _print_Mv(self, expr):
        if expr.obj == S(0):
            return '0'
        else:
            #print 'expr.obj =',expr.obj
            return expr.Mv_str()

    def _print_Pdop(self, expr):
        return expr.Pdop_str()

    def _print_Dop(self, expr):
        return expr.Dop_str()

    def _print_Sdop(self, expr):
        return expr.Sdop_str()

    def _print_Lt(self, expr):
        return expr.Lt_str()

    def _print_Mlt(self, expr):
        return expr.Mlt_str()


Basic.__str__ = lambda self: GaPrinter().doprint(self)
Matrix.__str__ = lambda self: GaPrinter().doprint(self)
Basic.__repr_ = lambda self: GaPrinter().doprint(self)

class GaLatexPrinter(LatexPrinter):
    """
    The latex printer is turned on with the function (in ga.py) -

        Format(Fmode=True,Dmode=True,ipy=False)

    where Fmode is the function printing mode that surpresses printing arguments,
    Dmode is the derivative printing mode that does not use fractions, and
    ipy=True is the Ipython notebook mode that does not redirect the print output.

    The latex output is post processed and displayed with the function (in GAPrint.py) -

        xpdf(filename='tmplatex.tex',debug=False)

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

        'grad' -> '\bm{\nabla} '
        '*'    -> ''
        '^'    -> '\W '
        '|'    -> '\cdot '
        '>'    -> '\lfloor '
        '<'    -> '\rfloor '

    In the case of a print statement of the form -

        print title,A

    everthing in the title processing still applies except that the multivector
    formatting is one multivector per line.

    For print statements of the form -

        print title

    where no program variables are printed if title contains '#' then title
    is printed as regular latex line.  If title does not contain '#' then
    title is printed in equation mode. '%' has the same effect in title as
    in the Fmt() member function.
    """
    fmt = 1
    latex_flg = False
    latex_str = ''
    ipy = False

    preamble = \
"""
\\pagestyle{empty}
\\usepackage[latin1]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{amsbsy}
\\usepackage{tensor}
\\usepackage{listings}
\\usepackage{color}
\\usepackage{xcolor}
\\usepackage{bm}
\\definecolor{gray}{rgb}{0.95,0.95,0.95}
\\setlength{\\parindent}{0pt}
\\DeclareMathOperator{\Tr}{Tr}
\\DeclareMathOperator{\Adj}{Adj}
\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}
\\newcommand{\\lp}{\\left (}
\\newcommand{\\rp}{\\right )}
\\newcommand{\\paren}[1]{\\lp {#1} \\rp}
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
    macros = '\\newcommand{\\f}[2]{{#1}\\left ({#2}\\right )}'

    Dmode = False  # True - Print derivative contracted
    Fmode = False  # True - Print function contracted
    latex_flg = False
    ipy = False

    # translate name, supers and subs to tex keywords
    greek = set(['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                 'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu',
                 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon',
                 'phi', 'chi', 'psi', 'omega'])

    greek_translated = {'lamda': 'lambda', 'Lamda': 'Lambda'}

    other = set(['aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth',
                 'hbar', 'hslash', 'mho'])

    special_alphabet = list(reversed(sorted(list(greek) + list(other), key=len)))

    @staticmethod
    def split_super_sub(text):
        """Split a symbol name into a name, superscripts and subscripts

           The first part of the symbol name is considered to be its actual
           'name', followed by super- and subscripts. Each superscript is
           preceded with a "^" character or by "__". Each subscript is preceded
           by a "_" character.  The three return values are the actual name, a
           list with superscripts and a list with subscripts.

           >>> from sympy.printing.conventions import split_super_sub
           >>> split_super_sub('a_x^1')
           ('a', ['1'], ['x'])
           >>> split_super_sub('var_sub1__sup_sub2')
           ('var', ['sup'], ['sub1', 'sub2'])

        """

        def sub_split_super_sub(text):

            pos = 0
            name = None
            supers = []
            subs = []
            while pos < len(text):
                start = pos + 1
                if text[pos:pos + 2] == "__":
                    start += 1
                pos_hat = text.find("^", start)
                if pos_hat < 0:
                    pos_hat = len(text)
                pos_usc = text.find("_", start)
                if pos_usc < 0:
                    pos_usc = len(text)
                pos_next = min(pos_hat, pos_usc)
                part = text[pos:pos_next]
                pos = pos_next
                if name is None:
                    name = part
                elif part.startswith("^"):
                    supers.append(part[1:])
                elif part.startswith("__"):
                    supers.append(part[2:])
                elif part.startswith("_"):
                    subs.append(part[1:])
                else:
                    raise RuntimeError("This should never happen.")

            # make a little exception when a name ends with digits, i.e. treat them
            # as a subscript too.
            m = re.match('(^[a-zA-Z]+)([0-9]+)$', name)
            if m is not None:
                name, sub = m.groups()
                subs.insert(0, sub)

            return name, supers, subs

        if '*' not in text and '^' not in text:
            name, supers, subs = sub_split_super_sub(text)
            return '*', [name], [supers], [subs]

        if '*' in text:
            basis = text.split('*')
            split_flg = '*'
        if '^' in text:
            basis = text.split('^')
            split_flg = '^'

        name_lst = []
        supers_lst = []
        subs_lst = []

        for base in basis:
            name, supers, subs = sub_split_super_sub(base)
            name_lst.append(name)
            supers_lst.append(supers)
            subs_lst.append(subs)

        return split_flg, name_lst, supers_lst, subs_lst

    @staticmethod
    def redirect():
        GaLatexPrinter.latex_flg = True
        GaLatexPrinter.Basic__str__ = Basic.__str__
        GaLatexPrinter.Matrix__str__ = Matrix.__str__
        Basic.__str__ = lambda self: GaLatexPrinter().doprint(self)
        Matrix.__str__ = lambda self: GaLatexPrinter().doprint(self)
        if GaLatexPrinter.ipy:
            pass
        else:
            GaLatexPrinter.stdout = sys.stdout
            sys.stdout = StringIO.StringIO()
        return

    @staticmethod
    def restore():
        if GaLatexPrinter.latex_flg:
            if not GaLatexPrinter.ipy:
                GaLatexPrinter.latex_str += sys.stdout.getvalue()
            GaLatexPrinter.latex_flg = False
            if not GaLatexPrinter.ipy:
                sys.stdout = GaLatexPrinter.stdout
            Basic.__str__ = GaLatexPrinter.Basic__str__
            Matrix.__str__ = GaLatexPrinter.Matrix__str__
        return

    def _print_Pow(self, expr):
        base = self._print(expr.base)
        if ('_' in base or '^' in base) and 'cdot' not in base:
            mode = True
        else:
            mode = False

        # Treat x**Rational(1,n) as special case
        if expr.exp.is_Rational and abs(expr.exp.p) == 1 and expr.exp.q != 1:
            #base = self._print(expr.base)
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
                return r"{\left ( %s \right )}^{%s/%s}" % (base, p, q)
            else:
                return r"%s^{%s/%s}" % (base, p, q)

        elif expr.exp.is_Rational and expr.exp.is_negative and expr.base.is_Function:
            # Things like 1/x
            return r"\frac{%s}{%s}" % \
                (1, self._print(C.Pow(expr.base, -expr.exp)))
        else:
            if expr.base.is_Function:
                return self._print(expr.base, self._print(expr.exp))
            else:
                if expr.is_commutative and expr.exp == -1:
                    #solves issue 1030
                    #As Mul always simplify 1/x to x**-1
                    #The objective is achieved with this hack
                    #first we get the latex for -1 * expr,
                    #which is a Mul expression
                    tex = self._print(S.NegativeOne * expr).strip()
                    #the result comes with a minus and a space, so we remove
                    if tex[:1] == "-":
                        return tex[1:].strip()
                if self._needs_brackets(expr.base):
                    tex = r"\left(%s\right)^{%s}"
                else:
                    if mode:
                        tex = r"{\left ( %s \right )}^{%s}"
                    else:
                        tex = r"%s^{%s}"

                return tex % (self._print(expr.base),
                              self._print(expr.exp))

    def _print_Symbol(self, expr):

        mode_dict = {'*': '', '^': '\\wedge '}

        def str_symbol(name_str):
            (mode, name_lst, supers_lst, subs_lst) = GaLatexPrinter.split_super_sub(name_str)

            def translate(s):
                #tmp = s.lower()
                tmp = s

                parse_dict = {}
                i_sub = 1

                for glyph in GaLatexPrinter.special_alphabet:
                    if glyph in tmp:
                        parse_sym = '????' + str(i_sub)
                        parse_dict[parse_sym] = '\\' + glyph + ' '
                        tmp = tmp.replace(glyph, parse_sym)

                for parse_sym in parse_dict:
                    tmp = tmp.replace(parse_sym, parse_dict[parse_sym])

                for glyph in GaLatexPrinter.greek_translated:
                    if glyph in tmp:
                        tmp = tmp.replace(glyph, GaLatexPrinter.greek_translated[glyph])

                return tmp

                """
                if tmp in GaLatexPrinter.greek or tmp in other:
                    return "\\" + s
                if s in GaLatexPrinter.greek_translated:
                    return "\\" + GaLatexPrinter.greek_translated[s]
                else:
                    return s
                """

            s = ''

            for (name, supers, subs) in zip(name_lst, supers_lst, subs_lst):

                name = translate(name)

                if supers != []:
                    supers = map(translate, supers)

                if subs != []:
                    subs = map(translate, subs)

                # glue all items together:
                if len(supers) > 0:
                    name += "^{%s}" % " ".join(supers)
                if len(subs) > 0:
                    name += "_{%s}" % " ".join(subs)

                s += name + mode_dict[mode]

            if mode == '^':
                s = s[:-7]

            if not expr.is_commutative and mv.Mv.latex_flg:
                s = '\\boldsymbol{' + s + '}'

            return s

        if expr in self._settings['symbol_names']:
            return self._settings['symbol_names'][expr]

        name_str = expr.name

        if '.' in name_str and name_str[0] == '(' and name_str[-1] == ')':
            name_str = name_str[1:-1]
            name_lst = name_str.split('.')
            name_str = r'\left ( ' + str_symbol(name_lst[0]) + r'\cdot ' + str_symbol(name_lst[1]) + r'\right ) '
            return(name_str)

        return(str_symbol(expr.name))

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
                len(args) == 1 and not self._needs_function_brackets(expr.args[0])

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
                    name = latex(Symbol(func)) + ' '
                    if '_' in func or '^' in func:
                        name = r'{\left ( ' + name + r'\right ) }^{' + exp + '}'
                    else:
                        name += '^{' + exp + '}'
            else:
                if func in accepted_latex_functions:
                    name = r"\%s" % func
                else:
                    name = latex(Symbol(func)) + ' '
                    if exp is not None:
                        if '_' in name or '^' in name:
                            name = r'\left ( ' + name + r'\right )^{' + exp + '}'
                        else:
                            name += '^{' + exp + '}'

            if can_fold_brackets:
                if func in accepted_latex_functions:
                    # Wrap argument safely to avoid parse-time conflicts
                    # with the function name itself
                    name += r" {%s}"
                else:
                    if not GaLatexPrinter.Fmode:
                        name += r"%s"
            else:
                if func in accepted_latex_functions or not GaLatexPrinter.Fmode:
                    name += r"{\left (%s \right )}"

            if inv_trig_power_case and exp is not None:
                name += r"^{%s}" % exp

            if func in accepted_latex_functions or not GaLatexPrinter.Fmode:
                if len(args) == 1:
                    name = name % args[0]
                else:
                    name = name % ",".join(args)

            if 'det(g)' in name:
                name = name.replace('det(g)',r'\det\left ( g \right )')

            return name

    def _print_Derivative(self, expr):
        dim = len(expr.variables)
        imax = 1
        if dim == 1:
            if GaLatexPrinter.Dmode:
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

            if GaLatexPrinter.Dmode:
                tex = ''
                for x, i in multiplicity:
                    if i == 1:
                        tex += r"\partial_{%s}" % (self._print(x),)
                    else:
                        tex += r"\partial^{%i}_{%s}" % (i, self._print(x))
            else:
                for x, i in multiplicity:
                    if i == 1:
                        tex += r"\partial %s" % self._print(x)
                    else:
                        tex += r"\partial^{%s} %s" % (i, self._print(x))
                tex = r"\frac{\partial^{%s}}{%s} " % (dim, tex)

        if isinstance(expr.expr, C.AssocOp):
            s = r"%s\left(%s\right)" % (tex, self._print(expr.expr))
        else:
            s = r"%s %s" % (tex, self._print(expr.expr))
        return s

    def _print_Mv(self, expr):
        if expr.obj == S(0):
            return('0 \n')
        else:
            ostr = expr.Mv_latex_str()
            return(ostr)

    def _print_Pdop(self, expr):
        return expr.Pdop_latex_str()

    def _print_Dop(self, expr):
        return expr.Dop_latex_str()

    def _print_Sdop(self, expr):
        return expr.Sdop_latex_str()

    def _print_Lt(self, expr):
        return expr.Lt_latex_str()

    def _print_Mlt(self, expr):
        return expr.Mlt_latex_str()

    def _print_Matrix(self, expr):
        lines = []
        for line in range(expr.rows):  # horrible, should be 'rows'
            lines.append(" & ".join([self._print(i) for i in expr[line:]]))

        ncols = 0
        for i in expr[line:]:
            ncols += 1

        out_str = ' \\left [ \\begin{array}{' + ncols * 'c' + '} '
        for line in lines[:-1]:
            out_str += line + ' \\\\ '
        out_str += lines[-1] + ' \\end{array}\\right ] '
        return out_str

    @staticmethod
    def latex(expr, **settings):

        if not isinstance(expr, list):
            return GaLatexPrinter(settings).doprint(expr)
        else:
            s = '\\begin{align*}'
            for x in expr:
                s += '\n & ' + latex(x) + ' \\\\'
            s += '\n\\end{align*}'
            return s


def latex(expr, **settings):

    if not isinstance(expr, list):
        return GaLatexPrinter(settings).doprint(expr)
    else:
        s = '\\begin{align*}'
        for x in expr:
            s += '\n & ' + latex(x) + ' \\\\'
        s += '\n\\end{align*}'
        return s


def print_latex(expr, **settings):
    """Prints LaTeX representation of the given expression."""
    print latex(expr, **settings)


def Format(Fmode=True, Dmode=True, ipy=False, dop=1):
    """
    Set modes for latex printer -

        Fmode:  Suppress function arguments (True)          Use sympy latex for functions (False)
        Dmode:  Use compact form of derivatives (True)      Use sympy latex for derivatives (False)
        ipy:    Do not redirect ouput for ipython (True)    Redirect output for latex post processing (False)

    and redirects printer output so that latex compiler can capture it.
    """

    GaLatexPrinter.Dmode = Dmode
    GaLatexPrinter.Fmode = Fmode
    GaLatexPrinter.ipy = ipy
    GaLatexPrinter.dop = dop
    if ipy:
        from IPython.core.display import display, Math, Latex

    GaLatexPrinter.latex_flg = True
    GaLatexPrinter.redirect()
    if ipy:
        return Latex(ip_cmds)
    return


def xpdf(filename=None, paper=(14, 11), crop=False, png=False, prog=False, debug=False, pt='10pt'):

    """
    Post processes LaTeX output (see comments below), adds preamble and
    postscript, generates tex file, inputs file to latex, displays resulting
    pdf file.

    Arg    Value    Result
    crop   True     Use "pdfcrop" to crop output file (pdfcrop must be installed, linux only)
    png    True     Use "convert" to produce png output (imagemagick must be installed, linux only)
    """
    if GaLatexPrinter.ipy:
        GaLatexPrinter.restore()
        return

    sys_cmd = SYS_CMD[sys.platform]

    latex_str = GaLatexPrinter.latex_str + sys.stdout.getvalue()
    GaLatexPrinter.latex_str = ''
    GaLatexPrinter.restore()
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
                if r'\begin{align*}' in latex_line:
                    latex_line = r'\begin{align*}' + latex_line.replace(r'\begin{align*}', '')
                else:
                    latex_line = '\\begin{equation*} ' + latex_line + ' \\end{equation*}\n'

        else:
            latex_line = latex_line.replace(r'\left.', r'@@') # Disabiguate '.' in '\left.'
            latex_line = latex_line.replace('.', r' \cdot ')  # For components of metric tensor
            latex_line = latex_line.replace(r'@@',r'\left.')  # Restore '\left.'
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
                    lhs = lhs.replace('rgrad', r'\bar{\boldsymbol{\nabla}} ')
                    lhs = lhs.replace('grad', r'\boldsymbol{\nabla} ')
                    lhs = lhs.replace('<', r'\rfloor ')
                    lhs = lhs.replace('>', r'\lfloor ')
                latex_line = lhs + latex_line

            if r'\begin{align*}' in latex_line:  # insert lhs of align environment
                latex_line = latex_line.replace(lhs, '')
                latex_line = latex_line.replace(r'\begin{align*}', r'\begin{align*} ' + lhs)
                lhs = ''
            else:  # normal sympy equation
                latex_line = latex_line.strip()
                if len(latex_line) > 0:
                    latex_line = '\\begin{equation*} ' + latex_line + ' \\end{equation*}'
            #print 'latex_line =', latex_line

        latex_str += latex_line + '\n'

    latex_str = latex_str.replace('\n\n', '\n')

    if prog:
        prog_file = open(sys.argv[0], 'r')
        prog_str = prog_file.read()
        prog_file.close()
        prog_str = '{\\Large \\bf Program:}\\begin{lstlisting}[language=Python,showspaces=false,' + \
                   'showstringspaces=false]\n' + \
                   prog_str + '\n\\end{lstlisting}\n {\\Large \\bf Code Output:} \n'
        latex_str = prog_str + latex_str

    if debug:
        print latex_str

    if paper == 'letter':
        paper_size = \
"""
\\documentclass[@10pt@,fleqn]{report}
"""
    else:
        paper_size = \
"""
\\documentclass[@10pt@,fleqn]{report}
\\usepackage[vcentering]{geometry}
"""
        if paper == 'landscape':
            paper = [11,8.5]
        paper_size += '\\geometry{papersize={' + str(paper[0]) + \
                      'in,' + str(paper[1]) + 'in},total={' + str(paper[0] - 1) + \
                      'in,' + str(paper[1] - 1) + 'in}}\n'

    paper_size = paper_size.replace('@10pt@',pt)
    latex_str = paper_size + GaLatexPrinter.preamble + latex_str + GaLatexPrinter.postscript

    if filename is None:
        pyfilename = sys.argv[0]
        rootfilename = pyfilename.replace('.py', '')
        filename = rootfilename + '.tex'

    print 'latex file =', filename

    latex_file = open(filename, 'w')
    latex_file.write(latex_str)
    latex_file.close()

    latex_str = None

    pdflatex = find_executable('pdflatex')

    print 'pdflatex path =', pdflatex

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
        print print_cmd

        os.system(print_cmd)
        raw_input('!!!!Return to continue!!!!\n')

        if debug:
            os.system(sys_cmd['rm'] + ' ' + filename[:-4] + '.aux ' + filename[:-4] + '.log')
        else:
            os.system(sys_cmd['rm'] + ' ' + filename[:-4] + '.aux ' + filename[:-4] + '.log ' + filename[:-4] + '.tex')
        if crop:
            os.system('pdfcrop ' + filename[:-4] + '.pdf')
            os.remove(filename[:-4] + '.pdf')
            os.rename(filename[:-4] + '-crop.pdf', filename[:-4] + '.pdf')
        if png:
            os.system('Pdf2Png ' + filename[:-4])
    return


def LatexFormat(Fmode=True, Dmode=True, ipy=False):
    GaLatexPrinter.Dmode = Dmode
    GaLatexPrinter.Fmode = Fmode
    GaLatexPrinter.ipy = ipy
    GaLatexPrinter.redirect(ipy)
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
    if GaLatexPrinter.latex_flg:
        #print '#Code for '+fct_name
        print '##\\begin{lstlisting}[language=Python,showspaces=false,' + \
              'showstringspaces=false,backgroundcolor=\color{gray},frame=single]'
        print tmp_str
        print '##\\end{lstlisting}'
        print '#Code Output:'
    else:
        print '\n' + 80 * '*'
        #print '\nCode for '+fct_name
        print tmp_str
        print 'Code output:\n'
    return


import re as regrep

op_cntrct = regrep.compile(r'(([A-Za-z0-9\_\#]+)(\||<|>)([A-Za-z0-9\_\#]+))')
op_wedge = regrep.compile(r'(([A-Za-z0-9\_\#]+)[\^]{1}([A-Za-z0-9\_\#]+)([\^]{1}([A-Za-z0-9\_\#]+))*)')
ops = r'[\^\|\<\>]+'
ops_search = regrep.compile(r'(\^|\||<|>)+')
parse_paren_calls = 0
global_dict = {}
op_dict = {}
op_lst = []

OPS = {'<>|': r'(([A-Za-z0-9\_\#]+)(\||<|>)([A-Za-z0-9\_\#]+))',
       '^': r'(([A-Za-z0-9\_\#]+)[\^]{1}([A-Za-z0-9\_\#]+)([\^]{1}([A-Za-z0-9\_\#]+))*)',
       '*': r'(([A-Za-z0-9\_\#]+)[\*]{1}([A-Za-z0-9\_\#]+)([\*]{1}([A-Za-z0-9\_\#]+))*)'}


def def_prec(gd, op_ord='<>|,^,*'):  # Default is Doran and Lasenby convention
    global global_dict, op_dict, op_lst
    global_dict = gd
    op_lst = op_ord.split(',')
    op_dict = {}
    for op in op_lst:
        op_dict[op] = regrep.compile(OPS[op])
    return


def contains_interval(interval1, interval2):  # interval1 inside interval2
    if interval1[0] > interval2[0] and interval1[1] < interval2[1]:
        return True
    else:
        return False


def parse_paren(line):
    global parse_paren_calls
    parse_paren_calls += 1

    if ('(' not in line) or (')' not in line):
        return [[[line]]]
    level = 0
    max_level = 0
    ich = 0
    paren_lst = []
    for ch in line:
        if ch == '(':
            level += 1
            paren_lst.append([level, ich])
        if ch == ')':
            if level < 1:
                raise ValueError('Mismathed Parenthesis in: ' + line + '\n')
            paren_lst.reverse()
            iparen = 0
            for elem in paren_lst:
                if elem[0] == level:
                    paren_lst[iparen].append(ich)
                    break
                iparen += 1
            paren_lst.reverse()
            level -= 1
        max_level = max(max_level, level)
        ich += 1
    if level != 0:
        raise ValueError('Mismatched Parenthesis in: ' + line + '\n')
    if max_level > 0:
        level_lst = []
        for x in range(max_level + 1):
            level_lst.append([])
        for group in paren_lst:
            level_lst[group[0]].append(group[1:])
        ilevel = max_level
        while ilevel > 1:
            level = level_lst[ilevel]
            level_down = level_lst[ilevel - 1]
            igroup = 0
            for group in level:
                igroup_down = 0
                for group_down in level_down:
                    if contains_interval(group, group_down):
                        level_lst[ilevel][igroup].append(igroup_down)
                    igroup_down += 1
                igroup += 1
            ilevel -= 1
        ilevel = 1
        for level in level_lst[1:]:
            igroup = 0
            for group in level:
                token = '#' + str(parse_paren_calls) + '_' + str(ilevel) + '_' + str(igroup) + '#'
                level_lst[ilevel][igroup].append(line[group[0]:group[1] + 1])
                level_lst[ilevel][igroup].append(token)
                igroup += 1
            ilevel += 1
        ilevel = 1
        for level in level_lst[1:]:
            igroup = 0
            for group in level:
                group.append(group[-2])
                level_lst[ilevel][igroup] = group
                igroup += 1
            ilevel += 1
        ilevel = max_level
        while ilevel > 1:
            igroup = 0
            for group in level_lst[ilevel]:
                group_down = level_lst[ilevel - 1][group[2]]
                replace_text = group_down[-1].replace(group[-3], group[-2])
                level_lst[ilevel - 1][group[2]][-1] = replace_text
                igroup += 1
            ilevel -= 1
        for group in level_lst[1]:
            line = line.replace(group[2], group[3])
        ilevel = 1
        level_lst[0] = [[line]]
    return level_lst


def unparse_paren(level_lst):
    line = level_lst[0][0][0]
    for level in level_lst[1:]:
        for group in level:
            new_string = group[-1]
            if new_string[:2] == '((' and new_string[-2:] == '))':
                new_string = new_string[1:-1]
            line = line.replace(group[-2], new_string)
    return line


def sub_paren(s):
    string = s.group(0)
    return '(%s)' % string


def add_paren(line, re_exprs):
    paren_flg = False
    if (line[0] == '(') and (line[-1] == ')'):
        paren_flg = True
        line = line[1:-1]
    if ('(' in line) or (')' in line):
        line_levels = parse_paren(line)
        ilevel = 0
        for level in line_levels:
            igroup = 0
            for group in level:
                group[-1] = regrep.sub(re_exprs, sub_paren, group[-1])
                line_levels[ilevel][igroup] = group
                igroup += 1
            ilevel += 1
        line = unparse_paren(line_levels)
    else:
        line = regrep.sub(re_exprs, sub_paren, line)
    if paren_flg:
        line = '(' + line + ')'
    return line


def parse_line(line):
    global op_lst, op_dict
    line = line.replace(' ', '')
    level_lst = parse_paren(line)
    ilevel = 0
    for level in level_lst:
        igroup = 0
        for group in level:
            string = group[-1]
            for op in op_lst:
                string = add_paren(string, op_dict[op])
            level_lst[ilevel][igroup][-1] = string
            igroup += 1
        ilevel += 1
    line = unparse_paren(level_lst)
    return line


def GAeval(s, pstr=False):
    """
    GAeval converts a string to a multivector expression where the
    user can control the precedence of the of the multivector operators so
    that one does not need to put parenthesis around every multivector
    operation.  The default precedence used (high to low) is <,>, and | have
    an have the highest precedence, then comes ^, and finally *.  The
    default precedence can be changed with the def_prec function.
    """

    seval = parse_line(s)
    if pstr:
        print s
        print seval
    return eval(seval, global_dict)

if __name__ == "__main__":
    pass
