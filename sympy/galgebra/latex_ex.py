#latex_ex.py

import sys
#if sys.version.find('Stackless') >= 0:
#    sys.path.append('/usr/lib/python2.5/site-packages')

import os,types,StringIO

from sympy.core import S, C, Basic, Symbol
from sympy.printing.printer import Printer
from sympy.simplify import fraction
import re as regrep

import sympy.galgebra.GA
#import sympy.galgebra.OGA
import numpy

def debug(txt):
    sys.stderr.write(txt+'\n')
    return

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

def debug(tstr):
    return

def len_cmp(str1,str2):
    return(len(str2)-len(str1))

def process_equals(xstr):
    eq1 = xstr.find('=')
    eq2 = xstr.rfind('=')
    if eq1 == eq2:
        return(xstr)
    xstr = xstr[:eq1]+xstr[eq2:]
    return(xstr)

class LatexPrinter(Printer):
    """
    A printer class which converts an expression into its LaTeX equivalent. This
    class extends the LatexPrinter class currently in sympy in the following ways:

        1. Variable and function names can now encode multiple greek symbols,
           number, greek, and roman super and sub scipts and accents plus bold
           math in an alphanumeric ascii string consisting of [A-Za-z0-9_]
           symbols
            1 - Accents and bold math are implemented in reverse notation. For
                example if you wished the LaTeX output to be '\bm{\hat{\sigma}}'
                you would give the variable the name sigmahatbm.
            2 - Subscripts are denoted by a single underscore and superscripts
                by a double underscore so that A_{\rho\beta}^{25} would be
                input as A_rhobeta__25.
        2. Some standard function names have been improved such as asin is now
           denoted by Sin^{-1} and log by ln.
        3. Several LaTeX formats for multivectors are available:
            1 - Print multivector on one line
            2 - Print each grade of multivector on one line
            3 - Print each base of multivector on one line
        4. A LaTeX output for numpy arrays containing sympy expressions is
           implemented for upto a three dimensional array.
        5. LaTeX formatting for raw LaTeX, eqnarray, and array is available
           in simple output strings.
            1 - The delimeter for raw LaTeX imput is '%'.  The raw input starts
                on the line where '%' is first encountered and continues untill
                the next line where '%' is encountered. It does not matter where
                '%' is in the line.
            2 - The delimeter for eqnarray input is '@'. The rules are the same
                as for raw input except that '=' in the first line is replaced
                be '&=&' and '\begin{eqnarray*}' is added before the first line
                and '\end{eqnarray*}' to after the last line in the group of
                lines.
            3 - The delimeter for array input is '#'. The rules are the same
                as for raw input except that '\begin{equation*}' is added before
                the first line and '\end{equation*}' to after the last line in
                the group of lines.
        6. Additonal formats for partial derivatives:
            0 - Same as sympy latex module
            1 - Use subscript notation with partial symbol to indicate which
                variable the differentiation is with respect to.  Symbol is of
                form \partial_{differentiation variable}
    """

    #printmethod ='_latex_ex_'
    sym_fmt = 0
    fct_fmt = 0
    pdiff_fmt = 0
    mv_fmt = 0
    str_fmt = 1
    LaTeX_flg = False

    mode = ('_','^')

    fmt_dict = {'sym':0,'fct':0,'pdiff':0,'mv':0,'str':1}

    fct_dict = {'sin':'sin','cos':'cos','tan':'tan','cot':'cot',\
                'asin':'Sin^{-1}','acos':'Cos^{-1}',\
                'atan':'Tan^{-1}','acot':'Cot^{-1}',\
                'sinh':'sinh','cosh':'cosh','tanh':'tanh','coth':'coth',\
                'asinh':'Sinh^{-1}','acosh':'Cosh^{-1}',
                'atanh':'Tanh^{-1}','acoth':'Coth^{-1}',\
                'sqrt':'sqrt','exp':'exp','log':'ln'}

    fct_dict_keys = fct_dict.keys()

    greek_keys = sorted(('alpha','beta','gamma','delta','varepsilon','epsilon','zeta',\
                         'vartheta','theta','iota','kappa','lambda','mu','nu','xi',\
                         'varpi','pi','rho','varrho','varsigma','sigma','tau','upsilon',\
                         'varphi','phi','chi','psi','omega','Gamma','Delta','Theta',\
                         'Lambda','Xi','Pi','Sigma','Upsilon','Phi','Psi','Omega','partial',\
                         'nabla','eta'),len_cmp)

    accent_keys = sorted(('hat','check','dot','breve','acute','ddot','grave','tilde',\
                          'mathring','bar','vec','bm','prm','abs'),len_cmp)

    greek_cnt = 0
    greek_dict = {}
    accent_cnt = 0
    accent_dict = {}

    preamble = '\\documentclass[10pt,letter,fleqn]{report}\n'+\
               '\\pagestyle{empty}\n'+\
               '\\usepackage[latin1]{inputenc}\n'+\
               '\\usepackage[dvips,landscape,top=1cm,nohead,nofoot]{geometry}\n'+\
               '\\usepackage{amsmath}\n'+\
               '\\usepackage{bm}\n'+\
               '\\usepackage{amsfonts}\n'+\
               '\\usepackage{amssymb}\n'+\
               '\\setlength{\\parindent}{0pt}\n'+\
               '\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}\n'+\
               '\\newcommand{\\lp}{\\left (}\n'+\
               '\\newcommand{\\rp}{\\right )}\n'+\
               '\\newcommand{\\half}{\\frac{1}{2}}\n'+\
               '\\newcommand{\\llt}{\\left <}\n'+\
               '\\newcommand{\\rgt}{\\right >}\n'+\
               '\\newcommand{\\abs}[1]{\\left |{#1}\\right | }\n'+\
               '\\newcommand{\\pdiff}[2]{\\bfrac{\\partial {#1}}{\\partial {#2}}}\n'+\
               '\\newcommand{\\lbrc}{\\left \\{}\n'+\
               '\\newcommand{\\rbrc}{\\right \\}}\n'+\
               '\\newcommand{\\W}{\\wedge}\n'+\
               "\\newcommand{\\prm}[1]{{#1}'}\n"+\
               '\\newcommand{\\ddt}[1]{\\bfrac{d{#1}}{dt}}\n'+\
               '\\newcommand{\\R}{\\dagger}\n'+\
               '\\begin{document}\n'
    postscript = '\\end{document}\n'

    @staticmethod
    def latex_bases():
        """
        Generate LaTeX strings for multivector bases
        """
        if type(sympy.galgebra.GA.MV.basislabel_lst) == types.IntType:
            sys.stderr.write('MV.setup() must be executed before LatexPrinter.format()!\n')
            sys.exit(1)
        LatexPrinter.latexbasis_lst = [['']]
        for grades in sympy.galgebra.GA.MV.basislabel_lst[1:]:
            grades_lst = []
            for grade in grades:
                grade_lst = []
                for base in grade:
                    latex_base = LatexPrinter.extended_symbol(base)
                    grade_lst.append(latex_base)
                grades_lst.append(grade_lst)
            LatexPrinter.latexbasis_lst.append(grades_lst)
        return

    @staticmethod
    def build_base(igrade,iblade,bld_flg):
        if igrade == 0:
            return('')
        base_lst = LatexPrinter.latexbasis_lst[igrade][iblade]
        if len(base_lst) == 1:
            return(base_lst[0])
        base_str = ''
        for base in base_lst[:-1]:
            if bld_flg:
                base_str += base+'\\W '
            else:
                base_str += base
        base_str += base_lst[-1]
        return(base_str)

    @staticmethod
    def format(sym=0,fct=0,pdiff=0,mv=0):
        LatexPrinter.LaTeX_flg = True
        LatexPrinter.fmt_dict['sym']   = sym
        LatexPrinter.fmt_dict['fct']   = fct
        LatexPrinter.fmt_dict['pdiff'] = pdiff
        LatexPrinter.fmt_dict['mv']    = mv
        LatexPrinter.fmt_dict['str'] = 1
        if sympy.galgebra.GA.MV.is_setup:
            LatexPrinter.latex_bases()
        LatexPrinter.redirect()
        return

    @staticmethod
    def str_basic(in_str):
        if not LatexPrinter.LaTeX_flg:
            return(str(in_str))
        Basic.__str__ = LatexPrinter.Basic__str__
        out_str = str(in_str)
        Basic.__str__ = LaTeX
        return(out_str)

    @staticmethod
    def redirect():
        LatexPrinter.Basic__str__ = Basic.__str__
        LatexPrinter.MV__str__    = sympy.galgebra.GA.MV.__str__
        LatexPrinter.stdout = sys.stdout
        sys.stdout = StringIO.StringIO()
        Basic.__str__ = LaTeX
        sympy.galgebra.GA.MV.__str__ = LaTeX
        return

    @staticmethod
    def restore():
        LatexPrinter_stdout = sys.stdout
        LatexPrinter_Basic__str__ = Basic.__str__
        LatexPrinter_MV__str__ = sympy.galgebra.GA.MV.__str__

        sys.stdout = LatexPrinter.stdout
        Basic.__str__ = LatexPrinter.Basic__str__
        sympy.galgebra.GA.MV.__str__ = LatexPrinter.MV__str__

        LatexPrinter.stdout = LatexPrinter_stdout
        LatexPrinter.Basic__str__ = LatexPrinter_Basic__str__
        LatexPrinter.MV__str__ = LatexPrinter_MV__str__
        return

    @staticmethod
    def format_str(fmt='0 0 0 0'):
        fmt_lst = fmt.split()
        if '=' not in fmt:
            LatexPrinter.fmt_dict['sym']   = int(fmt_lst[0])
            LatexPrinter.fmt_dict['fct']   = int(fmt_lst[1])
            LatexPrinter.fmt_dict['pdiff'] = int(fmt_lst[2])
            LatexPrinter.fmt_dict['mv']    = int(fmt_lst[3])
        else:
            for fmt in fmt_lst:
                x = fmt.split('=')
                LatexPrinter.fmt_dict[x[0]] = int(x[1])

        if LatexPrinter.LaTeX_flg == False:
            if sympy.galgebra.GA.MV.is_setup:
                LatexPrinter.latex_bases()
            LatexPrinter.redirect()
            LatexPrinter.LaTeX_flg = True
        return

    @staticmethod
    def append_body(xstr):
        if LatexPrinter.body_flg:
            LatexPrinter.body += xstr
            return('')
        else:
            return(xstr[:-1])

    @staticmethod
    def tokenize_greek(name_str):
        for sym in LatexPrinter.greek_keys:
            isym = name_str.find(sym)
            if isym > -1:
                keystr = '@'+str(LatexPrinter.greek_cnt)
                LatexPrinter.greek_cnt += 1
                LatexPrinter.greek_dict[keystr] = sym
                name_str = name_str.replace(sym,keystr)
        return(name_str)

    @staticmethod
    def tokenize_accents(name_str):
        for sym in LatexPrinter.accent_keys:
            if name_str.find(sym) > -1:
                keystr = '#'+str(LatexPrinter.accent_cnt)+'#'
                LatexPrinter.accent_cnt += 1
                LatexPrinter.accent_dict[keystr] = '\\'+sym
                name_str = name_str.replace(sym,keystr)
        return(name_str)

    @staticmethod
    def replace_greek_tokens(name_str):
        if name_str.find('@') == -1:
            return(name_str)
        for token in LatexPrinter.greek_dict.keys():
            name_str = name_str.replace(token,'{\\'+LatexPrinter.greek_dict[token]+'}')
        LatexPrinter.greek_cnt  = 0
        LatexPrinter.greek_dict = {}
        return(name_str)

    @staticmethod
    def replace_accent_tokens(name_str):
        tmp_lst = name_str.split('#')
        name_str = tmp_lst[0]
        if len(tmp_lst) == 1:
            return(name_str)
        for x in tmp_lst[1:]:
            if x != '':
                name_str = '{}'+LatexPrinter.accent_dict['#'+x+'#']+'{'+name_str+'}'
        LatexPrinter.accent_cnt  = 0
        LatexPrinter.accent_dict = {}
        return(name_str)

    @staticmethod
    def extended_symbol(name_str):
        name_str = LatexPrinter.tokenize_greek(name_str)
        tmp_lst = name_str.split('_')
        subsup_str = ''
        sym_str = tmp_lst[0]
        sym_str = LatexPrinter.tokenize_accents(sym_str)
        sym_str = LatexPrinter.replace_accent_tokens(sym_str)
        if len(tmp_lst) > 1:
            imode = 0
            for x in tmp_lst[1:]:
                if x == '':
                    imode = (imode+1)%2
                else:
                    subsup_str += LatexPrinter.mode[imode]+'{'+x+'}'
                    #subsup_str += LatexPrinter.mode[imode]+x+' '
                    imode = (imode+1)%2
        name_str = sym_str+subsup_str
        name_str = LatexPrinter.replace_greek_tokens(name_str)
        return(name_str)

    def coefficient(self,coef,first_flg):
        if isinstance(coef, C.AssocOp) and isinstance(-coef, C.AssocOp):
                coef_str =  r"\lp %s\rp " % self._print(coef)
        else:
            coef_str = self._print(coef)
        if first_flg:
            first_flg = False
            if coef_str[0] == '+':
                coef_str = coef_str[1:]
        else:
            if coef_str[0] != '-':
                if coef_str[0] != '+':
                    coef_str = '+'+coef_str
        if coef_str in ('1','+1','-1'):
            if coef_str == '1':
                coef_str = ''
            else:
                coef_str = coef_str[0]
        return(coef_str,first_flg)

    def __init__(self,inline=True):
        Printer.__init__(self)
        self._inline = inline

    def doprint(self, expr):
        tex = Printer.doprint(self, expr)
        xstr = ''

        if self._inline:
            if LatexPrinter.fmt_dict['fct'] == 1:
                xstr = r"%s" % tex
            else:
                xstr = r"$%s$" % tex
        else:
            xstr = r"\begin{equation*}%s\end{equation*}" % tex
        return(xstr)

    def _needs_brackets(self, expr):
        return not ((expr.is_Integer and expr.is_nonnegative) or expr.is_Atom)

    def _do_exponent(self, expr, exp):
        if exp is not None:
            return r"\left(%s\right)^{%s}" % (expr, exp)
        else:
            return expr

    def _print_Add(self, expr):
        tex = str(self._print(expr.args[0]))

        for term in expr.args[1:]:
            coeff = term.as_coeff_terms()[0]

            if coeff.is_negative:
                tex += r" %s" % self._print(term)
            else:
                tex += r" + %s" % self._print(term)

        return tex

    def _print_Mul(self, expr):
        coeff, terms = expr.as_coeff_terms()

        if not coeff.is_negative:
            tex = ""
        else:
            coeff = -coeff
            tex = "- "

        numer, denom = fraction(C.Mul(*terms))

        def convert(terms):
            product = []

            if not terms.is_Mul:
                return str(self._print(terms))
            else:
                for term in terms.args:
                    pretty = self._print(term)

                    if term.is_Add:
                        product.append(r"\left(%s\right)" % pretty)
                    else:
                        product.append(str(pretty))

                return r" ".join(product)

        if denom is S.One:
            if coeff is not S.One:
                tex += str(self._print(coeff)) + " "

            if numer.is_Add:
                tex += r"\left(%s\right)" % convert(numer)
            else:
                tex += r"%s" % convert(numer)
        else:
            if numer is S.One:
                if coeff.is_Integer:
                    numer *= coeff.p
                elif coeff.is_Rational:
                    if coeff.p != 1:
                        numer *= coeff.p

                    denom *= coeff.q
                elif coeff is not S.One:
                    tex += str(self._print(coeff)) + " "
            else:
                if coeff.is_Rational and coeff.p == 1:
                    denom *= coeff.q
                elif coeff is not S.One:
                    tex += str(self._print(coeff)) + " "

            tex += r"\frac{%s}{%s}" % \
                (convert(numer), convert(denom))

        return tex

    def _print_Pow(self, expr):
        if expr.exp.is_Rational and expr.exp.q == 2:
            base, exp = self._print(expr.base), abs(expr.exp.p)

            if exp == 1:
                tex = r"\sqrt{%s}" % base
            else:
                tex = r"\sqrt[%s]{%s}" % (exp, base)

            if expr.exp.is_negative:
                return r"\frac{1}{%s}" % tex
            else:
                return tex
        else:
            if expr.base.is_Function:
                return self._print(expr.base, self._print(expr.exp))
            else:
                if expr.exp == S.NegativeOne:
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
                    tex = r"{%s}^{%s}"

                return tex % (self._print(expr.base),
                              self._print(expr.exp))

    def _print_Derivative(self, expr):
        dim = len(expr.symbols)

        if dim == 1:
            if LatexPrinter.fmt_dict['pdiff'] == 1:
                tex = r'\partial_{%s}' % self._print(expr.symbols[0])
            else:
                tex = r"\frac{\partial}{\partial %s}" % self._print(expr.symbols[0])
        else:
            multiplicity, i, tex = [], 1, ""
            current = expr.symbols[0]
            for symbol in expr.symbols[1:]:
                if symbol == current:
                    i = i + 1
                else:
                    multiplicity.append((current, i))
                    current, i = symbol, 1
            else:
                multiplicity.append((current, i))

            if LatexPrinter.fmt_dict['pdiff'] == 1:
                for x, i in multiplicity:
                    if i == 1:
                        tex += r"\partial_{%s}" % self._print(x)
                    else:
                        tex += r"\partial^{%s}_{%s}" % (i, self._print(x))
            else:
                for x, i in multiplicity:
                    if i == 1:
                        tex += r"\partial %s" % self._print(x)
                    else:
                        tex += r"\partial^{%s} %s" % (i, self._print(x))

                tex = r"\frac{\partial^{%s}}{%s} " % (dim, tex)

        if isinstance(expr.expr, C.AssocOp):
            return r"%s\left(%s\right)" % (tex, self._print(expr.expr))
        else:
            return r"%s %s" % (tex, self._print(expr.expr))

    def _print_Integral(self, expr):
        tex, symbols = "", []

        for symbol, limits in reversed(expr.limits):
            tex += r"\int"

            if limits is not None:
                if not self._inline:
                    tex += r"\limits"

                tex += "_{%s}^{%s}" % (self._print(limits[0]),
                                       self._print(limits[1]))

            symbols.insert(0, "d%s" % self._print(symbol))

        return r"%s %s\,%s" % (tex,
            str(self._print(expr.function)), " ".join(symbols))

    def _print_Limit(self, expr):
        tex = r"\lim_{%s \to %s}" % (self._print(expr.var),
                                     self._print(expr.varlim))

        if isinstance(expr.expr, C.AssocOp):
            return r"%s\left(%s\right)" % (tex, self._print(expr.expr))
        else:
            return r"%s %s" % (tex, self._print(expr.expr))

    def _print_Function(self, expr, exp=None):
        func = expr.func.__name__

        if hasattr(self, '_print_' + func):
            return getattr(self, '_print_' + func)(expr, exp)
        else:
            args = [ str(self._print(arg)) for arg in expr.args ]

            if LatexPrinter.fmt_dict['fct'] == 1:
                if func in LatexPrinter.fct_dict_keys:
                    if exp is not None:
                        name = r"\operatorname{%s}^{%s}" % (LatexPrinter.fct_dict[func], exp)
                    else:
                        name = r"\operatorname{%s}" % LatexPrinter.fct_dict[func]
                    name += r"\left(%s\right)" % ",".join(args)
                    return name
                else:
                    func = self.print_Symbol_name(func)
                    if exp is not None:
                        name = r"{%s}^{%s}" % (func, exp)
                    else:
                        name = r"{%s}" % func
                    return name
            else:
                if exp is not None:
                    name = r"\operatorname{%s}^{%s}" % (func, exp)
                else:
                    name = r"\operatorname{%s}" % func
                return name + r"\left(%s\right)" % ",".join(args)

    def _print_floor(self, expr, exp=None):
        tex = r"\lfloor{%s}\rfloor" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_ceiling(self, expr, exp=None):
        tex = r"\lceil{%s}\rceil" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_abs(self, expr, exp=None):
        tex = r"\lvert{%s}\rvert" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_re(self, expr, exp=None):
        if self._needs_brackets(expr.args[0]):
            tex = r"\Re\left(%s\right)" % self._print(expr.args[0])
        else:
            tex = r"\Re{%s}" % self._print(expr.args[0])

        return self._do_exponent(tex, exp)

    def _print_im(self, expr, exp=None):
        if self._needs_brackets(expr.args[0]):
            tex = r"\Im\left(%s\right)" % self._print(expr.args[0])
        else:
            tex = r"\Im{%s}" % self._print(expr.args[0])

        return self._do_exponent(tex, exp)

    def _print_conjugate(self, expr, exp=None):
        tex = r"\overline{%s}" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_exp(self, expr, exp=None):
        tex = r"{e}^{%s}" % self._print(expr.args[0])
        return self._do_exponent(tex, exp)

    def _print_gamma(self, expr, exp=None):
        tex = r"\left(%s\right)" % self._print(expr.args[0])

        if exp is not None:
            return r"\operatorname{\Gamma}^{%s}%s" % (exp, tex)
        else:
            return r"\operatorname{\Gamma}%s" % tex

    def _print_Factorial(self, expr, exp=None):
        x = expr.args[0]
        if self._needs_brackets(x):
            tex = r"\left(%s\right)!" % self._print(x)
        else:
            tex = self._print(x) + "!"

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_Binomial(self, expr, exp=None):
        tex = r"{{%s}\choose{%s}}" % (self._print(expr[0]),
                                      self._print(expr[1]))

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_RisingFactorial(self, expr, exp=None):
        tex = r"{\left(%s\right)}^{\left(%s\right)}" % \
            (self._print(expr[0]), self._print(expr[1]))

        return self._do_exponent(tex, exp)

    def _print_FallingFactorial(self, expr, exp=None):
        tex = r"{\left(%s\right)}_{\left(%s\right)}" % \
            (self._print(expr[0]), self._print(expr[1]))

        return self._do_exponent(tex, exp)

    def _print_Rational(self, expr):
        if expr.q != 1:
            sign = ""
            p = expr.p
            if expr.p < 0:
                sign = "- "
                p = -p
            return r"%s\frac{%d}{%d}" % (sign, p, expr.q)
        else:
            return self._print(expr.p)

    def _print_Infinity(self, expr):
        return r"\infty"

    def _print_NegativeInfinity(self, expr):
        return r"-\infty"

    def _print_ComplexInfinity(self, expr):
        return r"\tilde{\infty}"

    def _print_ImaginaryUnit(self, expr):
        return r"\mathbf{\imath}"

    def _print_NaN(self, expr):
        return r"\bot"

    def _print_Pi(self, expr):
        return r"\pi"

    def _print_Exp1(self, expr):
        return r"e"

    def _print_EulerGamma(self, expr):
        return r"\gamma"

    def _print_Order(self, expr):
        return r"\operatorname{\mathcal{O}}\left(%s\right)" % \
            self._print(expr.args[0])

    @staticmethod
    def print_Symbol_name(name_str):
        if len(name_str) == 1:
            return (name_str)
        if LatexPrinter.fmt_dict['sym'] == 1:
            return LatexPrinter.extended_symbol(name_str)
        else:
            return(name_str)

        #convert trailing digits to subscript
            m = regrep.match('(^[a-zA-Z]+)([0-9]+)$',name_str)
            if m is not None:
                name, sub=m.groups()
                tex=self._print_Symbol(Symbol(name))
                tex="%s_{%s}" %(tex, sub)
                return tex

            # insert braces to expresions containing '_' or '^'
            m = regrep.match('(^[a-zA-Z0-9]+)([_\^]{1})([a-zA-Z0-9]+)$',name_str)
            if m is not None:
                name, sep, rest=m.groups()
                tex=self._print_Symbol(Symbol(name))
                tex="%s%s{%s}" %(tex, sep, rest)
                return tex

            greek = set([ 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                          'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu',
                          'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon',
                          'phi', 'chi', 'psi', 'omega' ])

            other = set( ['aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth',
                          'hbar', 'hslash', 'mho' ])

            if name_str.lower() in greek:
                return "\\" + name_str
            elif name_str in other:
                return "\\" + name_str
            else:
                return name_str

    def _print_Symbol(self, expr):
        return LatexPrinter.print_Symbol_name(expr.name)

    def _print_str(self,expr):
        if LatexPrinter.fmt_dict['str'] > 0:
            expr = expr.replace('^','{\\wedge}')
            expr = expr.replace('|','{\\cdot}')
            expr = expr.replace('__','^')
        return(expr)

    def _print_ndarray(self,expr):
        shape = numpy.shape(expr)
        ndim = len(shape)
        expr_str = ''

        if ndim == 1:
            expr_str += '#\\left [ \\begin{array}{'+shape[0]*'c'+'}  \n'
            for col in expr:
                expr_str += self._print(col)+' & '
            expr_str = expr_str[:-2]+'\n\\end{array}\\right ]#\n'
            return(expr_str)

        if ndim == 2:
            expr_str += '#\\left [ \\begin{array}{'+shape[1]*'c'+'}  \n'
            for row in expr[:-1]:
                for xij in row[:-1]:
                    expr_str += self._print(xij) + ' & '
                expr_str += self._print(row[-1]) + ' \\\\ \n'
            for xij in expr[-1][:-1]:
                expr_str += self._print(xij) + ' & '
            expr_str += self._print(expr[-1][-1]) + '\n \\end{array} \\right ] #\n'
            return(expr_str)

        if ndim == 3:
            expr_str = '#\\left \\{ \\begin{array}{'+shape[0]*'c'+'} \n'
            for x in expr[:-1]:
                xstr = self._print(x).replace('#','')
                expr_str += xstr + ' , & '
            xstr = self._print(expr[-1]).replace('#','')
            expr_str += xstr+'\n\\end{array} \\right \\}#\n'
        return(expr_str)

    def _print_MV(self,expr):
        igrade = 0
        MV_str = ''
        line_lst = []
        first_flg = True
        for grade in expr.mv:
            if type(grade) != types.IntType:
                if type(grade) != types.IntType:
                    ibase = 0
                    for base in grade:
                        if base != 0:
                            tmp = Symbol('XYZW')
                            base_str = str(base*tmp)
                            if base_str[0] != '-':
                                base_str = '+'+base_str
                            base_str = base_str.replace('- ','-')
                            if base_str[1:5] == 'XYZW':
                                base_str = base_str.replace('XYZW','')
                            else:
                                base_str = base_str.replace('XYZW','1')
                            MV_str += base_str+\
                                      LatexPrinter.build_base(igrade,ibase,expr.bladeflg)
                            if LatexPrinter.fmt_dict['mv'] == 3:
                                line_lst.append(MV_str)
                                MV_str = ''
                        ibase += 1
                if LatexPrinter.fmt_dict['mv'] == 2:
                    if MV_str != '':
                        line_lst.append(MV_str)
                        MV_str = ''
            igrade += 1
        n_lines = len(line_lst)
        if MV_str == '':
            if n_lines > 0 and line_lst[0][0] == '+':
                line_lst[0] = line_lst[0][1:]
        else:
            if MV_str[0] == '+':
                MV_str = MV_str[1:]
        if n_lines == 1:
            MV_str = line_lst[0]
            n_lines = 0
        if LatexPrinter.fmt_dict['mv'] >= 2:
            MV_str = '@'+line_lst[0]+' \\\\ \n'
            for line in line_lst[1:-1]:
                MV_str += '& '+line+' \\\\ \n'
            MV_str += '& '+line_lst[-1]+'@\n'
        if MV_str == '':
            MV_str = '0'
        if expr.name != '':
            MV_str = LatexPrinter.extended_symbol(expr.name)+' = '+MV_str
        return(MV_str)

    def _print_OMV(self,expr):
        igrade = 0
        MV_str = ''
        line_lst = []
        first_flg = True
        for grade in expr.mv:
            if type(grade) != None:
                if type(grade) != None:
                    ibase = 0
                    for base in grade:
                        if base != 0:
                            tmp = Symbol('XYZW')
                            base_str = str(base*tmp)
                            if base_str[0] != '-':
                                base_str = '+'+base_str
                            base_str = base_str.replace('- ','-')
                            if base_str[1:5] == 'XYZW':
                                base_str = base_str.replace('XYZW','')
                            else:
                                base_str = base_str.replace('XYZW','1')
                            MV_str += base_str+\
                                      LatexPrinter.build_base(igrade,ibase,expr.bladeflg)
                            if LatexPrinter.fmt_dict['mv'] == 3:
                                line_lst.append(MV_str)
                                MV_str = ''
                        ibase += 1
                if LatexPrinter.fmt_dict['mv'] == 2:
                    if MV_str != '':
                        line_lst.append(MV_str)
                        MV_str = ''
            igrade += 1
        n_lines = len(line_lst)
        if MV_str == '':
            if n_lines > 0 and line_lst[0][0] == '+':
                line_lst[0] = line_lst[0][1:]
        else:
            if MV_str[0] == '+':
                MV_str = MV_str[1:]
        if n_lines == 1:
            MV_str = line_lst[0]
            n_lines = 0
        if LatexPrinter.fmt_dict['mv'] >= 2:
            MV_str = '@'+line_lst[0]+' \\\\ \n'
            for line in line_lst[1:-1]:
                MV_str += '& '+line+' \\\\ \n'
            MV_str += '& '+line_lst[-1]+'@\n'
        if MV_str == '':
            MV_str = '0'
        if expr.name != '':
            MV_str = LatexPrinter.extended_symbol(expr.name)+' = '+MV_str
        return(MV_str)

    def _print_Relational(self, expr):
        charmap = {
            "==" : "=",
            "<"  : "<",
            "<=" : r"\leq",
            "!=" : r"\neq",
        }

        return "%s %s %s" % (self._print(expr.lhs),
            charmap[expr.rel_op], self._print(expr.rhs))

    def _print_Matrix(self, expr):
        lines = []

        for line in range(expr.lines): # horrible, should be 'rows'
            lines.append(" & ".join([ self._print(i) for i in expr[line,:] ]))

        if self._inline:
            tex = r"\left(\begin{smallmatrix}%s\end{smallmatrix}\right)"
        else:
            tex = r"\begin{pmatrix}%s\end{pmatrix}"

        return tex % r"\\".join(lines)

    def _print_tuple(self, expr):
        return r"\begin{pmatrix}%s\end{pmatrix}" % \
            r", & ".join([ self._print(i) for i in expr ])

    def _print_list(self, expr):
        return r"\begin{bmatrix}%s\end{bmatrix}" % \
            r", & ".join([ self._print(i) for i in expr ])

    def _print_dict(self, expr):
        items = []

        keys = expr.keys()
        keys.sort(Basic.compare_pretty)
        for key in keys:
            val = expr[key]
            items.append("%s : %s" % (self._print(key), self._print(val)))

        return r"\begin{Bmatrix}%s\end{Bmatrix}" % r", & ".join(items)

    def _print_DiracDelta(self, expr):
        if len(expr.args) == 1 or expr.args[1] == 0:
            tex = r"\delta\left(%s\right)" % self._print(expr.args[0])
        else:
            tex = r"\delta^{\left( %s \right)}\left( %s \right)" % (\
            self._print(expr.args[1]), self._print(expr.args[0]))
        return tex

def LaTeX(expr, inline=True):
    """
    Convert the given expression to LaTeX representation.

    You can specify how the generated code will be delimited.
    If the 'inline' keyword is set then inline LaTeX $ $ will
    be used. Otherwise the resulting code will be enclosed in
    'equation*' environment (remember to import 'amsmath').

    >>> from sympy import Rational
    >>> from sympy.abc import tau, mu

    >>> latex((2*tau)**Rational(7,2))
    '$8 \\\\sqrt{2} \\\\sqrt[7]{\\\\tau}$'

    >>> latex((2*mu)**Rational(7,2), inline=False)
    '\\\\begin{equation*}8 \\\\sqrt{2} \\\\sqrt[7]{\\\\mu}\\\\end{equation*}'

    Besides all Basic based expressions, you can recursively
    convert Pyhon containers (lists, tuples and dicts) and
    also SymPy matrices:

    >>> latex([2/x, y])
    '$\\\\begin{bmatrix}\\\\frac{2}{x}, & y\\\\end{bmatrix}$'

    The extended latex printer will also append the output to a
    string (LatexPrinter.body) that will be processed by xdvi()
    for immediate display one xdvi() is called.
    """
    xstr = LatexPrinter(inline).doprint(expr)
    return (xstr)

def print_LaTeX(expr):
    """Prints LaTeX representation of the given expression."""
    print LaTeX(expr)

def Format(fmt='1 1 1 1'):
    LatexPrinter.format_str(fmt)
    return

def xdvi(filename='tmplatex.tex',debug=False):
    """
    Post processes LaTeX output (see comments below), adds preamble and
    postscript, generates tex file, inputs file to latex, displays resulting
    dvi file with xdvi or yap.
    """
    if not LatexPrinter.LaTeX_flg:
        return
    body = sys.stdout.getvalue()
    LatexPrinter.restore()

    body_lst = body.split('\n')
    body = ''
    array_flg = False
    eqnarray_flg = False
    raw_flg = False
    nline = len(body_lst)
    iline = 0
    i = iter(body_lst)
    line = i.next()

    while True:
        if '$' in line: #Inline math expression(s)
            if len(line) > 0:
                line += '\\newline \n'
            body += line
            try:
                line = i.next()
            except StopIteration:
                break

        elif '%' in line: #Raw LaTeX input
            """
            If % in line assume line is beginning of raw LaTeX input and stop
            post processing
            """
            line = line.replace('%','')
            raw_flg = True
            while raw_flg:
                if '%' in line:
                    """
                    If % in line assume line is end of LaTeX input and begin
                    post processing
                    """
                    raw_flg = False
                    line = line.replace('%','')+'\n'
                else:
                    line += '\n'
                line = process_equals(line)
                body += line
                try:
                    line = i.next()
                except StopIteration:
                    break

        elif '#' in line: #Array input
            """
            If # in line assume line is beginning of array input and contains
            \begin{array} statement
            """
            line = line.replace('#','')
            array_flg = True
            line = '\\begin{equation*}\n'+line
            while array_flg:
                if '#' in line:
                    """
                    If # in line assume line is end of array input and contains
                    \end{array} statment
                    """
                    array_flg = False
                    line = line.replace('#','')
                    line += '\\end{equation*}\n'
                else:
                    line += '\n'
                line = process_equals(line)
                body += line
                try:
                    line = i.next()
                except StopIteration:
                    break

        elif '@' in line: #Align input
            """
            If @ in line assume line is beginning of align input
            """
            line = line.replace('@','')
            line = line.replace('=','& = ')
            eqnarray_flg = True
            line = '\\begin{align*}\n'+line
            line = process_equals(line)
            body += line
            try:
                line = i.next()
            except StopIteration:
                break
            while eqnarray_flg:
                if '@' in line:
                    """
                    If @ in line assume line is end of align input
                    """
                    eqnarray_flg = False
                    line = line.replace('@','')
                    line += '\\end{align*}\n'
                else:
                    line+'\n'
                line = process_equals(line)
                body += line
                try:
                    line = i.next()
                except StopIteration:
                    break

        else:
            if '=' in line: #Single line equation
                line = '\\begin{equation*}\n'+line+'\n\\end{equation*}'
            else: #Text with no math expression(s)unless \ or _ in line
                if '\\' in line or '_' in line or '^' in line:
                    line = '\\begin{equation*}\n'+line+'\n\\end{equation*}'
                else:
                    if len(line) > 0:
                        line += '\\newline \n'
            line = process_equals(line)
            body += line
            try:
                line = i.next()
            except StopIteration:
                break
    body = LatexPrinter.preamble+body+LatexPrinter.postscript
    latex_file = open(filename,'w')
    latex_file.write(body)
    latex_file.close()

    latex_str = None
    xdvi_str  = None

    if find_executable('latex') != None:
        latex_str = 'latex'

    if find_executable('xdvi') != None:
        xdvi_str = 'xdvi'

    if find_executable('yap') != None:
        xdvi_str = 'yap'

    if latex_str != None and xdvi_str != None:
        if debug: #Display latex excution output for debugging purposes
            os.system(latex_str+' '+filename[:-4])
        else: #Works for Linux don't know about Windows
            if sys.platform == 'linux2':
                os.system(latex_str+' '+filename[:-4]+' > /dev/null')
            else:
                os.system(latex_str+' '+filename[:-4]+' > NUL')
        os.system(xdvi_str+' '+filename[:-4]+' &')
    LatexPrinter.LaTeX_flg = False
    return

def MV_format(mv_fmt):
    """
    0 or 1 - Print multivector on one line
    2      - Print each multivector grade on one line
    3      - Print each multivector base on one line
    """
    if LatexPrinter.LaTeX_flg:
        LatexPrinter.fmt_dict['mv'] = mv_fmt
    return

def fct_format(fct_fmt):
    """
    0 - Default sympy latex format
    1 - Do not print arguments of arbitrary functions.
        Use symbol font for arbitrary functions.
        Use enhanced symbol nameing for arbitrary functions.
        Use new names for standard functions (acos -> Cos^{-1})
    """
    if LatexPrinter.LaTeX_flg:
        LatexPrinter.fct = fct_fmt
    return

def pdiff_format(pdiff_fmt):
    """
    0 - Use default sympy partial derivative format
    1 - Contracted derivative format (no fraction symbols)
    """
    if LatexPrinter.LaTeX_flg:
        LatexPrinter.fmt_dict['pdiff'] = pdiff_fmt
    return

def sym_format(sym_fmt):
    """
    0 - Use default sympy format
    1 - Ues extended symbol format including multiple greek letters in
        basic symbol (symbol precceeding sub and superscripts)and in
        sub and superscripts of basic symbol and accents in basic symbol
    """
    if LatexPrinter.LaTeX_flg:
        LatexPrinter.fmt_dict['sym'] = sym_fmt
    return

def str_format(str_fmt):
    """
    0 - Use default sympy format
    1 - Ues extended symbol format including multiple greek letters in
        basic symbol (symbol precceeding sub and superscripts)and in
        sub and superscripts of basic symbol and accents in basic symbol
    """
    if LatexPrinter.LaTeX_flg:
        LatexPrinter.fmt_dict['str'] = str_fmt
    return

def ext_str(xstr):
    return(LatexPrinter.extended_symbol(xstr))

