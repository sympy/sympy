"""
A Printer which converts an expression into its LaTeX equivalent.
"""

from sympy.core import S, C, Basic, Symbol, Wild, var
from printer import Printer
from conventions import split_super_sub
from sympy.simplify import fraction
from sympy import Interval

import sympy.mpmath.libmpf as mlib
from sympy.mpmath.settings import prec_to_dps

import re

class LatexPrinter(Printer):
    printmethod = "_latex_"

    def __init__(self, profile=None):
        Printer.__init__(self)

        mul_symbol_table = {
            None : r" ",
            "ldot" : r" \,.\, ",
            "dot" : r" \cdot ",
            "times" : r" \times "
        }

        self._settings = {
            "inline" : True,
            "descending" : False,
            "mainvar" : None,
            "fold_frac_powers" : False,
            "fold_func_brackets" : False,
            "mul_symbol" : None,
            "inv_trig_style" : "abbreviated",
            "mat_str" : "smallmatrix",
            "mat_delim" : "(",
        }

        self._delim_dict = {'(':')','[':']'}

        if profile is not None:
            if profile.has_key('inline'):
                if not profile['inline']:
                    #change to "good" defaults for inline=False before
                    #updating with the settings from profile
                    self._settings['mat_str'] = 'bmatrix'
                    self._settings['mat_delim'] = None
            self._settings.update(profile)

        self._settings['mul_symbol_latex'] = \
            mul_symbol_table[self._settings['mul_symbol']]

    def doprint(self, expr):
        tex = Printer.doprint(self, expr)

        if self._settings['inline']:
            return r"$%s$" % tex
        else:
            return r"\begin{equation*}%s\end{equation*}" % tex

    def _needs_brackets(self, expr):
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed, False otherwise. For example: a + b => True; a => False;
        10 => False; -10 => True.
        """
        return not ((expr.is_Integer and expr.is_nonnegative) or expr.is_Atom)

    def _needs_function_brackets(self, expr):
        """
        Returns True if the expression needs to be wrapped in brackets when
        passed as an argument to a function, False otherwise. This is a more
        liberal version of _needs_brackets, in that many expressions which need
        to be wrapped in brackets when added/substracted/raised to a power do
        not need them when passed to a function. Such an example is a*b.
        """
        if not self._needs_brackets(expr):
            return False
        else:
            # Muls of the form a*b*c... can be folded
            if expr.is_Mul and not self._mul_is_clean(expr):
                return True
            # Pows which don't need brackets can be folded
            elif expr.is_Pow and not self._pow_is_clean(expr):
                return True
            # Add and Function always need brackets
            elif expr.is_Add or expr.is_Function:
                return True
            else:
                return False

    def _mul_is_clean(self, expr):
        for arg in expr.args:
            if arg.is_Function:
                return False
        return True

    def _pow_is_clean(self, expr):
        return not self._needs_brackets(expr.base)

    def _do_exponent(self, expr, exp):
        if exp is not None:
            return r"\left(%s\right)^{%s}" % (expr, exp)
        else:
            return expr

    def _print_Add(self, expr):
        args = list(expr.args)

        if self._settings['mainvar'] is not None:
            mainvar = self._settings['mainvar']
            if type(mainvar) == str:
                mainvar = var(mainvar)

            def compare_exponents(a, b):
               p1, p2 = Wild("p1"), Wild("p2")
               r_a = a.match(p1 * mainvar**p2)
               r_b = b.match(p1 * mainvar**p2)
               if r_a is None and r_b is None:
                   c = Basic._compare_pretty(a,b)
                   return c
               elif r_a is not None:
                   if r_b is None:
                       return 1
                   else:
                       c = Basic.compare(r_a[p2], r_b[p2])
                       if c!=0:
                           return c
                       else:
                           c = Basic._compare_pretty(a,b)
                           return c
               elif r_b is not None and r_a is None:
                    return -1

            args.sort(compare_exponents)
        else:
            args.sort(Basic._compare_pretty)
        if self._settings['descending']:
            args.reverse()

        tex = str(self._print(args[0]))

        for term in args[1:]:
            coeff = term.as_coeff_terms()[0]

            if coeff.is_negative:
                tex += r" %s" % self._print(term)
            else:
                tex += r" + %s" % self._print(term)

        return tex

    def _print_Real(self, expr):
        # Based off of that in StrPrinter
        dps = prec_to_dps(expr._prec)
        str_real = mlib.to_str(expr._mpf_, dps, strip_zeros=True)

        # Must always have a mul symbol (as 2.5 10^{20} just looks odd)
        seperator = r" \times "

        if self._settings['mul_symbol'] is not None:
            seperator = self._settings['mul_symbol_latex']

        if 'e' in str_real:
            (mant, exp) = str_real.split('e')

            if exp[0] == '+':
                exp = exp[1:]

            return r"%s%s10^{%s}" % (mant, seperator, exp)
        elif str_real == "+inf":
            return r"\infty"
        elif str_real == "-inf":
            return r"- \intfy"
        else:
            return str_real

    def _print_Mul(self, expr):
        coeff, terms = expr.as_coeff_terms()

        if not coeff.is_negative:
            tex = ""
        else:
            coeff = -coeff
            tex = "- "

        numer, denom = fraction(C.Mul(*terms))
        seperator = self._settings['mul_symbol_latex']

        def convert(terms):
            if not terms.is_Mul:
                return str(self._print(terms))
            else:
                _tex = last_term_tex = ""
                for term in terms.args:
                    pretty = self._print(term)

                    if term.is_Add:
                        term_tex = (r"\left(%s\right)" % pretty)
                    else:
                        term_tex = str(pretty)

                    # between two digits, \times must always be used,
                    # to avoid confusion
                    if seperator == " " and \
                            re.search("[0-9][} ]*$", last_term_tex) and \
                            re.match("[{ ]*[-+0-9]", term_tex):
                        _tex += r" \times "
                    elif _tex:
                        _tex += seperator

                    _tex += term_tex
                    last_term_tex = term_tex
                return _tex

        if denom is S.One:
            if numer.is_Add:
                _tex = r"\left(%s\right)" % convert(numer)
            else:
                _tex = r"%s" % convert(numer)

            if coeff is not S.One:
                tex += str(self._print(coeff))

                # between two digits, \times must always be used, to avoid
                # confusion
                if seperator == " " and re.search("[0-9][} ]*$", tex) and \
                        re.match("[{ ]*[-+0-9]", _tex):
                    tex +=  r" \times " + _tex
                else:
                    tex += seperator + _tex
            else:
                tex += _tex

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
        elif self._settings['fold_frac_powers'] \
             and expr.exp.is_Rational \
             and expr.exp.q != 1:
            base, p, q = self._print(expr.base), expr.exp.p, expr.exp.q
            return r"%s^{%s/%s}" % (base, p, q)
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
                    tex = r"%s^{%s}"

                return tex % (self._print(expr.base),
                              self._print(expr.exp))

    def _print_Derivative(self, expr):
        dim = len(expr.symbols)

        if dim == 1:
            tex = r"\frac{\partial}{\partial %s}" % \
                self._print(expr.symbols[0])
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
                if not self._settings['inline']:
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
                name = r"\operatorname{%s}^{-1}" % func
            elif exp is not None:
                name = r"\operatorname{%s}^{%s}" % (func, exp)
            else:
                name = r"\operatorname{%s}" % func

            if can_fold_brackets:
                name += r"%s"
            else:
                name += r"\left(%s\right)"

            if inv_trig_power_case and exp is not None:
                name += r"^{%s}" % exp

            return name % ",".join(args)

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
        tex = r"e^{%s}" % self._print(expr.args[0])
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

    def _print_Symbol(self, expr):
        name, supers, subs = split_super_sub(expr.name)

        # translate name, supers and subs to tex keywords
        greek = set([ 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                      'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu',
                      'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon',
                      'phi', 'chi', 'psi', 'omega' ])

        other = set( ['aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth',
                      'hbar', 'hslash', 'mho' ])

        def translate(s):
            tmp = s.lower()
            if tmp in greek or tmp in other:
                return "\\" + s
            else:
                return s

        name = translate(name)
        supers = [translate(sup) for sup in supers]
        subs = [translate(sub) for sub in subs]

        # glue all items together:
        if len(supers) > 0:
            name += "^{%s}" % " ".join(supers)
        if len(subs) > 0:
            name += "_{%s}" % " ".join(subs)

        return name

    def _print_Relational(self, expr):
        charmap = {
            "==" : "=",
            "<"  : "<",
            "<=" : r"\leq",
            "!=" : r"\neq",
        }

        return "%s %s %s" % (self._print(expr.lhs),
            charmap[expr.rel_op], self._print(expr.rhs))

    def _print_Piecewise(self, expr):
        ecpairs = [r"%s & for %s" % (self._print(e), self._print(c)) \
                       for e, c in expr.args[:-1]]
        if expr.args[-1].cond == True:
            ecpairs.append(r"%s & \textrm{otherwise}" % \
                               self._print(expr.args[-1].expr))
        else:
            ecpairs.append(r"%s & for %s" % \
                           (self._print(expr.args[-1].cond),
                            self._print(expr.args[-1].expr)))
        tex = r"\left\{\begin{array}{cl} %s \end{array}\right."
        return tex % r" \\".join(ecpairs)

    def _print_Matrix(self, expr):
        lines = []

        for line in range(expr.rows): # horrible, should be 'rows'
            lines.append(" & ".join([ self._print(i) for i in expr[line,:] ]))

        out_str = r'\begin{%MATSTR%}%s\end{%MATSTR%}'
        out_str = out_str.replace('%MATSTR%', self._settings['mat_str'])
        if self._settings['mat_delim']:
            left_delim = self._settings['mat_delim']
            right_delim = self._delim_dict[left_delim]
            out_str = r'\left' + left_delim + out_str + \
                      r'\right' + right_delim
        return out_str % r"\\".join(lines)

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

    def _print_Interval(self, i):
        if i.start == i.end:
            return r"\left{%s\right}" % self._print(i.start)

        else:
            if i.left_open:
                left = '('
            else:
                left = '['

            if i.right_open:
                right = ')'
            else:
                right = ']'

            return r"\left%s%s, %s\right%s" % \
                   (left, self._print(i.start), self._print(i.end), right)

    def _print_Union(self, u):
        other_sets, singletons = [], []
        for set in u.args:
            if isinstance(set, Interval) and set.measure == 0:
                singletons.append(set.start)
            else:
                other_sets.append(set)

        S2 = r"%s" % \
             r" \cup ".join([ self._print_Interval(i) for i in other_sets ])

        if len(singletons) > 0:
            S1 = r"\left{%s\right}" % \
                 r", ".join([ self._print(i) for i in singletons ])

            S = r"%s \cup %s" % (S1, S2)
        else:
            S = S2

        return S

    def _print_EmptySet(self, e):
        return r"\emptyset"

def latex(expr, profile=None, **kargs):
    r"""Convert the given expression to LaTeX representation.

        You can specify how the generated code will be delimited.
        If the 'inline' keyword is set then inline LaTeX $ $ will
        be used. Otherwise the resulting code will be enclosed in
        'equation*' environment (remember to import 'amsmath').

        >>> from sympy import *
        >>> from sympy.abc import *

        >>> latex((2*tau)**Rational(7,2))
        '$8 \\sqrt{2} \\sqrt[7]{\\tau}$'

        >>> latex((2*mu)**Rational(7,2), inline=False)
        '\\begin{equation*}8 \\sqrt{2} \\sqrt[7]{\\mu}\\end{equation*}'

        Besides all Basic based expressions, you can recursively
        convert Pyhon containers (lists, tuples and dicts) and
        also SymPy matrices:

        >>> latex([2/x, y])
        '$\\begin{bmatrix}\\frac{2}{x}, & y\\end{bmatrix}$'

    """

    if profile is not None:
        profile.update(kargs)
    else:
        profile = kargs

    return LatexPrinter(profile).doprint(expr)

def print_latex(expr):
    """Prints LaTeX representation of the given expression."""
    print latex(expr)
