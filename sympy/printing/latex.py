"""
A Printer which converts an expression into its LaTeX equivalent.
"""

from sympy.core import S, C, Add, Symbol
from sympy.core.function import _coeff_isneg
from printer import Printer
from conventions import split_super_sub
from sympy.simplify import fraction
from sympy.core.sympify import SympifyError

import sympy.mpmath.libmp as mlib
from sympy.mpmath.libmp import prec_to_dps

from sympy.utilities.misc import default_sort_key
from sympy.utilities.iterables import has_variety

import re, warnings

# Hand-picked functions which can be used directly in both LaTeX and MathJax
# Complete list at http://www.mathjax.org/docs/1.1/tex.html#supported-latex-commands
# This variable only contains those functions which sympy uses.
accepted_latex_functions = ['arcsin','arccos','arctan','sin','cos','tan',
                    'theta','beta','alpha','gamma','sinh','cosh','tanh','sqrt',
                    'ln','log','sec','csc','cot','coth','re','im','frac','root',
                    'arg','zeta','psi']

class LatexPrinter(Printer):
    printmethod = "_latex"

    _default_settings = {
        "order": None,
        "mode": "plain",
        "itex": False,
        "fold_frac_powers": False,
        "fold_func_brackets": False,
        "mul_symbol": None,
        "inv_trig_style": "abbreviated",
        "mat_str": "smallmatrix",
        "mat_delim": "[",
        "symbol_names": {},
    }

    def __init__(self, settings=None):
        if settings is not None and 'inline' in settings and not settings['inline']:
            # Change to "good" defaults for inline=False
            settings['mat_str'] = 'bmatrix'
            settings['mat_delim'] = None
        Printer.__init__(self, settings)

        if 'mode' in self._settings:
            valid_modes = ['inline', 'plain', 'equation', \
                            'equation*']
            if self._settings['mode'] not in valid_modes:
                raise ValueError("'mode' must be one of 'inline', 'plain', " \
                    "'equation' or 'equation*'")

        mul_symbol_table = {
            None : r" ",
            "ldot" : r" \,.\, ",
            "dot" : r" \cdot ",
            "times" : r" \times "
        }

        self._settings['mul_symbol_latex'] = \
            mul_symbol_table[self._settings['mul_symbol']]

        self._delim_dict = {'(':')','[':']'}

    def doprint(self, expr):
        tex = Printer.doprint(self, expr)

        if self._settings['mode'] == 'plain':
            return tex
        elif self._settings['mode'] == 'inline':
            return r"$%s$" % tex
        elif self._settings['itex']:
            return r"$$%s$$" % tex
        else:
            env_str = self._settings['mode']
            return r"\begin{%s}%s\end{%s}" % (env_str, tex, env_str)

    def _needs_brackets(self, expr):
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed, False otherwise. For example: a + b => True; a => False;
        10 => False; -10 => True.
        """
        return not ((expr.is_Integer and expr.is_nonnegative)
                or (expr.is_Atom and expr is not S.NegativeOne))

    def _needs_function_brackets(self, expr):
        """
        Returns True if the expression needs to be wrapped in brackets when
        passed as an argument to a function, False otherwise. This is a more
        liberal version of _needs_brackets, in that many expressions which need
        to be wrapped in brackets when added/subtracted/raised to a power do
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

    def _print_Add(self, expr, order=None):
        if self.order == 'none':
            terms = list(expr.args)
        else:
            terms = self._as_ordered_terms(expr, order=order)
        tex = self._print(terms[0])

        for term in terms[1:]:
            if not _coeff_isneg(term):
                tex += " +"

            tex += " " + self._print(term)

        return tex

    def _print_Float(self, expr):
        # Based off of that in StrPrinter
        dps = prec_to_dps(expr._prec)
        str_real = mlib.to_str(expr._mpf_, dps, strip_zeros=True)

        # Must always have a mul symbol (as 2.5 10^{20} just looks odd)
        separator = r" \times "

        if self._settings['mul_symbol'] is not None:
            separator = self._settings['mul_symbol_latex']

        if 'e' in str_real:
            (mant, exp) = str_real.split('e')

            if exp[0] == '+':
                exp = exp[1:]

            return r"%s%s10^{%s}" % (mant, separator, exp)
        elif str_real == "+inf":
            return r"\infty"
        elif str_real == "-inf":
            return r"- \infty"
        else:
            return str_real

    def _print_Mul(self, expr):
        coeff, tail = expr.as_coeff_Mul()

        if not coeff.is_negative:
            tex = ""
        else:
            coeff = -coeff
            tex = "- "

        numer, denom = fraction(tail, exact=True)
        separator = self._settings['mul_symbol_latex']

        def convert(expr):
            if not expr.is_Mul:
                return str(self._print(expr))
            else:
                _tex = last_term_tex = ""

                if self.order not in ('old', 'none'):
                    args = expr.as_ordered_factors()
                else:
                    args = expr.args

                for term in args:
                    pretty = self._print(term)

                    if term.is_Add:
                        term_tex = (r"\left(%s\right)" % pretty)
                    else:
                        term_tex = str(pretty)

                    # between two digits, \times must always be used,
                    # to avoid confusion
                    if separator == " " and \
                            re.search("[0-9][} ]*$", last_term_tex) and \
                            re.match("[{ ]*[-+0-9]", term_tex):
                        _tex += r" \times "
                    elif _tex:
                        _tex += separator

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
                if separator == " " and re.search("[0-9][} ]*$", tex) and \
                        re.match("[{ ]*[-+0-9]", _tex):
                    tex +=  r" \times " + _tex
                else:
                    tex += separator + _tex
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
        # Treat x**Rational(1,n) as special case
        if expr.exp.is_Rational and abs(expr.exp.p) == 1 and expr.exp.q != 1:
            base = self._print(expr.base)
            expq = expr.exp.q

            if expq == 2:
                tex = r"\sqrt{%s}" % base
            elif self._settings['itex']:
                tex = r"\root{%d}{%s}" % (expq,base)
            else:
                tex = r"\sqrt[%d]{%s}" % (expq,base)

            if expr.exp.is_negative:
                return r"\frac{1}{%s}" % tex
            else:
                return tex
        elif self._settings['fold_frac_powers'] \
             and expr.exp.is_Rational \
             and expr.exp.q != 1:
            base, p, q = self._print(expr.base), expr.exp.p, expr.exp.q
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
                    tex = r"%s^{%s}"

                return tex % (self._print(expr.base),
                              self._print(expr.exp))

    def _print_Sum(self, expr):
        if len(expr.limits) == 1:
            tex = r"\sum_{%s=%s}^{%s} " % \
                tuple([ self._print(i) for i in expr.limits[0] ])
        else:
            def _format_ineq(l):
                return r"%s \leq %s \leq %s" % \
                    tuple([self._print(s) for s in l[1], l[0], l[2]])

            tex = r"\sum_{\substack{%s}} " % \
                str.join('\\\\', [ _format_ineq(l) for l in expr.limits ])

        if isinstance(expr.function, Add):
            tex += r"\left(%s\right)" % self._print(expr.function)
        else:
            tex += self._print(expr.function)

        return tex

    def _print_Product(self, expr):
        if len(expr.limits) == 1:
            tex = r"\prod_{%s=%s}^{%s} " % \
                tuple([ self._print(i) for i in expr.limits[0] ])
        else:
            def _format_ineq(l):
                return r"%s \leq %s \leq %s" % \
                    tuple([self._print(s) for s in l[1], l[0], l[2]])

            tex = r"\prod_{\substack{%s}} " % \
                str.join('\\\\', [ _format_ineq(l) for l in expr.limits ])

        if isinstance(expr.function, Add):
            tex += r"\left(%s\right)" % self._print(expr.function)
        else:
            tex += self._print(expr.function)

        return tex

    def _print_Derivative(self, expr):
        dim = len(expr.variables)

        if dim == 1:
            tex = r"\frac{\partial}{\partial %s}" % \
                self._print(expr.variables[0])
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

    def _print_Subs(self, subs):
        expr, old, new = subs.args
        latex_expr = self._print(expr)
        latex_old = (self._print(e) for e in old)
        latex_new = (self._print(e) for e in new)
        latex_subs = r'\\ '.join(e[0] + '=' + e[1] for e in zip(latex_old, latex_new))
        return r'\left. %s \right|_{\substack{ %s }}' % (latex_expr, latex_subs)

    def _print_Integral(self, expr):
        tex, symbols = "", []

        # Only up to \iiiint exists
        if len(expr.limits) <= 4 and all(len(lim) == 1 for lim in expr.limits):
            # Use len(expr.limits)-1 so that syntax highlighters don't think
            # \" is an escaped quote
            tex = r"\i" + "i"*(len(expr.limits)-1) + "nt"
            symbols = [r"\, d%s" % self._print(symbol[0]) for symbol in expr.limits]

        else:
            for lim in reversed(expr.limits):
                symbol = lim[0]
                tex += r"\int"

                if len(lim) > 1:
                    if self._settings['mode'] in ['equation','equation*'] \
                       and not self._settings['itex']:
                        tex += r"\limits"

                    if len(lim) == 3:
                        tex += "_{%s}^{%s}" % (self._print(lim[1]),
                                               self._print(lim[2]))
                    if len(lim) == 2:
                        tex += "^{%s}" % (self._print(lim[1]))

                symbols.insert(0, r"\, d%s" % self._print(symbol))

        return r"%s %s%s" % (tex,
            str(self._print(expr.function)), "".join(symbols))

    def _print_Limit(self, expr):
        e, z, z0, dir = expr.args

        tex = r"\lim_{%s \to %s}" % (self._print(z),
                                     self._print(z0))

        if isinstance(e, C.AssocOp):
            return r"%s\left(%s\right)" % (tex, self._print(e))
        else:
            return r"%s %s" % (tex, self._print(e))

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
                if func in accepted_latex_functions:
                    name = r"\%s^{-1}" % func
                else:
                    name = r"\operatorname{%s}^{-1}" % func
            elif exp is not None:
                if func in accepted_latex_functions:
                    name = r"\%s^{%s}" % (func,exp)
                else:
                    # If the generic function name contains an underscore, handle it
                    name = r"\operatorname{%s}^{%s}" % (func.replace("_", r"\_"), exp)
            else:
                if func in accepted_latex_functions:
                    name = r"\%s" % func
                else:
                    # If the generic function name contains an underscore, handle it
                    name = r"\operatorname{%s}" % func.replace("_", r"\_")

            if can_fold_brackets:
                if func in accepted_latex_functions:
                    # Wrap argument safely to avoid parse-time conflicts
                    # with the function name itself
                    name += r" {%s}"
                else:
                    name += r"%s"
            else:
                name += r"{\left (%s \right )}"

            if inv_trig_power_case and exp is not None:
                name += r"^{%s}" % exp

            return name % ",".join(args)

    def _print_Lambda(self, expr):
        symbols, expr = expr.args

        if len(symbols) == 1:
            symbols = self._print(symbols[0])
        else:
            symbols = self._print(tuple(symbols))

        args = (symbols, self._print(expr))
        tex = r"\Lambda {\left (%s \right )}" % ", ".join(args)

        return tex

    def _print_Min(self, expr, exp=None):
        args = sorted(expr.args, key=default_sort_key)
        texargs = [r"%s" % self._print(symbol) for symbol in args]
        return r"\min\left(%s\right)" % ", ".join(texargs)

    def _print_Max(self, expr, exp=None):
        args = sorted(expr.args, key=default_sort_key)
        texargs = [r"%s" % self._print(symbol) for symbol in args]
        return r"\max\left(%s\right)" % ", ".join(texargs)

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

    def _print_Abs(self, expr, exp=None):
        tex = r"\lvert{%s}\rvert" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_re(self, expr, exp=None):
        if self._needs_brackets(expr.args[0]):
            tex = r"\Re {\left (%s \right )}" % self._print(expr.args[0])
        else:
            tex = r"\Re{%s}" % self._print(expr.args[0])

        return self._do_exponent(tex, exp)

    def _print_im(self, expr, exp=None):
        if self._needs_brackets(expr.args[0]):
            tex = r"\Im {\left ( %s \right )}" % self._print(expr.args[0])
        else:
            tex = r"\Im{%s}" % self._print(expr.args[0])

        return self._do_exponent(tex, exp)

    def _print_Not(self, e):
        if (e.args[0].is_Boolean):
            return r"\neg (%s)" % self._print(e.args[0])
        else:
            return r"\neg %s" % self._print(e.args[0])

    def _print_And(self, e):
        args = sorted(e.args, key=default_sort_key)
        arg = args[0]
        if arg.is_Boolean and not arg.is_Not:
            tex = r"\left(%s\right)" % self._print(arg)
        else:
            tex = r"%s" % self._print(arg)

        for arg in args[1:]:
            if arg.is_Boolean and not arg.is_Not:
                tex += r" \wedge \left(%s\right)" % (self._print(arg))
            else:
                tex += r" \wedge %s" % (self._print(arg))

        return tex

    def _print_Or(self, e):
        args = sorted(e.args, key=default_sort_key)
        arg = args[0]
        if arg.is_Boolean and not arg.is_Not:
            tex = r"\left(%s\right)" % self._print(arg)
        else:
            tex = r"%s" % self._print(arg)

        for arg in args[1:]:
            if arg.is_Boolean and not arg.is_Not:
                tex += r" \vee \left(%s\right)" % (self._print(arg))
            else:
                tex += r" \vee %s" % (self._print(arg))

        return tex

    def _print_Implies(self, e):
        return r"%s \Rightarrow %s" % (self._print(e.args[0]), self._print(e.args[1]))

    def _print_Equivalent(self, e):
        return r"%s \Leftrightarrow %s" % (self._print(e.args[0]), self._print(e.args[1]))

    def _print_conjugate(self, expr, exp=None):
        tex = r"\overline{%s}" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_ExpBase(self, expr, exp=None):
        # TODO should exp_polar be printed differently?
        #      what about exp_polar(0), exp_polar(1)?
        tex = r"e^{%s}" % self._print(expr.args[0])
        return self._do_exponent(tex, exp)

    def _print_gamma(self, expr, exp=None):
        tex = r"\left(%s\right)" % self._print(expr.args[0])

        if exp is not None:
            return r"\Gamma^{%s}%s" % (exp, tex)
        else:
            return r"\Gamma%s" % tex

    def _print_uppergamma(self, expr, exp=None):
        tex = r"\left(%s, %s\right)" % (self._print(expr.args[0]),
                                        self._print(expr.args[1]))

        if exp is not None:
            return r"\Gamma^{%s}%s" % (exp, tex)
        else:
            return r"\Gamma%s" % tex

    def _print_lowergamma(self, expr, exp=None):
        tex = r"\left(%s, %s\right)" % (self._print(expr.args[0]),
                                        self._print(expr.args[1]))

        if exp is not None:
            return r"\gamma^{%s}%s" % (exp, tex)
        else:
            return r"\gamma%s" % tex

    def _print_expint(self, expr, exp=None):
        tex = r"\left(%s\right)" % self._print(expr.args[1])
        nu = self._print(expr.args[0])

        if exp is not None:
            return r"\operatorname{E}_{%s}^{%s}%s" % (nu, exp, tex)
        else:
            return r"\operatorname{E}_{%s}%s" % (nu, tex)

    def _print_factorial(self, expr, exp=None):
        x = expr.args[0]
        if self._needs_brackets(x):
            tex = r"\left(%s\right)!" % self._print(x)
        else:
            tex = self._print(x) + "!"

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_factorial2(self, expr, exp=None):
        x = expr.args[0]
        if self._needs_brackets(x):
            tex = r"\left(%s\right)!!" % self._print(x)
        else:
            tex = self._print(x) + "!!"

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_binomial(self, expr, exp=None):
        tex = r"{\binom{%s}{%s}}" % (self._print(expr.args[0]),
                                     self._print(expr.args[1]))

        if exp is not None:
            return r"%s^{%s}" % (tex, exp)
        else:
            return tex

    def _print_RisingFactorial(self, expr, exp=None):
        tex = r"{\left(%s\right)}^{\left(%s\right)}" % \
            (self._print(expr.args[0]), self._print(expr.args[1]))

        return self._do_exponent(tex, exp)

    def _print_FallingFactorial(self, expr, exp=None):
        tex = r"{\left(%s\right)}_{\left(%s\right)}" % \
            (self._print(expr.args[0]), self._print(expr.args[1]))

        return self._do_exponent(tex, exp)

    def _hprint_BesselBase(self, expr, exp, sym):
        tex = r"%s" % (sym)

        need_exp = False
        if exp is not None:
            if tex.find('^') == -1:
                tex = r"%s^{%s}" % (tex, self._print(exp))
            else:
                need_exp = True

        tex = r"%s_{%s}\left(%s\right)" % (tex, self._print(expr.order),
                                           self._print(expr.argument))

        if need_exp:
            tex = self._do_exponent(tex, exp)
        return tex

    def _hprint_vec(self, vec):
        if len(vec) == 0:
            return ""
        s = ""
        for i in vec[:-1]:
            s += "%s, " % self._print(i)
        s += self._print(vec[-1])
        return s

    def _print_besselj(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'J')

    def _print_besseli(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'I')

    def _print_besselk(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'K')

    def _print_bessely(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'Y')

    def _print_yn(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'y')

    def _print_jn(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'j')

    def _print_hankel1(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'H^{(1)}')

    def _print_hankel2(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'H^{(2)}')

    def _print_fresnels(self, expr, exp=None):
        tex = r"\left(%s\right)" % self._print(expr.args[0])

        if exp is not None:
            return r"S^{%s}%s" % (exp, tex)
        else:
            return r"S%s" % tex

    def _print_fresnelc(self, expr, exp=None):
        tex = r"\left(%s\right)" % self._print(expr.args[0])

        if exp is not None:
            return r"C^{%s}%s" % (exp, tex)
        else:
            return r"C%s" % tex

    def _print_hyper(self, expr, exp=None):
        tex = r"{{}_{%s}F_{%s}\left(\begin{matrix} %s \\ %s \end{matrix}" \
              r"\middle| {%s} \right)}" % \
             (self._print(len(expr.ap)), self._print(len(expr.bq)),
              self._hprint_vec(expr.ap), self._hprint_vec(expr.bq),
              self._print(expr.argument))

        if exp is not None:
            tex = r"{%s}^{%s}" % (tex, self._print(exp))
        return tex

    def _print_meijerg(self, expr, exp=None):
        tex = r"{G_{%s, %s}^{%s, %s}\left(\begin{matrix} %s & %s \\" \
              r"%s & %s \end{matrix} \middle| {%s} \right)}" % \
             (self._print(len(expr.ap)), self._print(len(expr.bq)),
              self._print(len(expr.bm)), self._print(len(expr.an)),
              self._hprint_vec(expr.an), self._hprint_vec(expr.aother),
              self._hprint_vec(expr.bm), self._hprint_vec(expr.bother),
              self._print(expr.argument))

        if exp is not None:
            tex = r"{%s}^{%s}" % (tex, self._print(exp))
        return tex

    def _print_dirichlet_eta(self, expr, exp=None):
        tex = r"\left(%s\right)" % self._print(expr.args[0])
        if exp is not None:
            return r"\eta^{%s}%s" % (self._print(exp), tex)
        return r"\eta%s" % tex

    def _print_zeta(self, expr, exp=None):
        if len(expr.args) == 2:
            tex = r"\left(%s, %s\right)" % tuple(map(self._print, expr.args))
        else:
            tex = r"\left(%s\right)" % self._print(expr.args[0])
        if exp is not None:
            return r"\zeta^{%s}%s" % (self._print(exp), tex)
        return r"\zeta%s" % tex

    def _print_lerchphi(self, expr, exp=None):
        tex = r"\left(%s, %s, %s\right)" % tuple(map(self._print, expr.args))
        if exp is None:
            return r"\Phi%s" % tex
        return r"\Phi^{%s}%s" % (self._print(exp), tex)

    def _print_polylog(self, expr, exp=None):
        s, z = map(self._print, expr.args)
        tex = r"\left(%s\right)" % z
        if exp is None:
            return r"\operatorname{Li}_{%s}%s" % (s, tex)
        return r"\operatorname{Li}_{%s}^{%s}%s" % (s, self._print(exp), tex)

    def _print_jacobi(self, expr, exp=None):
        n, a, b, x = map(self._print, expr.args)
        tex = r"P_{%s}^{\left(%s,%s\right)}\left(%s\right)" % (n, a, b, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

    def _print_gegenbauer(self, expr, exp=None):
        n, a, x = map(self._print, expr.args)
        tex = r"C_{%s}^{\left(%s\right)}\left(%s\right)" % (n, a, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

    def _print_chebyshevt(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r"T_{%s}\left(%s\right)" % (n, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

    def _print_chebyshevu(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r"U_{%s}\left(%s\right)" % (n, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

    def _print_legendre(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r"P_{%s}\left(%s\right)" % (n, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

    def _print_assoc_legendre(self, expr, exp=None):
        n, a, x = map(self._print, expr.args)
        tex = r"P_{%s}^{\left(%s\right)}\left(%s\right)" % (n, a, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

    def _print_hermite(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r"H_{%s}\left(%s\right)" % (n, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

    def _print_laguerre(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r"L_{%s}\left(%s\right)" % (n, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

    def _print_assoc_laguerre(self, expr, exp=None):
        n, a, x = map(self._print, expr.args)
        tex = r"L_{%s}^{\left(%s\right)}\left(%s\right)" % (n, a, x)
        if exp is not None:
            tex = r"\left(" + tex + r"\right)^{%s}" % (self._print(exp))
        return tex

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
        return r"\mathcal{O}\left(%s\right)" % \
            self._print(expr.args[0])

    def _print_Symbol(self, expr):
        if expr in self._settings['symbol_names']:
            return self._settings['symbol_names'][expr]

        name, supers, subs = split_super_sub(expr.name)

        # translate name, supers and subs to tex keywords
        greek = set([ 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                      'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu',
                      'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon',
                      'phi', 'chi', 'psi', 'omega' ])

        greek_translated = {'lamda': 'lambda', 'Lamda': 'Lambda'}

        other = set( ['aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth',
                      'hbar', 'hslash', 'mho' ])

        def translate(s):
            tmp = s.lower()
            if tmp in greek or tmp in other:
                return "\\" + s
            if s in greek_translated:
                return "\\" + greek_translated[s]
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
        if self._settings['itex']:
            gt = r"\gt"
            lt = r"\lt"
        else:
            gt = ">"
            lt = "<"

        charmap = {
            "==" : "=",
            ">"  : gt,
            "<"  : lt,
            ">=" : r"\geq",
            "<=" : r"\leq",
            "!=" : r"\neq",
        }

        return "%s %s %s" % (self._print(expr.lhs),
            charmap[expr.rel_op], self._print(expr.rhs))

    def _print_Piecewise(self, expr):
        ecpairs = [r"%s & \text{for}\: %s" % (self._print(e), self._print(c)) \
                       for e, c in expr.args[:-1]]
        if expr.args[-1].cond == True:
            ecpairs.append(r"%s & \text{otherwise}" % \
                               self._print(expr.args[-1].expr))
        else:
            ecpairs.append(r"%s & \text{for}\: %s" % \
                           (self._print(expr.args[-1].expr),
                            self._print(expr.args[-1].cond)))
        tex = r"\begin{cases} %s \end{cases}"
        return tex % r" \\".join(ecpairs)

    def _print_MatrixBase(self, expr):
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
    _print_ImmutableMatrix = _print_MatrixBase
    _print_MutableMatrix = _print_MatrixBase

    def _print_BlockMatrix(self, expr):
        return self._print(expr.mat)

    def _print_Transpose(self, expr):
        mat = expr.arg
        if mat.is_Add or mat.is_Mul:
            return r"\left(%s\right)^T"%self._print(mat)
        else:
            return "%s^T"%self._print(mat)

    def _print_MatAdd(self, expr):
        return self._print_Add(expr)

    def _print_MatMul(self, expr):
        return self._print_Mul(expr)

    def _print_MatPow(self, expr):
        base, exp = expr.base, expr.exp
        if base.is_Add or base.is_Mul:
            return r"\left(%s\right)^{%s}"%(self._print(base), self._print(exp))
        else:
            return "%s^{%s}"%(self._print(base), self._print(exp))

    def _print_ZeroMatrix(self, Z):
        return r"\bold{0}"

    def _print_Identity(self, I):
        return r"\mathbb{I}"

    def _print_tuple(self, expr):
        return r"\begin{pmatrix}%s\end{pmatrix}" % \
            r", & ".join([ self._print(i) for i in expr ])

    def _print_Tuple(self, expr):
        return self._print_tuple(expr)

    def _print_list(self, expr):
        return r"\begin{bmatrix}%s\end{bmatrix}" % \
            r", & ".join([ self._print(i) for i in expr ])

    def _print_dict(self, d):
        keys = sorted(d.keys(), key=default_sort_key)
        items = []

        for key in keys:
            val = d[key]
            items.append("%s : %s" % (self._print(key), self._print(val)))

        return r"\begin{Bmatrix}%s\end{Bmatrix}" % r", & ".join(items)

    def _print_Dict(self, expr):
        return self._print_dict(expr)

    def _print_DiracDelta(self, expr):
        if len(expr.args) == 1 or expr.args[1] == 0:
            tex = r"\delta\left(%s\right)" % self._print(expr.args[0])
        else:
            tex = r"\delta^{\left( %s \right)}\left( %s \right)" % (\
            self._print(expr.args[1]), self._print(expr.args[0]))
        return tex

    def _print_ProductSet(self, p):
        if len(p.sets) > 1 and not has_variety(p.sets):
            return self._print(p.sets[0]) + "^%d"%len(p.sets)
        else:
            return r" \times ".join(self._print(set) for set in p.sets)

    def _print_RandomDomain(self, d):
        try:
            return 'Domain: '+ self._print(d.as_boolean())
        except:
            try:
                return ('Domain: ' + self._print(d.symbols) + ' in ' +
                        self._print(d.set))
            except:
                return 'Domain on ' + self._print(d.symbols)

    def _print_FiniteSet(self, s):
        items = sorted(s.args, key=default_sort_key)
        return self._print_set(items)

    def _print_set(self, s):
        items = sorted(s, key=default_sort_key)
        items = ", ".join(map(self._print, items))
        return r"\left\{%s\right\}" % items

    _print_frozenset = _print_set

    def _print_Range(self, s):
        if len(s) > 4:
            it = iter(s)
            printset = it.next(), it.next(), '\ldots', s._last_element
        else:
            printset = tuple(s)

        return (r"\left\{"
              + r", ".join(self._print(el) for el in printset)
              + r"\right\}")

    def _print_Interval(self, i):
        if i.start == i.end:
            return r"\left\{%s\right\}" % self._print(i.start)

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
        return r" \cup ".join([self._print(i) for i in u.args])

    def _print_Intersection(self, u):
        return r" \cap ".join([self._print(i) for i in u.args])

    def _print_EmptySet(self, e):
        return r"\emptyset"

    def _print_Naturals(self, n):
        return r"\mathbb{N}"

    def _print_Integers(self, i):
        return r"\mathbb{Z}"

    def _print_Reals(self, i):
        return r"\mathbb{R}"

    def _print_TransformationSet(self, s):
        return r"\left\{%s\; |\; %s \in %s\right\}"%(
                self._print(s.lamda.expr),
                ', '.join([self._print(var) for var in s.lamda.variables]),
                self._print(s.base_set))

    def _print_FiniteField(self, expr):
        return r"\mathbb{F}_{%s}" % expr.mod

    def _print_IntegerRing(self, expr):
        return r"\mathbb{Z}"

    def _print_RationalField(self, expr):
        return r"\mathbb{Q}"

    def _print_RealDomain(self, expr):
        return r"\mathbb{R}"

    def _print_ComplexDomain(self, expr):
        return r"\mathbb{C}"

    def _print_PolynomialRingBase(self, expr):
        domain = self._print(expr.dom)
        gens = ", ".join(map(self._print, expr.gens))
        inv = ""
        if not expr.is_Poly:
            inv = r"S_<^{-1}"
        return r"%s%s\left[%s\right]" % (inv, domain, gens)

    def _print_FractionField(self, expr):
        domain = self._print(expr.dom)
        gens = ", ".join(map(self._print, expr.gens))
        return r"%s\left(%s\right)" % (domain, gens)

    def _print_Poly(self, poly):
        cls = poly.__class__.__name__
        expr = self._print(poly.as_expr())
        gens = map(self._print, poly.gens)
        domain = "domain=%s" % self._print(poly.get_domain())

        args = ", ".join([expr] + gens + [domain])
        if cls in accepted_latex_functions:
            tex = r"\%s {\left (%s \right )}" % (cls, args)
        else:
            tex = r"\operatorname{%s}{\left( %s \right)}" % (cls, args)

        return tex

    def _print_RootOf(self, root):
        cls = root.__class__.__name__
        expr = self._print(root.expr)
        index = root.index
        if cls in accepted_latex_functions:
            return r"\%s {\left(%s, %d\right)}" % (cls, expr, index)
        else:
            return r"\operatorname{%s} {\left(%s, %d\right)}" % (cls, expr, index)


    def _print_RootSum(self, expr):
        cls = expr.__class__.__name__
        args = [self._print(expr.expr)]

        if expr.fun is not S.IdentityFunction:
            args.append(self._print(expr.fun))

        if cls in accepted_latex_functions:
            return r"\%s {\left(%s\right)}" % (cls, ", ".join(args))
        else:
            return r"\operatorname{%s} {\left(%s\right)}" % (cls, ", ".join(args))

    def _print_euler(self, expr):
        return r"E_{%s}" % self._print(expr.args[0])

    def _print_catalan(self, expr):
        return r"C_{%s}" % self._print(expr.args[0])

    def _print_MellinTransform(self, expr):
        return r"\mathcal{M}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseMellinTransform(self, expr):
        return r"\mathcal{M}^{-1}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_LaplaceTransform(self, expr):
        return r"\mathcal{L}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseLaplaceTransform(self, expr):
        return r"\mathcal{L}^{-1}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_FourierTransform(self, expr):
        return r"\mathcal{F}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseFourierTransform(self, expr):
        return r"\mathcal{F}^{-1}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_SineTransform(self, expr):
        return r"\mathcal{SIN}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseSineTransform(self, expr):
        return r"\mathcal{SIN}^{-1}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_CosineTransform(self, expr):
        return r"\mathcal{COS}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseCosineTransform(self, expr):
        return r"\mathcal{COS}^{-1}_{%s}\left[%s\right]\left(%s\right)" % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_DMP(self, p):
        try:
            if p.ring is not None:
                # TODO incorporate order
                return self._print(p.ring.to_sympy(p))
        except SympifyError:
            pass
        return self._print(repr(p))

    def _print_DMF(self, p):
        return self._print_DMP(p)

    def _print_Object(self, object):
        return self._print(Symbol(object.name))

    def _print_Morphism(self, morphism):
        domain = self._print(morphism.domain)
        codomain = self._print(morphism.codomain)
        return "%s\\rightarrow %s" % (domain, codomain)

    def _print_NamedMorphism(self, morphism):
        pretty_name = self._print(Symbol(morphism.name))
        pretty_morphism = self._print_Morphism(morphism)
        return "%s:%s" % (pretty_name, pretty_morphism)

    def _print_IdentityMorphism(self, morphism):
        from sympy.categories import NamedMorphism
        return self._print_NamedMorphism(NamedMorphism(
            morphism.domain, morphism.codomain, "id"))

    def _print_CompositeMorphism(self, morphism):
        from sympy.categories import NamedMorphism

        # All components of the morphism have names and it is thus
        # possible to build the name of the composite.
        component_names_list = [self._print(Symbol(component.name)) for \
                                component in morphism.components]
        component_names_list.reverse()
        component_names = "\\circ ".join(component_names_list) + ":"

        pretty_morphism = self._print_Morphism(morphism)
        return component_names + pretty_morphism

    def _print_Category(self, morphism):
        return "\\mathbf{%s}" % self._print(Symbol(morphism.name))

    def _print_Diagram(self, diagram):
        if not diagram.premises:
            # This is an empty diagram.
            return self._print(S.EmptySet)

        latex_result = self._print(diagram.premises)
        if diagram.conclusions:
            latex_result += "\\Longrightarrow %s" % \
                            self._print(diagram.conclusions)

        return latex_result

    def _print_DiagramGrid(self, grid):
        latex_result = "\\begin{array}{%s}\n" % ("c" * grid.width)

        for i in xrange(grid.height):
            for j in xrange(grid.width):
                if grid[i, j]:
                    latex_result += latex(grid[i, j])
                latex_result += " "
                if j != grid.width - 1:
                    latex_result += "& "

            if i != grid.height - 1:
                latex_result += "\\\\"
            latex_result += "\n"

        latex_result += "\\end{array}\n"
        return latex_result

    def _print_FreeModule(self, M):
        return '{%s}^{%s}' % (self._print(M.ring), self._print(M.rank))

    def _print_FreeModuleElement(self, m):
        # Print as row vector for convenience, for now.
        return r"\left[ %s \right]" % ",".join(
                '{' + self._print(x) + '}' for x in m)

    def _print_SubModule(self, m):
        return r"\left< %s \right>" % ",".join(
                '{' + self._print(x) + '}' for x in m.gens)

    def _print_ModuleImplementedIdeal(self, m):
        return r"\left< %s \right>" % ",".join(
                '{' + self._print(x) + '}' for [x] in m._module.gens)

    def _print_QuotientRing(self, R):
        # TODO nicer fractions for few generators...
        return r"\frac{%s}{%s}" % (self._print(R.ring), self._print(R.base_ideal))

    def _print_QuotientRingElement(self, x):
        return r"{%s} + {%s}" % (self._print(x.data), self._print(x.ring.base_ideal))

    def _print_QuotientModuleElement(self, m):
        return r"{%s} + {%s}" % (self._print(m.data),
                                 self._print(m.module.killed_module))

    def _print_QuotientModule(self, M):
        # TODO nicer fractions for few generators...
        return r"\frac{%s}{%s}" % (self._print(M.base),
                                   self._print(M.killed_module))

    def _print_MatrixHomomorphism(self, h):
        return r"{%s} : {%s} \to {%s}" % (self._print(h._sympy_matrix()),
            self._print(h.domain), self._print(h.codomain))

    def _print_BaseScalarField(self, field):
        string = field._coord_sys._names[field._index]
        return r'\boldsymbol{\mathrm{%s}}' % self._print(Symbol(string))

    def _print_BaseVectorField(self, field):
        string = field._coord_sys._names[field._index]
        return r'\partial_{%s}' % self._print(Symbol(string))

    def _print_Differential(self, diff):
        field = diff._form_field
        if hasattr(field, '_coord_sys'):
            string = field._coord_sys._names[field._index]
            return r'\mathbb{d}%s' % self._print(Symbol(string))
        else:
            return 'd(%s)'%self._print(field)
            string = self._print(field)
            return r'\mathbb{d}\left(%s\right)' % string

    def _print_Tr(self, p):
        #Todo: Handle indices
        contents = self._print(p.args[0])
        return r'\mbox{Tr}\left(%s\right)' % (contents)


def latex(expr, **settings):
    r"""
    Convert the given expression to LaTeX representation.

    >>> from sympy import latex, sin, asin, Matrix, Rational
    >>> from sympy.abc import x, y, mu, tau

    >>> latex((2*tau)**Rational(7,2))
    '8 \\sqrt{2} \\tau^{\\frac{7}{2}}'

    order: Any of the supported monomial orderings (currently "lex", "grlex", or
    "grevlex"), "old", and "none". This parameter does nothing for Mul objects.
    Setting order to "old" uses the compatibility ordering for Add defined in
    Printer. For very large expressions, set the 'order' keyword to 'none' if
    speed is a concern.

    mode: Specifies how the generated code will be delimited. 'mode' can be one
    of 'plain', 'inline', 'equation' or 'equation*'.  If 'mode' is set to
    'plain', then the resulting code will not be delimited at all (this is the
    default). If 'mode' is set to 'inline' then inline LaTeX $ $ will be used.
    If 'mode' is set to 'equation' or 'equation*', the resulting code will be
    enclosed in the 'equation' or 'equation*' environment (remember to import
    'amsmath' for 'equation*'), unless the 'itex' option is set. In the latter
    case, the ``$$ $$`` syntax is used.

    >>> latex((2*mu)**Rational(7,2), mode='plain')
    '8 \\sqrt{2} \\mu^{\\frac{7}{2}}'

    >>> latex((2*tau)**Rational(7,2), mode='inline')
    '$8 \\sqrt{2} \\tau^{\\frac{7}{2}}$'

    >>> latex((2*mu)**Rational(7,2), mode='equation*')
    '\\begin{equation*}8 \\sqrt{2} \\mu^{\\frac{7}{2}}\\end{equation*}'

    >>> latex((2*mu)**Rational(7,2), mode='equation')
    '\\begin{equation}8 \\sqrt{2} \\mu^{\\frac{7}{2}}\\end{equation}'

    itex: Specifies if itex-specific syntax is used, including emitting ``$$ $$``.

    >>> latex((2*mu)**Rational(7,2), mode='equation', itex=True)
    '$$8 \\sqrt{2} \\mu^{\\frac{7}{2}}$$'

    fold_frac_powers: Emit "^{p/q}" instead of "^{\frac{p}{q}}" for fractional
    powers.

    >>> latex((2*tau)**Rational(7,2), fold_frac_powers=True)
    '8 \\sqrt{2} \\tau^{7/2}'

    fold_func_brackets: Fold function brackets where applicable.

    >>> latex((2*tau)**sin(Rational(7,2)))
    '\\left(2 \\tau\\right)^{\\sin{\\left (\\frac{7}{2} \\right )}}'
    >>> latex((2*tau)**sin(Rational(7,2)), fold_func_brackets = True)
    '\\left(2 \\tau\\right)^{\\sin {\\frac{7}{2}}}'

    mul_symbol: The symbol to use for multiplication. Can be one of None,
    "ldot", "dot", or "times".

    >>> latex((2*tau)**sin(Rational(7,2)), mul_symbol="times")
    '\\left(2 \\times \\tau\\right)^{\\sin{\\left (\\frac{7}{2} \\right )}}'

    inv_trig_style: How inverse trig functions should be displayed. Can be one
    of "abbreviated", "full", or "power". Defaults to "abbreviated".

    >>> latex(asin(Rational(7,2)))
    '\\operatorname{asin}{\\left (\\frac{7}{2} \\right )}'
    >>> latex(asin(Rational(7,2)), inv_trig_style="full")
    '\\arcsin{\\left (\\frac{7}{2} \\right )}'
    >>> latex(asin(Rational(7,2)), inv_trig_style="power")
    '\\sin^{-1}{\\left (\\frac{7}{2} \\right )}'

    mat_str: Which matrix environment string to emit. "smallmatrix", "bmatrix",
    etc. Defaults to "smallmatrix".

    >>> latex(Matrix(2, 1, [x, y]), mat_str = "array")
    '\\left[\\begin{array}x\\\\y\\end{array}\\right]'

    mat_delim: The delimiter to wrap around matrices. Can be one of "[", "(",
    or the empty string. Defaults to "[".

    >>> latex(Matrix(2, 1, [x, y]), mat_delim="(")
    '\\left(\\begin{smallmatrix}x\\\\y\\end{smallmatrix}\\right)'

    symbol_names: Dictionary of symbols and the custom strings they should be
    emitted as.

    >>> latex(x**2, symbol_names={x:'x_i'})
    'x_i^{2}'

    Besides all Basic based expressions, you can recursively
    convert Python containers (lists, tuples and dicts) and
    also SymPy matrices:

    >>> latex([2/x, y], mode='inline')
    '$\\begin{bmatrix}\\frac{2}{x}, & y\\end{bmatrix}$'

    """

    return LatexPrinter(settings).doprint(expr)

def print_latex(expr, **settings):
    """Prints LaTeX representation of the given expression."""
    print latex(expr, **settings)
