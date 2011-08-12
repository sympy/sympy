"""
A Printer which converts an expression into its LaTeX equivalent.
"""

from sympy.core import S, C, Basic, Add, Mul, Wild, var
from printer import Printer
from conventions import split_super_sub
from sympy.simplify import fraction
from sympy import Interval

import sympy.mpmath.libmp as mlib
from sympy.mpmath.libmp import prec_to_dps

from sympy.core.compatibility import cmp_to_key

import re, warnings

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
        "mat_delim": "(",
    }

    def __init__(self, settings=None):
        if settings is not None and 'inline' in settings and not settings['inline']:
            # Change to "good" defaults for inline=False
            settings['mat_str'] = 'bmatrix'
            settings['mat_delim'] = None
        Printer.__init__(self, settings)

        if ('inline') in self._settings:
            warnings.warn("'inline' is deprecated, please use 'mode'. "
                "'mode' can be one of 'inline', 'plain', 'equation', or "
                "'equation*'.")
            if self._settings['inline']:
                self._settings['mode'] = 'inline'
            else:
                self._settings['mode'] = 'equation*'
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
        terms = self._as_ordered_terms(expr, order=order)
        tex = self._print(terms[0])

        for term in terms[1:]:
            coeff = term.as_coeff_mul()[0]

            if coeff >= 0:
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

                if self.order != 'old':
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
        # Treat x**(Rational(1,n)) as special case
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
        else:
            if expr.base.is_Function:
                return self._print(expr.base, self._print(expr.exp))
            else:
                if expr.is_commutative and expr.exp is S.NegativeOne:
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

    def _print_Integral(self, expr):
        tex, symbols = "", []

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

            symbols.insert(0, "d%s" % self._print(symbol))

        return r"%s %s\,%s" % (tex,
            str(self._print(expr.function)), " ".join(symbols))

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

    def _print_Lambda(self, expr):
        symbols, expr = expr.args

        if len(symbols) == 1:
            symbols = self._print(symbols[0])
        else:
            symbols = self._print(tuple(symbols))

        args = (symbols, self._print(expr))
        tex = r"\operatorname{\Lambda}\left(%s\right)" % ", ".join(args)

        return tex

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

    def _print_uppergamma(self, expr, exp=None):
        tex = r"\left(%s, %s\right)" % (self._print(expr.args[0]),
                                        self._print(expr.args[1]))

        if exp is not None:
            return r"\operatorname{\Gamma}^{%s}%s" % (exp, tex)
        else:
            return r"\operatorname{\Gamma}%s" % tex

    def _print_lowergamma(self, expr, exp=None):
        tex = r"\left(%s, %s\right)" % (self._print(expr.args[0]),
                                        self._print(expr.args[1]))

        if exp is not None:
            return r"\operatorname{\gamma}^{%s}%s" % (exp, tex)
        else:
            return r"\operatorname{\gamma}%s" % tex

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

    def _print_binomial(self, expr, exp=None):
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

    def _print_hyper(self, expr, exp=None):
        tex = r"{{}_{%s}F_{%s}\left.\left(\begin{matrix} %s \\ %s \end{matrix}" \
              r"\right| {%s} \right)}" % \
             (self._print(len(expr.ap)), self._print(len(expr.bq)),
              self._hprint_vec(expr.ap), self._hprint_vec(expr.bq),
              self._print(expr.argument))

        if exp is not None:
            tex = r"{%s}^{%s}" % (tex, self._print(exp))
        return tex

    def _print_meijerg(self, expr, exp=None):
        tex = r"{G_{%s, %s}^{%s, %s}\left.\left(\begin{matrix} %s & %s \\" \
              r"%s & %s \end{matrix} \right| {%s} \right)}" % \
             (self._print(len(expr.ap)), self._print(len(expr.bq)),
              self._print(len(expr.bm)), self._print(len(expr.an)),
              self._hprint_vec(expr.an), self._hprint_vec(expr.aother),
              self._hprint_vec(expr.bm), self._hprint_vec(expr.bother),
              self._print(expr.argument))

        if exp is not None:
            tex = r"{%s}^{%s}" % (tex, self._print(exp))
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
        if self._settings['itex']:
            lt = r"\lt"
        else:
            lt = "<"

        charmap = {
            "==" : "=",
            "<"  : lt,
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
                           (self._print(expr.args[-1].cond),
                            self._print(expr.args[-1].expr)))
        tex = r"\begin{cases} %s \end{cases}"
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
        keys.sort(key=cmp_to_key(Basic.compare_pretty))
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
    def _print_ProductSet(self, p):
        return r" \cross ".join(self._print(set) for set in p.sets)
    def _print_FiniteSet(self, s):
        if len(s) > 10:
            #take ten elements from the set at random
            q = iter(s)
            printset = [q.next() for i in xrange(10)]
        else:
            printset = s
        try:
            printset.sort()
        except:
            pass
        return r"\left\{" + r", ".join(self._print(el) for el in printset) + r"\right\}"
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
        return r" \cup ".join([self._print(i) for i in u.args])

    def _print_EmptySet(self, e):
        return r"\emptyset"

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

    def _print_PolynomialRing(self, expr):
        domain = self._print(expr.dom)
        gens = ", ".join(map(self._print, expr.gens))
        return r"%s\left\[%s\right\]" % (domain, gens)

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
        tex = r"\operatorname{%s}\left(%s\right)" % (cls, args)

        return tex

    def _print_RootOf(self, root):
        cls = root.__class__.__name__
        expr = self._print(root.expr)
        index = root.index
        return r"\operatorname{%s}\left(%s, %d\right)" % (cls, expr, index)

    def _print_RootSum(self, expr):
        cls = expr.__class__.__name__
        args = [self._print(expr.expr)]

        if expr.fun is not S.IdentityFunction:
            args.append(self._print(expr.fun))

        return r"\operatorname{%s}\left(%s\right)" % (cls, ", ".join(args))

def latex(expr, **settings):
    r"""Convert the given expression to LaTeX representation.

        You can specify how the generated code will be delimited using
        the 'mode' keyword. 'mode' can be one of 'plain', 'inline',
        'equation' or 'equation*'.  If 'mode' is set to 'plain', then
        the resulting code will not be delimited at all (this is the
        default). If 'mode' is set to 'inline' then inline LaTeX $ $ will be
        used.  If 'mode' is set to 'equation' or 'equation*', the resulting
        code will be enclosed in the 'equation' or 'equation*' environment
        (remember to import 'amsmath' for 'equation*'), unless the 'itex'
        option is set. In the latter case, the $$ $$ syntax is used.

        >>> from sympy import latex, Rational
        >>> from sympy.abc import x, y, mu, tau

        >>> latex((2*tau)**Rational(7,2))
        '8 \\sqrt{2} \\tau^{\\frac{7}{2}}'

        >>> latex((2*mu)**Rational(7,2), mode='plain')
        '8 \\sqrt{2} \\mu^{\\frac{7}{2}}'

        >>> latex((2*tau)**Rational(7,2), mode='inline')
        '$8 \\sqrt{2} \\tau^{\\frac{7}{2}}$'

        >>> latex((2*mu)**Rational(7,2), mode='equation*')
        '\\begin{equation*}8 \\sqrt{2} \\mu^{\\frac{7}{2}}\\end{equation*}'

        >>> latex((2*mu)**Rational(7,2), mode='equation')
        '\\begin{equation}8 \\sqrt{2} \\mu^{\\frac{7}{2}}\\end{equation}'

        >>> latex((2*mu)**Rational(7,2), mode='equation', itex=True)
        '$$8 \\sqrt{2} \\mu^{\\frac{7}{2}}$$'

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

