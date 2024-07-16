"""
A Printer which converts an expression into its Typst equivalent.
"""
from __future__ import annotations
from typing import Any, Callable

from sympy.core import Add, Mod, Mul, Number, S, Symbol, Expr
from sympy.core.alphabets import greeks
from sympy.core.operations import AssocOp
from sympy.core.power import Pow

from sympy.printing.precedence import precedence_traditional
from sympy.printing.printer import Printer, print_function
from sympy.printing.conventions import split_super_sub
from sympy.printing.precedence import PRECEDENCE

from mpmath.libmp.libmpf import prec_to_dps, to_str as mlib_to_str

from sympy.utilities.iterables import sift

import re

typst_greek_dictionary = {
    'Alpha': 'Alpha',
    'Beta': 'Beta',
    'Gamma': 'Gamma',
    'Delta': 'Delta',
    'Epsilon': 'Epsilon',
    'Zeta': 'Zeta',
    'Eta': 'Eta',
    'Theta': 'Theta',
    'Iota': 'Iota',
    'Kappa': 'Kappa',
    'Lambda': 'Lambda',
    'Mu': 'Mu',
    'Nu': 'Nu',
    'Xi': 'Xi',
    'omicron': 'omicron',
    'Omicron': 'Omicron',
    'Pi': 'Pi',
    'Rho': 'Rho',
    'Sigma': 'Sigma',
    'Tau': 'Tau',
    'Upsilon': 'Upsilon',
    'Phi': 'Phi',
    'Chi': 'Chi',
    'Psi': 'Psi',
    'Omega': 'Omega',
    'lamda': 'lambda',
    'Lamda': 'Lambda',
    'khi': 'chi',
    'Khi': 'Chi',
    'varepsilon': 'epsilon',
    'epsilon': 'epsilon.alt',
    'varkappa': 'kappa',
    'varphi': 'phi',
    'phi': 'phi.alt',
    'varpi': 'pi.alt',
    'varrho': 'rho.alt',
    'varsigma': 'sigma.alt',
    'vartheta': 'theta.alt',

    'hbar': 'plank.reduce' # not greek, may need another dictionary
}

other_symbols = {'aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth', 'hslash',
                     'mho', 'wp'}

# Variable name modifiers
modifier_dict: dict[str, Callable[[str], str]] = {
    # Accents
    'mathring': lambda s: 'accent('+s+', circle)',
    'ddddot': lambda s: 'accent('+s+', dot.quad)',
    'dddot': lambda s: 'accent('+s+', dot.triple)',
    'ddot': lambda s: 'accent('+s+', dot.double)',
    'dot': lambda s: 'accent('+s+', dot)',
    'check': lambda s: 'accent('+s+', caron)',
    'breve': lambda s: 'accent('+s+', breve)',
    'acute': lambda s: 'accent('+s+', acute)',
    'grave': lambda s: 'accent('+s+', grave)',
    'tilde': lambda s: 'accent('+s+', tilde)',
    'hat': lambda s: 'accent('+s+', hat)',
    'bar': lambda s: 'accent('+s+', macron)',
    'vec': lambda s: 'accent('+s+', arrow)',
    'prime': lambda s: "("+s+")'",
    'prm': lambda s: "("+s+")'",
    # Faces
    'bold': lambda s: r'bold('+s+r')',
    'bm': lambda s: r'bold('+s+r')',
    'cal': lambda s: r'cal('+s+r')',
    'scr': lambda s: r'cal('+s+r')',
    'frak': lambda s: r'frak('+s+r')',
    # Brackets
    'norm': lambda s: r'norm('+s+r')',
    'avg': lambda s: r'lr(angle.l '+s+r' angle.r)', # "lr()" for size match
    'abs': lambda s: r'abs('+s+r')',
    'mag': lambda s: r'abs('+s+r')',
}

greek_letters_set = frozenset(greeks)

_between_two_numbers_p = (
    re.compile(r'[0-9]$'),  # search
    re.compile(r'(\d|\d+/\d+)'),  # match
)

class TypstPrinter(Printer):
    printmethod = "_typst"

    _default_settings: dict[str, Any] = {
        "full_prec": False,
        "fold_frac_powers": False,
        "fold_func_brackets": False,
        "fold_short_frac": None,
        "inv_trig_style": "abbreviated",
        "itypst": False,
        "ln_notation": False,
        "long_frac_ratio": None,
        "mat_delim": "[",
        "mat_str": None,
        "mode": "plain",
        "mul_symbol": None,
        "order": None,
        "symbol_names": {},
        "root_notation": True,
        "mat_symbol_style": "plain",
        "imaginary_unit": "i",
        "gothic_re_im": False,
        "decimal_separator": "period",
        "perm_cyclic": True,
        "parenthesize_super": True,
        "min": None,
        "max": None,
        "diff_operator": "d",
        "adjoint_style": "dagger",
    }

    def __init__(self, settings=None):
        Printer.__init__(self, settings)

        mul_symbol_table = {
            None: r" ",
            "dot": r" dot ",
            "times": r" times "
        }
        try:
            self._settings['mul_symbol_typst'] = \
                mul_symbol_table[self._settings['mul_symbol']]
        except KeyError:
            self._settings['mul_symbol_typst'] = \
                self._settings['mul_symbol']
        try:
            self._settings['mul_symbol_typst_numbers'] = \
                mul_symbol_table[self._settings['mul_symbol'] or 'dot']
        except KeyError:
            if (self._settings['mul_symbol'].strip() in
                    ['', ' ']):
                self._settings['mul_symbol_typst_numbers'] = \
                    mul_symbol_table['dot']
            else:
                self._settings['mul_symbol_typst_numbers'] = \
                    self._settings['mul_symbol']

        # TODO need to be checked when implementing settings
        imaginary_unit_table = {
            None: r"i",
            "i": r"i",
            "ri": r"\mathrm(i)",
            "ti": r"text(i)",
            "j": r"j",
            "rj": r"bold(j)",
            "tj": r"text(j)",
        }
        imag_unit = self._settings['imaginary_unit']
        self._settings['imaginary_unit_typst'] = imaginary_unit_table.get(imag_unit, imag_unit)


    def _add_parens(self, s) -> str:
        return r"({})".format(s)

    def parenthesize(self, item, level, is_neg=False, strict=False) -> str:
        prec_val = precedence_traditional(item)
        if is_neg and strict:
            return self._add_parens(self._print(item))

        if (prec_val < level) or ((not strict) and prec_val <= level):
            return self._add_parens(self._print(item))
        else:
            return self._print(item)

    def parenthesize_super(self, s):
        """
        Protect superscripts in s

        If the parenthesize_super option is set, protect with parentheses, else
        wrap in braces.
        """
        if "^" in s and self._settings['parenthesize_super']:
            return self._add_parens(s)
        return s

    def doprint(self, expr) -> str:
        typ = Printer.doprint(self, expr)

        return typ

    def _needs_add_brackets(self, expr) -> bool:
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed as part of an Add, False otherwise.  This is False for most
        things.
        """
        if expr.is_Relational:
            return True
        if any(expr.has(x) for x in (Mod,)):
            return True
        if expr.is_Add:
            return True
        return False

    def _needs_mul_brackets(self, expr, first=False, last=False) -> bool:
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed as part of a Mul, False otherwise. This is True for Add,
        but also for some container objects that would not need brackets
        when appearing last in a Mul, e.g. an Integral. ``last=True``
        specifies that this expr is the last to appear in a Mul.
        ``first=True`` specifies that this expr is the first to appear in
        a Mul.
        """
        from sympy.concrete.products import Product
        from sympy.concrete.summations import Sum
        from sympy.integrals.integrals import Integral

        if expr.is_Mul:
            if not first and expr.could_extract_minus_sign():
                return True
        elif precedence_traditional(expr) < PRECEDENCE["Mul"]:
            return True
        elif expr.is_Relational:
            return True
        if expr.is_Piecewise:
            return True
        if any(expr.has(x) for x in (Mod,)):
            return True
        if (not last and
                any(expr.has(x) for x in (Integral, Product, Sum))):
            return True

        return False

    def _print_Add(self, expr, order=None):
        terms = self._as_ordered_terms(expr, order=order)

        typ = ""
        for i, term in enumerate(terms):
            if i == 0:
                pass
            elif term.could_extract_minus_sign():
                typ += " - "
                term = -term
            else:
                typ += " + "
            term_typ = self._print(term)
            if self._needs_add_brackets(term):
                term_typ = r"(%s)" % term_typ
            typ += term_typ

        return typ

    def _print_Float(self, expr):
        # Based off of that in StrPrinter
        dps = prec_to_dps(expr._prec)
        strip = False if self._settings['full_prec'] else True
        low = self._settings["min"] if "min" in self._settings else None
        high = self._settings["max"] if "max" in self._settings else None
        str_real = mlib_to_str(expr._mpf_, dps, strip_zeros=strip, min_fixed=low, max_fixed=high)

        # Must always have a mul symbol (as 2.5 10^{20} just looks odd)
        # thus we use the number separator
        separator = self._settings['mul_symbol_typst_numbers']

        if 'e' in str_real:
            (mant, exp) = str_real.split('e')

            if exp[0] == '+':
                exp = exp[1:]
            if self._settings['decimal_separator'] == 'comma':
                mant = mant.replace('.','{,}')

            return r"%s%s10^(%s)" % (mant, separator, exp)
        elif str_real == "+inf":
            return r"infinity"
        elif str_real == "-inf":
            return r"- infinity"
        else:
            if self._settings['decimal_separator'] == 'comma':
                str_real = str_real.replace('.',',')
            return str_real

    def _print_Mul(self, expr: Expr):
        from sympy.simplify import fraction
        separator: str = self._settings['mul_symbol_typst']
        numbersep: str = self._settings['mul_symbol_typst_numbers']

        def convert(expr) -> str:
            if not expr.is_Mul:
                return str(self._print(expr))
            else:
                if self.order not in ('old', 'none'):
                    args = expr.as_ordered_factors()
                else:
                    args = list(expr.args)

                # If there are quantities or prefixes, append them at the back.
                units, nonunits = sift(args, lambda x: (hasattr(x, "_scale_factor") or hasattr(x, "is_physical_constant")) or
                              (isinstance(x, Pow) and
                               hasattr(x.base, "is_physical_constant")), binary=True)
                prefixes, units = sift(units, lambda x: hasattr(x, "_scale_factor"), binary=True)
                return convert_args(nonunits + prefixes + units)

        def convert_args(args) -> str:
            _typ = last_term_typ = ""

            for i, term in enumerate(args):
                term_typ = self._print(term)
                if not (hasattr(term, "_scale_factor") or hasattr(term, "is_physical_constant")):
                    if self._needs_mul_brackets(term, first=(i == 0),
                                                last=(i == len(args) - 1)):
                        term_typ = r"(%s)" % term_typ

                    if  _between_two_numbers_p[0].search(last_term_typ) and \
                        _between_two_numbers_p[1].match(term_typ):
                        # between two numbers
                        _typ += numbersep
                    elif _typ:
                        _typ += separator
                elif _typ:
                    _typ += separator

                _typ += term_typ
                last_term_typ = term_typ
            return _typ

        # Check for unevaluated Mul. In this case we need to make sure the
        # identities are visible, multiple Rational factors are not combined
        # etc so we display in a straight-forward form that fully preserves all
        # args and their order.
        # XXX: _print_Pow calls this routine with instances of Pow...
        if isinstance(expr, Mul):
            args = expr.args
            if args[0] is S.One or any(isinstance(arg, Number) for arg in args[1:]):
                return convert_args(args)

        include_parens = False
        if expr.could_extract_minus_sign():
            expr = -expr
            typ = "- "
            if expr.is_Add:
                typ += "("
                include_parens = True
        else:
            typ = ""

        numer, denom = fraction(expr, exact=True)

        if denom is S.One and Pow(1, -1, evaluate=False) not in expr.args:
            # use the original expression here, since fraction() may have
            # altered it when producing numer and denom
            typ += convert(expr)

        else:
            snumer = convert(numer)
            sdenom = convert(denom)
            ldenom = len(sdenom.split())
            ratio = self._settings['long_frac_ratio']
            if self._settings['fold_short_frac'] and ldenom <= 2 and \
                    "^" not in sdenom:
                # handle short fractions
                if self._needs_mul_brackets(numer, last=False):
                    typ += r"(%s)/%s" % (snumer, sdenom)
                else:
                    typ += r"%s/%s" % (snumer, sdenom)
            elif ratio is not None and \
                    len(snumer.split()) > ratio*ldenom:
                # handle long fractions
                if self._needs_mul_brackets(numer, last=True):
                    typ += r"1/%s %s(%s)" \
                        % (sdenom, separator, snumer)
                elif numer.is_Mul:
                    # split a long numerator
                    a = S.One
                    b = S.One
                    for x in numer.args:
                        if self._needs_mul_brackets(x, last=False) or \
                                len(convert(a*x).split()) > ratio*ldenom or \
                                (b.is_commutative is x.is_commutative is False):
                            b *= x
                        else:
                            a *= x
                    if self._needs_mul_brackets(b, last=True):
                        typ += r"%s/%s %s(%s)" \
                            % (convert(a), sdenom, separator, convert(b))
                    else:
                        typ += r"%s/%s %s%s" \
                            % (convert(a), sdenom, separator, convert(b))
                else:
                    typ += r"1/%s %s%s" % (sdenom, separator, snumer)
            else:
                if numer.is_number or numer.is_symbol:
                    typ += r"%s/%s" % (snumer, sdenom)
                else:
                    typ += r"(%s)/%s" % (snumer, sdenom)

        if include_parens:
            typ += ")"
        return typ

    def _print_Pow(self, expr: Pow):
        # Treat x**Rational(1,n) as special case
        if expr.exp.is_Rational:
            p: int = expr.exp.p  # type: ignore
            q: int = expr.exp.q  # type: ignore
            if abs(p) == 1 and q != 1 and self._settings['root_notation']:
                base = self._print(expr.base)
                if q == 2:
                    typ = r"sqrt(%s)" % base
                else:
                    typ = r"root(%d, %s)" % (q, base)
                if expr.exp.is_negative:
                    return r"1/%s" % typ
                else:
                    return typ
            elif self._settings['fold_frac_powers'] and q != 1:
                base = self.parenthesize(expr.base, PRECEDENCE['Pow'])
                # issue #12886: add parentheses for superscripts raised to powers
                if expr.base.is_Symbol:
                    base = self.parenthesize_super(base)
                if expr.base.is_Function:
                    return self._print(expr.base, exp="%s/%s" % (p, q))
                return r"%s^(%s/%s)" % (base, p, q)
            elif expr.exp.is_negative and expr.base.is_commutative:
                # special case for 1^(-x), issue 9216
                if expr.base == 1:
                    return r"%s^(%s)" % (expr.base, expr.exp)
                # special case for (1/x)^(-y) and (-1/-x)^(-y), issue 20252
                if expr.base.is_Rational:
                    base_p: int = expr.base.p  # type: ignore
                    base_q: int = expr.base.q  # type: ignore
                    if base_p * base_q == abs(base_q):
                        if expr.exp == -1:
                            return r"1/(%s/%s)" % (base_p, base_q)
                        else:
                            return r"1/(%s/%s)^(%s)" % (base_p, base_q, abs(expr.exp))
                # things like 1/x
                return self._print_Mul(expr)
        if expr.base.is_Function:
            return self._print(expr.base, exp=self._print(expr.exp))
        typ = r"%s^(%s)"
        return self._helper_print_standard_power(expr, typ)

    def _helper_print_standard_power(self, expr, template: str) -> str:
        exp = self._print(expr.exp)
        # issue #12886: add parentheses around superscripts raised
        # to powers
        base = self.parenthesize(expr.base, PRECEDENCE['Pow'])
        if expr.base.is_Symbol:
            base = self.parenthesize_super(base)
        elif expr.base.is_Float:
            base = r"(%s)" % base
        # elif (isinstance(expr.base, Derivative)
        #     and base.startswith(r'(')
        #     and re.match(r'(\d?d?dot', base)
        #     and base.endswith(r'\right)')):
        #     # don't use parentheses around dotted derivative
        #     base = base[6: -7]  # remove outermost added parens
        return template % (base, exp)

    def _print_Sum(self, expr):
        if len(expr.limits) == 1:
            typst = r"sum_(%s=%s)^(%s) " % \
                tuple([self._print(i) for i in expr.limits[0]])
        else:
            def _format_ineq(l):
                return r"%s <= %s <= %s" % \
                    tuple([self._print(s) for s in (l[1], l[0], l[2])])

            typst = r"sum_(%s) " % \
                str.join(" \\ ",[_format_ineq(l) for l in expr.limits])

        if isinstance(expr.function, Add):
            typst += r"(%s)" % self._print(expr.function)
        else:
            typst += self._print(expr.function)

        return typst

    def _print_Limit(self, expr):
        e, z, z0, dir = expr.args

        typst = r"lim_(%s -> " % self._print(z)
        if str(dir) == '+-' or z0 in (S.Infinity, S.NegativeInfinity):
            typst += r"%s)" % self._print(z0)
        else:
            typst += r"%s^%s)" % (self._print(z0), self._print(dir))

        if isinstance(e, AssocOp):
            return r"%s(%s)" % (typst, self._print(e))
        else:
            return r"%s %s" % (typst, self._print(e))

    def _print_log(self, expr, exp=None):
        if not self._settings["ln_notation"]:
            typst = r"log(%s)" % self._print(expr.args[0])
        else:
            typst = r"ln(%s)" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^(%s)" % (typst, exp)
        else:
            return typst

    def _print_Symbol(self, expr: Symbol, style='plain'):
        name: str = self._settings['symbol_names'].get(expr)
        if name is not None:
            return name

        return self._deal_with_super_sub(expr.name, style=style)

    _print_RandomSymbol = _print_Symbol

    def _deal_with_super_sub(self, string: str, style='plain') -> str:
        if '{' in string:
            name, supers, subs = string, [], []
        else:
            name, supers, subs = split_super_sub(string)

            name = translate(name)
            supers = [translate(sup) for sup in supers]
            subs = [translate(sub) for sub in subs]

        # apply the style only to the name
        if style == 'bold':
            name = "bold({})".format(name)

        # glue all items together:
        if supers:
            name += "^(%s)" % " ".join(supers)
        if subs:
            name += "_(%s)" % " ".join(subs)

        return name


def translate(s: str) -> str:
    r'''
    Check for a modifier ending the string.  If present, convert the
    modifier to typst and translate the rest recursively.

    Given a description of a Greek letter or other special character,
    return the appropriate typst.

    Let everything else pass as given.

    >>> from sympy.printing.typst import translate
    >>> translate('alphahatdotprime')
    "accent(accent(alpha, hat), dot)'"
    '''
    # Process the rest
    typst = typst_greek_dictionary.get(s)
    if typst:
        return typst
    elif s.lower() in greek_letters_set:
        return s.lower()
    elif s in other_symbols:
        return s
    else:
        # Process modifiers, if any, and recurse
        for key in sorted(modifier_dict.keys(), key=len, reverse=True):
            if s.lower().endswith(key) and len(s) > len(key):
                return modifier_dict[key](translate(s[:-len(key)]))
        return s


@print_function(TypstPrinter)
def typst(expr, **settings):
    r"""Convert the given expression to Typst string representation.
    """
    return TypstPrinter(settings).doprint(expr)
