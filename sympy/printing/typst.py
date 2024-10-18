"""
A Printer which converts an expression into its Typst equivalent.
"""
from __future__ import annotations
from typing import Any, Callable, TYPE_CHECKING

import itertools

from sympy.core import Add, Mod, Mul, Number, S, Symbol, Expr
from sympy.core.alphabets import greeks
from sympy.core.function import Function, AppliedUndef
from sympy.core.operations import AssocOp
from sympy.core.power import Pow
from sympy.core.sorting import default_sort_key
from sympy.logic.boolalg import true, BooleanTrue, BooleanFalse

from sympy.printing.precedence import precedence_traditional
from sympy.printing.printer import Printer, print_function
from sympy.printing.conventions import split_super_sub, requires_partial
from sympy.printing.precedence import precedence, PRECEDENCE

from mpmath.libmp.libmpf import prec_to_dps, to_str as mlib_to_str

from sympy.utilities.iterables import sift

import re

if TYPE_CHECKING:
    from sympy.tensor.array import NDimArray
    from sympy.vector.basisdependent import BasisDependent

# need further check
accepted_typst_functions = ['arcsin', 'arccos', 'arctan', 'sin', 'cos', 'tan',
                            'sinh', 'cosh', 'tanh', 'sqrt', 'ln', 'log', 'sec',
                            'csc', 'cot', 'coth', 're', 'im', 'frac', 'root',
                            'arg',
                            ]

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

        self._delim_dict = {'(': ')', '[': ']', '{': '}'}

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

        diff_operator_table = {
            None: r"d",
            "d": r"d",
            "rd": r"bold(d)",
            "td": r"upright(d)",
        }
        diff_operator = self._settings['diff_operator']
        self._settings["diff_operator_typst"] = diff_operator_table.get(diff_operator, diff_operator)


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

    def _needs_brackets(self, expr) -> bool:
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed, False otherwise. For example: a + b => True; a => False;
        10 => False; -10 => True.
        """
        return not ((expr.is_Integer and expr.is_nonnegative)
                    or (expr.is_Atom and (expr is not S.NegativeOne
                                          and expr.is_Rational is False)))

    def _needs_function_brackets(self, expr) -> bool:
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

    def _mul_is_clean(self, expr) -> bool:
        for arg in expr.args:
            if arg.is_Function:
                return False
        return True

    def _pow_is_clean(self, expr) -> bool:
        return not self._needs_brackets(expr.base)

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

    def _do_exponent(self, expr: str, exp):
        if exp is not None:
            return r"(%s)^(%s)" % (expr, exp)
        else:
            return expr

    def _print_Basic(self, expr):
        name = self._deal_with_super_sub(expr.__class__.__name__)
        if expr.args:
            ls = [self._print(o) for o in expr.args]
            s = r'upright("{}")("{}")'
            return s.format(name, ", ".join(ls))
        else:
            return r'upright("{}")'.format(name)

    def _print_bool(self, e: bool | BooleanTrue | BooleanFalse):
        return r'upright("%s")' % e

    _print_BooleanTrue = _print_bool
    _print_BooleanFalse = _print_bool

    def _print_NoneType(self, e):
        return r'upright("%s")' % e

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

    def _print_Cycle(self, expr):
        from sympy.combinatorics.permutations import Permutation
        if expr.size == 0:
            return r"()"
        expr = Permutation(expr)
        expr_perm = expr.cyclic_form
        siz = expr.size
        if expr.array_form[-1] == siz - 1:
            expr_perm = expr_perm + [[siz - 1]]
        term_tex = ''
        for i in expr_perm:
            term_tex += str(i).replace(',', r" space")
        term_tex = term_tex.replace('[', r"(")
        term_tex = term_tex.replace(']', r")")
        return term_tex

    def _print_Permutation(self, expr):
        from sympy.combinatorics.permutations import Permutation
        from sympy.utilities.exceptions import sympy_deprecation_warning

        perm_cyclic = Permutation.print_cyclic
        if perm_cyclic is not None:
            sympy_deprecation_warning(
                f"""
                Setting Permutation.print_cyclic is deprecated. Instead use
                init_printing(perm_cyclic={perm_cyclic}).
                """,
                deprecated_since_version="1.6",
                active_deprecations_target="deprecated-permutation-print_cyclic",
                stacklevel=8,
            )
        else:
            perm_cyclic = self._settings.get("perm_cyclic", True)

        if perm_cyclic:
            return self._print_Cycle(expr)

        if expr.size == 0:
            return r"()"

        lower = [self._print(arg) for arg in expr.array_form]
        upper = [self._print(arg) for arg in range(len(lower))]

        row1 = ", ".join(upper)
        row2 = ", ".join(lower)
        mat = r"; ".join((row1, row2))
        return r"mat(%s)" % mat

    def _print_AppliedPermutation(self, expr):
        perm, var = expr.args
        return r"sigma_(%s)(%s)" % (self._print(perm), self._print(var))

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

    def _print_Cross(self, expr):
        vec1 = expr._expr1
        vec2 = expr._expr2
        return r"%s times %s" % (self.parenthesize(vec1, PRECEDENCE['Mul']),
                                  self.parenthesize(vec2, PRECEDENCE['Mul']))

    def _print_Curl(self, expr):
        vec = expr._expr
        return r"nabla times %s" % self.parenthesize(vec, PRECEDENCE['Mul'])

    def _print_Divergence(self, expr):
        vec = expr._expr
        return r"nabla dot %s" % self.parenthesize(vec, PRECEDENCE['Mul'])

    def _print_Dot(self, expr):
        vec1 = expr._expr1
        vec2 = expr._expr2
        return r"%s dot %s" % (self.parenthesize(vec1, PRECEDENCE['Mul']),
                                 self.parenthesize(vec2, PRECEDENCE['Mul']))

    def _print_Gradient(self, expr):
        func = expr._expr
        return r"nabla %s" % self.parenthesize(func, PRECEDENCE['Mul'])

    def _print_Laplacian(self, expr):
        func = expr._expr
        return r"Delta %s" % self.parenthesize(func, PRECEDENCE['Mul'])

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
                if numer.is_number or numer.is_symbol or numer.is_Pow:
                    if denom.is_number or denom.is_symbol or denom.is_Pow:
                        typ += r"%s/%s" % (snumer, sdenom)
                    else:
                        typ += r"%s/(%s)" % (snumer, sdenom)
                else:
                    typ += r"(%s)/%s" % (snumer, sdenom)

        if include_parens:
            typ += ")"
        return typ

    def _print_AlgebraicNumber(self, expr):
        if expr.is_aliased:
            return self._print(expr.as_poly().as_expr())
        else:
            return self._print(expr.as_expr())

    def _print_PrimeIdeal(self, expr):
        p = self._print(expr.p)
        if expr.is_inert:
            return rf'({p})'
        alpha = self._print(expr.alpha.as_expr())
        return rf'({p}, {alpha})'

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
        #     and base.endswith(r')')):
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

    def _print_BasisDependent(self, expr: 'BasisDependent'):
        from sympy.vector import Vector

        o1: list[str] = []
        if expr == expr.zero:
            return expr.zero._typst_form
        if isinstance(expr, Vector):
            items = expr.separate().items()
        else:
            items = [(0, expr)]

        for system, vect in items:
            inneritems = list(vect.components.items())
            inneritems.sort(key=lambda x: x[0].__str__())
            for k, v in inneritems:
                if v == 1:
                    o1.append(' + ' + k._typst_form)
                elif v == -1:
                    o1.append(' - ' + k._typst_form)
                else:
                    arg_str = r'(' + self._print(v) + r')'
                    o1.append(' + ' + arg_str + k._typst_form)

        outstr = (''.join(o1))
        if outstr[1] != '-':
            outstr = outstr[3:]
        else:
            outstr = outstr[1:]
        return outstr

    def _print_CoordSystem(self, coordsys):
        return r'upright(%s)^(upright(%s))_(%s)' % (
            self._print(coordsys.name), self._print(coordsys.patch.name), self._print(coordsys.manifold)
        )

    def _print_CovarDerivativeOp(self, cvd):
        return r'nabla_(%s)' % self._print(cvd._wrt)

    def _print_BaseScalarField(self, field):
        string = field._coord_sys.symbols[field._index].name
        return r'bold({})'.format(self._print(Symbol(string)))

    def _print_BaseVectorField(self, field):
        string = field._coord_sys.symbols[field._index].name
        return r'partial_({})'.format(self._print(Symbol(string)))

    def _print_Derivative(self, expr):
        if requires_partial(expr.expr):
            diff_symbol = r'partial'
        else:
            diff_symbol = self._settings["diff_operator_typst"]

        typst = ""
        dim = 0
        for i, (x, num) in enumerate(reversed(expr.variable_count)):
            # handle mixed diff
            if i == 0:
                separator = ""
            else:
                separator = " "
            typst += separator

            dim += num
            if num == 1:
                typst += r"%s %s" % (diff_symbol, self._print(x))
            else:
                typst += r"%s %s^(%s)" % (diff_symbol,
                                        self.parenthesize_super(self._print(x)),
                                        self._print(num))

        if dim == 1:
            typst = r"%s / (%s)" % (diff_symbol, typst)
        else:
            typst = r"%s^(%s) / (%s)" % (diff_symbol, self._print(dim), typst)

        if any(i.could_extract_minus_sign() for i in expr.args):
            return r"%s %s" % (typst, self.parenthesize(expr.expr,
                                                  PRECEDENCE["Mul"],
                                                  is_neg=True,
                                                  strict=True))

        return r"%s %s" % (typst, self.parenthesize(expr.expr,
                                                  PRECEDENCE["Mul"],
                                                  is_neg=False,
                                                  strict=True))

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

    def _hprint_Function(self, func: str) -> str:
        r'''
        Logic to decide how to render a function to typst
          - if it is a recognized typst name, use the appropriate typst function
          - if it is a single letter, excluding sub- and superscripts, just use that letter
          - if it is a longer name, then put upright() around it and be
            mindful of undercores in the name
        '''
        func = self._deal_with_super_sub(func)
        superscriptidx = func.find("^")
        subscriptidx = func.find("_")
        if func in accepted_typst_functions:
            name = r"%s" % func
        elif len(func) == 1 or func.startswith('\\') or subscriptidx == 1 or superscriptidx == 1:
            name = func
        else:
            if superscriptidx > 0 and subscriptidx > 0:
                name = r"upright(%s)%s" %(
                    func[:min(subscriptidx,superscriptidx)],
                    func[min(subscriptidx,superscriptidx):])
            elif superscriptidx > 0:
                name = r"upright(%s)%s" %(
                    func[:superscriptidx],
                    func[superscriptidx:])
            elif subscriptidx > 0:
                name = r"upright(%s)%s" %(
                    func[:subscriptidx],
                    func[subscriptidx:])
            else:
                name = r"upright(%s)" % func
        return name

    def _print_Function(self, expr: Function, exp=None) -> str:
        r'''
        Render functions to Typst, handling functions that typst knows about
        e.g., sin, cos, ... by using the proper typst function (sin, cos, ...).
        For single-letter function names, render them as regular Typst math
        symbols. For multi-letter function names that Typst does not know
        about, (e.g., Li, sech) use upright() so that the function name
        is rendered in Roman font and Typst handles spacing properly.

        expr is the expression involving the function
        exp is an exponent
        '''
        func = expr.func.__name__
        if hasattr(self, '_print_' + func) and \
                not isinstance(expr, AppliedUndef):
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

            inv_trig_table = [
                "asin", "acos", "atan",
                "acsc", "asec", "acot",
                "asinh", "acosh", "atanh",
                "acsch", "asech", "acoth",
            ]

            # If the function is an inverse trig function, handle the style
            if func in inv_trig_table:
                if inv_trig_style == "abbreviated":
                    pass
                elif inv_trig_style == "full":
                    func = ("ar" if func[-1] == "h" else "arc") + func[1:]
                elif inv_trig_style == "power":
                    func = func[1:]
                    inv_trig_power_case = True

                    # Can never fold brackets if we're raised to a power
                    if exp is not None:
                        can_fold_brackets = False

            if inv_trig_power_case:
                if func in accepted_typst_functions:
                    name = r"%s^(-1)" % func
                else:
                    name = r"upright(%s)^(-1)" % func
            elif exp is not None:
                func_tex = self._hprint_Function(func)
                func_tex = self.parenthesize_super(func_tex)
                name = r'%s^(%s)' % (func_tex, exp)
            else:
                name = self._hprint_Function(func)

            if can_fold_brackets:
                if func in accepted_typst_functions:
                    # Wrap argument safely to avoid parse-time conflicts
                    # with the function name itself
                    name += r" (%s)"
                else:
                    name += r"%s"
            else:
                name += r"(%s)"

            if inv_trig_power_case and exp is not None:
                name += r"^(%s)" % exp

            return name % ",".join(args)

    def _print_UndefinedFunction(self, expr):
        return self._hprint_Function(str(expr))

    def _print_ElementwiseApplyFunction(self, expr):
        return r"%s_circle.stroked.small (%s)" % (
            self._print(expr.function),
            self._print(expr.expr),
        )

    @property
    def _special_function_classes(self):
        from sympy.functions.special.tensor_functions import KroneckerDelta
        from sympy.functions.special.gamma_functions import gamma, lowergamma
        from sympy.functions.special.beta_functions import beta
        from sympy.functions.special.delta_functions import DiracDelta
        from sympy.functions.special.error_functions import Chi
        return {KroneckerDelta: r'delta',
                gamma:  r'Gamma',
                lowergamma: r'gamma',
                beta: r'Beta',
                DiracDelta: r'delta',
                Chi: r'Chi'}

    def _print_FunctionClass(self, expr):
        for cls in self._special_function_classes:
            if issubclass(expr, cls) and expr.__name__ == cls.__name__:
                return self._special_function_classes[cls]
        return self._hprint_Function(str(expr))

    def _print_Lambda(self, expr):
        symbols, expr = expr.args

        if len(symbols) == 1:
            symbols = self._print(symbols[0])
        else:
            symbols = self._print(tuple(symbols))

        tex = r"( %s arrow.r.bar %s )" % (symbols, self._print(expr))

        return tex

    def _print_IdentityFunction(self, expr):
        return r"( x arrow.r.bar x )"

    def _hprint_variadic_function(self, expr, exp=None) -> str:
        args = sorted(expr.args, key=default_sort_key)
        texargs = [r"%s" % self._print(symbol) for symbol in args]
        tex = r"%s(%s)" % (str(expr.func).lower(),
                                       ", ".join(texargs))
        if exp is not None:
            return r"%s^(%s)" % (tex, exp)
        else:
            return tex

    _print_Min = _print_Max = _hprint_variadic_function

    def _print_floor(self, expr, exp=None):
        typst = r"floor(%s)" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^(%s)" % (typst, exp)
        else:
            return typst

    def _print_ceiling(self, expr, exp=None):
        typst = r"ceil(%s)" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^(%s)" % (typst, exp)
        else:
            return typst

    def _print_log(self, expr, exp=None):
        if not self._settings["ln_notation"]:
            typst = r"log(%s)" % self._print(expr.args[0])
        else:
            typst = r"ln(%s)" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^(%s)" % (typst, exp)
        else:
            return typst

    def _print_Abs(self, expr, exp=None):
        typst = r"abs(%s)" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^(%s)" % (typst, exp)
        else:
            return typst

    def _print_re(self, expr, exp=None):
        if self._settings['gothic_re_im']:
            tex = r"Re(%s)" % expr.args[0]
        else:
            tex = r'upright("re")({})'.format(expr.args[0])

        return self._do_exponent(tex, exp)

    def _print_im(self, expr, exp=None):
        if self._settings['gothic_re_im']:
            tex = r"Im(%s)" % expr.args[0]
        else:
            tex = r'upright("im")({})'.format(expr.args[0])

        return self._do_exponent(tex, exp)

    def _print_Not(self, e):
        from sympy.logic.boolalg import (Equivalent, Implies)
        if isinstance(e.args[0], Equivalent):
            return self._print_Equivalent(e.args[0], r"arrow.l.r.double.not")
        if isinstance(e.args[0], Implies):
            return self._print_Implies(e.args[0], r"arrow.r.double.not")
        if (e.args[0].is_Boolean):
            return r"not (%s)" % self._print(e.args[0])
        else:
            return r"not %s" % self._print(e.args[0])

    def _print_LogOp(self, args, char):
        arg = args[0]
        if arg.is_Boolean and not arg.is_Not:
            tex = r"(%s)" % self._print(arg)
        else:
            tex = r"%s" % self._print(arg)

        for arg in args[1:]:
            if arg.is_Boolean and not arg.is_Not:
                tex += r" %s (%s)" % (char, self._print(arg))
            else:
                tex += r" %s %s" % (char, self._print(arg))

        return tex

    def _print_And(self, e):
        args = sorted(e.args, key=default_sort_key)
        return self._print_LogOp(args, r"and")

    def _print_Or(self, e):
        args = sorted(e.args, key=default_sort_key)
        return self._print_LogOp(args, r"or")

    def _print_Xor(self, e):
        args = sorted(e.args, key=default_sort_key)
        return self._print_LogOp(args, r"\u{22BB}")

    def _print_Implies(self, e, altchar=None):
        return self._print_LogOp(e.args, altchar or r"=>")

    def _print_Equivalent(self, e, altchar=None):
        args = sorted(e.args, key=default_sort_key)
        return self._print_LogOp(args, altchar or r"<=>")

    def _print_conjugate(self, expr, exp=None):
        typ = r"overline(%s)" % self._print(expr.args[0])

        if exp is not None:
            return r"%s^(%s)" % (typ, exp)
        else:
            return typ

    def _print_polar_lift(self, expr, exp=None):
        func = r'upright("polar_lift")'
        arg = r'(%s)' % self._print(expr.args[0])

        if exp is not None:
            return r'%s^(%s)%s' % (func, exp, arg)
        else:
            return r'%s%s' % (func, arg)

    def _print_ExpBase(self, expr, exp=None):
        typst = r"e^(%s)" % self._print(expr.args[0])
        return self._do_exponent(typst, exp)

    def _print_Exp1(self, expr, exp=None):
        return "e"

    def _print_Rational(self, expr):
        if expr.q != 1:
            sign = ""
            p = expr.p
            if expr.p < 0:
                sign = "- "
                p = -p
            if self._settings['fold_short_frac']:
                return r"%s%d / %d" % (sign, p, expr.q)
            return r"%s%d/%d" % (sign, p, expr.q)
        else:
            return self._print(expr.p)

    def _print_Piecewise(self, expr):
        ecpairs = [r'%s & upright("for") %s' % (self._print(e), self._print(c))
                   for e, c in expr.args[:-1]]
        if expr.args[-1].cond == true:
            ecpairs.append(r'%s & upright("otherwise")' %
                           self._print(expr.args[-1].expr))
        else:
            ecpairs.append(r'%s & upright("for") %s' %
                           (self._print(expr.args[-1].expr),
                            self._print(expr.args[-1].cond)))
        typ = r'cases( %s )'
        return typ % r', '.join(ecpairs)

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

    def _print_Relational(self, expr):
        charmap = {
            "==": "=",
            ">": ">",
            "<": "<",
            ">=": ">=",
            "<=": "<=",
            "!=": "!=",
        }

        return "%s %s %s" % (self._print(expr.lhs),
                             charmap[expr.rel_op], self._print(expr.rhs))

    def _print_matrix_contents(self, expr):
        lines = []

        for line in range(expr.rows):  # horrible, should be 'rows'
            lines.append(", ".join([self._print(i) for i in expr[line, :]]))

        out_str = r'%s'
        return out_str % r"; ".join(lines)

    def _print_MatrixBase(self, expr):
        out_str = self._print_matrix_contents(expr)
        if self._settings['mat_delim']:
            left_delim = self._settings['mat_delim']
            out_str = r'mat(delim: "' + left_delim + '", ' + \
                        out_str + ")"
        else:
            out_str = r'mat(delim: #none, ' + out_str + ')'
        return out_str

    def _print_MatrixElement(self, expr):
        matrix_part = self.parenthesize(expr.parent, PRECEDENCE['Atom'], strict=True)
        index_part = f"{self._print(expr.i)},{self._print(expr.j)}"
        return f"{{{matrix_part}}}_{{{index_part}}}"

    def _print_MatrixSlice(self, expr):
        def latexslice(x, dim):
            x = list(x)
            if x[2] == 1:
                del x[2]
            if x[0] == 0:
                x[0] = None
            if x[1] == dim:
                x[1] = None
            return ':'.join(self._print(xi) if xi is not None else '' for xi in x)
        return (self.parenthesize(expr.parent, PRECEDENCE["Atom"], strict=True) + r'[' +
                latexslice(expr.rowslice, expr.parent.rows) + ', ' +
                latexslice(expr.colslice, expr.parent.cols) + r']')

    def _print_BlockMatrix(self, expr):
        return self._print(expr.blocks)

    def _print_Transpose(self, expr):
        mat = expr.arg
        from sympy.matrices import MatrixSymbol, BlockMatrix
        if (not isinstance(mat, MatrixSymbol) and
            not isinstance(mat, BlockMatrix) and mat.is_MatrixExpr):
            return r"(%s)^(T)" % self._print(mat)
        else:
            s = self.parenthesize(mat, precedence_traditional(expr), True)
            if '^' in s:
                return r"(%s)^(T)" % s
            else:
                return "%s^(T)" % s

    def _print_Trace(self, expr):
        mat = expr.arg
        return r'upright("tr")(%s )' % self._print(mat)

    def _print_Adjoint(self, expr):
        style_to_latex = {
            "dagger"   : r"dagger",
            "star"     : r"ast",
            "hermitian": r"sans(upright(H))"
        }
        adjoint_style = style_to_latex.get(self._settings["adjoint_style"], r"dagger")
        mat = expr.arg
        from sympy.matrices import MatrixSymbol, BlockMatrix
        if (not isinstance(mat, MatrixSymbol) and
            not isinstance(mat, BlockMatrix) and mat.is_MatrixExpr):
            return r"(%s)^(%s)" % (self._print(mat), adjoint_style)
        else:
            s = self.parenthesize(mat, precedence_traditional(expr), True)
            if '^' in s:
                return r"(%s)^(%s)" % (s, adjoint_style)
            else:
                return r"%s^(%s)" % (s, adjoint_style)

    def _print_MatMul(self, expr):
        from sympy import MatMul

        # Parenthesize nested MatMul but not other types of Mul objects:
        parens = lambda x: self._print(x) if isinstance(x, Mul) and not isinstance(x, MatMul) else \
            self.parenthesize(x, precedence_traditional(expr), False)

        args = list(expr.args)
        if expr.could_extract_minus_sign():
            if args[0] == -1:
                args = args[1:]
            else:
                args[0] = -args[0]
            return '- ' + ' '.join(map(parens, args))
        else:
            return ' '.join(map(parens, args))

    def _print_Determinant(self, expr):
        mat = expr.arg
        if mat.is_MatrixExpr:
            from sympy.matrices.expressions.blockmatrix import BlockMatrix
            if isinstance(mat, BlockMatrix):
                return r"|(%s)|" % self._print_matrix_contents(mat.blocks)
            return r"|(%s)|" % self._print(mat)
        return r"|(%s)|" % self._print_matrix_contents(mat)

    def _print_MatPow(self, expr):
        base, exp = expr.base, expr.exp
        from sympy.matrices import MatrixSymbol
        if not isinstance(base, MatrixSymbol) and base.is_MatrixExpr:
            return "(%s)^(%s)" % (self._print(base),
                                              self._print(exp))
        else:
            base_str = self._print(base)
            if '^' in base_str:
                return r"(%s)^(%s)" % (base_str, self._print(exp))
            else:
                return "%s^(%s)" % (base_str, self._print(exp))

    def _print_MatrixSymbol(self, expr):
        return self._print_Symbol(expr, style=self._settings[
            'mat_symbol_style'])

    def _print_ZeroMatrix(self, Z):
        return "0" if self._settings[
            'mat_symbol_style'] == 'plain' else r"bold(0)"

    def _print_OneMatrix(self, O):
        return "1" if self._settings[
            'mat_symbol_style'] == 'plain' else r"bold(1)"

    def _print_Identity(self, I):
        return r"II" if self._settings[
            'mat_symbol_style'] == 'plain' else r"bold(I)"

    def _print_PermutationMatrix(self, P):
        perm_str = self._print(P.args[0])
        return "P_(%s)" % perm_str

    def _print_NDimArray(self, expr: NDimArray):

        if expr.rank() == 0:
            return self._print(expr[()])


        if self._settings['mat_delim']:
            left_delim: str = self._settings['mat_delim']
            block_str = r'mat(delim: "' + left_delim + '", ' + \
                         "%s )"
        else:
            block_str = out_str = r'mat(delim: #none, ' + '%s )'

        if expr.rank() == 0:
            return block_str % ""

        level_str: list[list[str]] = [[] for i in range(expr.rank() + 1)]
        shape_ranges = [list(range(i)) for i in expr.shape]
        for outer_i in itertools.product(*shape_ranges):
            level_str[-1].append(self._print(expr[outer_i]))
            even = True
            for back_outer_i in range(expr.rank()-1, -1, -1):
                if len(level_str[back_outer_i+1]) < expr.shape[back_outer_i]:
                    break
                if even:
                    level_str[back_outer_i].append(
                        r", ".join(level_str[back_outer_i+1]))
                else:
                    level_str[back_outer_i].append(
                        block_str % (r"; ".join(level_str[back_outer_i+1])))
                    if len(level_str[back_outer_i+1]) == 1:
                        level_str[back_outer_i][-1] = r"[" + \
                            level_str[back_outer_i][-1] + r"]"
                even = not even
                level_str[back_outer_i+1] = []

        out_str = level_str[0][0]

        if expr.rank() % 2 == 1:
            out_str = block_str % out_str

        return out_str

    def _printer_tensor_indices(self, name, indices, index_map: dict):
        out_str = self._print(name)
        last_valence = None
        prev_map = None
        for index in indices:
            new_valence = index.is_up
            if index.is_up:
                out_str += ' scripts("")^('
            else:
                out_str += ' scripts("")_('

            # add "," if needed
            if ((index in index_map) or prev_map) and \
                    last_valence == new_valence:
                out_str += ","

            out_str += self._print(index.args[0])
            if index in index_map:
                out_str += "="
                out_str += self._print(index_map[index])
                prev_map = True
            else:
                prev_map = False
            last_valence = new_valence
            out_str += ')'
        return out_str

    def _print_Tensor(self, expr):
        name = expr.args[0].args[0]
        indices = expr.get_indices()
        return self._printer_tensor_indices(name, indices, {})

    def _print_TensorElement(self, expr):
        name = expr.expr.args[0].args[0]
        indices = expr.expr.get_indices()
        index_map = expr.index_map
        return self._printer_tensor_indices(name, indices, index_map)

    def _print_TensMul(self, expr):
        # prints expressions like "A(a)", "3*A(a)", "(1+x)*A(a)"
        sign, args = expr._get_args_for_traditional_printer()
        return sign + "".join(
            [self.parenthesize(arg, precedence(expr)) for arg in args]
        )

    def _print_TensAdd(self, expr):
        a = []
        args = expr.args
        for x in args:
            a.append(self.parenthesize(x, precedence(expr)))
        a.sort()
        s = ' + '.join(a)
        s = s.replace('+ -', '- ')
        return s

    def _print_TensorIndex(self, expr):
        return 'scripts("")%s(%s)' % (
            "^" if expr.is_up else "_",
            self._print(expr.args[0])
        )

    def _print_PartialDerivative(self, expr):
        if len(expr.variables) == 1:
            return r"partial / (partial %s) %s" % (
                self._print(expr.variables[0]),
                self.parenthesize(expr.expr, PRECEDENCE["Mul"], False)
            )
        else:
            return r"partial^(%s) / (%s) %s" % (
                len(expr.variables),
                " ".join([r"partial %s" % self._print(i) for i in expr.variables]),
                self.parenthesize(expr.expr, PRECEDENCE["Mul"], False)
            )

    def _print_ArraySymbol(self, expr):
        return self._print(expr.name)

    def _print_ArrayElement(self, expr):
        return r'%s_(%s)' % (
            self.parenthesize(expr.name, PRECEDENCE["Func"], True),
            ", ".join([f"{self._print(i)}" for i in expr.indices]))


def translate(s: str) -> str:
    r'''
    Check for a modifier ending the string.  If present, convert the
    modifier to typst and translate the rest recursively.

    Given a description of a Greek letter or other special character,
    return the appropriate typst.

    Let everything else pass as given.

    >>> from sympy.printing.typst import translate
    >>> translate('alphahatdotprime')
    "(accent(accent(alpha, hat), dot))'"
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
