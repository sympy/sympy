from sympy.core import Basic
from printer import Printer
from stringpict import *


class PrettyPrinter(Printer):
    """
    A class that prints a prettified expression, one that is not limited
    to one dimension like casting the expression to a string would return.
    """
    def __init__(self, use_unicode=False):
        Printer.__init__(self)
        self._use_unicode = use_unicode
        if self._use_unicode:
            self._str = unicode
        self.emptyPrinter = lambda x : prettyForm(self._str(x), prettyForm.ATOM)

    def doprint(self, expr):
        Printer.doprint.__doc__
        return self._print(expr).terminal_string()

    def _print_Symbol(self, e):
        if self._use_unicode:
            # letter: (uppercase, lowercase)
            greek = {
                'alpha': (u'\u0391', u'\u03b1'),
                'beta': (u'\u0392', u'\u03b2'),
                'gamma': (u'\u0393', u'\u03b3'),
                'delta': (u'\u0394', u'\u03b4'),
                'epsilon': (u'\u0395', u'\u03b5'),
                'zeta': (u'\u0396', u'\u03b6'),
                'eta': (u'\u0397', u'\u03b7'),
                'theta': (u'\u0398', u'\u03b8'),
                'iota': (u'\u0399', u'\u03b9'),
                'kappa': (u'\u039a', u'\u03ba'),
                'lambda': (u'\u039b', u'\u03bb'),
                'mu': (u'\u039c', u'\u03bc'),
                'nu': (u'\u039d', u'\u03bd'),
                'xi': (u'\u039e', u'\u03be'),
                'omicron': (u'\u039f', u'\u03bf'),
                'pi': (u'\u03a0', u'\u03c0'),
                'rho': (u'\u03a1', u'\u03c1'),
                'sigma': (u'\u03a3', u'\u03c3'),
                'tau': (u'\u03a4', u'\u03c4'),
                'upsilon': (u'\u03a5', u'\u03c5'),
                'phi': (u'\u03a6', u'\u03c6'),
                'chi': (u'\u03a7', u'\u03c7'),
                'psi': (u'\u03a8', u'\u03c8'),
                'omega': (u'\u03a9', u'\u03c9')
            }

            name = e.name.lower()
            if name in greek:
                # If first character lowercase, use lowercase greek letter
                if name[0] == e.name[0]:
                    return prettyForm(greek[name][1], binding=prettyForm.ATOM)
                else:
                    return prettyForm(greek[name][0], binding=prettyForm.ATOM)

    def _print_Factorial(self, e):
        x = e.args[0]
        if (isinstance(x, Basic.Integer) and x.is_nonnegative) or \
            isinstance(x, Basic.Symbol):
            s = self._print(x)
        else:
            s = "(" + self._print(x) + ")"
        return s + "!"

    def _print_Exp1(self, e):
        if self._use_unicode:
            return prettyForm(u'\u212f', binding=prettyForm.ATOM)

    def _print_Pi(self, e):
        if self._use_unicode:
            return prettyForm(u'\u03c0', binding=prettyForm.ATOM)

    def _print_Infinity(self, e):
        if self._use_unicode:
            return prettyForm(u'\u221e', binding=prettyForm.ATOM)

    def _print_NegativeInfinity(self, e):
        if self._use_unicode:
            return prettyForm(u'-\u221e', binding=prettyForm.ATOM)

    def _print_ImaginaryUnit(self, e):
        if self._use_unicode:
            return prettyForm(u'\u03b9', binding=prettyForm.ATOM)

    def _print_Relational(self, e):
        charmap = {
            '==': ('=', '='),
            '<':  ('<', '<'),
            '<=': ('<=', u'\u2264'),
            '!=': ('!=', u'\u2260')
        }
        if self._use_unicode:
            op = charmap[e.rel_op][1]
        else:
            op = charmap[e.rel_op][0]
        op = prettyForm(' ' + op + ' ')

        l = self._print(e.lhs)
        r = self._print(e.rhs)
        pform = prettyForm(*stringPict.next(l, op))
        pform = prettyForm(*stringPict.next(pform, r))
        return pform

    def _print_ApplyConjugate(self, e):
        pform = self._print(e.args[0])
        return prettyForm(*stringPict.top(pform, '_'*pform.width()))

    def _print_ApplyAbs(self, e):
        pform = self._print(e.args[0])
        pform.baseline = 0
        bars = '|' + ('\n|' * (pform.height()-1))
        pform = prettyForm(*stringPict.next(bars, pform))
        pform = prettyForm(*stringPict.next(pform, bars))
        return pform

    def _print_Derivative(self, deriv):
        syms = list(deriv.symbols)
        syms.reverse()
        x = None
        for sym in syms:
            if x is None:
                x = prettyForm('d' + str(sym))
            else:
                x = prettyForm(*stringPict.next(x, ' d' + str(sym)))

        f = prettyForm(binding=prettyForm.FUNC, *self._print(deriv.expr).parens())

        pform = prettyForm('d')
        if len(syms) > 1:
            pform = pform ** prettyForm(str(len(deriv.symbols)))

        pform = prettyForm(*pform.below(stringPict.LINE, x))
        pform.baseline = pform.baseline + 1
        pform = prettyForm(*stringPict.next(pform, f))
        return pform

    def _print_Integral(self, integral):
        f,x = integral.f, integral.x
        a,b = integral.a, integral.b

        # Add parentheses if a sum and create pretty form for argument
        prettyF = self._print(f)
        if isinstance(f, Basic.Add):
            prettyF = prettyForm(*prettyF.parens())

        arg = prettyForm( *stringPict.next(prettyF, " d" , self._print(x)) )
        arg.baseline = 0

        # Create pretty forms for endpoints, if definite integral
        if a is not None:
            prettyA = self._print(a)
            prettyB = self._print(b)

        # Create bar based on the height of the argument
        bar = '  |   ' + ('\r  |   ' * (arg.height()+1))

        # Construct the pretty form with the integral sign and the argument
        pform = prettyForm(bar)
        pform = prettyForm(*pform.below('/  '))
        pform = prettyForm(*pform.top(' /'))
        pform.baseline = (arg.height() + 3)/2

        if a is not None:
            pform = prettyForm(*stringPict.top(pform, prettyB))
            pform = prettyForm(*stringPict.below(pform, prettyA))
        pform = prettyForm(*stringPict.right(pform, arg))
        return pform

    def _print_ApplyExp(self, e):
        if self._use_unicode:
            base = prettyForm(u'\u212f', binding=prettyForm.ATOM)
        else:
            base = prettyForm('e', binding=prettyForm.ATOM)
        return base ** self._print(e.args[0])

    def _print_Apply(self, e):
        func = e.func
        args = e.args
        n = len(args)

        prettyFunc = self._print(Basic.Symbol(func.name));
        prettyArgs = self._print(args[0])
        for i in xrange(1, n):
            pform = self._print(args[i])
            prettyArgs = prettyForm(*stringPict.next(prettyArgs, ', '))
            prettyArgs = prettyForm(*stringPict.next(prettyArgs, pform))

        pform = prettyForm(*stringPict.next(prettyFunc, '('))
        pform = prettyForm(*stringPict.next(pform, prettyArgs))
        pform = stringPict.next(pform, ')')
        return prettyForm(binding=prettyForm.FUNC, *pform)

    def _print_Add(self, sum):
        pforms = []
        for x in sum:
            # Check for negative "things" so that this information can be enforce upon
            # the pretty form so that it can be made of use (such as in a sum).
            if isinstance(x, Basic.Mul) and x.as_coeff_terms()[0] < 0:
                pform1 = self._print(-x)
                if len(pforms) == 0:
                    if pform1.height() > 1:
                        pform2 = '- '
                    else:
                        pform2 = '-'
                else:
                    pform2 = ' - '
                pform = stringPict.next(pform2, pform1)
                pforms.append(prettyForm(binding=prettyForm.NEG, *pform))
            elif isinstance(x, Basic.Number) and x < 0:
                pform1 = self._print(-x)
                if len(pforms) == 0:
                    if pform1.height() > 1:
                        pform2 = '- '
                    else:
                        pform2 = '-'
                    pform = stringPict.next(pform2, pform1)
                else:
                    pform = stringPict.next(' - ', pform1)
                pforms.append(prettyForm(binding=prettyForm.NEG, *pform))
            else:
                pforms.append(self._print(x))
        return prettyForm.__add__(*pforms)

    def _print_Mul(self, product):
        a = [] # items in the numerator
        b = [] # items that are in the denominator (if any)

        # Gather terms for numerator/denominator
        for item in product:
            if isinstance(item, Basic.Pow) and item.exp == -1:
                b.append(item.base)
            elif isinstance(item, Basic.Rational):
                if item.p != 1:
                    a.append( Basic.Rational(item.p) )
                if item.q != 1:
                    b.append( Basic.Rational(item.q) )
            else:
                a.append(item)

        # Convert to pretty forms. Add parens to Add instances if there
        # is more than one term in the numer/denom
        for i in xrange(0, len(a)):
            if isinstance(a[i], Basic.Add) and len(a) > 1:
                a[i] = prettyForm(*self._print(a[i]).parens())
            else:
                a[i] = self._print(a[i])

        for i in xrange(0, len(b)):
            if isinstance(b[i], Basic.Add) and len(b) > 1:
                b[i] = prettyForm(*self._print(b[i]).parens())
            else:
                b[i] = self._print(b[i])

        # Construct a pretty form
        if len(b) == 0:
            return prettyForm.__mul__(*a)
        else:
            if len(a) == 0:
                a.append( self._print(Basic.One()) )
            return prettyForm.__mul__(*a) / prettyForm.__mul__(*b)

    def _print_Pow(self, power):
        if isinstance(power.exp, Basic.Half):
            # If it's a square root
            bpretty = self._print(power.base)
            bl = int((bpretty.height() / 2.0) + 0.5)

            s2 = stringPict("\\/")
            for x in xrange(1, bpretty.height()):
                s3 = stringPict(" " * (2*x+1) + "/")
                s2 = stringPict(*s2.top(s3))
            s2.baseline = -1

            s = prettyForm("__" + "_" * bpretty.width())
            s = prettyForm(*stringPict.below(s, bpretty))
            s = prettyForm(*stringPict.left(s, s2))
            s.baseline = bl
            return s
        elif power.exp == -1:
            # Things like 1/x
            return prettyForm("1") / self._print(power.base)

        # None of the above special forms, do a standard power
        b,e = power.as_base_exp()
        return self._print(b)**self._print(e)

    def _print_Rational(self, r):
        if r.q == 1:
            return prettyForm(str(r.p), prettyForm.ATOM)
        elif abs(r.p) > 10 and abs(r.q) > 10:
            # If more than one digit in numer and denom, print larger fraction
            if r.is_negative:
                pform = prettyForm(str(-r.p))/prettyForm(str(r.q))
                return prettyForm(binding=prettyForm.NEG, *pform.left('- '))
            else:
                return prettyForm(str(r.p))/prettyForm(str(r.q))

def pretty(expr, use_unicode=False):
    """
    Returns a string containing the prettified form of expr. If use_unicode
    is set to True then certain expressions will use unicode characters,
    such as the greek letter pi for Basic.Pi instances.
    """
    pp = PrettyPrinter(use_unicode)
    return pp.doprint(expr)

def pretty_print(expr, use_unicode=False):
    """
    Prints expr in pretty form.

    pprint is just a shortcut for this function
    """
    print pretty(expr, use_unicode)

pprint = pretty_print
