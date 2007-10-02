from sympy.core import Basic
from printer import Printer
from stringpict import *

from pretty_symbology import xstr, hobj, vobj, xobj, xrel, pretty_symbol, pretty_atom, pretty_use_unicode 

pprint_use_unicode = pretty_use_unicode


def prettyAtom(s):
    return prettyForm(s, binding=prettyForm.ATOM)

class PrettyPrinter(Printer):
    """
    A class that prints a prettified expression, one that is not limited
    to one dimension like casting the expression to a string would return.
    """
    def __init__(self, use_unicode=None):
        Printer.__init__(self)
        self.emptyPrinter = lambda x : prettyAtom(xstr(x))

    def doprint(self, expr):
        Printer.doprint.__doc__
        return self._print(expr).terminal_string()

    def _print_Symbol(self, e):
        symb = pretty_symbol(e.name)
        return prettyAtom(symb)

    def _print_Atom(self, e):
        try:
            # print atoms like Exp1 or Pi
            return prettyAtom(pretty_atom(e.__class__.__name__))
        except KeyError:
            pass

    # Infinity inherits from Rational, so we have to override _print_XXX order
    _print_Infinity         = _print_Atom
    _print_NegativeInfinity = _print_Atom


    def _print_Factorial(self, e):
        x = e[0]
        if (isinstance(x, Basic.Integer) and x.is_nonnegative) or \
            isinstance(x, Basic.Symbol):
            s = self._print(x)
        else:
            # XXX parens
            s = "(" + self._print(x) + ")"
        return s + "!"


    def _print_Relational(self, e):
        op = prettyForm(' ' + xrel(e.rel_op) + ' ')

        l = self._print(e.lhs)
        r = self._print(e.rhs)
        pform = prettyForm(*stringPict.next(l, op, r))
        return pform

    def _print_conjugate(self, e):
        pform = self._print(e[0])
        return prettyForm( *pform.above( hobj('_',pform.width())) )

    def _print_abs(self, e):
        pform = self._print(e[0])

        vbar = vobj('|', pform.height())
        vbar = stringPict(vbar, baseline=pform.baseline)

        pform  = prettyForm(*pform.left (vbar))
        pform  = prettyForm(*pform.right(vbar))
        return pform


    def _print_Derivative(self, deriv):
        # XXX use U('PARTIAL DIFFERENTIAL') here ?
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
        f   = integral.function

        # Add parentheses if a sum and create pretty form for argument
        prettyF = self._print(f)
        # XXX generalize parents
        if isinstance(f, Basic.Add):
            prettyF = prettyForm(*prettyF.parens())

        # dx dy dz ...
        arg = prettyF
        for x,ab in integral.limits:
            prettyArg = self._print(x)
            # XXX qparens   (parens if needs-parens)
            if prettyArg.width() > 1:
                prettyArg = prettyForm(*prettyArg.parens())

            arg = prettyForm(*arg.right(' d', prettyArg))


        # \int \int \int ...
        firstterm = True
        S = None
        for x,ab in integral.limits:
            # Create bar based on the height of the argument
            h = arg.height()
            H = h+2

            # XXX hack!
            ascii_mode = not pretty_use_unicode()
            if ascii_mode:
                H += 2

            vint= vobj('int', H)

            # Construct the pretty form with the integral sign and the argument
            pform = prettyForm(vint)
            #pform.baseline = pform.height()//2  # vcenter
            pform.baseline = arg.baseline + (H-h)//2    # covering the whole argument


            if ab is not None:
                # Create pretty forms for endpoints, if definite integral
                prettyA = self._print(ab[0])
                prettyB = self._print(ab[1])

                if ascii_mode:  # XXX hack
                    # Add spacing so that endpoint can more easily be
                    # identified with the correct integral sign
                    spc = max(1, 3 - prettyB.width())
                    prettyB = prettyForm(*prettyB.left(' ' * spc))

                    spc = max(1, 4 - prettyA.width())
                    prettyA = prettyForm(*prettyA.right(' ' * spc))

                pform = prettyForm(*pform.above(prettyB))
                pform = prettyForm(*pform.below(prettyA))

                #if ascii_mode:  # XXX hack
                #    # too much vspace beetween \int and argument
                #    # but I left it as is
                #    pform = prettyForm(*pform.right(' '))

            if not ascii_mode:  # XXX hack
                pform = prettyForm(*pform.right(' '))

            if firstterm:
                S = pform   # first term
                firstterm = False
            else:
                S = prettyForm(*S.left(pform))

        pform = prettyForm(*arg.left(S))
        return pform

    def _print_exp(self, e):
        base = prettyAtom(pretty_atom('Exp1', 'e'))
        return base ** self._print(e[0])

    def _print_Function(self, e):
        # XXX works only for applied functions
        func = e.func
        args = e[:]
        n = len(args)

        func_name = func.__name__

        prettyFunc = self._print(Basic.Symbol(func_name));
        prettyArgs = self._print(args[0])
        for i in xrange(1, n):
            pform = self._print(args[i])
            prettyArgs = prettyForm(*stringPict.next(prettyArgs, ', '))
            prettyArgs = prettyForm(*stringPict.next(prettyArgs, pform))

        prettyArgs = prettyForm(*prettyArgs.parens(ifascii_nougly=True))

        pform = prettyForm(binding=prettyForm.FUNC, *stringPict.next(prettyFunc, prettyArgs))

        # store pform parts so it can be reassembled e.g. when powered
        pform.prettyFunc = prettyFunc
        pform.prettyArgs = prettyArgs

        return pform

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
            H = bpretty.height()

            _zZ= xobj('/',1)
            s2 = stringPict(xobj('\\',1)+_zZ+' '*(H-1))
            for x in xrange(1, H):
                s3 = stringPict(' '*(x+1) + _zZ + ' '*(H-(x+1)))
                s2 = stringPict(*s2.above(s3))

            s2.baseline = bpretty.baseline  # vertical: each-to-each

            s = prettyForm(hobj('_', 2+ bpretty.width()))
            s = prettyForm(*bpretty.above(s))
            s = prettyForm(*s.left(s2))
            return s
        elif power.exp == -1:
            # Things like 1/x
            return prettyForm("1") / self._print(power.base)

        # None of the above special forms, do a standard power
        b,e = power.as_base_exp()
        return self._print(b)**self._print(e)

    def _print_Rational(self, r):
        if r.q == 1:
            return prettyAtom(str(r.p))
        elif abs(r.p) >= 10 and abs(r.q) >= 10:
            # If more than one digit in numer and denom, print larger fraction
            if r.is_negative:
                pform = prettyForm(str(-r.p))/prettyForm(str(r.q))
                return prettyForm(binding=prettyForm.NEG, *pform.left('- '))
            else:
                return prettyForm(str(r.p))/prettyForm(str(r.q))

def pretty(expr, use_unicode=None):
    """
    Returns a string containing the prettified form of expr. If use_unicode
    is set to True then certain expressions will use unicode characters,
    such as the greek letter pi for Basic.Pi instances.
    """
    uflag = pretty_use_unicode(use_unicode)
    try:
        pp = PrettyPrinter()
        return pp.doprint(expr)
    finally:
        pretty_use_unicode(uflag)

def pretty_print(expr, use_unicode=None):
    """
    Prints expr in pretty form.

    pprint is just a shortcut for this function
    """
    print pretty(expr, use_unicode)

pprint = pretty_print
