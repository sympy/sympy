from sympy.core import S, C, Basic
from sympy.printing.printer import Printer
from sympy.printing.str import sstr
from stringpict import prettyForm, stringPict
from sympy import Interval

from pretty_symbology import xstr, hobj, vobj, xobj, xsym, pretty_symbol,\
        pretty_atom, pretty_use_unicode, pretty_try_use_unicode, greek

# rename for usage from outside
pprint_use_unicode = pretty_use_unicode
pprint_try_use_unicode = pretty_try_use_unicode


class PrettyPrinter(Printer):
    """Printer, which converts an expression into 2D ascii-art figure."""

    def __init__(self, profile=None):
        Printer.__init__(self)
        self.emptyPrinter = lambda x : prettyForm(xstr(x))

        self._settings = {
                "full_prec" : "auto",
                "use_unicode" : True,
                "wrap_line" : True,
        }

        if profile is not None:
            self._settings.update(profile)

    def doprint(self, expr):
        return self._print(expr).render(**self._settings)

    # empty op so _print(stringPict) returns the same
    def _print_stringPict(self, e):
        return e

    def _print_basestring(self, e):
        return prettyForm(e)

    def _print_Symbol(self, e):
        symb = pretty_symbol(e.name)
        return prettyForm(symb)

    def _print_Real(self, e):
        # we will use StrPrinter's Real printer, but we need to handle the
        # full_prec ourselves, according to the self._print_level
        full_prec = self._settings["full_prec"]
        if  full_prec == "auto":
            full_prec = self._print_level == 1
        return prettyForm(sstr(e, full_prec=full_prec))

    def _print_Atom(self, e):
        try:
            # print atoms like Exp1 or Pi
            return prettyForm(pretty_atom(e.__class__.__name__))
        except KeyError:
            return self.emptyPrinter(e)

    # Infinity inherits from Rational, so we have to override _print_XXX order
    _print_Infinity         = _print_Atom
    _print_NegativeInfinity = _print_Atom
    _print_EmptySet         = _print_Atom

    def _print_Factorial(self, e):
        x = e.args[0]
        pform = self._print(x)
        # Add parentheses if needed
        if not ((x.is_Integer and x.is_nonnegative) or x.is_Symbol):
            pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.right('!'))
        return pform

    def _print_Relational(self, e):
        op = prettyForm(' ' + xsym(e.rel_op) + ' ')

        l = self._print(e.lhs)
        r = self._print(e.rhs)
        pform = prettyForm(*stringPict.next(l, op, r))
        return pform

    def _print_conjugate(self, e):
        pform = self._print(e.args[0])
        return prettyForm( *pform.above( hobj('_',pform.width())) )

    def _print_abs(self, e):
        pform = self._print(e.args[0])
        pform = prettyForm(*pform.parens('|', '|'))
        return pform

    def _print_floor(self, e):
        if pretty_use_unicode():
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens('lfloor', 'rfloor'))
            return pform
        else:
            return self._print_Function(e)

    def _print_ceiling(self, e):
        if pretty_use_unicode():
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens('lceil', 'rceil'))
            return pform
        else:
            return self._print_Function(e)

    def _print_Derivative(self, deriv):
        # XXX use U('PARTIAL DIFFERENTIAL') here ?
        syms = list(deriv.symbols)
        syms.reverse()
        x = None
        for sym in syms:
            S = self._print(sym)
            dS= prettyForm(*S.left('d'))

            if x is None:
                x = dS
            else:
                x = prettyForm(*x.right(' '))
                x = prettyForm(*x.right(dS))

        f = prettyForm(binding=prettyForm.FUNC, *self._print(deriv.expr).parens())

        pform = prettyForm('d')
        if len(syms) > 1:
            pform = pform ** prettyForm(str(len(deriv.symbols)))

        pform = prettyForm(*pform.below(stringPict.LINE, x))
        pform.baseline = pform.baseline + 1
        pform = prettyForm(*stringPict.next(pform, f))
        return pform

    def _print_PDF(self, pdf):
        lim = self._print(pdf.pdf.args[0])
        lim = prettyForm(*lim.right(', '))
        lim = prettyForm(*lim.right(self._print(pdf.domain[0])))
        lim = prettyForm(*lim.right(', '))
        lim = prettyForm(*lim.right(self._print(pdf.domain[1])))
        lim = prettyForm(*lim.parens())

        f = self._print(pdf.pdf.args[1])
        f = prettyForm(*f.right(', '))
        f = prettyForm(*f.right(lim))
        f = prettyForm(*f.parens())

        pform = prettyForm('PDF')
        pform = prettyForm(*pform.right(f))
        return pform

    def _print_Integral(self, integral):
        f   = integral.function

        # Add parentheses if a sum and create pretty form for argument
        prettyF = self._print(f)
        # XXX generalize parents
        if f.is_Add:
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
                # Create pretty forms for endpoints, if definite integral.
                # Do not print empty endpoints.
                if ab[0] is None:
                    prettyA = prettyForm("")
                else:
                    prettyA = self._print(ab[0])
                if ab[1] is None:
                    prettyB = prettyForm("")
                else:
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


    def _print_Limit(self, l):
        # XXX we do not print dir ...
        e, z, z0, dir = l.args

        E       = self._print(e)
        Lim     = prettyForm('lim')

        LimArg  = self._print(z)
        LimArg  = prettyForm(*LimArg.right('->'))
        LimArg  = prettyForm(*LimArg.right(self._print(z0)))

        Lim     = prettyForm(*Lim.below(LimArg))
        Lim     = prettyForm(*Lim.right(E))


        return Lim

    # Matrix is special:
    #
    # it can exist in SymPy in two forms:
    # - as Matrix
    # - as _MatrixAsBasic
    #
    # see _MatrixAsBasic docstring, and #420
    def _print__MatrixAsBasic(self, e):
        return self._print_Matrix(e.m)

    def _print_Matrix(self, e):
        M = e   # matrix
        S = {}  # i,j -> pretty(M[i,j])
        for i in range(M.rows):
            for j in range(M.cols):
                S[i,j] = self._print(M[i,j])

        # h- and v- spacers
        hsep = 2
        vsep = 1

        # max width for columns
        maxw = [-1] * M.cols

        for j in range(M.cols):
            maxw[j] = max([S[i,j].width()  for i in range(M.rows)])


        # drawing result
        D = None

        for i in range(M.rows):

            D_row = None
            for j in range(M.cols):
                s = S[i,j]

                # reshape s to maxw
                # XXX this should be generalized, and go to stringPict.reshape ?
                assert s.width()  <= maxw[j]

                # hcenter it, +0.5 to the right                        2
                # ( it's better to align formula starts for say 0 and r )
                # XXX this is not good in all cases -- maybe introduce vbaseline?
                wdelta = maxw[j] - s.width()
                wleft  = wdelta // 2
                wright = wdelta - wleft

                s = prettyForm(*s.right(' '*wright))
                s = prettyForm(*s.left (' '*wleft))

                # we don't need vcenter cells -- this is automatically done in
                # a pretty way because when their baselines are taking into
                # account in .right()

                if D_row is None:
                    D_row = s   # first box in a row
                    continue

                D_row = prettyForm(*D_row.right(' '*hsep))  # h-spacer
                D_row = prettyForm(*D_row.right(s))

            if D is None:
                D = D_row       # first row in a picture
                continue

            # v-spacer
            for _ in range(vsep):
                D = prettyForm(*D.below(' '))

            D = prettyForm(*D.below(D_row))

        D = prettyForm(*D.parens('[',']'))
        return D

    def _print_Piecewise(self, pexpr):

        P = {}
        for n, ec in enumerate(pexpr.args):
            P[n,0] = self._print(ec.expr)
            if ec.cond == True:
                P[n,1] = prettyForm('otherwise')
            else:
                P[n,1] = prettyForm(*prettyForm('for ').right(self._print(ec.cond)))
        hsep = 2
        vsep = 1
        len_args = len(pexpr.args)

        # max widths
        maxw = [max([P[i,j].width() for i in xrange(len_args)]) \
                    for j in xrange(2)]

        # FIXME: Refactor this code and matrix into some tabular environment.
        # drawing result
        D = None

        for i in xrange(len_args):
            D_row = None
            for j in xrange(2):
                p = P[i,j]
                assert p.width() <= maxw[j]

                wdelta = maxw[j] - p.width()
                wleft  = wdelta // 2
                wright = wdelta - wleft

                p = prettyForm(*p.right(' '*wright))
                p = prettyForm(*p.left (' '*wleft))

                if D_row is None:
                    D_row = p
                    continue

                D_row = prettyForm(*D_row.right(' '*hsep))  # h-spacer
                D_row = prettyForm(*D_row.right(p))
            if D is None:
                D = D_row       # first row in a picture
                continue

            # v-spacer
            for _ in range(vsep):
                D = prettyForm(*D.below(' '))

            D = prettyForm(*D.below(D_row))

        D = prettyForm(*D.parens('{',''))
        return D

    def _print_exp(self, e):
        base = prettyForm(pretty_atom('Exp1', 'e'))
        return base ** self._print(e.args[0])

    def _print_Function(self, e):
        # XXX works only for applied functions
        func = e.func
        args = e.args
        n = len(args)

        func_name = func.__name__

        prettyFunc = self._print(C.Symbol(func_name));
        prettyArgs = prettyForm(*self._print_seq(args).parens())

        pform = prettyForm(binding=prettyForm.FUNC, *stringPict.next(prettyFunc, prettyArgs))

        # store pform parts so it can be reassembled e.g. when powered
        pform.prettyFunc = prettyFunc
        pform.prettyArgs = prettyArgs

        return pform

    def _print_gamma(self, e):
        if pretty_use_unicode():
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens())
            pform = prettyForm(*pform.left(greek['gamma'][1]))
            return pform
        else:
            return self._print_Function(e)

    def _print_Add(self, sum):
        args = list(sum.args)
        args.sort(Basic._compare_pretty)
        pforms = []
        for x in args:
            # Check for negative "things" so that this information can be enforce upon
            # the pretty form so that it can be made of use (such as in a sum).
            if x.is_Mul and x.as_coeff_terms()[0] < 0:
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
            elif x.is_Number and x < 0:
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
        for item in product.args:
            if item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                b.append(C.Pow(item.base, -item.exp))
            elif item.is_Rational and item is not S.Infinity:
                if item.p != 1:
                    a.append( C.Rational(item.p) )
                if item.q != 1:
                    b.append( C.Rational(item.q) )
            else:
                a.append(item)

        # Convert to pretty forms. Add parens to Add instances if there
        # is more than one term in the numer/denom
        for i in xrange(0, len(a)):
            if a[i].is_Add and len(a) > 1:
                a[i] = prettyForm(*self._print(a[i]).parens())
            else:
                a[i] = self._print(a[i])

        for i in xrange(0, len(b)):
            if b[i].is_Add and len(b) > 1:
                b[i] = prettyForm(*self._print(b[i]).parens())
            else:
                b[i] = self._print(b[i])

        # Construct a pretty form
        if len(b) == 0:
            return prettyForm.__mul__(*a)
        else:
            if len(a) == 0:
                a.append( self._print(S.One) )
            return prettyForm.__mul__(*a) / prettyForm.__mul__(*b)

    def _print_Pow(self, power):
        # square roots, other roots or n-th roots
        #test for fraction 1/n or power x**-1
        if (isinstance(power.exp, C.Rational) and power.exp.p==1) or \
           (   isinstance(power.exp, C.Pow) and
               isinstance(power.exp.args[0], C.Symbol) and
               power.exp.args[1]==S.NegativeOne):
            bpretty = self._print(power.base)

            #construct root sign, start with the \/ shape
            _zZ= xobj('/',1)
            rootsign = xobj('\\',1)+_zZ
            #make exponent number to put above it
            if isinstance(power.exp, C.Rational):
                exp = str(power.exp.q)
                if exp=='2': exp = ''
            else: exp = str(power.exp.args[0])
            exp = exp.ljust(2)
            if len(exp)>2: rootsign = ' '*(len(exp)-2)+rootsign
            #stack the exponent
            rootsign = stringPict(exp+'\n'+rootsign)
            rootsign.baseline = 0
            #diagonal: length is one less than height of base
            linelength = bpretty.height()-1
            diagonal = stringPict('\n'.join(
                ' '*(linelength-i-1)+_zZ+' '*i
                for i in range(linelength)
                ))
            #put baseline just below lowest line: next to exp
            diagonal.baseline = linelength-1
            #make the root symbol
            rootsign = prettyForm(*rootsign.right(diagonal))
            #set the baseline to match contents to fix the height
            #but if the height of bpretty is one, the rootsign must be one higher
            rootsign.baseline = max(1, bpretty.baseline)
            #build result
            s = prettyForm(hobj('_', 2+ bpretty.width()))
            s = prettyForm(*bpretty.above(s))
            s = prettyForm(*s.left(rootsign))
            return s
        elif power.exp.is_Rational and power.exp.is_negative:
            # Things like 1/x
            return prettyForm("1") / self._print(C.Pow(power.base, -power.exp))

        # None of the above special forms, do a standard power
        b,e = power.as_base_exp()
        return self._print(b)**self._print(e)

    def _print_Rational(self, r):
        if r.q == 1:
            if r.is_negative:
                return prettyForm(str(r.p),binding=prettyForm.NEG)
            else:
                return prettyForm(str(r.p))
        elif abs(r.p) >= 10 and abs(r.q) >= 10:
            # If more than one digit in numer and denom, print larger fraction
            if r.is_negative:
                pform = prettyForm(str(-r.p))/prettyForm(str(r.q))
                return prettyForm(binding=prettyForm.NEG, *pform.left('- '))
            else:
                return prettyForm(str(r.p))/prettyForm(str(r.q))
        else:
            return self.emptyPrinter(r)

    def _print_Interval(self, i):
        if i.start == i.end:
            return self._print_seq(i.args[:1], '{', '}')

        else:
            if i.left_open:
                left = '('
            else:
                left = '['

            if i.right_open:
                right = ')'
            else:
                right = ']'

            return self._print_seq(i.args[:2], left, right)

    def _print_Union(self, u):
        other_sets, singletons = [], []
        for set in u.args:
            if isinstance(set, Interval) and set.measure == 0:
                singletons.append(set.start)
            else:
                other_sets.append(set)

        union_delimiter = ' %s ' % pretty_atom('Union')

        S2 = self._print_seq(other_sets, None, None, union_delimiter)

        if len(singletons) > 0:
            S1 = self._print_seq(singletons, '{', '}')

            S = prettyForm(*stringPict.next(S1, union_delimiter))
            S = prettyForm(*stringPict.next(S, S2))
        else:
            S = S2

        return S

    def _print_seq(self, seq, left=None, right=None, delimiter=', '):
        S = None

        for item in seq:
            pform = self._print(item)

            if S is None:
                # first element
                S = pform
            else:
                S = prettyForm(*stringPict.next(S, delimiter))
                S = prettyForm(*stringPict.next(S, pform))

        if S is None:
            S = stringPict('')

        S = prettyForm(*S.parens(left, right, ifascii_nougly=True))
        return S

    def _print_list(self, l):
        return self._print_seq(l, '[', ']')

    def _print_tuple(self, t):
        if len(t) == 1:
            ptuple = prettyForm(*stringPict.next(self._print(t[0]), ','))
            return prettyForm(*ptuple.parens('(', ')', ifascii_nougly=True))
        else:
            return self._print_seq(t, '(', ')')

    def _print_dict(self, d):
        items = []

        keys = d.keys()
        keys.sort( Basic.compare_pretty )

        for k in keys:
            K = self._print(k)
            V = self._print(d[k])
            S = prettyForm(*stringPict.next(K, ': ', V))

            items.append(S)

        return self._print_seq(items, '{', '}')

    def __print_set(self, s):
        items = list(s)
        items.sort( Basic.compare_pretty )

        S = self._print_seq(items, '(', ')')
        S = prettyForm(*stringPict.next(type(s).__name__, S))
        return S

    _print_set       = __print_set
    _print_frozenset = __print_set


def pretty(expr, profile=None, **kargs):
    """
    Returns a string containing the prettified form of expr.

    Arguments
    ---------
    expr: the expression to print
    wrap_line: line wrapping enabled/disabled, should be a boolean value (default to True)
    use_unicode: use unicode characters, such as the greek letter pi instead of
        the string pi. Values should be boolean or None
    full_prec: use full precission. Default to "auto"
    """
    if profile is not None:
        profile.update(kargs)
    else:
        profile = kargs
    uflag = pretty_use_unicode(kargs.get("use_unicode", None))
    try:
        pp = PrettyPrinter(profile)
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
