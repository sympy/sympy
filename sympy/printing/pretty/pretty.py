from sympy.core import S, C, Basic, Interval
from sympy.utilities import group

from sympy.printing.printer import Printer
from sympy.printing.str import sstr

from stringpict import prettyForm, stringPict
from pretty_symbology import xstr, hobj, vobj, xobj, xsym, pretty_symbol,\
        pretty_atom, pretty_use_unicode, pretty_try_use_unicode, greek, U

from sympy.core.compatibility import cmp_to_key

# rename for usage from outside
pprint_use_unicode = pretty_use_unicode
pprint_try_use_unicode = pretty_try_use_unicode

class PrettyPrinter(Printer):
    """Printer, which converts an expression into 2D ASCII-art figure."""
    printmethod = "_pretty"

    _default_settings = {
        "order": None,
        "full_prec": "auto",
        "use_unicode": None,
        "wrap_line": True,
    }

    def __init__(self, settings=None):
        Printer.__init__(self, settings)
        self.emptyPrinter = lambda x: prettyForm(xstr(x))

    @property
    def _use_unicode(self):
        if self._settings['use_unicode']:
            return True
        else:
            return pretty_use_unicode()

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

    def _print_Float(self, e):
        # we will use StrPrinter's Float printer, but we need to handle the
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

    def _print_factorial(self, e):
        x = e.args[0]
        pform = self._print(x)
        # Add parentheses if needed
        if not ((x.is_Integer and x.is_nonnegative) or x.is_Symbol):
            pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.right('!'))
        return pform

    def _print_binomial(self, e):
        n, k = e.args

        n_pform = self._print(n)
        k_pform = self._print(k)

        bar = ' '*max(n_pform.width(), k_pform.width())

        pform = prettyForm(*k_pform.above(bar))
        pform = prettyForm(*pform.above(n_pform))
        pform = prettyForm(*pform.parens('(', ')'))

        pform.baseline = (pform.baseline + 1)//2

        return pform

    def _print_Relational(self, e):
        op = prettyForm(' ' + xsym(e.rel_op) + ' ')

        l = self._print(e.lhs)
        r = self._print(e.rhs)
        pform = prettyForm(*stringPict.next(l, op, r))
        return pform

    def _print_Not(self, e):
        if self._use_unicode:
            arg = e.args[0]
            pform = self._print(arg)

            if arg.is_Boolean and not arg.is_Not:
                pform = prettyForm(*pform.parens())

            return prettyForm(*pform.left(u"\u00ac "))
        else:
            return self._print_Function(e)

    def __print_Boolean(self, e, char):
        arg = e.args[0]
        pform = self._print(arg)

        if arg.is_Boolean and not arg.is_Not:
            pform = prettyForm(*pform.parens())

        for arg in e.args[1:]:
            pform_arg = self._print(arg)

            if arg.is_Boolean and not arg.is_Not:
                pform_arg = prettyForm(*pform_arg.parens())

            pform = prettyForm(*pform.right(u' %s ' % char))
            pform = prettyForm(*pform.right(pform_arg))

        return pform

    def _print_And(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, u"\u2227")
        else:
            return self._print_Function(e)

    def _print_Or(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, u"\u2228")
        else:
            return self._print_Function(e)

    def _print_Xor(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, u"\u22bb")
        else:
            return self._print_Function(e)

    def _print_Nand(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, u"\u22bc")
        else:
            return self._print_Function(e)

    def _print_Nor(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, u"\u22bd")
        else:
            return self._print_Function(e)

    def _print_Implies(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, u"\u2192")
        else:
            return self._print_Function(e)

    def _print_Equivalent(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, u"\u2261")
        else:
            return self._print_Function(e)

    def _print_conjugate(self, e):
        pform = self._print(e.args[0])
        return prettyForm( *pform.above( hobj('_',pform.width())) )

    def _print_Abs(self, e):
        pform = self._print(e.args[0])
        pform = prettyForm(*pform.parens('|', '|'))
        return pform

    def _print_floor(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens('lfloor', 'rfloor'))
            return pform
        else:
            return self._print_Function(e)

    def _print_ceiling(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens('lceil', 'rceil'))
            return pform
        else:
            return self._print_Function(e)

    def _print_Derivative(self, deriv):
        # XXX use U('PARTIAL DIFFERENTIAL') here ?
        syms = list(reversed(deriv.variables))
        x = None

        for sym, num in group(syms, multiple=False):
            s = self._print(sym)
            ds = prettyForm(*s.left('d'))

            if num > 1:
                ds = ds**prettyForm(str(num))

            if x is None:
                x = ds
            else:
                x = prettyForm(*x.right(' '))
                x = prettyForm(*x.right(ds))

        f = prettyForm(binding=prettyForm.FUNC, *self._print(deriv.expr).parens())

        pform = prettyForm('d')

        if len(syms) > 1:
            pform = pform**prettyForm(str(len(syms)))

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

        # Add parentheses if arg involves addition of terms and
        # create a pretty form for the argument
        prettyF = self._print(f)
        # XXX generalize parens
        if f.is_Add:
            prettyF = prettyForm(*prettyF.parens())

        # dx dy dz ...
        arg = prettyF
        for x in integral.limits:
            prettyArg = self._print(x[0])
            # XXX qparens   (parens if needs-parens)
            if prettyArg.width() > 1:
                prettyArg = prettyForm(*prettyArg.parens())

            arg = prettyForm(*arg.right(' d', prettyArg))


        # \int \int \int ...
        firstterm = True
        s = None
        for lim in integral.limits:
            x = lim[0]
            # Create bar based on the height of the argument
            h = arg.height()
            H = h+2

            # XXX hack!
            ascii_mode = not self._use_unicode
            if ascii_mode:
                H += 2

            vint= vobj('int', H)

            # Construct the pretty form with the integral sign and the argument
            pform = prettyForm(vint)
            #pform.baseline = pform.height()//2  # vcenter
            pform.baseline = arg.baseline + (H-h)//2    # covering the whole argument


            if len(lim) > 1:
                # Create pretty forms for endpoints, if definite integral.
                # Do not print empty endpoints.
                if len(lim) == 2:
                    prettyA = prettyForm("")
                    prettyB = self._print(lim[1])
                if len(lim) == 3:
                    prettyA = self._print(lim[1])
                    prettyB = self._print(lim[2])

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
                s = pform   # first term
                firstterm = False
            else:
                s = prettyForm(*s.left(pform))

        pform = prettyForm(*arg.left(s))
        return pform

    def _print_Sum(self, expr):
        def asum(hrequired, lower, upper):
            def adjust(s, wid=None, how='<^>'):
                if not wid or len(s)>wid:
                    return s
                need = wid - len(s)
                if how == '<^>' or how == "<" or how not in list('<^>'):
                    return s + ' '*need
                half = need//2
                lead = ' '*half
                if how == ">":
                    return " "*need + s
                return lead + s + ' '*(need - len(lead))

            h = max(hrequired, 2)
            d = h//2
            wrequired = max(lower, upper)
            w = d + 1
            more = hrequired % 2

            lines = []
            lines.append("_"*(w) + ' ')
            lines.append("\%s`" % (' '*(w - 1)))
            for i in range(1, d):
              lines.append('%s\\%s' % (' '*i, ' '*(w - i)))
            if more:
                lines.append('%s)%s' % (' '*(d), ' '*(w - d)))
            for i in reversed(range(1, d)):
              lines.append('%s/%s' % (' '*i, ' '*(w - i)))
            lines.append("/" + "_"*(w - 1) + ',')
            return d, h + more, lines

        f = expr.function

        prettyF = self._print(f)

        if f.is_Add: # add parens
            prettyF = prettyForm(*prettyF.parens())

        H = prettyF.height() + 2

        # \sum \sum \sum ...
        first = True
        max_upper = 0
        sign_height = 0

        for lim in expr.limits:
            if len(lim) == 3:
                prettyUpper = self._print(lim[2])
                prettyLower = self._print(C.Equality(lim[0], lim[1]))
            elif len(lim) == 2:
                prettyUpper = self._print("")
                prettyLower = self._print(C.Equality(lim[0], lim[1]))
            elif len(lim) == 1:
                prettyUpper = self._print("")
                prettyLower = self._print(lim[0])

            max_upper = max(max_upper, prettyUpper.height())

            # Create sum sign based on the height of the argument
            d, h, slines = asum(H, prettyLower.width(), prettyUpper.width())
            prettySign = stringPict('')
            prettySign = prettyForm(*prettySign.stack(*slines))

            if first:
                sign_height = prettySign.height()

            prettySign = prettyForm(*prettySign.above(prettyUpper))
            prettySign = prettyForm(*prettySign.below(prettyLower))

            if first:
                # change F baseline so it centers on the sign
                prettyF.baseline -= d - (prettyF.height()//2 -
                                         prettyF.baseline)
                first = False

            # put padding to the right
            pad = stringPict('')
            pad = prettyForm(*pad.stack(*[' ']*h))
            prettySign = prettyForm(*prettySign.right(pad))
            # put the present prettyF to the right
            prettyF = prettyForm(*prettySign.right(prettyF))

        prettyF.baseline = max_upper + sign_height//2
        return prettyF

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
        Ms = {}  # i,j -> pretty(M[i,j])
        for i in range(M.rows):
            for j in range(M.cols):
                Ms[i,j] = self._print(M[i,j])

        # h- and v- spacers
        hsep = 2
        vsep = 1

        # max width for columns
        maxw = [-1] * M.cols

        for j in range(M.cols):
            maxw[j] = max([Ms[i,j].width()  for i in range(M.rows)])


        # drawing result
        D = None

        for i in range(M.rows):

            D_row = None
            for j in range(M.cols):
                s = Ms[i,j]

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

        if D is None:
            D = prettyForm('') # Empty Matrix

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

    def _hprint_vec(self, v):
        D = None

        for a in v:
            p = a
            if D is None:
                D = p
            else:
                D = prettyForm(*D.right(', '))
                D = prettyForm(*D.right(p))
        if D is None:
            D = stringPict(' ')

        return D

    def _hprint_vseparator(self, p1, p2):
        tmp = prettyForm(*p1.right(p2))
        sep = stringPict(vobj('|', tmp.height()), baseline=tmp.baseline)
        return prettyForm(*p1.right(sep, p2))

    def _print_hyper(self, e):
        # FIXME refactor Matrix, Piecewise, and this into a tabular environment
        ap = [self._print(a) for a in e.ap]
        bq = [self._print(b) for b in e.bq]

        P = self._print(e.argument)

        # Drawing result - first create the ap, bq vectors
        D = None
        for v in [ap, bq]:
            D_row = self._hprint_vec(v)
            if D is None:
                D = D_row       # first row in a picture
            else:
                D = prettyForm(*D.below(' '))
                D = prettyForm(*D.below(D_row))

        # make sure that the argument `z' is centred vertically
        D.baseline = D.height()//2

        # insert horizontal separator
        P = prettyForm(*P.left(' '))
        D = prettyForm(*D.right(' '))

        # insert separating `|`
        D = self._hprint_vseparator(D, P)

        # add parens
        D = prettyForm(*D.parens('(', ')'))

        # create the F symbol
        above = D.height()//2 - 1
        below = D.height() - above - 1

        if self._use_unicode:
            pic = (2, 0, 2, u'\u250c\u2500\n\u251c\u2500\n\u2575')
        else:
            pic = ((3, 0, 3, ' _\n|_\n|\n'))

        add = 0
        sz, t, b, img = pic
        F = prettyForm('\n' * (above - t) + img + '\n' * (below - b),
                       baseline = above + sz)
        add = (sz+1)//2

        F = prettyForm(*F.left(self._print(len(e.ap))))
        F = prettyForm(*F.right(self._print(len(e.bq))))
        F.baseline = above + add

        D = prettyForm(*F.right(' ', D))

        return D

    def _print_meijerg(self, e):
        # FIXME refactor Matrix, Piecewise, and this into a tabular environment

        v = {}
        v[(0, 0)] = [self._print(a) for a in e.an]
        v[(0, 1)] = [self._print(a) for a in e.aother]
        v[(1, 0)] = [self._print(b) for b in e.bm]
        v[(1, 1)] = [self._print(b) for b in e.bother]

        P = self._print(e.argument)

        vp = {}
        for idx in v:
            vp[idx] = self._hprint_vec(v[idx])

        for i in range(2):
            maxw = max(vp[(0, i)].width(), vp[(1, i)].width())
            for j in range(2):
                s = vp[(j, i)]
                left = (maxw - s.width()) // 2
                right = maxw - left - s.width()
                s = prettyForm(*s.left(' ' * left))
                s = prettyForm(*s.right(' ' * right))
                vp[(j, i)] = s

        D1 = prettyForm(*vp[(0, 0)].right('  ', vp[(0, 1)]))
        D1 = prettyForm(*D1.below(' '))
        D2 = prettyForm(*vp[(1, 0)].right('  ', vp[(1, 1)]))
        D  = prettyForm(*D1.below(D2))

        # make sure that the argument `z' is centred vertically
        D.baseline = D.height()//2

        # insert horizontal separator
        P = prettyForm(*P.left(' '))
        D = prettyForm(*D.right(' '))

        # insert separating `|`
        D = self._hprint_vseparator(D, P)

        # add parens
        D = prettyForm(*D.parens('(', ')'))

        # create the G symbol
        above = D.height()//2 - 1
        below = D.height() - above - 1

        if self._use_unicode:
            pic = (3, 0, 3, 1,
                   u'\u256d\u2500\u256e\n\u2502\u2576\u2510\n\u2570\u2500\u256f')
        else:
            pic = (3, 0, 3, 1, ' __\n/__\n\_|')

        add = 0
        sz, t, b, add, img = pic
        F = prettyForm('\n' * (above - t) + img + '\n' * (below - b),
                       baseline = above + sz)

        pp = self._print(len(e.ap))
        pq = self._print(len(e.bq))
        pm = self._print(len(e.bm))
        pn = self._print(len(e.an))
        def adjust(p1, p2):
            diff = p1.width() - p2.width()
            if diff == 0:
                return p1, p2
            elif diff > 0:
                return p1, prettyForm(*p2.left(' '*diff))
            else:
                return prettyForm(*p1.left(' '*-diff)), p2
        pp, pm = adjust(pp, pm)
        pq, pn = adjust(pq, pn)
        pu = prettyForm(*pm.right(', ', pn))
        pl = prettyForm(*pp.right(', ', pq))

        ht = F.baseline - above - 2
        if ht > 0:
            pu = prettyForm(*pu.below('\n'*ht))
        p = prettyForm(*pu.below(pl))

        F.baseline = above
        F = prettyForm(*F.right(p))

        F.baseline = above + add

        D = prettyForm(*F.right(' ', D))

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

        prettyFunc = self._print(C.Symbol(func_name))
        prettyArgs = prettyForm(*self._print_seq(args).parens())

        pform = prettyForm(binding=prettyForm.FUNC, *stringPict.next(prettyFunc, prettyArgs))

        # store pform parts so it can be reassembled e.g. when powered
        pform.prettyFunc = prettyFunc
        pform.prettyArgs = prettyArgs

        return pform

    def _print_Lambda(self, e):
        symbols, expr = e.args

        if len(symbols) == 1:
            symbols = self._print(symbols[0])
        else:
            symbols = self._print(tuple(symbols))

        args = (symbols, self._print(expr))

        prettyFunc = self._print(C.Symbol("Lambda"))
        prettyArgs = prettyForm(*self._print_seq(args).parens())

        pform = prettyForm(binding=prettyForm.FUNC, *stringPict.next(prettyFunc, prettyArgs))

        # store pform parts so it can be reassembled e.g. when powered
        pform.prettyFunc = prettyFunc
        pform.prettyArgs = prettyArgs

        return pform

    def _print_Order(self, e):
        pform = self._print(e.expr)
        pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.left('O'))
        return pform

    def _print_gamma(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens())
            pform = prettyForm(*pform.left(greek['gamma'][1]))
            return pform
        else:
            return self._print_Function(e)

    def _print_uppergamma(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.right(', ', self._print(e.args[1])))
            pform = prettyForm(*pform.parens())
            pform = prettyForm(*pform.left(greek['gamma'][1]))
            return pform
        else:
            return self._print_Function(e)

    def _print_lowergamma(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.right(', ', self._print(e.args[1])))
            pform = prettyForm(*pform.parens())
            pform = prettyForm(*pform.left(greek['gamma'][0]))
            return pform
        else:
            return self._print_Function(e)

    def _print_Add(self, expr, order=None):
        terms = self._as_ordered_terms(expr, order=order)
        pforms, indices = [], []

        def pretty_negative(pform, index):
            """Prepend a minus sign to a pretty form. """
            if index == 0:
                if pform.height() > 1:
                    pform_neg = '- '
                else:
                    pform_neg = '-'
            else:
                pform_neg = ' - '

            pform = stringPict.next(pform_neg, pform)
            return prettyForm(binding=prettyForm.NEG, *pform)

        for i, term in enumerate(terms):
            if term.is_Mul and term.as_coeff_mul()[0] < 0:
                pform = self._print(-term)
                pforms.append(pretty_negative(pform, i))
            elif term.is_Rational and term.q > 1:
                pforms.append(None)
                indices.append(i)
            elif term.is_Number and term < 0:
                pform = self._print(-term)
                pforms.append(pretty_negative(pform, i))
            else:
                pforms.append(self._print(term))

        if indices:
            large = True

            for pform in pforms:
                if pform is not None and pform.height() > 1:
                    break
            else:
                large = False

            for i in indices:
                term, negative = terms[i], False

                if term < 0:
                    term, negative = -term, True

                if large:
                    pform = prettyForm(str(term.p))/prettyForm(str(term.q))
                else:
                    pform = self._print(term)

                if negative:
                    pform = pretty_negative(pform, i)

                pforms[i] = pform

        return prettyForm.__add__(*pforms)

    def _print_Mul(self, product):
        a = [] # items in the numerator
        b = [] # items that are in the denominator (if any)

        if self.order != 'old':
            args = product.as_ordered_factors()
        else:
            args = product.args

        # Gather terms for numerator/denominator
        for item in args:
            if item.is_commutative and item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
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
            return prettyForm.__mul__(*a)/prettyForm.__mul__(*b)

    def _print_Pow(self, power):
        # square roots, other roots or n-th roots
        #test for fraction 1/n or power x**-1
        if power.is_commutative:
            if (isinstance(power.exp, C.Rational) and power.exp.p==1 and power.exp.q !=1) or \
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

    def __print_numer_denom(self, p, q):
        if q == 1:
            if p < 0:
                return prettyForm(str(p),binding=prettyForm.NEG)
            else:
                return prettyForm(str(p))
        elif abs(p) >= 10 and abs(q) >= 10:
            # If more than one digit in numer and denom, print larger fraction
            if p < 0:
                pform = prettyForm(str(-p))/prettyForm(str(q))
                return prettyForm(binding=prettyForm.NEG, *pform.left('- '))
            else:
                return prettyForm(str(p))/prettyForm(str(q))
        else:
            return None

    def _print_Rational(self, expr):
        result = self.__print_numer_denom(expr.p, expr.q)

        if result is not None:
            return result
        else:
            return self.emptyPrinter(expr)

    def _print_Fraction(self, expr):
        result = self.__print_numer_denom(expr.numerator, expr.denominator)

        if result is not None:
            return result
        else:
            return self.emptyPrinter(expr)

    def _print_ProductSet(self, p):
        prod_char = u'\xd7'
        return self._print_seq(p.sets, None, None, ' %s '%prod_char,
                parenthesize = lambda set:set.is_Union )

    def _print_FiniteSet(self, s):
        if len(s) > 10:
            # Take ten elements from the set at random
            q = iter(s)
            printset = [q.next() for i in xrange(10)]
            printset.append('...')
        else:
            printset = s
        try:
            printset = sorted(printset)
        except:  pass

        return self._print_seq(printset, '{', '}', ', ' )


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

        union_delimiter = ' %s ' % pretty_atom('Union')

        return self._print_seq(u.args, None, None, union_delimiter,
                parenthesize = lambda set:set.is_ProductSet)

    def _print_seq(self, seq, left=None, right=None, delimiter=', ',
            parenthesize = lambda x:False):
        s = None

        for item in seq:
            pform = self._print(item)
            if parenthesize(item):
                pform = prettyForm(*pform.parens())
            if s is None:
                # first element
                s = pform
            else:
                s = prettyForm(*stringPict.next(s, delimiter))
                s = prettyForm(*stringPict.next(s, pform))

        if s is None:
            s = stringPict('')

        s = prettyForm(*s.parens(left, right, ifascii_nougly=True))
        return s

    def _print_list(self, l):
        return self._print_seq(l, '[', ']')

    def _print_tuple(self, t):
        if len(t) == 1:
            ptuple = prettyForm(*stringPict.next(self._print(t[0]), ','))
            return prettyForm(*ptuple.parens('(', ')', ifascii_nougly=True))
        else:
            return self._print_seq(t, '(', ')')

    def _print_Tuple(self, expr):
        return self._print_tuple(expr)

    def _print_dict(self, d):
        items = []

        keys = d.keys()
        keys.sort( key=cmp_to_key(Basic.compare_pretty) )

        for k in keys:
            K = self._print(k)
            V = self._print(d[k])
            s = prettyForm(*stringPict.next(K, ': ', V))

            items.append(s)

        return self._print_seq(items, '{', '}')

    def __print_set(self, set_):
        items = list(set_)
        items.sort( key=cmp_to_key(Basic.compare_pretty) )

        s = self._print_seq(items, '(', ')')
        s = prettyForm(*stringPict.next(type(set_).__name__, s))
        return s

    _print_set       = __print_set
    _print_frozenset = __print_set

    def _print_AlgebraicNumber(self, expr):
        if expr.is_aliased:
            return self._print(expr.as_poly().as_expr())
        else:
            return self._print(expr.as_expr())

    def _print_RootOf(self, expr):
        args = [self._print_Add(expr.expr, order='lex'), expr.index]
        pform = prettyForm(*self._print_seq(args).parens())
        pform = prettyForm(*pform.left('RootOf'))
        return pform

    def _print_RootSum(self, expr):
        args = [self._print_Add(expr.expr, order='lex')]

        if expr.fun is not S.IdentityFunction:
            args.append(self._print(expr.fun))

        pform = prettyForm(*self._print_seq(args).parens())
        pform = prettyForm(*pform.left('RootSum'))

        return pform

    def _print_FiniteField(self, expr):
        if self._use_unicode:
            form = u'\u2124_%d'
        else:
            form = 'GF(%d)'

        return prettyForm(pretty_symbol(form % expr.mod))

    def _print_IntegerRing(self, expr):
        if self._use_unicode:
            return prettyForm(u'\u2124')
        else:
            return prettyForm('ZZ')

    def _print_RationalField(self, expr):
        if self._use_unicode:
            return prettyForm(u'\u211A')
        else:
            return prettyForm('QQ')

    def _print_RealDomain(self, expr):
        if self._use_unicode:
            return prettyForm(u'\u211D')
        else:
            return prettyForm('RR')

    def _print_ComplexDomain(self, expr):
        if self._use_unicode:
            return prettyForm(u'\u2102')
        else:
            return prettyForm('CC')

    def _print_PolynomialRing(self, expr):
        pform = self._print_seq(expr.gens, '[', ']')
        pform = prettyForm(*pform.left(self._print(expr.dom)))

        return pform

    def _print_FractionField(self, expr):
        pform = self._print_seq(expr.gens, '(', ')')
        pform = prettyForm(*pform.left(self._print(expr.dom)))

        return pform

    def _print_Subs(self, e):
        pform = self._print(e.expr)
        pform = prettyForm(*pform.parens())

        h = pform.height() if pform.height() > 1 else 2
        rvert = stringPict(vobj('|', h), baseline=pform.baseline)
        pform = prettyForm(*pform.right(rvert))

        b = pform.baseline
        pform.baseline = pform.height() - 1
        pform = prettyForm(*pform.right(self._print_seq([
            self._print_seq((self._print(v[0]), xsym('=='), self._print(v[1])),
                delimiter='') for v in zip(e.variables, e.point) ])))

        pform.baseline = b
        return pform

def pretty(expr, **settings):
    """
    Returns a string containing the prettified form of expr.

    Arguments
    ---------
    expr: the expression to print
    wrap_line: line wrapping enabled/disabled, should be a boolean value (default to True)
    use_unicode: use unicode characters, such as the Greek letter pi instead of
        the string pi. Values should be boolean or None
    full_prec: use full precision. Default to "auto"
    """
    pp = PrettyPrinter(settings)

    # XXX: this is an ugly hack, but at least it works
    use_unicode = pp._settings['use_unicode']
    uflag = pretty_use_unicode(use_unicode)

    try:
        return pp.doprint(expr)
    finally:
        pretty_use_unicode(uflag)

def pretty_print(expr, **settings):
    """
    Prints expr in pretty form.

    pprint is just a shortcut for this function
    """
    print pretty(expr, **settings)

pprint = pretty_print

