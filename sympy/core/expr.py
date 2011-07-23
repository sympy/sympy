from core import C
from basic import Basic, Atom
from singleton import S
from evalf import EvalfMixin
from decorators import _sympifyit, call_highest_priority
from cache import cacheit
from compatibility import reduce

class Expr(Basic, EvalfMixin):
    __slots__ = []

    def sort_key(self, order=None):
        # XXX: The order argument does not actually work
        from sympy.core import S

        def key_inner(arg):
            if isinstance(arg, Basic):
                return arg.sort_key(order=order)
            elif hasattr(arg, '__iter__'):
                return tuple(key_inner(arg) for arg in args)
            else:
                return arg

        coeff, expr = self.as_coeff_Mul()
        if expr.is_Pow:
            expr, exp = expr.args
        else:
            expr, exp = expr, S.One

        if expr.is_Atom:
            if expr.is_Symbol:
                args = (str(expr),)
            else:
                args = (expr,)
        else:
            if expr.is_Add:
                args = expr.as_ordered_terms(order=order)
            else:
                args = expr.args

            args = tuple(key_inner(arg) for arg in args)

            if expr.is_Mul:
                args = sorted(args)

        args = (len(args), args)
        exp = exp.sort_key(order=order)

        return expr.class_key(), args, exp, coeff


    # ***************
    # * Arithmetics *
    # ***************

    # Expr and its sublcasses use _op_priority to determine which object
    # passed to a binary special method (__mul__, etc.) will handle the
    # operation. In general, the 'call_highest_priority' decorator will choose
    # the object with the highest _op_priority to handle the call.
    # Custom subclasses that want to define their own binary special methods
    # should set an _op_priority value that is higher than the default.
    #
    # **NOTE**:
    # This is a temporary fix, and will eventually be replaced with
    # something better and more powerful.  See issue 2411.
    _op_priority = 10.0

    def __pos__(self):
        return self
    def __neg__(self):
        return Mul(S.NegativeOne, self)
    def __abs__(self):
        return C.Abs(self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return Add(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return Add(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return Add(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return Add(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return Mul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return Mul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return Pow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        return Pow(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return Mul(self, Pow(other, S.NegativeOne))
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        return Mul(other, Pow(self, S.NegativeOne))

    __truediv__ = __div__
    __rtruediv__ = __rdiv__


    def __float__(self):
        result = self.evalf()
        if result.is_Number:
            return float(result)
        else:
            raise ValueError("Symbolic value, can't compute")

    def __complex__(self):
        result = self.evalf()
        re, im = result.as_real_imag()
        return complex(float(re), float(im))

    @_sympifyit('other', False) # sympy >  other
    def __lt__(self, other):
        dif = self - other
        if dif.is_negative != dif.is_nonnegative:
            return dif.is_negative
        return C.StrictInequality(self, other)

    @_sympifyit('other', True)  # sympy >  other
    def __gt__(self, other):
        dif = self - other
        if dif.is_positive !=  dif.is_nonpositive:
            return dif.is_positive
        return C.StrictInequality(other, self)

    @_sympifyit('other', False) # sympy >  other
    def __le__(self, other):
        dif = self - other
        if dif.is_nonpositive != dif.is_positive:
            return dif.is_nonpositive
        return C.Inequality(self, other)

    @_sympifyit('other', True)  # sympy >  other
    def __ge__(self, other):
        dif = self - other
        if dif.is_nonnegative != dif.is_negative:
            return dif.is_nonnegative
        return C.Inequality(other, self)

    @staticmethod
    def _from_mpmath(x, prec):
        if hasattr(x, "_mpf_"):
            return C.Float._new(x._mpf_, prec)
        elif hasattr(x, "_mpc_"):
            re, im = x._mpc_
            re = C.Float._new(re, prec)
            im = C.Float._new(im, prec)*S.ImaginaryUnit
            return re+im
        else:
            raise TypeError("expected mpmath number (mpf or mpc)")

    @property
    def is_number(self):
        """Returns True if 'self' is a number.

           >>> from sympy import log, Integral
           >>> from sympy.abc import x, y

           >>> x.is_number
           False
           >>> (2*x).is_number
           False
           >>> (2 + log(2)).is_number
           True
           >>> (2 + Integral(2, x)).is_number
           False
           >>> (2 + Integral(2, (x, 1, 2))).is_number
           True

        """
        if not self.args:
            return False
        return all(obj.is_number for obj in self.iter_basic_args())

    def _eval_interval(self, x, a, b):
        """
        Returns evaluation over an interval.  For most functions this is:

        self.subs(x, b) - self.subs(x, a),

        possibly using limit() if NaN is returned from subs.

        If b or a is None, it only evaluates -self.subs(x, a) or self.subs(b, x),
        respectively.

        """
        from sympy.series import limit
        if (a is None and b is None):
            raise ValueError('Both interval ends cannot be None.')

        if a is None:
            A = 0
        else:
            A = self.subs(x, a)
            if A is S.NaN:
                A = limit(self, x, a)
                if A is S.NaN:
                    return A

        if b is None:
            B = 0
        else:
            B = self.subs(x, b)
            if B is S.NaN:
                B = limit(self, x, b)

        return B - A

    def _eval_power(self, other):
        return None

    def _eval_conjugate(self):
        if self.is_real:
            return self

    def conjugate(self):
        from sympy.functions.elementary.complexes import conjugate as c
        return c(self)


    @classmethod
    def _parse_order(cls, order):
        """Parse and configure the ordering of terms. """
        from sympy.polys.monomialtools import monomial_key

        try:
            reverse = order.startswith('rev-')
        except AttributeError:
            reverse = False
        else:
            if reverse:
                order = order[4:]

        monom_key = monomial_key(order)

        def neg(monom):
            result = []

            for m in monom:
                if isinstance(m, tuple):
                    result.append(neg(m))
                else:
                    result.append(-m)

            return tuple(result)

        def key(term):
            _, ((re, im), monom, ncpart) = term

            monom = neg(monom_key(monom))
            ncpart = tuple([ e.sort_key(order=order) for e in ncpart ])
            coeff = ((bool(im), im), (re, im))

            return monom, ncpart, coeff

        return key, reverse

    def as_ordered_factors(self, order=None):
        """
        Transform an expression to an ordered list of factors.

        **Examples**

        >>> from sympy import sin, cos
        >>> from sympy.abc import x, y

        >>> (2*x*y*sin(x)*cos(x)).as_ordered_factors()
        [2, x, y, sin(x), cos(x)]

        """
        if not self.is_Mul:
            return [self]

        cpart = []
        ncpart = []

        for arg in self.args:
            if arg.is_commutative:
                cpart.append(arg)
            else:
                ncpart.append(arg)

        return sorted(cpart, key=lambda expr: expr.sort_key(order=order)) + ncpart

    def as_ordered_terms(self, order=None, data=False):
        """
        Transform an expression to an ordered list of terms.

        **Examples**

        >>> from sympy import sin, cos
        >>> from sympy.abc import x, y

        >>> (sin(x)**2*cos(x) + sin(x)**2 + 1).as_ordered_terms()
        [sin(x)**2*cos(x), sin(x)**2, 1]

        """
        key, reverse = self._parse_order(order)
        terms, gens = self.as_terms()

        if not any(term.is_Order for term, _ in terms):
            ordered = sorted(terms, key=key, reverse=reverse)
        else:
            _terms, _order = [], []

            for term, repr in terms:
                if not term.is_Order:
                    _terms.append((term, repr))
                else:
                    _order.append((term, repr))

            ordered = sorted(_terms, key=key, reverse=True) \
                    + sorted(_order, key=key, reverse=True)

        if data:
            return ordered, gens
        else:
            return [ term for term, _ in ordered ]

    def as_terms(self):
        """Transform an expression to a list of terms. """
        from sympy.core import Add, Mul, S
        from sympy.core.exprtools import decompose_power
        from sympy.utilities import default_sort_key

        gens, terms = set([]), []

        for term in Add.make_args(self):
            coeff, _term = term.as_coeff_Mul()

            coeff = complex(coeff)
            cpart, ncpart = {}, []

            if _term is not S.One:
                for factor in Mul.make_args(_term):
                    if factor.is_number:
                        try:
                            coeff *= complex(factor)
                        except ValueError:
                            pass
                        else:
                            continue

                    if factor.is_commutative:
                        base, exp = decompose_power(factor)

                        cpart[base] = exp
                        gens.add(base)
                    else:
                        ncpart.append(factor)

            coeff = coeff.real, coeff.imag
            ncpart = tuple(ncpart)

            terms.append((term, (coeff, cpart, ncpart)))

        gens = sorted(gens, key=default_sort_key)

        k, indices = len(gens), {}

        for i, g in enumerate(gens):
            indices[g] = i

        result = []

        for term, (coeff, cpart, ncpart) in terms:
            monom = [0]*k

            for base, exp in cpart.iteritems():
                monom[indices[base]] = exp

            result.append((term, (coeff, tuple(monom), ncpart)))

        return result, gens


    def removeO(self):
        """Removes the additive O(..) symbol if there is one"""
        return self

    def getO(self):
        """Returns the additive O(..) symbol if there is one, else None."""
        return None

    def getn(self):
        """
        Returns the order of the expression.

        The order is determined either from the O(...) term. If there
        is no O(...) term, it returns None.

        Example:

        >>> from sympy import O
        >>> from sympy.abc import x
        >>> (1 + x + O(x**2)).getn()
        2
        >>> (1 + x).getn()

        """
        o = self.getO()
        if o is None:
            return None
        elif o.is_Order:
            o = o.expr
            if o is S.One:
                return S.Zero
            if o.is_Symbol:
                return S.One
            if o.is_Pow:
                return o.args[1]
            if o.is_Mul: # x**n*log(x)**n or x**n/log(x)**n
                for oi in o.args:
                    if oi.is_Symbol:
                        return S.One
                    if oi.is_Pow:
                        syms = oi.atoms(C.Symbol)
                        if len(syms) == 1:
                            x = syms.pop()
                            oi = oi.subs(x, C.Dummy('x', positive=True))
                            if oi.base.is_Symbol and oi.exp.is_Rational:
                                return abs(oi.exp)

        raise NotImplementedError('not sure of order of %s' % o)

    def count_ops(self, visual=None):
        """wrapper for count_ops that returns the operation count."""
        from sympy import count_ops
        return count_ops(self, visual)

    def args_cnc(self):
        """treat self as Mul and split it into tuple (set, list)
        where ``set`` contains the commutative parts and ``list`` contains
        the ordered non-commutative args.

        A special treatment is that -1 is separated from a Rational:

        >>> from sympy import symbols
        >>> A, B = symbols('A B', commutative=0)
        >>> x, y = symbols('x y')
        >>> (-2*x*y).args_cnc()
        [set([-1, 2, x, y]), []]
        >>> (-2*x*A*B*y).args_cnc()
        [set([-1, 2, x, y]), [A, B]]

        The arg is treated as a Mul:

        >>> (-2 + x + A).args_cnc()
        [set(), [x - 2 + A]]
        """

        if self.is_Mul:
            args = list(self.args)
        else:
            args = [self]
        for i, mi in enumerate(args):
            if not mi.is_commutative:
                c = args[:i]
                nc = args[i:]
                break
        else:
            c = args
            nc = []

        if c and c[0].is_Rational and c[0].is_negative and c[0] != S.NegativeOne:
            c[:1] = [S.NegativeOne, -c[0]]

        return [set(c), nc]

    def coeff(self, x, right=False):
        """
        Returns the coefficient of the exact term "x" or None if there is no "x".

        When x is noncommutative, the coeff to the left (default) or right of x
        can be returned. The keyword 'right' is ignored when x is commutative.

        Examples::

        >>> from sympy import symbols
        >>> from sympy.abc import x, y, z

        You can select terms that have an explicit negative in front of them:

        >>> (-x+2*y).coeff(-1)
        x
        >>> (x-2*y).coeff(-1)
        2*y

        You can select terms with no rational coefficient:

        >>> (x+2*y).coeff(1)
        x
        >>> (3+2*x+4*x**2).coeff(1)

        You can select terms that have a numerical term in front of them:

        >>> (-x-2*y).coeff(2)
        -y
        >>> from sympy import sqrt
        >>> (x+sqrt(2)*x).coeff(sqrt(2))
        x

        The matching is exact:

        >>> (3+2*x+4*x**2).coeff(x)
        2
        >>> (3+2*x+4*x**2).coeff(x**2)
        4
        >>> (3+2*x+4*x**2).coeff(x**3)
        >>> (z*(x+y)**2).coeff((x+y)**2)
        z
        >>> (z*(x+y)**2).coeff(x+y)

        In addition, no factoring is done, so 2 + y is not obtained from the
        following:

        >>> (2*x+2+(x+1)*y).coeff(x+1)
        y

        >>> n, m, o = symbols('n m o', commutative=False)
        >>> n.coeff(n)
        1
        >>> (3*n).coeff(n)
        3
        >>> (n*m + m*n*m).coeff(n) # = (1 + m)*n*m
        1 + m
        >>> (n*m + m*n*m).coeff(n, right=True) # = (1 + m)*n*m
        m

        If there is more than one possible coefficient None is returned:

        >>> (n*m + m*n).coeff(n)

        If there is only one possible coefficient, it is returned:

        >>> (n*m + o*m*n).coeff(m*n)
        o
        >>> (n*m + o*m*n).coeff(m*n, right=1)
        1
        """
        x = sympify(x)
        if not x: # 0 or None
            return None
        if x == self:
            return S.One
        if x is S.One:
            try:
                assert Add.make_args(S.Zero) and Mul.make_args(S.One)
                # replace try/except with this
                co = [a for a in Add.make_args(self)
                      if not any(ai.is_number for ai in Mul.make_args(a))]
            except AssertionError:
                co = [a for a in (Add.make_args(self) or [S.Zero])
                      if not any(ai.is_number for ai in (Mul.make_args(a) or [S.One]))]
            if not co:
                return None
            return Add(*co)

        def incommon(l1, l2):
            if not l1 or not l2:
                return []
            n = min(len(l1), len(l2))
            for i in xrange(n):
                if l1[i] != l2[i]:
                    return l1[:i]
            return l1[:]

        def arglist(x):
            """ Return list of x's args when treated as a Mul after checking
            to see if a negative Rational is present (in which case it is made
            positive and a -1 is added to the list).
            """

            margs = list(Mul.make_args(x))
            try:
                assert Mul.make_args(S.One)
                # replace try/except with the following
                if margs[0].is_Rational and margs[0].is_negative and margs[0] != S.NegativeOne:
                    margs.append(S.NegativeOne)
                    margs[0] *= -1
            except AssertionError:
                if margs and margs[0].is_Rational and margs[0].is_negative and margs[0] != S.NegativeOne:
                    margs.append(S.NegativeOne)
                    margs[0] *= -1
            return margs

        def find(l, sub, first=True):
            """ Find where list sub appears in list l. When ``first`` is True
            the first occurance from the left is returned, else the last
            occurance is returned. Return None if sub is not in l.

            >> l = range(5)*2
            >> find(l, [2, 3])
            2
            >> find(l, [2, 3], first=0)
            7
            >> find(l, [2, 4])
            None

            """
            if not sub or not l or len(sub) > len(l):
                return None
            n = len(sub)
            if not first:
                l.reverse()
                sub.reverse()
            for i in xrange(0, len(l) - n + 1):
                if all(l[i + j] == sub[j] for j in range(n)):
                    break
            else:
                i = None
            if not first:
                l.reverse()
                sub.reverse()
            if i is not None and not first:
                i = len(l) - (i + n)
            return i

        co = []
        try:
            assert Add.make_args(S.Zero)
            # replace try/except with this
            args = Add.make_args(self)
        except AssertionError:
            args = Add.make_args(self) or [S.Zero]
        self_c = self.is_commutative
        x_c = x.is_commutative
        if self_c and not x_c:
            return None

        if self_c:
            xargs = set(arglist(x))
            for a in Add.make_args(self):
                margs = set(arglist(a))
                if len(xargs) > len(margs):
                    continue
                resid = margs.difference(xargs)
                if len(resid) + len(xargs) == len(margs):
                    co.append(Mul(*resid))
            if co == []:
                return None
            elif co:
                return Add(*co)
        elif x_c:
            xargs = set(arglist(x))
            for a in Add.make_args(self):
                margs, nc = a.args_cnc()
                if len(xargs) > len(margs):
                    continue
                resid = margs.difference(xargs)
                if len(resid) + len(xargs) == len(margs):
                    co.append(Mul(*(list(resid) + nc)))
            if co == []:
                return None
            elif co:
                return Add(*co)
        else: # both nc
            xargs, nx = x.args_cnc()
            # find the parts that pass the commutative terms
            for a in Add.make_args(self):
                margs, nc = a.args_cnc()
                if len(xargs) > len(margs):
                    continue
                resid = margs.difference(xargs)
                if len(resid) + len(xargs) == len(margs):
                    co.append((resid, nc))
            # now check the non-comm parts
            if not co:
                return None
            if all(n == co[0][1] for r, n in co):
                ii = find(co[0][1], nx, right)
                if not ii is None:
                    if not right:
                        return Mul(Add(*[Mul(*r) for r, c in co]), Mul(*co[0][1][:ii]))
                    else:
                        return Mul(*co[0][1][ii+len(nx):])
            beg = reduce(incommon, (n[1] for n in co))
            if beg:
                ii = find(beg, nx, right)
                if not ii is None:
                    if not right:
                        gcdc = co[0][0]
                        for i in xrange(1, len(co)):
                            gcdc = gcdc.intersection(co[i][0])
                            if not gcdc:
                                break
                        return Mul(*(list(gcdc) + beg[:ii]))
                    else:
                        m = ii + len(nx)
                        return Add(*[Mul(*(list(r) + n[m:])) for r, n in co])
            end = list(reversed(reduce(incommon, (list(reversed(n[1])) for n in co))))
            if end:
                ii = find(end, nx, right)
                if not ii is None:
                    if not right:
                        return Add(*[Mul(*(list(r) + n[:-len(end)+ii])) for r, n in co])
                    else:
                        return Mul(*end[ii+len(nx):])
            # look for single match
            hit = None
            for i, (r, n) in enumerate(co):
                ii = find(n, nx, right)
                if not ii is None:
                    if not hit:
                        hit = ii, r, n
                    else:
                        break
            else:
                if hit:
                    ii, r, n = hit
                    if not right:
                        return Mul(*(list(r) + n[:ii]))
                    else:
                        return Mul(*n[ii+len(nx):])

            return None

    def as_expr(self, *gens):
        """
        Convert a polynomial to a SymPy expression.

        **Examples**

        >>> from sympy import sin
        >>> from sympy.abc import x, y

        >>> f = (x**2 + x*y).as_poly(x, y)
        >>> f.as_expr()
        x**2 + x*y

        >>> sin(x).as_expr()
        sin(x)

        """
        return self

    def as_coefficient(self, expr):
        """Extracts symbolic coefficient at the given expression. In
           other words, this functions separates 'self' into product
           of 'expr' and 'expr'-free coefficient. If such separation
           is not possible it will return None.

           >>> from sympy import E, pi, sin, I, symbols
           >>> from sympy.abc import x, y

           >>> E.as_coefficient(E)
           1
           >>> (2*E).as_coefficient(E)
           2
           >>> (2*sin(E)*E).as_coefficient(E)

           >>> (2*E + x*E).as_coefficient(E)
           x + 2
           >>> (2*E*x + x).as_coefficient(E)

           >>> (E*(x + 1) + x).as_coefficient(E)

           >>> (2*pi*I).as_coefficient(pi*I)
           2
           >>> (2*I).as_coefficient(pi*I)

        """

        r = self.extract_multiplicatively(expr)
        if r and not r.has(expr):
            return r

    def as_independent(self, *deps, **hint):
        """
        A mostly naive separation of a Mul or Add into arguments that are not
        are dependent on deps. To obtain as complete a separation of variables
        as possible, use a separation method first, e.g.:

        * separatevars() to change Mul, Add and Pow (including exp) into Mul
        * .expand(mul=True) to change Add or Mul into Add
        * .expand(log=True) to change log expr into an Add

        The only non-naive thing that is done here is to respect noncommutative
        ordering of variables.

        The returned tuple (i, d) has the following interpretation:

        * i will has no variable that appears in deps
        * d will be 1 or else have terms that contain variables that are in deps
        * if self is an Add then self = i + d
        * if self is a Mul then self = i*d
        * if self is anything else, either tuple (self, S.One) or (S.One, self)
          is returned.

        To force the expression to be treated as an Add, use the hint as_Add=True

        Examples:

        -- self is an Add

        >>> from sympy import sin, cos, exp
        >>> from sympy.abc import x, y, z

        >>> (x + x*y).as_independent(x)
        (0, x*y + x)
        >>> (x + x*y).as_independent(y)
        (x, x*y)
        >>> (2*x*sin(x) + y + x + z).as_independent(x)
        (y + z, 2*x*sin(x) + x)
        >>> (2*x*sin(x) + y + x + z).as_independent(x, y)
        (z, 2*x*sin(x) + x + y)

        -- self is a Mul

        >>> (x*sin(x)*cos(y)).as_independent(x)
        (cos(y), x*sin(x))

        non-commutative terms cannot always be separated out when self is a Mul

        >>> from sympy import symbols
        >>> n1, n2, n3 = symbols('n1 n2 n3', commutative=False)
        >>> (n1 + n1*n2).as_independent(n2)
        (n1, n1*n2)
        >>> (n2*n1 + n1*n2).as_independent(n2)
        (0, n1*n2 + n2*n1)
        >>> (n1*n2*n3).as_independent(n1)
        (1, n1*n2*n3)
        >>> (n1*n2*n3).as_independent(n2)
        (n1, n2*n3)
        >>> ((x-n1)*(x-y)).as_independent(x)
        (1, (x - y)*(x - n1))

        -- self is anything else:

        >>> (sin(x)).as_independent(x)
        (1, sin(x))
        >>> (sin(x)).as_independent(y)
        (sin(x), 1)
        >>> exp(x+y).as_independent(x)
        (1, exp(x + y))

        -- force self to be treated as an Add:

        >>> (3*x).as_independent(x, as_Add=1)
        (0, 3*x)

        -- force self to be treated as a Mul:

        >>> (3+x).as_independent(x, as_Add=0)
        (1, x + 3)
        >>> (-3+x).as_independent(x, as_Add=0)
        (1, x - 3)

        Note how the below differs from the above in making the
        constant on the dep term positive.

        >>> (y*(-3+x)).as_independent(x)
        (y, x - 3)

        Note: when trying to get independent terms, a separation method
        might need to be used first. In this case, it is important to keep
        track of what you send to this routine so you know how to interpret
        the returned values

        >>> from sympy import separatevars, log
        >>> separatevars(exp(x+y)).as_independent(x)
        (exp(y), exp(x))
        >>> (x + x*y).as_independent(y)
        (x, x*y)
        >>> separatevars(x + x*y).as_independent(y)
        (x, y + 1)
        >>> (x*(1 + y)).as_independent(y)
        (x, y + 1)
        >>> (x*(1 + y)).expand(mul=True).as_independent(y)
        (x, x*y)
        >>> a, b=symbols('a b',positive=True)
        >>> (log(a*b).expand(log=True)).as_independent(b)
        (log(a), log(b))

        See also: .separatevars(), .expand(log=True),
                  .as_two_terms(), .as_coeff_add(), .as_coeff_mul()
        """
        from sympy.utilities.iterables import sift

        func = self.func
        if hint.get('as_Add', func is Add):
            want = Add
        else:
            want = Mul
        if (want is not func or
            func is not Add and func is not Mul):
            if self.has(*deps):
                return (want.identity, self)
            else:
                return (self, want.identity)
        else:
            if func is Add:
                args = list(self.args)
            else:
                args, nc = self.args_cnc()

        d = sift(args, lambda x: x.has(*deps))
        depend = d.pop(True, [])
        indep = d.pop(False, [])
        if func is Add: # all terms were treated as commutative
            return (Add(*indep),
                    Add(*depend))
        else: # handle noncommutative by stopping at first dependent term
            for i, n in enumerate(nc):
                if n.has(*deps):
                    depend.extend(nc[i:])
                    break
                indep.append(n)
            return Mul(*indep), Mul(*depend)

    def as_real_imag(self, deep=True):
        """Performs complex expansion on 'self' and returns a tuple
           containing collected both real and imaginary parts. This
           method can't be confused with re() and im() functions,
           which does not perform complex expansion at evaluation.

           However it is possible to expand both re() and im()
           functions and get exactly the same results as with
           a single call to this function.

           >>> from sympy import symbols, I

           >>> x, y = symbols('x,y', real=True)

           >>> (x + y*I).as_real_imag()
           (x, y)

           >>> from sympy.abc import z, w

           >>> (z + w*I).as_real_imag()
           (-im(w) + re(z), im(z) + re(w))

        """
        return (C.re(self), C.im(self))

    def as_powers_dict(self):
        return dict([self.as_base_exp()])

    def as_base_exp(self):
        # a -> b ** e
        return self, S.One

    def as_coeff_terms(self, *deps):
        import warnings
        warnings.warn("\nuse as_coeff_mul() instead of as_coeff_terms().",
                      DeprecationWarning)

    def as_coeff_factors(self, *deps):
        import warnings
        warnings.warn("\nuse as_coeff_add() instead of as_coeff_factors().",
                      DeprecationWarning)

    def as_coeff_mul(self, *deps):
        """Return the tuple (c, args) where self is written as a Mul, ``m``.

        c should be a Rational multiplied by any terms of the Mul that are
        independent of deps.

        args should be a tuple of all other terms of m; args is empty
        if self is a Number or if self is independent of deps (when given).

        This should be used when you don't know if self is a Mul or not but
        you want to treat self as a Mul or if you want to process the
        individual arguments of the tail of self as a Mul.

        - if you know self is a Mul and want only the head, use self.args[0];
        - if you don't want to process the arguments of the tail but need the
          tail then use self.as_two_terms() which gives the head and tail;
        - if you want to split self into an independent and dependent parts
          use self.as_independent(\*deps)

        >>> from sympy import S
        >>> from sympy.abc import x, y
        >>> (S(3)).as_coeff_mul()
        (3, ())
        >>> (3*x*y).as_coeff_mul()
        (3, (x, y))
        >>> (3*x*y).as_coeff_mul(x)
        (3*y, (x,))
        >>> (3*y).as_coeff_mul(x)
        (3*y, ())
        """
        if deps:
            if not self.has(*deps):
                return self, tuple()
        return S.One, (self,)

    def as_coeff_add(self, *deps):
        """Return the tuple (c, args) where self is written as an Add, ``a``.

        c should be a Rational added to any terms of the Add that are
        independent of deps.

        args should be a tuple of all other terms of ``a``; args is empty
        if self is a Number or if self is independent of deps (when given).

        This should be used when you don't know if self is an Add or not but
        you want to treat self as an Add or if you want to process the
        individual arguments of the tail of self as an Add.

        - if you know self is an Add and want only the head, use self.args[0];
        - if you don't want to process the arguments of the tail but need the
          tail then use self.as_two_terms() which gives the head and tail.
        - if you want to split self into an independent and dependent parts
          use self.as_independent(\*deps)

        >>> from sympy import S
        >>> from sympy.abc import x, y
        >>> (S(3)).as_coeff_add()
        (3, ())
        >>> (3 + x + y).as_coeff_add()
        (3, (y, x))
        >>> (3 + x +y).as_coeff_add(x)
        (y + 3, (x,))
        >>> (3 + y).as_coeff_add(x)
        (y + 3, ())

        """
        if deps:
            if not self.has(*deps):
                return self, tuple()
        return S.Zero, (self,)

    def as_numer_denom(self):
        """ a/b -> a,b

        This is just a stub that should be defined by
        an object's class methods to get anything else."""

        return self, S.One

    def normal(self):
        n, d = self.as_numer_denom()
        if d is S.One:
            return n
        return n/d

    def extract_multiplicatively(self, c):
        """Return None if it's not possible to make self in the form
           c * something in a nice way, i.e. preserving the properties
           of arguments of self.

           >>> from sympy import symbols, Rational

           >>> x, y = symbols('x,y', real=True)

           >>> ((x*y)**3).extract_multiplicatively(x**2 * y)
           x*y**2

           >>> ((x*y)**3).extract_multiplicatively(x**4 * y)

           >>> (2*x).extract_multiplicatively(2)
           x

           >>> (2*x).extract_multiplicatively(3)

           >>> (Rational(1,2)*x).extract_multiplicatively(3)
           x/6

        """
        c = sympify(c)
        if c is S.One:
            return self
        elif c == self:
            return S.One
        elif c.is_Mul:
            x = self.extract_multiplicatively(c.as_two_terms()[0])
            if x != None:
                return x.extract_multiplicatively(c.as_two_terms()[1])
        quotient = self / c
        if self.is_Number:
            if self is S.Infinity:
                if c.is_positive:
                    return S.Infinity
            elif self is S.NegativeInfinity:
                if c.is_negative:
                    return S.Infinity
                elif c.is_positive:
                    return S.NegativeInfinity
            elif self is S.ComplexInfinity:
                if not c.is_zero:
                    return S.ComplexInfinity
            elif self is S.NaN:
                return S.NaN
            elif self.is_Integer:
                if not quotient.is_Integer:
                    return None
                elif self.is_positive and quotient.is_negative:
                    return None
                else:
                    return quotient
            elif self.is_Rational:
                if not quotient.is_Rational:
                    return None
                elif self.is_positive and quotient.is_negative:
                    return None
                else:
                    return quotient
            elif self.is_Float:
                if not quotient.is_Float:
                    return None
                elif self.is_positive and quotient.is_negative:
                    return None
                else:
                    return quotient
        elif self.is_NumberSymbol or self.is_Symbol or self is S.ImaginaryUnit:
            if quotient.is_Mul and len(quotient.args) == 2:
                if quotient.args[0].is_Integer and quotient.args[0].is_positive and quotient.args[1] == self:
                    return quotient
            elif quotient.is_Integer:
                return quotient
        elif self.is_Add:
            newargs = []
            for arg in self.args:
                newarg = arg.extract_multiplicatively(c)
                if newarg != None:
                    newargs.append(newarg)
                else:
                    return None
            return Add(*newargs)
        elif self.is_Mul:
            for i in xrange(len(self.args)):
                newargs = list(self.args)
                del(newargs[i])
                tmp = self._new_rawargs(*newargs).extract_multiplicatively(c)
                if tmp != None:
                    return tmp * self.args[i]
        elif self.is_Pow:
            if c.is_Pow and c.base == self.base:
                new_exp = self.exp.extract_additively(c.exp)
                if new_exp != None:
                    return self.base ** (new_exp)
            elif c == self.base:
                new_exp = self.exp.extract_additively(1)
                if new_exp != None:
                    return self.base ** (new_exp)

    def extract_additively(self, c):
        """Return None if it's not possible to make self in the form
           something + c in a nice way, i.e. preserving the properties
           of arguments of self.

           >>> from sympy import symbols

           >>> x, y = symbols('x,y', real=True)

           >>> ((x*y)**3).extract_additively(1)

           >>> (x+1).extract_additively(x)
           1

           >>> (x+1).extract_additively(2*x)

           >>> (x+1).extract_additively(-x)
           2*x + 1

           >>> (-x+1).extract_additively(2*x)
           -3*x + 1

        """
        c = sympify(c)
        if c is S.Zero:
            return self
        elif c == self:
            return S.Zero
        elif self is S.Zero:
            return None
        elif c.is_Add:
            x = self.extract_additively(c.as_two_terms()[0])
            if x != None:
                return x.extract_additively(c.as_two_terms()[1])
        sub = self - c
        if self.is_Number:
            if self.is_Integer:
                if not sub.is_Integer:
                    return None
                elif self.is_positive and sub.is_negative:
                    return None
                else:
                    return sub
            elif self.is_Rational:
                if not sub.is_Rational:
                    return None
                elif self.is_positive and sub.is_negative:
                    return None
                else:
                    return sub
            elif self.is_Float:
                if not sub.is_Float:
                    return None
                elif self.is_positive and sub.is_negative:
                    return None
                else:
                    return sub
        elif self.is_NumberSymbol or self.is_Symbol or self is S.ImaginaryUnit:
            if sub.is_Mul and len(sub.args) == 2:
                if sub.args[0].is_Integer and sub.args[0].is_positive and sub.args[1] == self:
                    return sub
            elif sub.is_Integer:
                return sub
        elif self.is_Add:
            terms = self.as_two_terms()
            subs0 = terms[0].extract_additively(c)
            if subs0 != None:
                return subs0 + terms[1]
            else:
                subs1 = terms[1].extract_additively(c)
                if subs1 != None:
                    return subs1 + terms[0]
        elif self.is_Mul:
            self_coeff, self_terms = self.as_coeff_mul()
            if c.is_Mul:
                c_coeff, c_terms = c.as_coeff_mul()
                if c_terms == self_terms:
                    new_coeff = self_coeff.extract_additively(c_coeff)
                    if new_coeff != None:
                        return new_coeff * c._new_rawargs(*c_terms)
            elif c == self_terms:
                new_coeff = self_coeff.extract_additively(1)
                if new_coeff != None:
                    return new_coeff * c

    def could_extract_minus_sign(self):
        """Canonical way to choose an element in the set {e, -e} where
           e is any expression. If the canonical element is e, we have
           e.could_extract_minus_sign() == True, else
           e.could_extract_minus_sign() == False.

           For any expression, the set ``{e.could_extract_minus_sign(),
           (-e).could_extract_minus_sign()}`` must be ``{True, False}``.

           >>> from sympy.abc import x, y
           >>> (x-y).could_extract_minus_sign() != (y-x).could_extract_minus_sign()
           True

        """
        negative_self = -self
        self_has_minus = (self.extract_multiplicatively(-1) != None)
        negative_self_has_minus = ((negative_self).extract_multiplicatively(-1) != None)
        if self_has_minus != negative_self_has_minus:
            return self_has_minus
        else:
            if self.is_Add:
                # We choose the one with less arguments with minus signs
                all_args = len(self.args)
                negative_args = len([False for arg in self.args if arg.could_extract_minus_sign()])
                positive_args = all_args - negative_args
                if positive_args > negative_args:
                    return False
                elif positive_args < negative_args:
                    return True
            elif self.is_Mul:
                # We choose the one with an odd number of minus signs
                num, den = self.as_numer_denom()
                args = Mul.make_args(num) + Mul.make_args(den)
                arg_signs = [arg.could_extract_minus_sign() for arg in args]
                negative_args = filter(None, arg_signs)
                return len(negative_args) % 2 == 1

            # As a last resort, we choose the one with greater hash
            return hash(self) < hash(negative_self)

    def _eval_is_polynomial(self, syms):
        if self.free_symbols.intersection(syms) == set([]):
            return True
        return False

    def is_polynomial(self, *syms):
        """
        Return True if self is a polynomial in syms and False otherwise.

        This checks if self is an exact polynomial in syms.  This function
        returns False for expressions that are "polynomials" with symbolic
        exponents.  Thus, you should be able to apply polynomial algorithms to
        expressions for which this returns True, and Poly(expr, \*syms) should
        work only if and only if expr.is_polynomial(\*syms) returns True. The
        polynomial does not have to be in expanded form.  If no symbols are
        given, all free symbols in the expression will be used.

        This is not part of the assumptions system.  You cannot do
        Symbol('z', polynomial=True).

        **Examples**

        >>> from sympy import Symbol
        >>> x = Symbol('x')
        >>> ((x**2 + 1)**4).is_polynomial(x)
        True
        >>> ((x**2 + 1)**4).is_polynomial()
        True
        >>> (2**x + 1).is_polynomial(x)
        False


        >>> n = Symbol('n', nonnegative=True, integer=True)
        >>> (x**n + 1).is_polynomial(x)
        False

        This function does not attempt any nontrivial simplifications that may
        result in an expression that does not appear to be a polynomial to
        become one.

        >>> from sympy import sqrt, factor, cancel
        >>> y = Symbol('y', positive=True)
        >>> a = sqrt(y**2 + 2*y + 1)
        >>> a.is_polynomial(y)
        False
        >>> factor(a)
        y + 1
        >>> factor(a).is_polynomial(y)
        True

        >>> b = (y**2 + 2*y + 1)/(y + 1)
        >>> b.is_polynomial(y)
        False
        >>> cancel(b)
        y + 1
        >>> cancel(b).is_polynomial(y)
        True

        See also .is_rational_function()

        """
        if syms:
            syms = set(map(sympify, syms))
        else:
            syms = self.free_symbols

        if syms.intersection(self.free_symbols) == set([]):
            # constant polynomial
            return True
        else:
            return self._eval_is_polynomial(syms)

    def _eval_is_rational_function(self, syms):
        if self.free_symbols.intersection(syms) == set([]):
            return True
        return False

    def is_rational_function(self, *syms):
        """
        Test whether function is a ratio of two polynomials in the given
        symbols, syms. When syms is not given, all free symbols will be used.
        The rational function does not have to be in expanded or in any kind of
        canonical form.

        This function returns False for expressions that are "rational
        functions" with symbolic exponents.  Thus, you should be able to call
        .as_numer_denom() and apply polynomial algorithms to the result for
        expressions for which this returns True.

        This is not part of the assumptions system.  You cannot do
        Symbol('z', rational_function=True).

        Example:

        >>> from sympy import Symbol, sin
        >>> from sympy.abc import x, y

        >>> (x/y).is_rational_function()
        True

        >>> (x**2).is_rational_function()
        True

        >>> (x/sin(y)).is_rational_function(y)
        False

        >>> n = Symbol('n', integer=True)
        >>> (x**n + 1).is_rational_function(x)
        False

        This function does not attempt any nontrivial simplifications that may
        result in an expression that does not appear to be a rational function
        to become one.

        >>> from sympy import sqrt, factor, cancel
        >>> y = Symbol('y', positive=True)
        >>> a = sqrt(y**2 + 2*y + 1)/y
        >>> a.is_rational_function(y)
        False
        >>> factor(a)
        (y + 1)/y
        >>> factor(a).is_rational_function(y)
        True

        See also is_rational_function().

        """
        if syms:
            syms = set(map(sympify, syms))
        else:
            syms = self.free_symbols

        if syms.intersection(self.free_symbols) == set([]):
            # constant rational function
            return True
        else:
            return self._eval_is_rational_function(syms)

    ###################################################################################
    ##################### SERIES, LEADING TERM, LIMIT, ORDER METHODS ##################
    ###################################################################################

    def series(self, x=None, x0=0, n=6, dir="+"):
        """
        Series expansion of "self" around ``x = x0`` yielding either terms of
        the series one by one (the lazy series given when n=None), else
        all the terms at once when n != None.

        Note: when n != None, if an O() term is returned then the x in the
        in it and the entire expression represents x - x0, the displacement
        from x0. (If there is no O() term then the series was exact and x has
        it's normal meaning.) This is currently necessary since sympy's O()
        can only represent terms at x0=0. So instead of::

          cos(x).series(x0=1, n=2) --> (1 - x)*sin(1) + cos(1) + O((x - 1)**2)

        which graphically looks like this::

               |
              .|.         . .
             . | \      .     .
            ---+----------------------
               |   . .          . .
               |    \
              x=0

        the following is returned instead::

        -x*sin(1) + cos(1) + O(x**2)

        whose graph is this::

               \ |
              . .|        . .
             .   \      .     .
            -----+\------------------.
                 | . .          . .
                 |  \
                x=0

        which is identical to ``cos(x + 1).series(n=2)``.

        Usage:
            Returns the series expansion of "self" around the point ``x = x0``
            with respect to ``x`` up to O(x**n) (default n is 6).

            If ``x=None`` and ``self`` is univariate, the univariate symbol will
            be supplied, otherwise an error will be raised.

            >>> from sympy import cos, exp
            >>> from sympy.abc import x, y
            >>> cos(x).series()
            1 - x**2/2 + x**4/24 + O(x**6)
            >>> cos(x).series(n=4)
            1 - x**2/2 + O(x**4)
            >>> e = cos(x + exp(y))
            >>> e.series(y, n=2)
            cos(x + 1) - y*sin(x + 1) + O(y**2)
            >>> e.series(x, n=2)
            cos(exp(y)) - x*sin(exp(y)) + O(x**2)

            If ``n=None`` then an iterator of the series terms will be returned.

            >>> term=cos(x).series(n=None)
            >>> [term.next() for i in range(2)]
            [1, -x**2/2]

            For ``dir=+`` (default) the series is calculated from the right and
            for ``dir=-`` the series from the left. For smooth functions this
            flag will not alter the results.

            >>> abs(x).series(dir="+")
            x
            >>> abs(x).series(dir="-")
            -x

        """
        if x is None:
            syms = self.atoms(C.Symbol)
            if len(syms) > 1:
                raise ValueError('x must be given for multivariate functions.')
            x = syms.pop()

        if not self.has(x):
            if n is None:
                return (s for s in [self])
            else:
                return self

        ## it seems like the following should be doable, but several failures
        ## then occur. Is this related to issue 1747 et al See also XPOS below.
        #if x.is_positive is x.is_negative is None:
        #    # replace x with an x that has a positive assumption
        #    xpos = C.Dummy('x', positive=True)
        #    rv = self.subs(x, xpos).series(xpos, x0, n, dir)
        #    if n is None:
        #        return (s.subs(xpos, x) for s in rv)
        #    else:
        #        return rv.subs(xpos, x)

        if len(dir) != 1 or dir not in '+-':
            raise ValueError("Dir must be '+' or '-'")

        if x0 in [S.Infinity, S.NegativeInfinity]:
            dir = {S.Infinity: '+', S.NegativeInfinity: '-'}[x0]
            s = self.subs(x, 1/x).series(x, n=n, dir=dir)
            if n is None:
                return (si.subs(x, 1/x) for si in s)
            # don't include the order term since it will eat the larger terms
            return s.removeO().subs(x, 1/x)

        # use rep to shift origin to x0 and change sign (if dir is negative)
        # and undo the process with rep2
        if x0 or dir == '-':
            if dir == '-':
                rep = -x + x0
                rep2 = -x
                rep2b = x0
            else:
                rep = x + x0
                rep2 = x
                rep2b = -x0
            s = self.subs(x, rep).series(x, x0=0, n=n, dir='+')
            if n is None: # lseries...
                return (si.subs(x, rep2 + rep2b) for si in s)
            # nseries...
            o = s.getO() or S.Zero
            s = s.removeO()
            if o and x0:
                rep2b = 0 # when O() can handle x0 != 0 this can be removed
            return s.subs(x, rep2 + rep2b) + o

        # from here on it's x0=0 and dir='+' handling

        if n != None: # nseries handling
            s1 = self._eval_nseries(x, n=n, logx=None)
            o = s1.getO() or S.Zero
            if o:
                # make sure the requested order is returned
                ngot = o.getn()
                if ngot > n:
                    # leave o in its current form (e.g. with x*log(x)) so
                    # it eats terms properly, then replace it below
                    s1 += o.subs(x, x**C.Rational(n, ngot))
                elif ngot < n:
                    # increase the requested number of terms to get the desired
                    # number keep increasing (up to 9) until the received order
                    # is different than the original order and then predict how
                    # many additional terms are needed
                    for more in range(1, 9):
                        s1 = self._eval_nseries(x, n=n + more, logx=None)
                        newn = s1.getn()
                        if newn != ngot:
                            ndo = n + (n - ngot)*more/(newn - ngot)
                            s1 = self._eval_nseries(x, n=ndo, logx=None)
                            # if this assertion fails then our ndo calculation
                            # needs modification
                            assert s1.getn() == n
                            break
                    else:
                        raise ValueError('Could not calculate %s terms for %s'
                                         % (str(n), self))
                o = s1.getO()
                s1 = s1.removeO()
            else:
                o = C.Order(x**n)
                if (s1 + o).removeO() == s1:
                    o = S.Zero

            return s1 + o

        else: # lseries handling
            def yield_lseries(s):
                """Return terms of lseries one at a time."""
                for si in s:
                    if not si.is_Add:
                        yield si
                        continue
                    # yield terms 1 at a time if possible
                    # by increasing order until all the
                    # terms have been returned
                    yielded = 0
                    o = C.Order(si)*x
                    ndid = 0
                    ndo = len(si.args)
                    while 1:
                        do = (si - yielded + o).removeO()
                        o *= x
                        if not do or do.is_Order:
                            continue
                        if do.is_Add:
                            ndid += len(do.args)
                        else:
                            ndid += 1
                        yield do
                        if ndid == ndo:
                            raise StopIteration
                        yielded += do

            return yield_lseries(self.removeO()._eval_lseries(x))

    def lseries(self, x=None, x0=0, dir='+'):
        """
        Wrapper for series yielding an iterator of the terms of the series.

        Note: an infinite series will yield an infinite iterator. The following,
        for exaxmple, will never terminate. It will just keep printing terms
        of the sin(x) series::

          for term in sin(x).lseries(x):
              print term

        The advantage of lseries() over nseries() is that many times you are
        just interested in the next term in the series (i.e. the first term for
        example), but you don't know how many you should ask for in nseries()
        using the "n" parameter.

        See also nseries().
        """
        return self.series(x, x0, n=None, dir=dir)

    def _eval_lseries(self, x):
        # default implementation of lseries is using nseries(), and adaptively
        # increasing the "n". As you can see, it is not very efficient, because
        # we are calculating the series over and over again. Subclasses should
        # override this method and implement much more efficient yielding of
        # terms.
        n = 0
        series = self._eval_nseries(x, n=n, logx=None)
        if not series.is_Order:
            if series.is_Add:
                yield series.removeO()
            else:
                yield series
            raise StopIteration

        while series.is_Order:
            n += 1
            series = self._eval_nseries(x, n=n, logx=None)
        e = series.removeO()
        yield e
        while 1:
            while 1:
                n += 1
                series = self._eval_nseries(x, n=n, logx=None).removeO()
                if e != series:
                    break
            yield series - e
            e = series

    def nseries(self, x=None, x0=0, n=6, dir='+',logx=None):
        """
        Wrapper to _eval_nseries if assumptions allow, else to series.

        If x is given, x0 is 0, dir='+', and self has x, then _eval_nseries is
        called. This calculates "n" terms in the innermost expressions and
        then builds up the final series just by "cross-multiplying" everything
        out.

        Advantage -- it's fast, because we don't have to determine how many
        terms we need to calculate in advance.

        Disadvantage -- you may end up with less terms than you may have
        expected, but the O(x**n) term appended will always be correct and
        so the result, though perhaps shorter, will also be correct.

        If any of those assumptions is not met, this is treated like a
        wrapper to series which will try harder to return the correct
        number of terms.

        See also lseries().
        """
        if x and not self.has(x):
            return self
        if x is None or x0 or dir != '+':#{see XPOS above} or (x.is_positive == x.is_negative == None):
            assert logx == None
            return self.series(x, x0, n, dir)
        else:
            return self._eval_nseries(x, n=n, logx=logx)

    def _eval_nseries(self, x, n, logx):
        """
        Return terms of series for self up to O(x**n) at x=0
        from the positive direction.

        This is a method that should be overridden in subclasses. Users should
        never call this method directly (use .nseries() instead), so you don't
        have to write docstrings for _eval_nseries().
        """
        raise NotImplementedError("(%s).nseries(%s, %s, %s)" % (self, x, x0, n))

    def limit(self, x, xlim, dir='+'):
        """ Compute limit x->xlim.
        """
        from sympy.series.limits import limit
        return limit(self, x, xlim, dir)

    def compute_leading_term(self, x, skip_abs=False, logx=None):
        """ as_leading_term is only allowed for results of .series()
            This is a wrapper to compute a series first.
            If skip_abs is true, the absolute term is assumed to be zero.
            (This is necessary because sometimes it cannot be simplified
             to zero without a lot of work, but is still known to be zero.
             See log._eval_nseries for an example.)
            If skip_log is true, log(x) is treated as an independent symbol.
            (This is needed for the gruntz algorithm.)
        """
        from sympy.series.gruntz import calculate_series
        from sympy import cancel, expand_mul
        if self.removeO() == 0:
            return self
        if logx is None:
            d = C.Dummy('logx')
            s = calculate_series(self, x, skip_abs, d).subs(d, C.log(x))
        else:
            s = calculate_series(self, x, skip_abs, logx)
        s = cancel(s)
        if skip_abs:
            s = expand_mul(s).as_independent(x)[1]
        return s.as_leading_term(x)

    @cacheit
    def as_leading_term(self, *symbols):
        """
        Returns the leading term.

        Example:

        >>> from sympy.abc import x
        >>> (1+x+x**2).as_leading_term(x)
        1
        >>> (1/x**2+x+x**2).as_leading_term(x)
        x**(-2)

        Note:

        self is assumed to be the result returned by Basic.series().
        """
        from sympy import powsimp
        if len(symbols)>1:
            c = self
            for x in symbols:
                c = c.as_leading_term(x)
            return c
        elif not symbols:
            return self
        x = sympify(symbols[0])
        assert x.is_Symbol, repr(x)
        if not self.has(x):
            return self
        obj = self._eval_as_leading_term(x)
        if obj is not None:
            return powsimp(obj, deep=True, combine='exp')
        raise NotImplementedError('as_leading_term(%s, %s)' % (self, x))

    def _eval_as_leading_term(self, x):
        return self

    def as_coeff_exponent(self, x):
        """ c*x**e -> c,e where x can be any symbolic expression.
        """
        x = sympify(x)
        wc = Wild('wc')
        we = Wild('we')
        p  = wc*x**we
        from sympy import collect
        s = collect(self, x)
        d = s.match(p)
        if d is not None and we in d:
            return d[wc], d[we]
        return s, S.Zero

    def leadterm(self, x):
        """
        Returns the leading term a*x**b as a tuple (a, b).

        Example:

        >>> from sympy.abc import x
        >>> (1+x+x**2).leadterm(x)
        (1, 0)
        >>> (1/x**2+x+x**2).leadterm(x)
        (1, -2)

        Note:

        self is assumed to be the result returned by Basic.series().
        """
        from sympy import powsimp
        x = sympify(x)
        c, e = self.as_leading_term(x).as_coeff_exponent(x)
        c = powsimp(c, deep=True, combine='exp')
        if not c.has(x):
            return c, e
        raise ValueError("cannot compute leadterm(%s, %s), got c=%s" % (self, x, c))

    def as_coeff_Mul(self):
        """Efficiently extract the coefficient of a product. """
        return S.One, self

    ###################################################################################
    ##################### DERIVATIVE, INTEGRAL, FUNCTIONAL METHODS ####################
    ###################################################################################

    def diff(self, *symbols, **assumptions):
        new_symbols = map(sympify, symbols) # e.g. x, 2, y, z
        assumptions.setdefault("evaluate", True)
        return Derivative(self, *new_symbols, **assumptions)

    ###########################################################################
    ###################### EXPRESSION EXPANSION METHODS #######################
    ###########################################################################

    # These should be overridden in subclasses

    def _eval_expand_basic(self, deep=True, **hints):
        return self

    def _eval_expand_power_exp(self, deep=True, **hints):
        return self

    def _eval_expand_power_base(self, deep=True, **hints):
        return self

    def _eval_expand_mul(self, deep=True, **hints):
        return self

    def _eval_expand_multinomial(self, deep=True, **hints):
        return self

    def _eval_expand_log(self, deep=True, **hints):
        return self

    def _eval_expand_complex(self, deep=True, **hints):
        return self

    def _eval_expand_trig(self, deep=True, **hints):
        return self

    def _eval_expand_func(self, deep=True, **hints):
        return self

    def expand(self, deep=True, modulus=None, power_base=True, power_exp=True, \
            mul=True, log=True, multinomial=True, basic=True, **hints):
        """
        Expand an expression using hints.

        See the docstring in function.expand for more information.
        """
        hints.update(power_base=power_base, power_exp=power_exp, mul=mul, \
           log=log, multinomial=multinomial, basic=basic)

        expr = self
        for hint, use_hint in hints.iteritems():
            if use_hint:
                func = getattr(expr, '_eval_expand_'+hint, None)
                if func is not None:
                    expr = func(deep=deep, **hints)

        if modulus is not None:
            modulus = sympify(modulus)

            if not modulus.is_Integer or modulus <= 0:
                raise ValueError("modulus must be a positive integer, got %s" % modulus)

            terms = []

            for term in Add.make_args(expr):
                coeff, tail = term.as_coeff_mul()

                coeff %= modulus

                if coeff:
                    terms.append(Mul(*((coeff,) + tail)))

            expr = Add(*terms)

        return expr

    ###########################################################################
    ################### GLOBAL ACTION VERB WRAPPER METHODS ####################
    ###########################################################################

    def integrate(self, *args, **kwargs):
        """See the integrate function in sympy.integrals"""
        from sympy.integrals import integrate
        return integrate(self, *args, **kwargs)

    def simplify(self):
        """See the simplify function in sympy.simplify"""
        from sympy.simplify import simplify
        return simplify(self)

    def nsimplify(self, constants=[], tolerance=None, full=False):
        """See the nsimplify function in sympy.simplify"""
        from sympy.simplify import nsimplify
        return nsimplify(self, constants, tolerance, full)

    def separate(self, deep=False, force=False):
        """See the separate function in sympy.simplify"""
        from sympy.simplify import separate
        return separate(self, deep)

    def collect(self, syms, evaluate=True, exact=False):
        """See the collect function in sympy.simplify"""
        from sympy.simplify import collect
        return collect(self, syms, evaluate, exact)

    def together(self, *args, **kwargs):
        """See the together function in sympy.polys"""
        from sympy.polys import together
        return together(self, *args, **kwargs)

    def apart(self, x=None, **args):
        """See the apart function in sympy.polys"""
        from sympy.polys import apart
        return apart(self, x, **args)

    def ratsimp(self):
        """See the ratsimp function in sympy.simplify"""
        from sympy.simplify import ratsimp
        return ratsimp(self)

    def trigsimp(self, deep=False, recursive=False):
        """See the trigsimp function in sympy.simplify"""
        from sympy.simplify import trigsimp
        return trigsimp(self, deep, recursive)

    def radsimp(self):
        """See the radsimp function in sympy.simplify"""
        from sympy.simplify import radsimp
        return radsimp(self)

    def powsimp(self, deep=False, combine='all'):
        """See the powsimp function in sympy.simplify"""
        from sympy.simplify import powsimp
        return powsimp(self, deep, combine)

    def combsimp(self):
        """See the combsimp function in sympy.simplify"""
        from sympy.simplify import combsimp
        return combsimp(self)

    def factor(self, *gens, **args):
        """See the factor() function in sympy.polys.polytools"""
        from sympy.polys import factor
        return factor(self, *gens, **args)

    def refine(self, assumption=True):
        """See the refine function in sympy.assumptions"""
        from sympy.assumptions import refine
        return refine(self, assumption)

    def cancel(self, *gens, **args):
        """See the cancel function in sympy.polys"""
        from sympy.polys import cancel
        return cancel(self, *gens, **args)

    def invert(self, g):
        """See the invert function in sympy.polys"""
        from sympy.polys import invert
        return invert(self, g)


class AtomicExpr(Atom, Expr):
    """
    A parent class for object which are both atoms and Exprs.

    Examples: Symbol, Number, Rational, Integer, ...
    But not: Add, Mul, Pow, ...
    """

    is_Atom = True

    __slots__ = []

    def _eval_derivative(self, s):
        if self == s:
            return S.One
        return S.Zero

    def as_numer_denom(self):
        return self, S.One

    def _eval_is_polynomial(self, syms):
        return True

    def _eval_is_rational_function(self, syms):
        return True

    def _eval_nseries(self, x, n, logx):
        return self

from mul import Mul
from add import Add
from power import Pow
from function import Derivative
from sympify import sympify
from symbol import Wild
