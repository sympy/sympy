
from basic import Basic, S, C
from evalf import EvalfMixin
from decorators import _sympifyit
from cache import cacheit


class Expr(Basic, EvalfMixin):
    __slots__ = []

    # ***************
    # * Arithmetics *
    # ***************

    def __pos__(self):
        return self
    def __neg__(self):
        return Mul(S.NegativeOne, self)
    def __abs__(self):
        return C.abs(self)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        return Add(self, other)
    @_sympifyit('other', NotImplemented)
    def __radd__(self, other):
        return Add(other, self)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        return Add(self, -other)
    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        return Add(other, -self)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        return Mul(self, other)
    @_sympifyit('other', NotImplemented)
    def __rmul__(self, other):
        return Mul(other, self)

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        return Pow(self, other)
    @_sympifyit('other', NotImplemented)
    def __rpow__(self, other):
        return Pow(other, self)

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        return Mul(self, Pow(other, S.NegativeOne))
    @_sympifyit('other', NotImplemented)
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
            return C.Real._new(x._mpf_, prec)
        elif hasattr(x, "_mpc_"):
            re, im = x._mpc_
            re = C.Real._new(re, prec)
            im = C.Real._new(im, prec)*S.ImaginaryUnit
            return re+im
        else:
            raise TypeError("expected mpmath number (mpf or mpc)")


    def _eval_interval(self, x, a, b):
        """
        Returns evaluation over an interval.  For most functions this is:

        self.subs(x, b) - self.subs(x, a),

        possibly using limit() if NaN is returned from subs.

        If b or a is None, it only evaluates -self.subs(x, a) or self.subs(b, x),
        respectively.

        """
        from sympy.series import limit
        if a is None:
            A = 0
        else:
            A = self.subs(x, a)

        if A is S.NaN:
            A = limit(self, x, a)
            if A is S.NaN:
                return self
        if b is None:
            B = 0
        else:
            B = self.subs(x, b)

        if B is S.NaN:
            B = limit(self, x, b)
        if B is S.NaN:
            return self
        return B - A

    def _eval_power(self, other):
        return None

    def _eval_derivative(self, s):
        return

    def _eval_conjugate(self):
        if self.is_real:
            return self

    def conjugate(self):
        from sympy.functions.elementary.complexes import conjugate as c
        return c(self)

    def removeO(self):
        "Removes the O(..) symbol if there is one"
        if self.is_Order:
            return Integer(0)
        for i,x in enumerate(self.args):
            if x.is_Order:
                return Add(*(self.args[:i]+self.args[i+1:]))
        return self

    def getO(e):
        "Returns the O(..) symbol, or None if there is none."
        if e.is_Order:
            return e
        if e.is_Add:
            for x in e.args:
                if x.is_Order:
                    return x

    def coeff(self, x, expand=True):
        """
        Returns the coefficient of the term "x" or None if there is no "x".

        Optional expand keyword argument allows one to control whether the
        expression is expanded before terms are collected, which can be useful
        if the term "x" isn't nested inside of terms and you don't want the
        returned coefficient to be expanded.

        Example:

        >>> from sympy import symbols
        >>> from sympy.abc import x, y, z
        >>> (3+2*x+4*x**2).coeff(1)
        >>> (3+2*x+4*x**2).coeff(x)
        2
        >>> (3+2*x+4*x**2).coeff(x**2)
        4
        >>> (3+2*x+4*x**2).coeff(x**3)
        >>> (z*(x+y)**2).coeff(z)
        2*x*y + x**2 + y**2
        >>> (z*(x+y)**2).coeff(z, expand=False)
        (x + y)**2
        >>>

        """
        from sympy import collect
        x = sympify(x)
        const = x.as_coeff_terms()[0] # constant multiplying x
        if const != S.One: # get rid of constants
            result = self.coeff(x/const)
            if result is not None:
                return(result/const)
            else:
                return None
        if x.is_Integer:
            return

        if expand:
            self = self.expand() # collect expects its arguments in expanded form
        result = collect(self, x, evaluate=False, exact=True)
        if x in result:
            return result[x]
        else:
            return None

    def as_coefficient(self, expr):
        """Extracts symbolic coefficient at the given expression. In
           other words, this functions separates 'self' into product
           of 'expr' and 'expr'-free coefficient. If such separation
           is not possible it will return None.

           >>> from sympy import E, pi, sin, I
           >>> from sympy.abc import x, y

           >>> E.as_coefficient(E)
           1
           >>> (2*E).as_coefficient(E)
           2

           >>> (2*E + x).as_coefficient(E)
           >>> (2*sin(E)*E).as_coefficient(E)

           >>> (2*pi*I).as_coefficient(pi*I)
           2

           >>> (2*I).as_coefficient(pi*I)

        """
        if expr.is_Add:
            return None
        else:
            w = Wild('w')

            coeff = self.match(w * expr)

            if coeff is not None:
                if expr.is_Mul:
                    args = expr.args
                else:
                    args = [expr]

                if coeff[w].has(*args):
                    return None
                else:
                    return coeff[w]
            else:
                return None

    def as_independent(self, *deps):
        """Returns a pair with separated parts of a given expression
           independent of specified symbols in the first place and
           dependent on them in the other. Both parts are valid
           SymPy expressions.

           >>> from sympy import sin, cos
           >>> from sympy.abc import x, y

           >>> (2*x*sin(x)+y+x).as_independent(x)
           (y, x + 2*x*sin(x))

           >>> (x*sin(x)*cos(y)).as_independent(x)
           (cos(y), x*sin(x))

           All other expressions are multiplicative:

           >>> (sin(x)).as_independent(x)
           (1, sin(x))

           >>> (sin(x)).as_independent(y)
           (sin(x), 1)

        """
        indeps, depend = [], []

        if self.is_Add or self.is_Mul:
            for term in self.args:
                if term.has(*deps):
                    depend.append(term)
                else:
                    indeps.append(term)

            return (self.__class__(*indeps),
                    self.__class__(*depend))
        else:
            if self.has(*deps):
                return (S.One, self)
            else:
                return (self, S.One)

    def as_real_imag(self):
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
        expr = self.expand(complex=True, deep=False)


        re_part, im_part = [], []
        if not expr.is_Add:
            args = [expr]
        else:
            args = expr.args

        for term in args:
            coeff = term.as_coefficient(S.ImaginaryUnit)

            if coeff is None:
                re_part.append(term)
            else:
                im_part.append(coeff)

        return (Add(*re_part), Add(*im_part))

    def as_powers_dict(self):
        return { self : S.One }

    def as_base_exp(self):
        # a -> b ** e
        return self, S.One

    def as_coeff_terms(self, x=None):
        # a -> c * t
        if x is not None:
            if not self.has(x):
                return self, tuple()
        return S.One, (self,)

    def as_coeff_factors(self, x=None):
        # a -> c + f
        if x is not None:
            if not self.has(x):
                return self, tuple()
        return S.Zero, (self,)

    def as_numer_denom(self):
        """ a/b -> a,b

        The following is a possible way to modify Eq which are now
        just returned as (Eq(), 1). It is not a trivial change,
        however, and it causes many failures.

        from sympy.core.relational import Equality
        from sympy import Eq
        if isinstance(self, Equality):
            l = Symbol('l', dummy=True)
            r = Symbol('r', dummy=True)
            n, d = (l*self.lhs - r*self.rhs).as_numer_denom()
            return Eq(n.subs({l: 1, r: 0}),
                      n.subs({l: 0, r: -1})), d.subs({l: 1, r: 1})
        """

        base, exp = self.as_base_exp()
        coeff, terms = exp.as_coeff_terms()
        if coeff.is_negative:
            # b**-e -> 1, b**e
            return S.One, base ** (-exp)
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
            elif self.is_Real:
                if not quotient.is_Real:
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
            return C.Add(*newargs)
        elif self.is_Mul:
            for i in xrange(len(self.args)):
                newargs = list(self.args)
                del(newargs[i])
                tmp = C.Mul(*newargs).extract_multiplicatively(c)
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
           1 + 2*x

           >>> (-x+1).extract_additively(2*x)
           1 - 3*x

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
            elif self.is_Real:
                if not sub.is_Real:
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
            self_coeff, self_terms = self.as_coeff_terms()
            if c.is_Mul:
                c_coeff, c_terms = c.as_coeff_terms()
                if c_terms == self_terms:
                    new_coeff = self_coeff.extract_additively(c_coeff)
                    if new_coeff != None:
                        return new_coeff * C.Mul(*self_terms)
            elif c == self_terms:
                new_coeff = self_coeff.extract_additively(1)
                if new_coeff != None:
                    return new_coeff * C.Mul(*self_terms)

    def could_extract_minus_sign(self):
        """Canonical way to choose an element in the set {e, -e} where
           e is any expression. If the canonical element is e, we have
           e.could_extract_minus_sign() == True, else
           e.could_extract_minus_sign() == False.

           For any expression, the set {e.could_extract_minus_sign(),
           (-e).could_extract_minus_sign()} must be {True, False}.

           >>> from sympy.abc import x, y
           >>> (x-y).could_extract_minus_sign() != (y-x).could_extract_minus_sign()
           True

        """
        from sympy.utilities.iterables import make_list
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
                args = (make_list(num, Mul)) + (make_list(den, Mul))
                arg_signs = [arg.could_extract_minus_sign() for arg in args]
                negative_args = filter(None, arg_signs)
                return len(negative_args) % 2 == 1

            # As a last resort, we choose the one with greater hash
            return hash(self) < hash(negative_self)


    ###################################################################################
    ##################### SERIES, LEADING TERM, LIMIT, ORDER METHODS ##################
    ###################################################################################

    def series(self, x, point=0, n=6, dir="+"):
        """
        Series expansion of "self" around "point".

        Usage:
            Returns the Taylor (Laurent or generalized) series of "self" around
            the point "point" (default 0) with respect to "x" until the n-th
            term (default n is 6).

            For dir="+" (default) it calculates the series from the right
            and for dir="-" the series from the left.
            For smooth functions this argument doesn't matter.

        Notes:
            This method is the most high level method and it returns the
            series including the O(x**n) term.

            Internally, it executes a method nseries(), see nseries() docstring
            for more information.
        """
        x = sympify(x)
        point = sympify(point)
        if dir == "+":
            return self.nseries(x, point, n)
        elif dir == "-":
            return self.subs(x, -x).nseries(x, -point, n).subs(x, -x)
        else:
            raise ValueError("Dir has to be '+' or '-'")


    def lseries(self, x, x0):
        """
        lseries is a generator yielding terms in the series.

        Example: if you do:

        for term in sin(x).lseries(x, 0):
            print term

        It will print all terms of the sin(x) series (i.e. it never
        terminates).

        The advantage of lseries() over nseries() is that many times you are
        just interested in the next term in the series (i.e. the first term for
        example), but you don't know how many you should ask for in nseries()
        using the "n" parameter.

        See also nseries().
        """
        return self._eval_lseries(x, x0)

    def _eval_lseries(self, x, x0):
        # default implementation of lseries is using nseries(), and adaptively
        # increasing the "n". As you can see, it is not very efficient, because
        # we are calculating the series over and over again. Subclasses should
        # override this method and implement much more efficient yielding of
        # terms.
        n = 0
        e = self.nseries(x, x0, n)
        while e.is_Order:
            n += 1
            e = self.nseries(x, x0, n)
        series = e.removeO()
        yield series
        while 1:
            n += 1
            e = self.nseries(x, x0, n).removeO()
            while series == e:
                n += 1
                e = self.nseries(x, x0, n).removeO()
            term = e - series
            series = e
            yield term

    def nseries(self, x, x0, n):
        """
        Calculates a generalized series expansion.

        nseries calculates "n" terms in the innermost expressions and then
        builds up the final series just by "cross-multiplying" everything out.

        Advantage -- it's fast, because we don't have to determine how many
        terms we need to calculate in advance.

        Disadvantage -- you may end up with less terms than you may have
        expected, but the O(x**n) term appended will always be correct and
        so the result, though perhaps shorter, will also be correct.

        See also lseries().
        """
        return self._eval_nseries(x, x0, n)

    def _eval_nseries(self, x, x0, n):
        """
        This is a method that should be overridden in subclasses. Users should
        never call this method directly (use .nseries() instead), so you don't
        have to write docstrings for _eval_nseries().
        """
        raise NotImplementedError("(%s).nseries(%s, %s, %s)" % (self, x, x0, n))

    def limit(self, x, xlim, direction='+'):
        """ Compute limit x->xlim.
        """
        from sympy.series.limits import limit
        return limit(self, x, xlim, direction)

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
        assert x.is_Symbol, `x`
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
        self = collect(self, x)
        d = self.match(p)
        if d is not None and we in d:
            return d[wc], d[we]
        return self, S.Zero

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
        c,e = self.as_leading_term(x).as_coeff_exponent(x)
        c = powsimp(c, deep=True, combine='exp')
        if not c.has(x):
            return c,e
        raise ValueError("cannot compute leadterm(%s, %s), got c=%s" % (self, x, c))

    def as_Add(self):
        """Returns `self` as it was `Add` instance. """
        return [self]

    def as_Mul(self):
        """Returns `self` as it was `Mul` instance. """
        return [self]

    def as_Pow(self):
        """Returns `self` as it was `Pow` instance. """
        return (self, S.One)

    ###################################################################################
    ##################### DERIVATIVE, INTEGRAL, FUNCTIONAL METHODS ####################
    ###################################################################################

    def diff(self, *symbols, **assumptions):
        new_symbols = map(sympify, symbols)
        if not "evaluate" in assumptions:
            assumptions["evaluate"] = True
        ret = Derivative(self, *new_symbols, **assumptions)
        return ret

    def fdiff(self, *indices):
        # FIXME FApply -> ?
        return C.FApply(C.FDerivative(*indices), self)

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

            for term in expr.as_Add():
                coeff, tail = term.as_coeff_terms()

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

    def together(self, *args, **kwargs):
        """See the together function in sympy.simplify"""
        from sympy.simplify import together
        return together(self, *args, **kwargs)

    def nsimplify(self, constants=[], tolerance=None, full=False):
        """See the nsimplify function in sympy.simplify"""
        from sympy.simplify import nsimplify
        return nsimplify(self, constants, tolerance, full)

    def separate(self, deep=False):
        """See the seperate function in sympy.simplify"""
        from sympy.simplify import separate
        return separate(self, deep)

    def collect(self, syms, evaluate=True, exact=False):
        """See the collect function in sympy.simplify"""
        from sympy.simplify import collect
        return collect(self, syms, evaluate, exact)

    def apart(self, x=None, **args):
        """See the apart function in sympy.simplify"""
        from sympy.polys import apart
        return apart(self, x=None, **args)

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
        """See the factor function in sympy.simplify"""
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

from mul import Mul
from power import Pow
from add import Add
from relational import Inequality, StrictInequality
from function import FunctionClass, Derivative
from numbers import Rational, Integer
from sympify import _sympify, sympify, SympifyError
from symbol import Wild

