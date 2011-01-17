from core import C
from basic import Basic, Atom
from singleton import S
from evalf import EvalfMixin
from decorators import _sympifyit, call_highest_priority
from cache import cacheit

class Expr(Basic, EvalfMixin):
    __slots__ = []

    # ***************
    # * Arithmetics *
    # ***************

    # Expr and its sublcasses use _op_priority to determine which object
    # passed to a binary special method (__mul__, etc.) will handle the
    # operation. In general, the 'call_highest_priority' decorator will choose
    # the object with the highest _op_priority to handle the call.
    # Custom subclasses that want to define their own binary special methods
    # should set an _op_priority value that is higher than the default.
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

    def _eval_conjugate(self):
        if self.is_real:
            return self

    def conjugate(self):
        from sympy.functions.elementary.complexes import conjugate as c
        return c(self)

    def removeO(self):
        """Removes the additive O(..) symbol if there is one"""
        if self.is_Order:
            return S.Zero
        if not self.is_Add:
            return self
        args = []
        for a in self.args:
            if a.is_Order:
                continue
            args.append(a)
        return self._new_rawargs(*args)

    def getO(self):
        """Returns the additive O(..) symbol if there is one, else None."""
        if self.is_Order:
            return self
        if not self.is_Add:
            return None
        args = []
        for a in self.args:
            if not a.is_Order:
                continue
            args.append(a)
        if args:
            return self._new_rawargs(*args)

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
        >>>

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
        const = x.as_coeff_mul()[0] # constant multiplying x
        if const != S.One: # get rid of constants
            result = self.coeff(x/const)
            if result is not None:
                return(result/const)
            else:
                return None
        if x.is_Integer:
            return

        result = self
        if expand:
            result = result.expand() # collect expects its arguments in expanded form
        result = collect(result, x, evaluate=False, exact=True)
        if x in result:
            return result[x]
        else:
            return None

    def as_expr(self, *gens):
        """
        Convert a polynomial to a SymPy expression.

        **Examples**

        >>> from sympy import sin
        >>> from sympy.abc import x, y

        >>> f = (x**2 + x*y).as_poly(x, y)
        >>> f.as_expr()
        x*y + x**2

        >>> sin(x).as_expr()
        sin(x)

        """
        return self

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

                if coeff[w].has(*Mul.make_args(expr)):
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

        See also:
         - .as_two_terms() to split expr into a head and tail
         - .as_coeff_add(*deps) to split expr like as_independent(),
                                but return the dependent pieces as
                                a tuple of arguments when treating
                                expr as an Add
         - .as_coeff_add(*deps) to split expr like as_independent(),
                                but return the dependent pieces as
                                a tuple of arguments when treating
                                expr as a Mul
        """
        indeps, depend = [], []

        if self.is_Add or self.is_Mul:
            for term in self.args:
                if term.has(*deps):
                    depend.append(term)
                else:
                    indeps.append(term)

            return (self._new_rawargs(*indeps),
                    self._new_rawargs(*depend))
        else:
            if self.has(*deps):
                return (S.One, self)
            else:
                return (self, S.One)

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
        """Return the tuple (c, args) where self is written as a Mul, `m`.

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
          use self.as_independent(*deps)

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
        """Return the tuple (c, args) where self is written as an Add, `a`.

        c should be a Rational added to any terms of the Add that are
        independent of deps.

        args should be a tuple of all other terms of `a`; args is empty
        if self is a Number or if self is independent of deps (when given).

        This should be used when you don't know if self is an Add or not but
        you want to treat self as an Add or if you want to process the
        individual arguments of the tail of self as an Add.

        - if you know self is an Add and want only the head, use self.args[0];
        - if you don't want to process the arguments of the tail but need the
          tail then use self.as_two_terms() which gives the head and tail.
        - if you want to split self into an independent and dependent parts
          use self.as_independent(*deps)

        >>> from sympy import S
        >>> from sympy.abc import x, y
        >>> (S(3)).as_coeff_add()
        (3, ())
        >>> (3 + x + y).as_coeff_add()
        (3, (y, x))
        >>> (3 + x +y).as_coeff_add(x)
        (3 + y, (x,))
        >>> (3 + y).as_coeff_add(x)
        (3 + y, ())
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

           For any expression, the set {e.could_extract_minus_sign(),
           (-e).could_extract_minus_sign()} must be {True, False}.

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


    ###################################################################################
    ##################### SERIES, LEADING TERM, LIMIT, ORDER METHODS ##################
    ###################################################################################

    def series(self, x=None, x0=0, n=6, dir="+"):
        """
        Series expansion of "self" around `x = x0` yielding either terms of
        the series one by one (the lazy series given when n=None), else
        all the terms at once when n != None.

        Note: when n != None, if an O() term is returned then the x in the
        in it and the entire expression represents x - x0, the displacement
        from x0. (If there is no O() term then the series was exact and x has
        it's normal meaning.) This is currently necessary since sympy's O()
        can only represent terms at x0=0. So instead of

            >> cos(x).series(x0=1, n=2)
            (1 - x)*sin(1) + cos(1) + O((x - 1)**2)

        which graphically looks like this:

               \
              .|.         . .
             . | \      .     .
            ---+----------------------
               |   . .          . .
               |    \
              x=0

        the following is returned instead

            -x*sin(1) + cos(1) + O(x**2)

        whose graph is this

               \ |
              . .|        . .
             .   \      .     .
            -----+\------------------.
                 | . .          . .
                 |  \
                x=0

        which is identical to cos(x + 1).series(n=2).

        Usage:
            Returns the series expansion of "self" around the point `x = x0`
            with respect to `x` up to O(x**n) (default n is 6).

            If `x=None` and `self` is univariate, the univariate symbol will
            be supplied, otherwise an error will be raised.

            >>> from sympy import cos, exp
            >>> from sympy.abc import x, y
            >>> cos(x).series()
            1 - x**2/2 + x**4/24 + O(x**6)
            >>> cos(x).series(n=4)
            1 - x**2/2 + O(x**4)
            >>> e = cos(x + exp(y))
            >>> e.series(y, n=2)
            -y*sin(1 + x) + cos(1 + x) + O(y**2)
            >>> e.series(x, n=2)
            -x*sin(exp(y)) + cos(exp(y)) + O(x**2)

            If `n=None` then an iterator of the series terms will be returned.

            >>> term=cos(x).series(n=None)
            >>> [term.next() for i in range(2)]
            [1, -x**2/2]

            For `dir=+` (default) the series is calculated from the right and
            for `dir=-` the series from the left. For smooth functions this
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
        ## then occur. Is this related to issue 1747 et al? See also XPOS below.
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
            s1 = self._eval_nseries(x, n=n)
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
                        s1 = self._eval_nseries(x, n=n + more)
                        newn = s1.getn()
                        if newn != ngot:
                            ndo = n + (n - ngot)*more/(newn - ngot)
                            s1 = self._eval_nseries(x, n=ndo)
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
        of the sin(x) series:

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
        series = self._eval_nseries(x, n=n)
        if not series.is_Order:
            if series.is_Add:
                yield series.removeO()
            else:
                yield series
            raise StopIteration

        while series.is_Order:
            n += 1
            series = self._eval_nseries(x, n=n)
        e = series.removeO()
        yield e
        while 1:
            while 1:
                n += 1
                series = self._eval_nseries(x, n=n).removeO()
                if e != series:
                    break
            yield series - e
            e = series

    def nseries(self, x=None, x0=0, n=6, dir='+'):
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
            return self.series(x, x0, n, dir)
        else:
            return self._eval_nseries(x, n=n)

    def _eval_nseries(self, x, n):
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

    def separate(self, deep=False):
        """See the seperate function in sympy.simplify"""
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

    @property
    def is_number(self):
        return True

    def _eval_nseries(self, x, n):
        return self

from mul import Mul
from add import Add
from power import Pow
from relational import Inequality, StrictInequality
from function import Derivative
from sympify import _sympify, sympify, SympifyError
from symbol import Wild
