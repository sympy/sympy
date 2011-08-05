from sympy.core import Basic, S, C, sympify, Expr, oo, Rational, Symbol, Dummy
from sympy.core import Add, Mul
from sympy.core.cache import cacheit
from sympy.core.compatibility import cmp_to_key

class Order(Expr):
    """
    Represents O(f(x)) at the point x = 0.

    Definition
    ==========

    g(x) = O(f(x)) as x->0  if and only if
    |g(x)|<=M|f(x)| near x=0                     (1)

    for some positive but finite M. An equivalent way of saying (1) is:

    lim_{x->0}  |g(x)/f(x)|  < oo

    Let's illustrate it on the following example:

    sin x = x - x**3/3! + O(x**5)

    where in this case O(x**5) = x**5/5! - x**7/7! + .... and the definition
    of O means:

    |x**5/5! - x**7/7! + ....| <= M|x**5|      near x=0

    or equivalently:

    lim_{x->0} | (x**5/5! - x**7/7! + ....) / x**5| < oo

    which surely is true, because

    lim_{x->0} | (x**5/5! - x**7/7! + ....) / x**5| = 1/5!


    So intuitively O(x**3) means: all terms x**3, x**4 and
    higher. But not x**2, x or 1.

    Examples:
    =========
    >>> from sympy import O
    >>> from sympy.abc import x
    >>> O(x)
    O(x)
    >>> O(x)*x
    O(x**2)
    >>> O(x)-O(x)
    O(x)

       External links
       --------------

         U{Big O notation<http://en.wikipedia.org/wiki/Big_O_notation>}

    Properties:
    ===========

      g(x) = O(f(x)) as x->0  <->  |g(x)| <= M|f(x)| near x=0
                              <->  lim_{x->0}  |g(x)/f(x)| < oo

      g(x,y) = O(f(x,y))  <->  lim_{x,y->0}  |g(x,y)/f(x,y)|  < oo;
                               it is assumed that limits commute.

    Notes:
    ======

    In O(f(x), x) the expression f(x) is assumed to have a leading term.
    O(f(x), x) is automatically transformed to O(f(x).as_leading_term(x),x).

        O(expr*f(x), x) is O(f(x), x)
        O(expr, x) is O(1)
        O(0, x) is 0.

    Multivariate O is also supported:

        O(f(x, y), x, y) is transformed to
        O(f(x, y).as_leading_term(x,y).as_leading_term(y), x, y)

    If no symbols are passed then all symbols in the expression are used:

    """

    is_Order = True

    __slots__ = []

    @cacheit
    def __new__(cls, expr, *symbols, **assumptions):

        expr = sympify(expr).expand()
        if expr is S.NaN:
            return S.NaN

        if symbols:
            symbols = map(sympify, symbols)
            if not all(isinstance(s, Symbol) for s in symbols):
                raise NotImplementedError('Order at points other than 0 not supported.')
        else:
            symbols = list(expr.free_symbols)

        if expr.is_Order:

            new_symbols = list(expr.variables)
            for s in symbols:
                if s not in new_symbols:
                    new_symbols.append(s)
            if len(new_symbols) == len(expr.variables):
                return expr
            symbols = new_symbols

        elif symbols:

            if expr.is_Add:
                lst = expr.extract_leading_order(*symbols)
                expr = Add(*[f.expr for (e,f) in lst])
            elif expr:
                if len(symbols) > 1:
                    # TODO
                    # We cannot use compute_leading_term because that only
                    # works in one symbol.
                    expr = expr.as_leading_term(*symbols)
                else:
                    expr = expr.compute_leading_term(symbols[0])
                coeff, terms = expr.as_coeff_mul()
                expr = Mul(*[t for t in terms if t.has(*symbols)])

        if expr is S.Zero:
            return expr
        elif not expr.has(*symbols):
            expr = S.One

        # create Order instance:
        symbols.sort(key=cmp_to_key(Basic.compare))
        obj = Expr.__new__(cls, expr, *symbols, **assumptions)

        return obj

    def _hashable_content(self):
        return self.args

    def oseries(self, order):
        return self

    def _eval_nseries(self, x, n, logx):
        return self

    @property
    def expr(self):
        return self._args[0]

    @property
    def variables(self):
        return self._args[1:]

    @property
    def free_symbols(self):
        return self.expr.free_symbols

    def _eval_power(b, e):
        if e.is_Number:
            return Order(b.expr ** e, *b.variables)
        return

    def as_expr_variables(self, order_symbols):
        if order_symbols is None:
            order_symbols = self.variables
        else:
            for s in self.variables:
                if s not in order_symbols:
                    order_symbols = order_symbols + (s,)
        return self.expr, order_symbols

    def removeO(self):
        return S.Zero

    def getO(self):
        return self


    @cacheit
    def contains(self, expr):
        """
        Return True if expr belongs to Order(self.expr, *self.variables).
        Return False if self belongs to expr.
        Return None if the inclusion relation cannot be determined
        (e.g. when self and expr have different symbols).
        """
        # NOTE: when multiplying out series a lot of queries like
        #       O(...).contains(a*x**b) with many a and few b are made.
        #       Separating out the independent part allows for better caching.
        c, m = expr.as_coeff_mul(*self.variables)
        if m != ():
            return self._contains(Mul(*m))
        else:
            # Mul(*m) == 1, and O(1) treatment is somewhat peculiar ...
            # some day this else should not be necessary
            return self._contains(expr)

    @cacheit
    def _contains(self, expr):
        from sympy import powsimp, limit
        if expr is S.Zero:
            return True
        if expr is S.NaN:
            return False
        if expr.is_Order:
            if self.variables and expr.variables:
                common_symbols = tuple([s for s in self.variables if s in expr.variables])
            elif self.variables:
                common_symbols = self.variables
            else:
                common_symbols = expr.variables
            if not common_symbols:
                if not (self.variables or expr.variables): # O(1),O(1)
                    return True
                return None
            r = None
            for s in common_symbols:
                l = limit(powsimp(self.expr/expr.expr, deep=True,\
                combine='exp'), s, 0) != 0
                if r is None:
                    r = l
                else:
                    if r != l:
                        return
            return r
        obj = Order(expr, *self.variables)
        return self.contains(obj)

    def _eval_subs(self, old, new):
        if self == old:
            return new
        if isinstance(old, Symbol) and old in self.variables:
            i = list(self.variables).index(old)
            if isinstance(new, Symbol):
                return Order(self.expr._eval_subs(old, new), *(self.variables[:i]+(new,)+self.variables[i+1:]))
            return Order(self.expr._eval_subs(old, new), *(self.variables[:i]+self.variables[i+1:]))
        return Order(self.expr._eval_subs(old, new), *self.variables)

    def _eval_derivative(self, x):
        return self.func(self.expr.diff(x), *self.variables) or self

    def _sage_(self):
        #XXX: SAGE doesn't have Order yet. Let's return 0 instead.
        return Rational(0)._sage_()

O = Order
