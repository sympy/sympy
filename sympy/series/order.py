from sympy.core.basic import Basic, S, C, sympify
from sympy.core import oo, Rational, Pow
from sympy.core.cache import cacheit

class Order(Basic):
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
    >>> from sympy import *
    >>> x = Symbol("x")
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

      g(x) = O(f(x)) as x->0  <->  |g(x)|<=M|f(x)| near x=0  <->  lim_{x->0}  |g(x)/f(x)|  < oo

      g(x,y) = O(f(x,y))  <->  lim_{x,y->0}  |g(x,y)/f(x,y)|  < oo, we'll assume that limits commute.

    Notes:
    ======

      In O(f(x),x) the expression f(x) is assumed to have a leading term.
      O(f(x),x) is automatically transformed to O(f(x).as_leading_term(x),x).
      O(expr*f(x),x) is O(f(x),x)
      O(expr,x) is O(1)
      O(0, x) is 0.

      Multivariate O is also supported:
        O(f(x,y),x,y) is transformed to O(f(x,y).as_leading_term(x,y).as_leading_term(y), x, y)

      If O is used with only expression argument then the symbols are
      all symbols in the expression.
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
        else:
            symbols = list(expr.atoms(C.Symbol))

        symbols.sort(Basic.compare)

        if expr.is_Order:

            new_symbols = list(expr.symbols)
            for s in symbols:
                if s not in new_symbols:
                    new_symbols.append(s)
            if len(new_symbols)==len(expr.symbols):
                return expr
            symbols = new_symbols

        elif symbols:

            symbol_map = {}
            new_symbols = []
            for s in symbols:
                if isinstance(s, C.Symbol):
                    new_symbols.append(s)
                    continue
                z = C.Symbol('z',dummy=True)
                x1,s1 = s.solve4linearsymbol(z)
                expr = expr.subs(x1,s1)
                symbol_map[z] = s
                new_symbols.append(z)

            if symbol_map:
                r = Order(expr, *new_symbols, **assumptions)
                expr = r.expr.subs(symbol_map)
                symbols = []
                for s in r.symbols:
                    if symbol_map.has_key(s):
                        symbols.append(symbol_map[s])
                    else:
                        symbols.append(s)
            else:
                if expr.is_Add:
                    lst = expr.extract_leading_order(*symbols)
                    expr = C.Add(*[f.expr for (e,f) in lst])
                else:
                    expr = expr.as_leading_term(*symbols)
                    coeff, terms = expr.as_coeff_terms()
                    if coeff is S.Zero:
                        return coeff
                    expr = C.Mul(*[t for t in terms if t.has(*symbols)])

        elif expr is not S.Zero:
            expr = S.One

        if expr is S.Zero:
            return expr

        # create Order instance:
        obj = Basic.__new__(cls, expr, *symbols, **assumptions)

        return obj

    def _hashable_content(self):
        if self.args[0].is_number:
            return (self.args[0],)
        return self.args

    def oseries(self, order):
        return self

    def _eval_nseries(self, x, x0, n):
        return self

    @classmethod
    def find_limit(cls, f, x):
        """Basically identical to:

        return limit(f, x, 0, dir="+")

        but first trying some easy cases (like x**2) using heuristics, to avoid
        infinite recursion. This is only needed in the Order class and series
        expansion (that shouldn't rely on the Gruntz algorithm too much),
        that's why find_limit() is defined here.
        """

        from sympy import limit, Wild, log

        if f.is_Pow:
            if f.args[0] == x:
                if f.args[1].is_Rational:
                    if f.args[1] > 0:
                        return S.Zero
                    else:
                        return oo
                if f.args[1].is_number:
                    if f.args[1].evalf() > 0:
                        return S.Zero
                    else:
                        return oo
        if f == x:
            return S.Zero
        p, q = Wild("p"), Wild("q")
        r = f.match(x**p * log(x)**q)
        if r:
            p, q = r[p], r[q]
            if q.is_number and p.is_number:
                if q > 0:
                    if p > 0:
                        return S.Zero
                    else:
                        return -oo
                elif q < 0:
                    if p >= 0:
                        return S.Zero
                    else:
                        return -oo

        return limit(f, x, 0, dir="+")

    @property
    def expr(self):
        return self._args[0]

    @property
    def symbols(self):
        return self._args[1:]

    def _eval_power(b, e):
        if e.is_Number:
            return Order(b.expr ** e, *b.symbols)
        return

    def as_expr_symbols(self, order_symbols):
        if order_symbols is None:
            order_symbols = self.symbols
        else:
            for s in self.symbols:
                if s not in order_symbols:
                    order_symbols = order_symbols + (s,)
        return self.expr, order_symbols

    @cacheit
    def contains(self, expr):
        """
        Return True if expr belongs to Order(self.expr, *self.symbols).
        Return False if self belongs to expr.
        Return None if the inclusion relation cannot be determined (e.g. when self and
        expr have different symbols).
        """
        from sympy import powsimp
        if expr is S.Zero:
            return True
        if expr is S.NaN:
            return False
        if expr.is_Order:
            if self.symbols and expr.symbols:
                common_symbols = tuple([s for s in self.symbols if s in expr.symbols])
            elif self.symbols:
                common_symbols = self.symbols
            else:
                common_symbols = expr.symbols
            if not common_symbols:
                if not (self.symbols or expr.symbols): # O(1),O(1)
                    return True
                return None
            r = None
            for s in common_symbols:
                l = Order.find_limit(powsimp(self.expr/expr.expr, deep=True,\
                combine='exp'), s) != 0
                if r is None:
                    r = l
                else:
                    if r != l:
                        return
            return r
        obj = Order(expr, *self.symbols)
        return self.contains(obj)

    def _eval_subs(self, old, new):
        if self==old:
            return new
        if isinstance(old, C.Symbol) and old in self.symbols:
            i = list(self.symbols).index(old)
            if isinstance(new, C.Symbol):
                return Order(self.expr._eval_subs(old, new), *(self.symbols[:i]+(new,)+self.symbols[i+1:]))
            return Order(self.expr._eval_subs(old, new), *(self.symbols[:i]+self.symbols[i+1:]))
        return Order(self.expr._eval_subs(old, new), *self.symbols)

    def _sage_(self):
        #XXX: SAGE doesn't have Order yet. Let's return 0 instead.
        return Rational(0)._sage_()

Basic.singleton['O'] = lambda : Order
