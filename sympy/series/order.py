
from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.methods import ArithMeths, RelMeths

class Order(Basic, ArithMeths, RelMeths):
    """
    Represents O(f(x)) at the point x = 0.

    Definition
    ==========

    g(x) = O(f(x)) as x->0  if and only if
    |g(x)|<=M|f(x)| near x=0                     (1)

    An equivalent way of saying (1) is:

    lim_{x->0}  |g(x)/f(x)|  < oo
    
    Let's illustrate it on the following example:

    sin x = x - x**3/3! + O(x**5)

    where in this case O(x**5) = x**5/5! - x**7/7! + .... and the definition
    of O means:

    |x**5/5! - x**7/7! + ....| <= M|x**5|      near x=0

    or equivalently:

    lim_{x->0} |x**5/5! - x**7/7! + .... / x**5| < oo

    which surely is true, because 
    
    lim_{x->0} |x**5/5! - x**7/7! + .... / x**5| = 1/5!


    So intuitively O(x**3) means (in our case): all terms x**3, x**4 and
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

      g(x,y) = O(f(x,y))  <->  lim_{x,y->0}  |g(x,y)/f(x,y)|  < oo, will assume that limits commute.

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

    precedence = Basic.Apply_precedence

    _cache = {}

    @cache_it_immutable
    def __new__(cls, expr, *symbols, **assumptions):
        expr = Basic.sympify(expr).expand(trig=True)
        if isinstance(expr, Basic.NaN):
            return S.NaN
        
        if symbols:
            symbols = map(Basic.sympify, symbols)
        else:
            symbols = list(expr.atoms(Basic.Symbol))

        symbols.sort(Basic.compare)

        if isinstance(expr, Order):

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
                if isinstance(s, Basic.Symbol):
                    new_symbols.append(s)
                    continue
                z = Basic.Symbol('z',dummy=True)
                x1,s1 = s.solve4linearsymbol(z)
                expr = expr.subs(x1,s1)
                symbol_map[z] = s
                new_symbols.append(z)
                
            if symbol_map:
                r = Order(expr, *new_symbols, **assumptions)
                expr = r.expr.subs_dict(symbol_map)
                symbols = []
                for s in r.symbols:
                    if symbol_map.has_key(s):
                        symbols.append(symbol_map[s])
                    else:
                        symbols.append(s)
            else:
                if isinstance(expr, Basic.Add):
                    lst = expr.extract_leading_order(*symbols)
                    expr = Basic.Add(*[f.expr for (e,f) in lst])
                else:
                    expr = expr.as_leading_term(*symbols)
                    coeff, terms = expr.as_coeff_terms()
                    if isinstance(coeff, Basic.Zero):
                        return coeff
                    expr = Basic.Mul(*[t for t in terms if t.has(*symbols)]) 

        elif not isinstance(expr, Basic.Zero):
            expr = Basic.One()

        if isinstance(expr, Basic.Zero):
            return expr

        # remove unused symbols
        #symbols = tuple([s for s in symbols if expr.has(s)])
        symbols = tuple(symbols)

        # look Order symbols from cache, TODO: make cache a dictionary
        cache = Order._cache.get(symbols,[])
        for o in cache:
            if o.expr==expr:
                return o

        # Order symbols are assumed to be close to 0 from right:
        for s in symbols:
            assume_dict = {}
            if not s.is_infinitesimal:
                assume_dict['infinitesimal'] = True
            #XXX This causes problems, that it changes the assumption in the
            #   symbol, outside the scope of Order and breaks code. Don't know
            #   why
            #   But sometimes it's necessary for simplifications...
            #   well, how to solve that? I don't know...
            #   ok - so the problem is in caching - in core/function.py:63
            # see the issue 369
            if s.is_positive is None:
                assume_dict['positive'] = True

            #
            if assume_dict:
                s.assume(**assume_dict)

        # create Order instance:
        obj = Basic.__new__(cls, expr, *symbols, **assumptions)

        # cache univariate Order symbols:
        if len(symbols)>1:
            for s in symbols:
                Order(expr, s)._get_cache_index(s)
        elif symbols:
            obj._get_cache_index(symbols[0])

        # cache multivariate Order symbols:
        cache.append(obj)
        Order._cache[symbols] = cache
        
        return obj

    def _get_cache_index(obj, symbol):
        if len(obj.symbols)>1:
            obj = Order(obj.expr, symbol)
        elif not obj.symbols:
            obj = Order(obj.expr, symbol)
        cache = Order._cache.get(symbol,[])
        try: return cache.index(obj)
        except ValueError: pass
        i = -1
        for o in cache:
            i += 1
            l = (obj.expr/o.expr).limit(symbol, 0, direction='<')
            if l.is_unbounded:
                cache.insert(i,obj)
                break
            if l.is_bounded:
                continue
            print obj.expr/o.expr,l
            raise NotImplementedError("failed to determine the inclusion relation between %s and %s (got lim=%s)" % (o, obj, l))
        else:
            cache.append(obj)
        Order._cache[symbol] = cache
        return cache.index(obj)

    @property
    def expr(self):
        return self._args[0]

    @property
    def symbols(self):
        return self._args[1:]

    def tostr(self, level = 0):
        if len(self.symbols) <= 1:
            r = 'O(%s)' % self.expr.tostr()
        else:
            r = 'O(%s)' % (', '.join([s.tostr() for s in self]))
        if self.precedence <= level:
            r = '(%s)' % (r)
        return r

    def _eval_power(b, e):
        if isinstance(e, Basic.Number):
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

    @cache_it_immutable
    def contains(self, expr):
        """
        Return True if expr belongs to Order(self.expr, *self.symbols).
        Return False if self belongs to expr.
        Return None if the inclusion relation cannot be determined (e.g. when self and
        expr have different symbols).
        """
        if isinstance(expr, Basic.Zero):
            return True
        if isinstance(expr, Order):
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
                i1 = self._get_cache_index(s)
                i2 = expr._get_cache_index(s)
                if r is None:
                    r = (i1<=i2)
                else:
                    if r != (i1<=i2):
                        return None
            return r
        obj = Order(expr, *self.symbols)
        return self.contains(obj)

    def _eval_subs(self, old, new):
        if self==old:
            return new
        if isinstance(old, Basic.Symbol) and old in self.symbols:
            i = list(self.symbols).index(old)
            if isinstance(new, Basic.Symbol):
                return Order(self.expr.subs(old, new), *(self.symbols[:i]+(new,)+self.symbols[i+1:]))
            return Order(self.expr.subs(old, new), *(self.symbols[:i]+self.symbols[i+1:]))
        return Order(self.expr.subs(old, new), *self.symbols)

    def _calc_splitter(self, d):
        return Basic.Zero()

Basic.singleton['O'] = lambda : Order
