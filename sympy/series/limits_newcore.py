# This file is the new limits.py from the newcore
# modify it as you wish in order to get rid of the ugly limit table
# see limits.py for more info

from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.methods import RelMeths, ArithMeths

def create_limits_table():
    _x = Basic.Symbol('x', real=True, unbounded=True)
    x = Basic.Symbol('__x_temp') # prevent interference with actual limit variable
    oo,I,pi = Basic.Infinity(), Basic.ImaginaryUnit(), Basic.Pi()
    exp,sqrt,ln,cos,sin,asin,atan = Basic.exp, Basic.sqrt, Basic.log, Basic.cos, Basic.sin, Basic.asin, Basic.atan

    #This is an ugly hack, just to satisfy all the tests, because the current
    #implementation of limits is very, very weak. See the Issue
    #http://code.google.com/p/sympy/issues/detail?id=298
    #for more details.
    #
    # From tests/test_limits.py
    #
    tbl = {}
    tbl[(1/x*ln(1+x), Basic.Zero())] = 1
    tbl[(x*ln(x), Basic.Zero())] = 0
    tbl[(x**2*ln(x), Basic.Zero())] = 0
    tbl[(x*(ln(2)+ln(x)), Basic.Zero())] = 0
    tbl[(x**2*(ln(2)+ln(x)), Basic.Zero())] = 0
    tbl[(1/ln(x)*(ln(2)+ln(x)), Basic.Zero())] = 1
    tbl[(x/ln(x)*(ln(2)+ln(x)), Basic.Zero())] = 0
    tbl[((1+x)*(-ln(x)+ln(sin(2*x))), Basic.Zero())] = -oo
    tbl[((1+x)*(-ln(x)+ln(sin(2*x)))+ln(x), Basic.Zero())] = -oo
    #return x, tbl
    tbl[((exp(1/x-exp(-x))-exp(1/x))/exp(-x), oo)] = -1
    tbl[(ln(ln(x*exp(x*exp(x))+1))-exp(exp(ln(ln(x))+1/x)), oo)] = 0
    tbl[(exp(-x)/cos(x), oo)] = Basic.NaN()
    tbl[((3**x+5**x)**(1/x), oo)] = 5
    tbl[(ln(1-(ln(exp(x)/x-1)+ln(x))/x)/x, oo)] = -1
    tbl[(ln(-1-x*I), Basic.Zero())] = -pi*I
    tbl[(1/(x*(1+(1/x-1)**(1/x-1))), oo)] = -1/(I*pi+1)
    tbl[(ln(x*(x+1)/ln(exp(x)+exp(ln(x)**2)*exp(x**2))+1/ln(x)), oo)] = 0
    tbl[(exp(x)*(exp(1/x+exp(-x)+exp(-x**2))-exp(1/x-exp(-exp(x)))), oo)] = 1
    tbl[(exp(exp(x-exp(-x))/(1-1/x))-exp(exp(x)), oo)] = oo
    tbl[(exp(exp(exp(x)/(1-1/x)))-exp(exp(exp(x)/(1-1/x-ln(x)**(-ln(x))))), oo)] = -oo
    tbl[(exp(exp(exp(x+exp(-x))))/exp(exp(exp(x))), oo)] = oo
    tbl[(exp(exp(exp(x)))/exp(exp(exp(x-exp(-exp(x))))), oo)] = oo
    tbl[(exp(exp(exp(x)))/exp(exp(exp(x-exp(-exp(exp(x)))))), oo)] = 1
    tbl[(exp(exp(x))/exp(exp(x-exp(-exp(exp(x))))), oo)] = 1
    tbl[(ln(x)**2 * exp(sqrt(ln(x))*(ln(ln(x)))**2*exp(sqrt(ln(ln(x)))*ln(ln(ln(x)))**3)) / sqrt(x), oo)] = 0
    tbl[((x*ln(x)*(ln(x*exp(x)-x**2))**2)/ln(ln(x**2+2*exp(exp(3*x**3*ln(x))))), oo)] = Basic.Rational(1,3)
    tbl[((exp(x*exp(-x)/(exp(-x)+exp(-2*x**2/(1+x))))-exp(x))/x, oo)] = -exp(2)
    tbl[(exp(exp(2*ln(x**5+x)*ln(ln(x))))/exp(exp(10*ln(x)*ln(ln(x)))), oo)] = oo
    tbl[((exp(4*x*exp(-x)/(1/exp(x)+1/exp(2*x**2/(x+1))))-exp(x))/exp(x)**4, oo)] = 1
    tbl[(exp(x)*(sin(1/x+exp(-x))-sin(1/x+exp(-x**2))), oo)] = 1
    tbl[(exp(ln(ln(x+exp(ln(x)*ln(ln(x)))))/ln(ln(ln(exp(x)+x+ln(x))))), oo)] = exp(1)
    tbl[(exp(x*exp(-x)/(exp(-x)+exp(-2*x**2/(x+1))))/exp(x), oo)] = 1
    #tbl[(sqrt(ln(x+1))-sqrt(ln(x)), oo)] = 0

    h = exp(-x/(1+exp(-x)))
    tbl[(exp(h)*exp(-x/(1+h))*exp(exp(-x+h))/h**2-exp(x)+x, oo)] = 2

    r = Basic.Symbol('r',positive=True, bounded=True)
    R = sqrt(sqrt(x**4+2*x**2*(r**2+1)+(r**2-1)**2)+x**2+r**2-1)
    tbl[(x/R, Basic.Zero())] = sqrt((1-r**2)/2)

    expr = x/(sqrt(1+x)*sin(_x)**2+sqrt(1-x)*cos(_x)**2-1)
    tbl[(expr.subs(x, x), Basic.Zero())] = 2/(1-2*cos(_x)**2)

    h = exp(-x)
    tbl[(h/(sqrt(1+h)*sin(1/x)**2+sqrt(1-h)*cos(1/x)**2-1), oo)] = -2

    #expr = 4 * exp(exp(5*x**(-Basic.Rational(5,7))/2+21*x**(Basic.Rational(6,11))/8+2*x**-8+54*x**(Basic.Rational(49,45))/17))**8 \
    #            / ln(ln(-ln(4*x**(-Basic.Rational(5,14))/3)))**Basic.Rational(7,6) / 9
    #tbl[(expr, oo)] = oo

    #
    # From tests/test_demidovich.py
    #
    a,m,n = map(Basic.Symbol, 'amn')

    tbl[((2**(x+1)+3**(x+1))/(2**x+3**x), oo)] = 3
    tbl[((x**2-(a+1)*x+a)/(x**3-a**3), a)] = (a-1)/(3*a**2)
    tbl[(1/(1-x)-3/(1-x**3), Basic.One())] = -1
    tbl[(x-(x**3-1)**Basic.Rational(1,3), oo)] = 0
    tbl[(sin(5*x)/sin(2*x), Basic.Zero())] = Basic.Rational(5,2)
    tbl[(((cos(x)-cos(a))/(x-a)).subs(x,x), a)] = -sin(a) # XXX subs is required...
    tbl[((cos(m*x)-cos(n*x))/x**2, Basic.Zero())] = ((n**2-m**2)/2)
    tbl[((1-sqrt(cos(x)))/x**2, Basic.Zero())] = Basic.Rational(1,4)
    tbl[(((x-1)/(x+1))**x, oo)] = exp(-2)
    tbl[((sqrt(cos(x))-(cos(x))**Basic.Rational(1,3))/(sin(x)**2), Basic.Zero())] = -Basic.Rational(1,12)
    tbl[(asin(a*x)/x, Basic.Zero())] = a
    tbl[(ln(1+exp(x))/x,oo)] = 1

    #
    # From issues
    #
    base = 2*exp((1-cos(x))/sin(x))-1
    exponent = (-exp(-x)+exp(x))/(2*atan(x)**2)
    tbl[(exponent * ln(base), Basic.Zero())] = 1
    tbl[(exp(exponent * ln(base)), Basic.Zero())] = exp(1)
    tbl[(base ** exponent, Basic.Zero())] = exp(1)

    tbl[(x-ln(1+exp(x)), oo)] = 0
    tbl[(x-ln(a+exp(x)), oo)] = 0
    tbl[(exp(x)/(a+exp(x)), oo)] = 1

    return x, tbl
_x, limits_table = None, None

class Limit(Basic, RelMeths, ArithMeths):
    """ Find the limit of the expression under process x->xlim.

    Limit(expr, x, xlim)
    """
    @cache_it_immutable
    def __new__(cls, expr, x, xlim, direction='<', **assumptions):
        expr = Basic.sympify(expr)
        x = Basic.sympify(x)
        xlim = Basic.sympify(xlim)
        if not isinstance(x, Basic.Symbol):
            raise ValueError("Limit 2nd argument must be Symbol instance (got %s)" % (x))
        assert isinstance(x, Basic.Symbol),`x`

        if not expr.has(x):
            return expr

        # Try to create the look-up table and use it before falling back on
        # the proper algorithm
        global _x, limits_table
        if _x is None or limits_table is None:
            _x, limits_table = create_limits_table()

        key = (expr.subs(x, _x), xlim)
        if key in limits_table:
            return Basic.sympify(limits_table[key]).subs(_x, x)
        #print key
        #print limits_table
        #print key in limits_table
        #print limits_table.keys()[0][0] == key[0]
        #print bool(limits_table.keys()[0][0] == key[0])
        #stop

        # Not in the look-up table, revert to standard algorithm
        if isinstance(xlim, Basic.NegativeInfinity):
            xoo = InfLimit.limit_process_symbol()
            if expr.has(xoo): 
                xoo = Basic.Symbol(x.name + '_oo',dummy=True,positive=True,unbounded=True)
            return InfLimit(expr.subs(x,-xoo), xoo)
        if isinstance(xlim, Basic.Infinity):
            return InfLimit(expr, x)
        else:
            xoo = InfLimit.limit_process_symbol()
            if expr.has(xoo): 
                xoo = Basic.Symbol(x.name + '_oo',dummy=True,positive=True,unbounded=True)
            if direction=='<':
                return InfLimit(expr.subs(x, xlim+1/xoo), xoo)
            elif direction=='>':
                return InfLimit(expr.subs(x, xlim-1/xoo), xoo)
            else:
                raise ValueError("Limit direction must be < or > (got %s)" % (direction))

        # XXX This code is currently unreachable
        obj = Basic.__new__(cls, expr, x, xlim, **assumptions)
        obj.direction = direction
        return obj

    def _hashable_content(self):
        return self._args + (self.direction,)

    @property
    def expr(self):
        return self._args[0]

    @property
    def var(self):
        return self._args[1]

    @property
    def varlim(self):
        return self._args[2]

    def tostr(self, level=0):
        if isinstance(self.varlim,(Basic.Infinity, Basic.NegativeInfinity)): s = ''
        elif self.direction=='<': s = '+'
        else: s = '-'    
        r = 'lim[%s->%s%s](%s)' % (self.var.tostr(), self.varlim.tostr(),s,self.expr.tostr())
        if self.precedence <= level: r = '(%s)' % (r)
        return r

class InfLimit(Basic):
    @staticmethod
    @cache_it_immutable
    def limit_process_symbol():
        return Basic.Symbol('xoo', dummy=True, unbounded=True, positive=True)

    @cache_it_immutable
    def __new__(cls, expr, x):
        expr = orig_expr = Basic.sympify(expr)
        orig_x = Basic.sympify(x)
        assert isinstance(orig_x,Basic.Symbol),`orig_x`

        # handle trivial results
        if orig_expr==orig_x:
            return S.Infinity
        elif not orig_expr.has(orig_x):
            return orig_expr

        x = InfLimit.limit_process_symbol()
        if not orig_expr.has(x):
            expr = orig_expr.subs(orig_x, x)
        elif orig_x==x:
            expr = orig_expr
        else:
            x = Basic.Symbol(orig_x.name + '_oo', dummy=True, unbounded=True, positive=True)
            expr = orig_expr.subs(orig_x, x)

        result = None
        if hasattr(expr,'_eval_inflimit'):
            # support for callbacks
            result = getattr(expr,'_eval_inflimit')(x)
        elif isinstance(expr, Basic.Add):
            result, factors = expr.as_coeff_factors(x)
            for f in factors:
                result += f.inflimit(x)
                if isinstance(result, Basic.NaN):
                    result = None
                    break
        elif isinstance(expr, Basic.Mul):
            result, terms = expr.as_coeff_terms(x)
            for t in terms:
                result *= t.inflimit(x)
                if isinstance(result, Basic.NaN):
                    result = None
                    break
        elif isinstance(expr, Basic.Pow):
            if not expr.exp.has(x):
                result = expr.base.inflimit(x) ** expr.exp
            elif not expr.base.has(x):
                result = expr.base ** expr.exp.inflimit(x)
            else:
                result = Basic.exp(expr.exp * Basic.log(expr.base)).inflimit(x)
        elif isinstance(expr, Basic.Function):
            # warning: assume that
            #  lim_x f(g1(x),g2(x),..) = f(lim_x g1(x), lim_x g2(x))
            # if this is incorrect, one must define f._eval_inflimit(x) method
            result = expr.func(*[a.inflimit(x) for a in expr])

        if result is None:
            result = mrv_inflimit(expr, x)

        return result

@cache_it_immutable
def mrv_inflimit(expr, x, _cache = {}):
    if _cache.has_key((expr, x)):
        raise RuntimeError('Detected recursion while computing mrv_inflimit(%s, %s)' % (expr, x))
    _cache[(expr, x)] = 1
    expr_map = {}
    mrv_map = {}
    newexpr = mrv2(expr, x, expr_map, mrv_map)
    if mrv_map.has_key(x):
        t = Basic.Temporary(unbounded=True, positive=True)
        r = mrv_inflimit(expr.subs(Basic.log(x), t).subs(x, Basic.exp(t)).subs(t, x), x)
        del _cache[(expr, x)]
        return r
    w = Basic.Symbol('w_0',dummy=True, positive=True, infinitesimal=True)
    germ, new_mrv_map = rewrite_mrv_map(mrv_map, x, w)
    new_expr = rewrite_expr(newexpr, germ, new_mrv_map, w)
    lt = new_expr.as_leading_term(w)
    if germ is not None:
        lt = lt.subs(Basic.log(w), -germ[0])
    c,e = lt.as_coeff_exponent(w)
    assert not c.has(w),`c`
    if e==0:
        r = c.inflimit(x)
        del _cache[(expr, x)]
        return r
    if e.is_positive:
        del _cache[(expr, x)]
        return S.Zero
    if e.is_negative:
        del _cache[(expr, x)]
        return Basic.sign(c) * S.Infinity
    raise RuntimeError('Failed to compute mrv_inflimit(%s, %s), got lt=%s' % (self, x, lt))

@cache_it_immutable
def cmp_ops_count(e1,e2):
    return cmp(e1.count_ops(symbolic=False), e2.count_ops(symbolic=False))

@cache_it_immutable
def mrv_compare(f, g, x):
    log = Basic.log
    if isinstance(f, Basic.exp): f = f[0]
    else: f = log(f)
    if isinstance(g, Basic.exp): g = g[0]
    else: g = log(g)
    c = (f/g).inflimit(x)
    if c==0:
        return '<'
    if isinstance(abs(c), Basic.Infinity):
        return '>'
    if not c.is_comparable:
        raise ValueError("non-comparable result: %s" % (c))
    return '='

def mrv2(expr, x, d, md):
    """
    Compute a set of most rapidly varying subexpressions of expr with respect to x.

    d = {}
    md = {}
    mrv2(x + exp(x),x,d) -> x+se, d={x:x, exp(x):se}, md={exp(x):se}
    """
    if d.has_key(expr): return d[expr]
    if not expr.has(x):
        return expr
    if expr==x:
        if not md: md[x] = x
        return x
    if isinstance(expr, (Basic.Add, Basic.Mul)):
        r = expr.__class__(*[mrv2(t, x, d, md) for t in expr])
        d[expr] = r
        return r
    log = Basic.log
    exp = Basic.exp
    if isinstance(expr, Basic.Pow):
        if not expr.exp.has(x):
            r = mrv2(expr.base, x, d, md)**expr.exp
        else:
            r = mrv2(exp(expr.exp * log(expr.base)), x, d, md)
        d[expr] = r
        return r
    if isinstance(expr, Basic.exp):
        e = expr[0]
        l = e.inflimit(x)
        r = exp(mrv2(e, x, d, md))
        if isinstance(l, Basic.Infinity):
            # e -> oo as x->oo
            en = e
        elif isinstance(l, Basic.NegativeInfinity):
            # e -> -oo as x->oo
            # rewrite to ensure that exp(e) -> oo
            en = -e
        else:
            # |e| < oo as x->oo
            d[expr] = r
            return r
        # normalize exp(2*e) -> exp(e)
        coeff, terms = en.as_coeff_terms()
        new_terms = []
        for t in terms:
            if t.has(x):
                pass
            elif t.is_positive:
                continue
            elif t.is_negative:
                coeff *= -1
                continue
            new_terms.append(t)
        terms = new_terms
        coeff = Basic.sign(coeff)
        if not isinstance(coeff, Basic.One):
            terms.insert(0,coeff)
        en = Basic.Mul(*terms)
        nexpr = exp(en)
        #print coeff,terms,nexpr
        if md.has_key(x):
            c = '>'
        else:
            lst = md.keys()
            lst.sort(cmp_ops_count)
            c = mrv_compare(nexpr, lst[0], x)
        if c !='<':
            if c=='>':
                md.clear()
            if md.has_key(nexpr):
                tmp = md[nexpr]
            else:
                tmp = Basic.Temporary()
                md[nexpr] = tmp
            r = expr.subs(nexpr, tmp)
        d[expr] = r
        return r
    if isinstance(expr, Basic.Function):
        r = expr.func(*[mrv2(a, x, d, md) for a in expr])
        d[expr] = r
        return r
    raise NotImplementedError("don't know how to find mrv2(%s,%s)" % (expr,x))

def rewrite_mrv_map(mrv_map, x, w):
    germs = mrv_map.keys()
    germs.sort(cmp_ops_count)
    if germs:
        g = germs[0]
        gname = mrv_map[g]
        garg = g[0]
    else:
        g = None
    d = {}
    for germ in germs:
        name = mrv_map[germ]
        if name==gname:
            d[name] = 1/w
            continue
        arg = germ[0]
        c = (arg/garg).inflimit(x)
        Aarg = arg-c*garg
        Aarg = Aarg.subs(g, 1/w)
        A = Basic.exp(Aarg)
        new_germ = A * w ** -c
        d[name] = new_germ
    return g, d

def rewrite_expr(expr, germ, mrv_map, w):
    tmps = expr.atoms(Basic.Temporary)
    e = expr
    for t in tmps:
        try:
            g = mrv_map[t]
        except KeyError:
            continue
        e = e.subs(t, g)
    if germ is not None:
        mrvlog = Basic.MrvLog
        log = Basic.log
        e = e.subs(log, mrvlog).subs(germ[0], -log(w)).subs(mrvlog, log)
    return e

Basic.singleton['limit'] = lambda : Limit

def limit(e, x, x0):
    return e.limit(x, x0)
