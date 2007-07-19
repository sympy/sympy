
from basic import Basic, S, cache_it, cache_it_immutable
from methods import RelMeths, ArithMeths

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
        if not expr.has(x): return expr
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
                result = S.Exp(expr.exp * S.Log(expr.base)).inflimit(x)
        elif isinstance(expr, Basic.Apply):
            # warning: assume that
            #  lim_x f(g1(x),g2(x),..) = f(lim_x g1(x), lim_x g2(x))
            # if this is incorrect, one must define f._eval_inflimit(x) method
            result = expr.func(*[a.inflimit(x) for a in expr.args])

        if result is None:
            result = mrv_inflimit(expr, x)

        return result

@cache_it_immutable
def mrv_inflimit(expr, x):
    expr_map = {}
    mrv_map = {}
    newexpr = mrv2(expr, x, expr_map, mrv_map)
    if mrv_map.has_key(x):
        t = Basic.Temporary(unbounded=True, positive=True)
        return mrv_inflimit(expr.subs(S.Log(x), t).subs(x, S.Exp(t)).subs(t, x), x)
    w = Basic.Symbol('w_0',dummy=True, positive=True, infinitesimal=True)
    germ, new_mrv_map = rewrite_mrv_map(mrv_map, x, w)
    new_expr = rewrite_expr(newexpr, germ, new_mrv_map, w)
    lt = new_expr.as_leading_term(w)
    if germ is not None:
        lt = lt.subs(S.Log(w), -germ.args[0])
    c,e = lt.as_coeff_exponent(w)
    assert not c.has(w),`c`
    if e==0:
        return c.inflimit(x)
    if e.is_positive:
        return S.Zero
    if e.is_negative:
        return S.Sign(c) * S.Infinity
    raise RuntimeError('Failed to compute mrv_inflimit(%s, %s), got lt=%s' % (self, x, lt))

@cache_it_immutable
def cmp_ops_count(e1,e2):
    return cmp(e1.count_ops(symbolic=False), e2.count_ops(symbolic=False))

@cache_it_immutable
def mrv_compare(f, g, x):
    log = S.Log
    if isinstance(f, Basic.ApplyExp): f = f.args[0]
    else: f = log(f)
    if isinstance(g, Basic.ApplyExp): g = g.args[0]
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
    log = S.Log
    exp = S.Exp
    if isinstance(expr, Basic.Pow):
        if not expr.exp.has(x):
            r = mrv2(expr.base, x, d, md)**expr.exp
        else:
            r = mrv2(exp(expr.exp * log(expr.base)), x, d, md)
        d[expr] = r
        return r
    if isinstance(expr, Basic.ApplyExp):
        e = expr.args[0]
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
        coeff = Basic.Sign()(coeff)
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
    if isinstance(expr, Basic.Apply):
        r = expr.func(*[mrv2(a, x, d, md) for a in expr.args])
        d[expr] = r
        return r
    raise NotImplementedError("don't know how to find mrv2(%s,%s)" % (expr,x))

def rewrite_mrv_map(mrv_map, x, w):
    germs = mrv_map.keys()
    germs.sort(cmp_ops_count)
    if germs:
        g = germs[0]
        gname = mrv_map[g]
        garg = g.args[0]
    else:
        g = None
    d = {}
    for germ in germs:
        name = mrv_map[germ]
        if name==gname:
            d[name] = 1/w
            continue
        arg = germ.args[0]
        c = (arg/garg).inflimit(x)
        Aarg = arg-c*garg
        Aarg = Aarg.subs(g, 1/w)
        A = S.Exp(Aarg)
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
        mrvlog = S.MrvLog
        log = S.Log
        e = e.subs(log, mrvlog).subs(germ.args[0], -log(w)).subs(mrvlog, log)
    return e

Basic.singleton['limit'] = lambda : Limit
