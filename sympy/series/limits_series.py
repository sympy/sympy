# This is a very simplified (and not that reliable) version of the limit
# algorithm that is only used for series expansion. Use this file only in the
# series expansion code. Everywere else, use the general and robust limit
# algorithm in limits.py.

# The idea is to use the general and robust algorithm in limits.py even in the
# series expansion, but currently the limits.py is slower and there were some
# other bugs (mainly recursion), if it were used in the series expansion. So
# currently we use limits_series, until we move to limits.py completely.

from sympy.core.basic import Basic, S, C, sympify
from sympy.core.methods import RelMeths, ArithMeths
from sympy.core.cache import cacheit

class Limit_series(Basic, RelMeths, ArithMeths):
    """ Find the limit of the expression under process x->xlim.

    Limit(expr, x, xlim)
    """

    __slots__ = []

    @cacheit
    def __new__(cls, expr, x, xlim, direction='<', **assumptions):
        expr = sympify(expr)
        x = sympify(x)
        xlim = sympify(xlim)
        if not x.is_Symbol:
            raise ValueError("Limit 2nd argument must be Symbol instance (got %s)" % (x))

        if not expr.has(x):
            return expr

        if xlim is S.NegativeInfinity:
            xoo = InfLimit.limit_process_symbol()
            if expr.has(xoo):
                xoo = C.Symbol(x.name + '_oo',dummy=True,positive=True,unbounded=True)
            return InfLimit(expr.subs(x,-xoo), xoo)
        if xlim is S.Infinity:
            return InfLimit(expr, x)
        else:
            xoo = InfLimit.limit_process_symbol()
            if expr.has(xoo):
                xoo = C.Symbol(x.name + '_oo',dummy=True,positive=True,unbounded=True)
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
        if (self.varlim is S.Infinity) or (self.varlim is S.NegativeInfinity): s = ''
        elif self.direction=='<': s = '+'
        else: s = '-'
        r = 'lim[%s->%s%s](%s)' % (self.var.tostr(), self.varlim.tostr(),s,self.expr.tostr())
        if self.precedence <= level: r = '(%s)' % (r)
        return r

class InfLimit(Basic):
    _xoo = C.Symbol('xoo', dummy=True, unbounded=True, positive=True)

    __slots__ = []

    @staticmethod
    def limit_process_symbol():
        return InfLimit._xoo

    @cacheit
    def __new__(cls, expr, x):
        expr = orig_expr = sympify(expr)
        orig_x = sympify(x)
        assert orig_x.is_Symbol, `orig_x`

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
            x = C.Symbol(orig_x.name + '_oo', dummy=True, unbounded=True, positive=True)
            expr = orig_expr.subs(orig_x, x)

        result = None
        if hasattr(expr,'_eval_inflimit'):
            # support for callbacks
            result = getattr(expr,'_eval_inflimit')(x)
        elif expr.is_Add:
            result, factors = expr.as_coeff_factors(x)
            for f in factors:
                result += f.inflimit(x)
                if result is S.NaN:
                    result = None
                    break
        elif expr.is_Mul:
            result, terms = expr.as_coeff_terms(x)
            for t in terms:
                result *= t.inflimit(x)
                if result is S.NaN:
                    result = None
                    break
        elif expr.is_Pow:
            if not expr.exp.has(x):
                result = expr.base.inflimit(x) ** expr.exp
            elif not expr.base.has(x):
                result = expr.base ** expr.exp.inflimit(x)
            else:
                result = C.exp(expr.exp * C.log(expr.base)).inflimit(x)
        elif expr.is_Function:
            # warning: assume that
            #  lim_x f(g1(x),g2(x),..) = f(lim_x g1(x), lim_x g2(x))
            # if this is incorrect, one must define f._eval_inflimit(x) method
            result = expr.func(*[a.inflimit(x) for a in expr.args])

        if result is None:
            result = mrv_inflimit(expr, x)

        return result

@cacheit
def mrv_inflimit(expr, x, _cache = {}):
    if _cache.has_key((expr, x)):
        raise RuntimeError('Detected recursion while computing mrv_inflimit(%s, %s)' % (expr, x))
    _cache[(expr, x)] = 1
    expr_map = {}
    mrv_map = {}
    newexpr = mrv2(expr, x, expr_map, mrv_map)
    if mrv_map.has_key(x):
        t = C.Temporary(unbounded=True, positive=True)
        r = mrv_inflimit(expr.subs(C.log(x), t).subs(x, C.exp(t)).subs(t, x), x)
        del _cache[(expr, x)]
        return r
    w = C.Symbol('w_0',dummy=True, positive=True, infinitesimal=True)
    germ, new_mrv_map = rewrite_mrv_map(mrv_map, x, w)
    new_expr = rewrite_expr(newexpr, germ, new_mrv_map, w)
    lt = new_expr.as_leading_term(w)
    if germ is not None:
        lt = lt.subs(C.log(w), -germ.args[0])
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
        return C.sign(c) * S.Infinity
    raise RuntimeError('Failed to compute mrv_inflimit(%s, %s), got lt=%s' % (self, x, lt))

@cacheit
def cmp_ops_count(e1,e2):
    return cmp(e1.count_ops(symbolic=False), e2.count_ops(symbolic=False))

@cacheit
def mrv_compare(f, g, x):
    log = C.log
    if f.func is C.exp: f = f.args[0]
    else: f = log(f)
    if g.func is C.exp: g = g.args[0]
    else: g = log(g)
    c = (f/g).inflimit(x)
    if c==0:
        return '<'
    if abs(c) is S.Infinity:
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
    if expr.is_Add or expr.is_Mul:
        r = expr.__class__(*[mrv2(t, x, d, md) for t in expr.args])
        d[expr] = r
        return r
    log = C.log
    exp = C.exp
    if expr.is_Pow:
        if not expr.exp.has(x):
            r = mrv2(expr.base, x, d, md)**expr.exp
        else:
            r = mrv2(exp(expr.exp * log(expr.base)), x, d, md)
        d[expr] = r
        return r
    if expr.func is C.exp:
        e = expr.args[0]
        l = e.inflimit(x)
        r = exp(mrv2(e, x, d, md))
        if l is S.Infinity:
            # e -> oo as x->oo
            en = e
        elif l is S.NegativeInfinity:
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
        coeff = C.sign(coeff)
        if coeff is not S.One:
            terms.insert(0,coeff)
        en = C.Mul(*terms)
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
                tmp = C.Temporary()
                md[nexpr] = tmp
            r = expr.subs(nexpr, tmp)
        d[expr] = r
        return r
    if expr.is_Function:
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
        A = C.exp(Aarg)
        new_germ = A * w ** -c
        d[name] = new_germ
    return g, d

def rewrite_expr(expr, germ, mrv_map, w):
    tmps = expr.atoms(C.Temporary)
    e = expr
    for t in tmps:
        try:
            g = mrv_map[t]
        except KeyError:
            continue
        e = e.subs(t, g)
    if germ is not None:
        mrvlog = C.MrvLog
        log = C.log
        e = e.subs(log, mrvlog).subs(germ.args[0], -log(w)).subs(mrvlog, log)
    return e

#Basic.singleton['limit'] = lambda : Limit

def limit(e, x, x0):
    return e.limit(x, x0)
