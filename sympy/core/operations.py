
from basic import Basic, cache_it, cache_it_immutable

class AssocOp(Basic):
    """ Associative operations, can separate noncommutative and
    commutative parts.
    
    (a op b) op c == a op (b op c) == a op b op c.
    
    Base class for Add and Mul.
    """

    @cache_it_immutable
    def __new__(cls, *args, **assumptions):
        if len(args)==0:
            return cls.identity()
        if len(args)==1:
            return Basic.sympify(args[0])
        c_part, nc_part, lambda_args, order_symbols = cls.flatten(map(Basic.sympify, args))
        if len(c_part) + len(nc_part) <= 1:
            if c_part: obj = c_part[0]
            elif nc_part: obj = nc_part[0]
            else: obj = cls.identity()
        else:
            assumptions['commutative'] = not nc_part
            obj = Basic.__new__(cls, *(c_part + nc_part), **assumptions)
        if order_symbols is not None:
            obj = Basic.Order(obj, *order_symbols)
        if lambda_args is not None:
            obj = Basic.Lambda(obj, *lambda_args)
        return obj

    @classmethod
    def identity(cls):
        if cls is Basic.Mul: return Basic.One()
        if cls is Basic.Add: return Basic.Zero()
        if cls is Basic.Composition:
            s = Basic.Symbol('x',dummy=True)
            return Basic.Lambda(s,s)
        raise NotImplementedError,"identity not defined for class %r" % (cls.__name__)

    @classmethod
    def flatten(cls, seq):
        # apply associativity, no commutativity property is used
        new_seq = []
        while seq:
            o = seq.pop(0)
            if o.__class__ is cls: # classes must match exactly
                seq = list(o[:]) + seq
                continue
            new_seq.append(o)
        return [], new_seq, None, None

    _eval_subs = Basic._seq_subs

    def _matches_commutative(pattern, expr, repl_dict={}, evaluate=False):
        # apply repl_dict to pattern to eliminate fixed wild parts
        if evaluate:
            pat = pattern
            for old,new in repl_dict.items():
                pat = pat.subs(old, new)
            if pat != pattern:
                return pat.matches(expr, repl_dict)

        # handle simple patterns
        d = pattern._matches_simple(expr, repl_dict)
        if d is not None:
            return d

        # eliminate exact part from pattern: (2+a+w1+w2).matches(expr) -> (w1+w2).matches(expr-a-2)
        wild_part = []
        exact_part = []
        for p in pattern:
            if p.atoms(type=(Basic.Wild, Basic.WildFunction)):
                wild_part.append(p)
            else:
                exact_part.append(p)
        if exact_part:
            newpattern = pattern.__class__(*wild_part)
            newexpr = pattern.__class__._combine_inverse(expr, pattern.__class__(*exact_part))
            return newpattern.matches(newexpr, repl_dict)

        # now to real work ;)
        if isinstance(expr, pattern.__class__):
            expr_list = list(expr)
        else:
            expr_list = [expr]

        while expr_list:
            last_op = expr_list.pop()
            tmp = wild_part[:]
            while tmp:
                w = tmp.pop()
                d1 = w.matches(last_op, repl_dict)
                if d1 is not None:
                    d2 = pattern.matches(expr, d1, evaluate=True)
                    if d2 is not None:
                        return d2
        return

    def _eval_template_is_attr(self, is_attr):
        # return True if all elements have the property
        r = True
        for t in self:
            a = getattr(t, is_attr)
            if a is None: return
            if r and not a: r = False
        return r

    _eval_evalf = Basic._seq_eval_evalf
    
