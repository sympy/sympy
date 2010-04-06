from basic import S, C
from expr import Expr
from sympify import _sympify
from cache import cacheit

# from add import Add   /cyclic/
# from mul import Mul   /cyclic/
# from function import Lambda, WildFunction /cyclic/

class AssocOp(Expr):
    """ Associative operations, can separate noncommutative and
    commutative parts.

    (a op b) op c == a op (b op c) == a op b op c.

    Base class for Add and Mul.
    """

    # for performance reason, we don't let is_commutative go to assumptions,
    # and keep it right here
    __slots__ = ['is_commutative']

    @cacheit
    def __new__(cls, *args, **assumptions):
        if assumptions.get('evaluate') is False:
            return Expr.__new__(cls, *map(_sympify, args), **assumptions)
        if len(args)==0:
            return cls.identity()
        if len(args)==1:
            return _sympify(args[0])
        c_part, nc_part, order_symbols = cls.flatten(map(_sympify, args))
        if len(c_part) + len(nc_part) <= 1:
            if c_part: obj = c_part[0]
            elif nc_part: obj = nc_part[0]
            else: obj = cls.identity()
        else:
            obj = Expr.__new__(cls, *(c_part + nc_part), **assumptions)
            obj.is_commutative = not nc_part

        if order_symbols is not None:
            obj = C.Order(obj, *order_symbols)
        return obj


    def _new_rawargs(self, *args):
        """create new instance of own class with args exactly as provided by caller

           This is handy when we want to optimize things, e.g.

           >>> from sympy import Mul, symbols
           >>> from sympy.abc import x, y
           >>> e = Mul(3,x,y)
           >>> e.args
           (3, x, y)
           >>> Mul(*e.args[1:])
           x*y
           >>> e._new_rawargs(*e.args[1:])  # the same as above, but faster
           x*y

        """
        obj = Expr.__new__(type(self), *args)  # NB no assumptions for Add/Mul
        obj.is_commutative = self.is_commutative

        return obj

    @classmethod
    def identity(cls):
        from mul import Mul
        from add import Add
        if cls is Mul: return S.One
        if cls is Add: return S.Zero
        if cls is C.Composition:
            from symbol import Symbol
            s = Symbol('x',dummy=True)
            return Lambda(s,s)
        raise NotImplementedError("identity not defined for class %r" % (cls.__name__))

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
        # c_part, nc_part, order_symbols
        return [], new_seq, None

    _eval_subs = Expr._seq_subs

    def _matches_commutative(self, expr, repl_dict={}, evaluate=False):
        """
        Matches Add/Mul "pattern" to an expression "expr".

        repl_dict ... a dictionary of (wild: expression) pairs, that get
                      returned with the results
        evaluate .... if True, then repl_dict is first substituted into the
                      pattern, and then _matches_commutative is run

        This function is the main workhorse for Add/Mul.

        For instance:

        >> from sympy import symbols, Wild, sin
        >> a = Wild("a")
        >> b = Wild("b")
        >> c = Wild("c")
        >> x, y, z = symbols("x y z")
        >> (a+b*c)._matches_commutative(x+y*z)
        {a_: x, b_: y, c_: z}

        In the example above, "a+b*c" is the pattern, and "x+y*z" is the
        expression. Some more examples:

        >> (a+b*c)._matches_commutative(sin(x)+y*z)
        {a_: sin(x), b_: y, c_: z}
        >> (a+sin(b)*c)._matches_commutative(x+sin(y)*z)
        {a_: x, b_: y, c_: z}

        The repl_dict contains parts, that were already matched, and the
        "evaluate=True" kwarg tells _matches_commutative to substitute this
        repl_dict into pattern. For example here:

        >> (a+b*c)._matches_commutative(x+y*z, repl_dict={a: x}, evaluate=True)
        {a_: x, b_: y, c_: z}

        _matches_commutative substitutes "x" for "a" in the pattern and calls
        itself again with the new pattern "x+b*c" and evaluate=False (default):

        >> (x+b*c)._matches_commutative(x+y*z, repl_dict={a: x})
        {a_: x, b_: y, c_: z}

        the only function of the repl_dict now is just to return it in the
        result, e.g. if you omit it:

        >> (x+b*c)._matches_commutative(x+y*z)
        {b_: y, c_: z}

        the "a: x" is not returned in the result, but otherwise it is
        equivalent.

        """
        # apply repl_dict to pattern to eliminate fixed wild parts
        if evaluate:
            return self.subs(repl_dict.items()).matches(expr, repl_dict)

        # handle simple patterns
        d = self._matches_simple(expr, repl_dict)
        if d is not None:
            return d

        # eliminate exact part from pattern: (2+a+w1+w2).matches(expr) -> (w1+w2).matches(expr-a-2)
        wild_part = []
        exact_part = []
        from function import WildFunction
        from symbol import Wild
        for p in self.args:
            if p.atoms(Wild, WildFunction):
                # not all Wild should stay Wilds, for example:
                # (w2+w3).matches(w1) -> (w1+w3).matches(w1) -> w3.matches(0)
                if (not p in repl_dict) and (not p in expr):
                    wild_part.append(p)
                    continue

            exact_part.append(p)

        if exact_part:
            newpattern = self.__class__(*wild_part)
            newexpr = self.__class__._combine_inverse(expr, self.__class__(*exact_part))
            return newpattern.matches(newexpr, repl_dict)

        # now to real work ;)
        if isinstance(expr, self.__class__):
            expr_list = list(expr.args)
        else:
            expr_list = [expr]

        while expr_list:
            last_op = expr_list.pop()
            tmp = wild_part[:]
            while tmp:
                w = tmp.pop()
                d1 = w.matches(last_op, repl_dict)
                if d1 is not None:
                    d2 = self.subs(d1.items()).matches(expr, d1)
                    if d2 is not None:
                        return d2

        return

    def _eval_template_is_attr(self, is_attr):
        # return True if all elements have the property
        r = True
        for t in self.args:
            a = getattr(t, is_attr)
            if a is None: return
            if r and not a: r = False
        return r

    _eval_evalf = Expr._seq_eval_evalf

class ShortCircuit(Exception):
    pass

class LatticeOp(AssocOp):
    """
    Join/meet operations of an algebraic lattice[1].

    These binary operations are associative (op(op(a, b), c) = op(a, op(b, c))),
    commutative (op(a, b) = op(b, a)) and idempotent (op(a, a) = op(a) = a).
    Common examples are AND, OR, Union, Intersection, max or min. They have an
    identity element (op(identity, a) = a) and an absorbing element
    conventionally called zero (op(zero, a) = zero).

    This is an abstract base class, concrete derived classes must declare
    attributes zero and identity. All defining properties are then respected.

    >>> from sympy import Integer
    >>> from sympy.core.operations import LatticeOp
    >>> class my_join(LatticeOp):
    ...     zero = Integer(0)
    ...     identity = Integer(1)
    >>> my_join(2, 3) == my_join(3, 2)
    True
    >>> my_join(2, my_join(3, 4)) == my_join(2, 3, 4)
    True
    >>> my_join(0, 1, 4, 2, 3, 4)
    0
    >>> my_join(1, 2)
    2

    References:

    [1] - http://en.wikipedia.org/wiki/Lattice_(order)
    """

    is_commutative = True

    def __new__(cls, *args, **assumptions):
        args = (_sympify(arg) for arg in args)
        try:
            _args = frozenset(cls._new_args_filter(args))
        except ShortCircuit:
            return cls.zero
        if not _args:
            return cls.identity
        elif len(_args) == 1:
            return set(_args).pop()
        else:
            obj = Expr.__new__(cls, _args, **assumptions)
            obj._argset = _args
            return obj

    @classmethod
    def _new_args_filter(cls, arg_sequence):
        """Generator filtering args"""
        for arg in arg_sequence:
            if arg == cls.zero:
                raise ShortCircuit(arg)
            elif arg == cls.identity:
                continue
            elif arg.func == cls:
                for x in arg.iter_basic_args():
                    yield x
            else:
                yield arg

    @property
    def args(self):
        return tuple(self._argset)

    @staticmethod
    def _compare_pretty(a, b):
        return cmp(str(a), str(b))
