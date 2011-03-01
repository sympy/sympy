from core import C
from singleton import S
from expr import Expr
from sympify import _sympify, sympify
from cache import cacheit
from compatibility import all

# from add import Add   /cyclic/
# from mul import Mul   /cyclic/
# from function import Lambda, WildFunction /cyclic/

class AssocOp(Expr):
    """ Associative operations, can separate noncommutative and
    commutative parts.

    (a op b) op c == a op (b op c) == a op b op c.

    Base class for Add and Mul.

    This is an abstract base class, concrete derived classes must define
    the attribute `identity`.
    """

    # for performance reason, we don't let is_commutative go to assumptions,
    # and keep it right here
    __slots__ = ['is_commutative']

    @cacheit
    def __new__(cls, *args, **assumptions):
        if len(args) == 0:
            return cls.identity

        args = map(_sympify, args)
        if len(args) == 1:
            return args[0]

        if not assumptions.pop('evaluate', True):
            obj = Expr.__new__(cls, *args, **assumptions)
            obj.is_commutative = all(a.is_commutative for a in args)
            return obj

        c_part, nc_part, order_symbols = cls.flatten(args)
        if len(c_part) + len(nc_part) <= 1:
            if c_part:
                obj = c_part[0]
            elif nc_part:
                obj = nc_part[0]
            else:
                obj = cls.identity
        else:
            obj = Expr.__new__(cls, *(c_part + nc_part), **assumptions)
            obj.is_commutative = not nc_part

        if order_symbols is not None:
            obj = C.Order(obj, *order_symbols)
        return obj

    def _new_rawargs(self, *args, **kwargs):
        """create new instance of own class with args exactly as provided by caller
           but returning the self class identity if args is empty.

           This is handy when we want to optimize things, e.g.

               >>> from sympy import Mul, symbols, S
               >>> from sympy.abc import x, y
               >>> e = Mul(3, x, y)
               >>> e.args
               (3, x, y)
               >>> Mul(*e.args[1:])
               x*y
               >>> e._new_rawargs(*e.args[1:])  # the same as above, but faster
               x*y

           Note: use this with caution. There is no checking of arguments at
           all. This is best used when you are rebuilding an Add or Mul after
           simply removing one or more terms. If modification which result,
           for example, in extra 1s being inserted (as when collecting an
           expression's numerators and denominators) they will not show up in
           the result but a Mul will be returned nonetheless:

               >>> m = (x*y)._new_rawargs(S.One, x); m
               x
               >>> m == x
               False
               >>> m.is_Mul
               True

           Another issue to be aware of is that the commutativity of the result
           is based on the commutativity of self. If you are rebuilding the
           terms that came from a commutative object then there will be no
           problem, but if self was non-commutative then what you are
           rebuilding may now be commutative.

           Although this routine tries to do as little as possible with the
           input, getting the commutativity right is important, so this level
           of safety is enforced: commutativity will always be recomputed if
           either a) self has no is_commutate attribute or b) self is
           non-commutative and kwarg `reeval=False` has not been passed.

           If you don't have an existing Add or Mul and need one quickly, try
           the following.

               >>> m = object.__new__(Mul)
               >>> m._new_rawargs(x, y)
               x*y

           Note that the commutativity is always computed in this case since
           m doesn't have an is_commutative attribute; reeval is ignored:

               >>> _.is_commutative
               True
               >>> hasattr(m, 'is_commutative')
               False
               >>> m._new_rawargs(x, y, reeval=False).is_commutative
               True

           It is possible to define the commutativity of m. If it's False then
           the new Mul's commutivity will be re-evaluated:

               >>> m.is_commutative = False
               >>> m._new_rawargs(x, y).is_commutative
               True

           But if reeval=False then a non-commutative self can pass along
           its non-commutativity to the result (but at least you have to *work*
           to get this wrong):

               >>> m._new_rawargs(x, y, reeval=False).is_commutative
               False

        """
        if len(args) > 1:
            obj = Expr.__new__(type(self), *args)  # NB no assumptions for Add/Mul

            if (hasattr(self, 'is_commutative') and
                (self.is_commutative or
                not kwargs.pop('reeval', True))):
                obj.is_commutative = self.is_commutative
            else:
                obj.is_commutative = all(a.is_commutative for a in args)

        elif len(args) == 1:
            obj = args[0]
        else:
            obj = self.identity

        return obj


    @classmethod
    def flatten(cls, seq):
        # apply associativity, no commutivity property is used
        new_seq = []
        while seq:
            o = seq.pop(0)
            if o.__class__ is cls: # classes must match exactly
                seq = list(o[:]) + seq
                continue
            new_seq.append(o)
        # c_part, nc_part, order_symbols
        return [], new_seq, None

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
        >> x, y, z = symbols("x,y,z")
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
        if self == expr:
            return repl_dict

        d = self._matches_simple(expr, repl_dict)
        if d is not None:
            return d

        # eliminate exact part from pattern: (2+a+w1+w2).matches(expr) -> (w1+w2).matches(expr-a-2)
        wild_part = []
        exact_part = []
        from function import WildFunction
        from symbol import Wild
        for p in self.args:
            if p.has(Wild, WildFunction):
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

    @classmethod
    def make_args(cls, expr):
        """
        Return a sequence of elements `args` such that cls(*args) == expr

        >>> from sympy import Symbol, Mul, Add
        >>> x, y = map(Symbol, 'xy')

        >>> Mul.make_args(x*y)
        (x, y)
        >>> Add.make_args(x*y)
        (x*y,)
        >>> set(Add.make_args(x*y + y)) == set([y, x*y])
        True

        """
        if isinstance(expr, cls):
            return expr.args
        else:
            return (expr,)

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
        args = (sympify(arg) for arg in args)
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

    @classmethod
    def make_args(cls, expr):
        """
        Return a sequence of elements `args` such that cls(*args) == expr

        >>> from sympy import Symbol, Mul, Add
        >>> x, y = map(Symbol, 'xy')

        >>> Mul.make_args(x*y)
        (x, y)
        >>> Add.make_args(x*y)
        (x*y,)
        >>> set(Add.make_args(x*y + y)) == set([y, x*y])
        True

        """
        if isinstance(expr, cls):
            return expr._argset
        else:
            return frozenset([expr])

    @property
    def args(self):
        return tuple(self._argset)

    @staticmethod
    def _compare_pretty(a, b):
        return cmp(str(a), str(b))
