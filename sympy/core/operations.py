from sympy.core.core import C
from sympy.core.expr import Expr
from sympy.core.sympify import _sympify, sympify
from sympy.core.basic import Basic
from sympy.core.cache import cacheit
from sympy.core.compatibility import cmp, quick_sort
from sympy.core.logic import fuzzy_and

# from add import Add /cyclic/
# from mul import Mul /cyclic/
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
    def __new__(cls, *args, **options):
        args = map(_sympify, args)
        args = [a for a in args if a is not cls.identity]

        if not options.pop('evaluate', True):
            return cls._from_args(args)

        if len(args) == 0:
            return cls.identity
        if len(args) == 1:
            return args[0]

        c_part, nc_part, order_symbols = cls.flatten(args)
        is_commutative = not nc_part
        obj = cls._from_args(c_part + nc_part, is_commutative)

        if order_symbols is not None:
            return C.Order(obj, *order_symbols)
        return obj

    @classmethod
    def _from_args(cls, args, is_commutative=None):
        """Create new instance with already-processed args"""
        if len(args) == 0:
            return cls.identity
        elif len(args) == 1:
            return args[0]

        obj = Expr.__new__(cls, *args)
        if is_commutative is None:
            is_commutative = fuzzy_and(a.is_commutative for a in args)
        obj.is_commutative = is_commutative
        return obj

    def _new_rawargs(self, *args, **kwargs):
        """Create new instance of own class with args exactly as provided by
        caller but returning the self class identity if args is empty.

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
           self is non-commutative and kwarg `reeval=False` has not been
           passed.
        """
        if kwargs.pop('reeval', True) and self.is_commutative is False:
            is_commutative = None
        else:
            is_commutative = self.is_commutative
        return self._from_args(args, is_commutative)

    @classmethod
    def flatten(cls, seq):
        """Return seq so that none of the elements are of type `cls`. This is
        the vanilla routine that will be used if a class derived from AssocOp
        does not define its own flatten routine."""
        # apply associativity, no commutativity property is used
        new_seq = []
        while seq:
            o = seq.pop()
            if o.__class__ is cls: # classes must match exactly
                seq.extend(o.args)
            else:
                new_seq.append(o)
        # c_part, nc_part, order_symbols
        return [], new_seq, None

    def _matches_commutative(self, expr, repl_dict={}):
        """
        Matches Add/Mul "pattern" to an expression "expr".

        repl_dict ... a dictionary of (wild: expression) pairs, that get
                      returned with the results

        This function is the main workhorse for Add/Mul.

        For instance:

        >>> from sympy import symbols, Wild, sin
        >>> a = Wild("a")
        >>> b = Wild("b")
        >>> c = Wild("c")
        >>> x, y, z = symbols("x y z")
        >>> (a+sin(b)*c)._matches_commutative(x+sin(y)*z)
        {a_: x, b_: y, c_: z}

        In the example above, "a+sin(b)*c" is the pattern, and "x+sin(y)*z" is
        the expression.

        The repl_dict contains parts that were already matched. For example
        here:

        >>> (x+sin(b)*c)._matches_commutative(x+sin(y)*z, repl_dict={a: x})
        {a_: x, b_: y, c_: z}

        the only function of the repl_dict is to return it in the
        result, e.g. if you omit it:

        >>> (x+sin(b)*c)._matches_commutative(x+sin(y)*z)
        {b_: y, c_: z}

        the "a: x" is not returned in the result, but otherwise it is
        equivalent.

        """
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
            if p.has(Wild, WildFunction) and (not expr.has(p)):
                # not all Wild should stay Wilds, for example:
                # (w2+w3).matches(w1) -> (w1+w3).matches(w1) -> w3.matches(0)
                wild_part.append(p)
            else:
                exact_part.append(p)

        if exact_part:
            newpattern = self.func(*wild_part)
            newexpr = self._combine_inverse(expr, self.func(*exact_part))
            return newpattern.matches(newexpr, repl_dict)

        # now to real work ;)
        if expr.is_Add:
            i, d = expr.as_independent(C.Symbol)
            expr_list = (i,) + self.make_args(expr)
        else:
            expr_list = self.make_args(expr)
        for last_op in reversed(expr_list):
            for w in reversed(wild_part):
                d1 = w.matches(last_op, repl_dict)
                if d1 is not None:
                    d2 = self.xreplace(d1).matches(expr, d1)
                    if d2 is not None:
                        return d2

        return

    def _has_matcher(self):
        """Helper for .has()"""
        def _ncsplit(expr):
            # this is not the same as args_cnc because here
            # we don't assume expr is a Mul -- hence deal with args --
            # and always return a set.
            cpart, ncpart = [], []
            for arg in expr.args:
                if arg.is_commutative:
                    cpart.append(arg)
                else:
                    ncpart.append(arg)
            return set(cpart), ncpart

        c, nc = _ncsplit(self)
        cls = self.__class__
        def is_in(expr):
            if expr == self:
                return True
            elif not isinstance(expr, Basic):
                return False
            elif isinstance(expr, cls):
                _c, _nc = _ncsplit(expr)
                if (c & _c) == c:
                    if not nc:
                        return True
                    elif len(nc) <= len(_nc):
                        for i in xrange(len(_nc) - len(nc)):
                            if _nc[i:i+len(nc)] == nc:
                                return True
            return False
        return is_in

    def _eval_template_is_attr(self, is_attr, when_multiple=False):
        # return True if all elements have the property;
        # False if one doesn't have the property; and
        # if more than one doesn't have property, return
        #    False if when_multiple = False
        #    None if when_multiple is not False
        quick = when_multiple is None
        multi = False
        for t in self.args:
            a = getattr(t, is_attr)
            if a is True:
                continue
            if a is None:
                return
            if quick and multi:
                return None
            multi = True
        return not multi

    def _eval_evalf(self, prec):
        return self.func(*[s._evalf(prec) for s in self.args])

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

    def __new__(cls, *args, **options):
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
            obj = Expr.__new__(cls, _args)
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
    @cacheit
    def args(self):
        return quick_sort(self._argset)

    @staticmethod
    def _compare_pretty(a, b):
        return cmp(str(a), str(b))
