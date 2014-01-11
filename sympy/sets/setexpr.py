from sympy.core import Expr, Function, Add, Mul, Pow, Dummy
from sympy import sift, latex, Min, Max, Set
from sympy.core.sets import imageset, Interval, FiniteSet, Union
from sympy.core.compatibility import u

from sympy.core.decorators import call_highest_priority, _sympifyit


x = Dummy('x')


def setexpr(inp):
    if isinstance(inp, set):
        return SetExpr(FiniteSet(inp))
    if isinstance(inp, tuple) and len(inp) == 2:
        return SetExpr(Interval(inp[0], inp[1], True, True))
    if isinstance(inp, list) and len(inp) == 2:
        return SetExpr(Interval(inp[0], inp[1], False, False))
    return SetExpr(inp)


class SetExpr(Expr):
    """ An expression that can take on values of a set

    >>> from sympy import Interval, FiniteSet
    >>> from sympy.sets.setexpr import SetExpr, simplify

    >>> a = SetExpr(Interval(0, 5))
    >>> b = SetExpr(FiniteSet(1, 10))
    >>> simplify(a + b).set
    [1, 6] U [10, 15]

    >>> simplify(2*a + b).set
    [1, 20]
    """
    _op_priority = 11.0
    set = property(lambda self: self.args[0])

    def _latex(self, printer):
        return printer._print(self.set)

    def _pretty(self, printer):
        return printer._print(self.set)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return simplify(Add(self, other))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return simplify(Add(self, other))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return simplify(Mul(self, other))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return simplify(Mul(other, self))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return simplify(Add(self, -other))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return simplify(Add(other, -self))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return simplify(Pow(self, other))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        return simplify(Pow(other, self))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return simplify(Mul(self, 1/other))

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        return simplify(Mul(other, Pow(self, -1)))


    def _eval_func(self, func):
        return SetExpr(imageset(func, self.set))
    _eval_exp = _eval_log = _eval_sin = _eval_cos = _eval_tan =_eval_func


def simplify(inp):
    """ Collapse expression containing SetExprs to single SetExpr """
    if not inp.has(SetExpr):
        return inp

    # Recurse downwards
    # e.g. simplify(x + y) -> simplify(simplify(x) + simplify(y))
    inp = inp.func(*map(simplify, inp.args))

    op, args = inp.func, inp.args

    # If operation is commutative and we have non-SetExprs then
    # 1. Simplify op(*SetExprs)
    # 2. Turn all non-SetExprs into a function, f = x -> op(x, *non-set-exprs)
    # 3. Apply and simplify, return simplify(f(op(*SetExprs)))
    if (isinstance(inp, (Add, Mul))
            and not all(isinstance(arg, SetExpr) for arg in inp.args)):
        groups = sift(inp.args, lambda x: isinstance(x, SetExpr))
        setexprs, others = groups[True], groups[False]
        se = simplify(op(*setexprs))  # call out to many-setexpr function
        return SetExpr(imageset(x, op(x, *others), se.set))

    # Multiple dispatch to `_simplify_foo` functions defined below
    # Dispatched on operator, and then type of argument or type of set if
    # argument is SetExpr
    args2 = [arg.set if isinstance(arg, SetExpr) else arg for arg in args]
    for key, func in join_list:  # Multiple Dispatch
        if (issubclass(op, key[0])
            and len(key[1:]) == len(args2)
            and all(isinstance(se, typ)
                    for se, typ in zip(args2, key[1:]))):
            return func(op, *args)

    # Two args is a common case for _simplify_foo functions.
    # Lets try a reduction with join.
    if len(args) > 2:
        result = args[0]
        for se in args[1:]:
            result = simplify(op(result, se))
        return result

    return op(*args)


def _simplify_Function(func, a):
    """ f(SetExpr(set)) -> SetExpr(imageset(f, set)) """
    return SetExpr(imageset(func, a.set))


def _simplify_Pow(_, base, exp):
    """ SetExpr(set)**expr -> SetExpr(imageset(x, x**expr, set)) """
    return SetExpr(imageset(x, Pow(x, exp), base.set))


def _simplify_Add_Intervals(_, a, b):
    """ SetExpr([x, y]) + SetExpr([m, n]) -> SetExpr([x + m, y + n]) """
    return SetExpr(Interval(a.set.start + b.set.start,
                            a.set.end + b.set.end,
                            a.set.left_open or b.set.left_open,
                            a.set.right_open or b.set.right_open))


def _simplify_Mul_Intervals(_, a, b):
    """ SetExpr([x, y]) * SetExpr([m, n]) -> SetExpr([...]) """
    start = Min(*[x * y for x in (a.set.start, a.set.end)
                        for y in (b.set.start, b.set.end)])
    end =   Max(*[x * y for x in (a.set.start, a.set.end)
                        for y in (b.set.start, b.set.end)])
    # TODO: Handle left_open right_open
    return SetExpr(Interval(start, end))


def _simplify_FiniteSet(op, a, b):
    if isinstance(b.set, FiniteSet):
        a, b = b, a
    return SetExpr(Union(*[simplify(op(x, b)).set for x in a.set]))


join_list = [[(Function, Set),              _simplify_Function],
             [(Pow, Set, Expr),             _simplify_Pow],
             [(Add, Interval, Interval),    _simplify_Add_Intervals],
             [(Mul, Interval, Interval),    _simplify_Mul_Intervals],
             [((Add, Mul), FiniteSet, Set), _simplify_FiniteSet],
             [((Add, Mul), Set, FiniteSet), _simplify_FiniteSet]]
