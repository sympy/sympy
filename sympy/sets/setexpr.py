from sympy.core import Expr, Function, Add, Mul, Pow, Dummy
from sympy import sift, latex, Min, Max, Set
from sympy.core.sets import imageset, Interval, FiniteSet, Union
from sympy.core.compatibility import u

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
    set = property(lambda self: self.args[0])

    def _latex(self, printer):
        return r"x \mid x \in " + printer._print(self.set)

    def _pretty(self, printer):
        if printer._use_unicode:
            inn = u("\u220a")
        else:
            inn = 'in'
        bar = printer._print('|')
        xx = printer._print(x)
        set = printer._print(self.set)

        return printer._print_seq((x, bar, xx, inn, set), '', '', ' ')


def simplify(inp):
    """ Collapse expression containing SetExprs to single SetExpr """
    if not inp.has(SetExpr):
        return inp

    # Recurse downwards
    # e.g. simplify(x + y) -> simplify(simplify(x) + simplify(y))
    inp = inp.func(*map(simplify, inp.args))

    # Handle cases like exp(...) or sin(...)
    op, args = inp.func, inp.args

    # If operation is commutative and we have non-SetExprs then
    if (isinstance(inp, (Add, Mul))
            and not all(isinstance(arg, SetExpr) for arg in inp.args)):
        groups = sift(inp.args, lambda x: isinstance(x, SetExpr))
        setexprs, others = groups[True], groups[False]
        se = simplify(op(*setexprs))  # call out to many-setexpr function
        return SetExpr(imageset(x, op(x, *others), se.set))

    args2 = [arg.set if isinstance(arg, SetExpr) else arg for arg in args]
    for key, func in join_list:  # Multiple Dispatch
        if (issubclass(op, key[0])
            and len(key[1:]) == len(args2)
            and all(isinstance(se, typ)
                    for se, typ in zip(args2, key[1:]))):
            return func(op, *args)

    # Two args is a common case for _simplify_foo functions.  Lets reduce with join.
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
