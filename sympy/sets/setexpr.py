from sympy.core import Expr, Function, Add, Mul, Pow, Dummy
from sympy import sift, latex, Min, Max, Set
from sympy.core.sets import imageset, Interval, FiniteSet, Union
from sympy.core.compatibility import u

x = Dummy('x')


class SetExpr(Expr):
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


def simplify(setexpr):
    """ Collapse expression containing SetExprs to single SetExpr """
    if not setexpr.has(SetExpr):
        return setexpr

    # Recurse downwards
    # e.g. simplify(x + y) -> simplify(simplify(x) + simplify(y))
    setexpr = setexpr.func(*map(simplify, setexpr.args))

    # Handle cases like exp(...) or sin(...)
    if isinstance(setexpr, Function):
        return SetExpr(imageset(type(setexpr), setexpr.args[0].set))

    elif isinstance(setexpr, (Add, Mul)):
        groups = sift(setexpr.args, lambda x: isinstance(x, SetExpr))
        setexprs, others = groups[True], groups[False]
        op = type(setexpr)
        se = join(op, setexprs)  # call out to many-setexpr function
        return SetExpr(imageset(x, op(x, *others), se.set))

    elif isinstance(setexpr, Pow) and isinstance(setexpr.base, SetExpr):
        return SetExpr(imageset(x, Pow(x, setexpr.exp), setexpr.base.set))

    return setexpr


def join(op, setexprs):
    ''' Join many setexprs into one

    Relies on many functions named join_foo below
    '''
    assert all(isinstance(se, SetExpr) for se in setexprs)

    for key, func in join_list:  # Multiple Dispatch
        if (issubclass(op, key[0])
            and len(key[1:]) == len(setexprs)
            and all(isinstance(se.set, func)
                    for se, func in zip(setexprs, key[1:]))):
            return func(op, *setexprs)

    # Two args is a common case for join_foo functions.  Lets reduce with join.
    if len(setexprs) > 2:
        result = setexprs[0]
        for se in setexprs[1:]:
            result = join(op, [result, se])
        return result

    return op(*setexprs)


def join_Add_Intervals(_, a, b):
    return SetExpr(Interval(a.set.start + b.set.start,
                            a.set.end + b.set.end,
                            a.set.left_open or b.set.left_open,
                            a.set.right_open or b.set.right_open))


def join_Mul_Intervals(_, a, b):
    bounds = {(a.set.start, a.set.end), (b.set.start, b.set.end)}
    start = Min(*[x * y for x in (a.set.start, a.set.end)
                        for y in (b.set.start, b.set.end)])
    end =   Max(*[x * y for x in (a.set.start, a.set.end)
                        for y in (b.set.start, b.set.end)])
    # TODO: Handle left_open right_open
    return SetExpr(Interval(start, end))


def join_Add_FiniteSet(op, a, b):
    if isinstance(b.set, FiniteSet):
        a, b = b, a
    return SetExpr(Union(*[simplify(op(x, b)).set for x in a.set]))


join_list = [[(Add, Interval, Interval),    join_Add_Intervals],
             [(Mul, Interval, Interval),    join_Mul_Intervals],
             [((Add, Mul), FiniteSet, Set), join_Add_FiniteSet],
             [((Add, Mul), Set, FiniteSet), join_Add_FiniteSet]]
