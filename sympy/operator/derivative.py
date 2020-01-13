"""Provides derivative operator, i.e. f'(x)."""

from .operator import OperatorExpr, Operator
from sympy.core.compatibility import iterable
from sympy.core.symbol import Dummy, symbols

class DerivatedOperator(OperatorExpr):
    """
    Derivated operator.

    Examples
    ========

    >>> from sympy import sin
    >>> from sympy.operator import DerivatedOperator
    >>> from sympy.abc import x

    >>> DerivatedOperator(sin)
    sin'
    >>> DerivatedOperator(sin)(x, evaluate=False)
    (sin')(x)
    >>> _.doit()
    cos(x)
    >>> DerivatedOperator(sin, 2)(1, evaluate=False)
    (sin'')(1)
    >>> _.doit()
    -sin(1)

    Parameters
    ==========
    operator : instance of Expr or FunctionClass

    argidxs : tuple of tuples of two numbers
            First number indicates the index of derivating argument, and second
            indicates how many times the operator is derivated with respect to
            the argument.
    """
    def __new__(cls, operator, *argidxs, **kwargs):
        operator = Operator(operator)

        if not argidxs:
            argidxs = ((0, 1),)
        elif len(argidxs) == 1 and not iterable(argidxs[0]):
            argidxs = ((0, argidxs[0]),)
        elif len(argidxs) > 1:
            new_argidxs = []
            pass_step = False
            for i,a in enumerate(argidxs):
                if iterable(a):
                    new_argidxs.append(a)
                elif not iterable(a) and not pass_step:
                    new_argidxs.append((argidxs[i], argidxs[i+1]))
                    pass_step = True
                else:
                    pass_step = False
                    continue
            argidxs = tuple(new_argidxs)


        argidxs_count = {}
        for (argidx,n) in argidxs:
            if argidx not in argidxs_count:
                argidxs_count[argidx] = n
            else:
                argidxs_count[argidx] += n
        argidxs_count = sorted(list(argidxs_count.items()))

        obj = super(cls, DerivatedOperator).__new__(cls, operator, *argidxs_count, **kwargs)
        obj.operator = operator
        obj.argidxs = argidxs_count

        return obj

    def _eval_operation(self, *args, **kwargs):
        dummies = symbols('xi_:%s' % len(args), cls=Dummy)
        operated = self.operator(*dummies)
        for argidx,n in self.argidxs:
            dummy_s, arg = dummies[argidx], args[argidx]
            operated = operated.diff(dummy_s, n)
            operated = operated.subs(dummy_s, arg)
        return operated
