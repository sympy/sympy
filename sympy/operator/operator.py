"""Provides the basic classes for operators."""

from sympy.core.basic import Atom
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.core.function import FunctionClass
from sympy.core.expr import Expr
from sympy.core.singleton import S
from sympy.core.sympify import sympify

class OperatorExpr(Expr):
    """Basic class for every operator classes.
    """

    _op_priority = 20 # This is some random big value

    def __call__(self, *args, evaluate=True, **kwargs):
        args = sympify(args)
        return AppliedOperator(self, args, evaluate=evaluate, **kwargs)

    def _eval_operation(self, *args, **kwargs):
        """Every subclasses of OperatorExpr should override this.
        """
        raise NotImplementedError()

    def __pos__(self):
        return self

    def __neg__(self):
        from .operatorop import OperatorMul
        return OperatorMul(-1, self)

    def __add__(self, other):
        from .operatorop import OperatorAdd
        return OperatorAdd(self, other)

    def __radd__(self, other):
        from .operatorop import OperatorAdd
        return OperatorAdd(other, self)

    def __sub__(self, other):
        from .operatorop import OperatorAdd
        return OperatorAdd(self, -other)

    def __rsub__(self, other):
        from .operatorop import OperatorAdd
        return OperatorAdd(-self, other)

    def __mul__(self, other):
        from .operatorop import OperatorMul
        return OperatorMul(self, other)

    def __rmul__(self, other):
        from .operatorop import OperatorMul
        return OperatorMul(other, self)

    def __pow__(self, other):
        from .operatorop import OperatorPow
        return OperatorPow(self, other)

    def __rpow__(self, other):
        from .operatorop import OperatorPow
        return OperatorPow(other, self)

    def __div__(self, other):
        from .operatorop import OperatorMul, OperatorPow
        return OperatorMul(self, OperatorPow(other, S.NegativeOne))

    def __rdiv__(self, other):
        from .operatorop import OperatorMul, OperatorPow
        return OperatorMul(other, OperatorPow(self, S.NegativeOne))

    def __matmul__(self, other):
        from .composite import CompositeOperator
        return CompositeOperator(self, other)

    def __rmatmul__(self, other):
        from .composite import CompositeOperator
        return CompositeOperator(other, self)

    def deriv(self, *argidxs, **kwargs):
        from .derivative import DerivatedOperator
        return DerivatedOperator(self, *argidxs, **kwargs)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__


    @classmethod
    def _process_operator(cls, *args):
        """Unnest the Operator and wrap the arguments with Operator.
        """
        args = [sympify(a) for a in args]
        newargs = []
        for a in args:
            if isinstance(a, FunctionClass):
                newargs.append(Operator(a))
            elif isinstance(a, Expr):
                if isinstance(a, Operator):
                    if isinstance(a.operator, FunctionClass):
                        newargs.append(a)
                    elif isinstance(a.operator, Expr):
                        newargs.append(a.operator)
                else:
                    newargs.append(a)
            else:
                newargs.append(a)
        return newargs

class Operator(OperatorExpr, Atom):
    def __new__(cls, operator, argidxs=None, **kwargs):
        """
        Atomic operator.

        Explanation
        ===========

        This class wraps a single ``FunctionClass`` or ``Expr`` instance.
        Operator wrapping ``Expr`` instance is considered as a constant function.
        Calling it with arguments will always return the original ``Expr`` instance.
        Operator wrapping ``FunctionClass`` instance is considered as a operator.
        Calling it with arguments return result of calling original ``FunctionClass``
        instance.
        One may choose which arguments among the passed will be used with
        ``argidxs`` parameter.

        Examples
        ========

        >>> from sympy import sin, cos
        >>> from sympy.operator import Operator
        >>> from sympy.abc import x,y,z

        >>> Operator(1)
        1
        >>> Operator(sin)
        sin

        >>> Operator(1)(x,y)
        1
        >>> Operator(cos)(x,y,z)
        cos(x)

        >>> Operator(sin, (1,))(x,y,z)
        sin(y)

        Basic operation between operators are supported.

        >>> (Operator(sin) + Operator(cos))(x)
        sin(x) + cos(x)
        >>> (sin+cos)(x)    # This can be used instead of the example above.
        sin(x) + cos(x)

        >>> (Operator(sin) + Operator(cos, (1,)))(x,y)
        sin(x) + cos(y)


        Parameters
        ==========

        operator : Expr or FunctionClass instance

        argidxs : tuple of numbers or tuples, optional
                number represents the index, and tuple represent the slicing of
                the arguments that will be actually passed.
                If not given, automatially set according to the arity of operator.
        """
        if isinstance(operator, OperatorExpr):
            return operator

        operator = sympify(operator)
        if argidxs is None and isinstance(operator, FunctionClass):
            nargs = operator.nargs
            if nargs is S.Naturals0:
                argidxs = Tuple(Tuple(0,))
            else:
                argidxs = Tuple(Tuple(0, max(nargs)))

        argidxs = cls._canonicalize_argidxs(argidxs)

        obj = super(cls, Operator).__new__(cls, operator, argidxs, **kwargs)
        obj.operator = operator
        obj.argidxs = argidxs
        return obj

    @classmethod
    def _canonicalize_argidxs(cls, argidxs):
        if argidxs is None:
            argidxs = Tuple(Tuple(0,))

        new_argidxs = []
        for a in argidxs:
            if iterable(a):
                if len(a) == 2 and a[1] == a[0]+1:
                    new_argidxs.append(sympify(a[0]))
                else:
                    new_argidxs.append(Tuple(*a))
            else:
                new_argidxs.append(sympify(a))
        result = Tuple(*new_argidxs)
        return result

    def _eval_operation(self, *args, **kwargs):
        if isinstance(self.operator, Expr):
            result = self.operator
        elif isinstance(self.operator, FunctionClass):
            picked_args = []
            for a in self.argidxs:
                if not iterable(a):
                    picked_args.append(args[a])
                else:
                    if len(a) == 1:
                        i, = a
                        picked_args.extend(args[i:])
                    if len(a) == 2:
                        i,j = a
                        picked_args.extend(args[i:j])
            result = self.operator(*picked_args, **kwargs)
        return result

class AppliedOperator(OperatorExpr):
    """
    Applied, but unevaluated operator.

    Explanation
    ===========

    If ``operator`` parameter is combined operator, instance of this class is
    returned. Else, the evaluated result is always returned.

    Examples
    ========

    >>> from sympy import sin, cos, S
    >>> from sympy.operator import Operator, AppliedOperator
    >>> from sympy.abc import x,y,z

    >>> Operator(1)(x, evaluate=False)
    1
    >>> _ is S.One
    True

    >>> Operator(sin)(x, evaluate=False)
    sin(x)
    >>> isinstance(_, sin)
    True

    >>> (1+sin)(x, evaluate=False)
    (1 + sin)(x)
    >>> isinstance(_, AppliedOperator)
    True

    Parameters
    ==========

    operator : operator instance

    args : tuple of arguments
    """
    def __new__(cls, operator, args, evaluate=False, **kwargs):
        if not isinstance(operator, OperatorExpr):
            raise TypeError("%s must be instance of OperatorExpr." % operator)
        if isinstance(operator, Operator):
            if isinstance(operator.operator, FunctionClass):
                return operator._eval_operation(*args, evaluate=evaluate, **kwargs)
            elif isinstance(operator.operator, Expr):
                return operator.operator

        if evaluate:
            obj = operator._eval_operation(*args, **kwargs)
        else:
            obj = super(cls, AppliedOperator).__new__(cls, operator, args, **kwargs)
            obj.operator = obj.args[0]
            obj.arguments = obj.args[1]

        return obj

    def doit(self, **kwargs):
        return self.func(*self.args, evaluate=True)
