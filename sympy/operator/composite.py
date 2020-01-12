"""Operated operators."""
from .operator import OperatorExpr, Operator
from sympy.core.add import Add
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.singleton import S

class CompositeOperator(OperatorExpr):
    def __new__(cls, f, *gs, **kwargs):
        """
        gs는 ((0,g1,(2,)), (1,g2,(0,1)), ...)와 같은 형식을 가짐.
        첫 번째 원소는 g가 f의 몇 번째 argument가 되는가에 대한 정보.
        두 번째 원소는 f에 합성되는 g.
        세 번째 원소는 CompositeOperator가 call되었을 때 몇 번째 argument들을 g가 받는지에 대한 정보.
        """
        f, = cls._sep_basic_basicmeta(f)
        gs = cls._canonicalize_gs(*gs)
        obj = super(cls, CompositeOperator).__new__(cls, f, *gs, **kwargs)
        obj.f = obj.args[0]
        obj.gs = obj.args[1:]
        return obj

    @classmethod
    def _canonicalize_gs(cls, *gs):
        new_gs = []
        i = 0
        for g in gs:
            i += 1
            if not iterable(g):
                g, = cls._sep_basic_basicmeta(g)
                new_gs.append(Tuple(S(i-1), g, Tuple(S.Zero)))
            else:
                argidx1, g, argidx2 = g
                g, = cls._sep_basic_basicmeta(g)
                new_gs.append(Tuple(S(argidx1), g, Tuple(*argidx2)))

        new_gs.sort()
        return Tuple(*new_gs)

    def _eval_operation(self, *args, **kwargs):
        kwargs.pop('evaluate', None)
        f_args = []
        for idx1, g, idx2s in self.gs:
            g_args = [args[i] for i in idx2s]
            f_args.append(g(*g_args, evaluate=True, **kwargs))
        return self.f(*f_args, evaluate=True, **kwargs)
