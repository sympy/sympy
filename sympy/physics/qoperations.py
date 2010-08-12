from sympy.core.expr import Expr
from sympy.core.decorators import call_highest_priority
from sympy.core.sympify import sympify
from sympy.physics.quantum import Operator, KetBase, BraBase, OuterProduct, InnerProduct, StateBase
from sympy.core.numbers import Number
from sympy import Symbol
from sympy.core.mul import Mul
from sympy.printing.str import sstr
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.physics.quantumbasic import QuantumError, QuantumBasic
from sympy.printing.pretty.stringpict import prettyForm

class QAssocOp(QuantumBasic):
    """
        QAssocOp is Quantum version of Sympy's AssocOp
        This is the Base class of QAdd and QMul
    """

    def __new__(cls, *args, **assumptions):
        if len(args) == 1:
            return sympify(args[0])
        if len(args) == 0:
            return sympify(1)
        #try to instantiate an object of cls type
        return cls.instantiate(map(sympify, args))

    @classmethod
    def instantiate(cls, seq):
        """
            Apply rules to each element in the input sequence
            The rules check to see if the input is valid and then instantiate the objects if they can
        """
        rules = getattr(cls, '_rules_%s' % cls.__name__)
        result = seq[0]
        for i in range(len(seq)-1):
            result = rules(result, seq[i+1])
        return result

    def __getitem__(self, number):
        """
            This retrieves arguments from self.args
        """
        return self.args[number]

    def _sympystr(self, printer, *args):
        from sympy.physics.qmul import QMul
        from sympy.physics.qadd import QAdd
        length = len(self.args)
        string = ''
        for i in range(length):
            if isinstance(self.args[i], (Add, Pow, QAdd)):
                string = string + '('
            string = string + sstr(self.args[i])
            if isinstance(self.args[i], (Add, Pow, QAdd)):
                string = string + ')'
            if i != length-1:
                string = string + self.__class__.binop
        return string

    def _pretty(self, printer, *args):
        from sympy.physics.qmul import QMul
        from sympy.physics.qadd import QAdd
        length = len(self.args)
        pform = printer._print('', *args)
        for i in range(length):
            if isinstance(self.args[i], (Add, Pow, QAdd)):
                pform = prettyForm(*pform.right(prettyForm(u'\u0028')))
            if hasattr(self.args[i], '_pretty'):
                pform = prettyForm(*pform.right(self.args[i]._pretty(printer, *args)))
            else:
                pform = prettyForm(*pform.right(sstr(self.args[i])))
            if isinstance(self.args[i], (Add, Pow, QAdd)):
                pform = prettyForm(*pform.right(prettyForm(u'\u0029')))
            if i != length-1:
                pform = prettyForm(*pform.right(self.binopPretty))
        return pform


    def _new_rawargs(self, evaluates, hilbert_space, *args):
        """ Create new instance of own class with args exactly as provided by caller;
            Sets evaluates and hilbert_space attributes to those specified in inputs

            This is handy when we want to optimize things, e.g.
        """
        if len(args) == 1:
            return args[0]
        obj = Expr.__new__(type(self), *args)  # NB no assumptions for Add/Mul
        obj.evaluates = evaluates
        obj.hilbert_space = hilbert_space
        return obj

