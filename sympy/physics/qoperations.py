from sympy.core.expr import Expr
from sympy.core.sympify import sympify
from sympy.printing.str import sstr
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.physics.quantumbasic import QuantumBasic
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
            The rules check to see if the input is valid and then instantiate
            the objects if they can
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
        from sympy.physics.qadd import QAdd
        length = len(self.args)
        pform = printer._print('', *args)
        for i in range(length):
            if isinstance(self.args[i], (Add, Pow, QAdd)):
                pform = prettyForm(*pform.right(prettyForm(u'\u0028')))
            if hasattr(self.args[i], '_pretty'):
                pform = prettyForm(*pform.right(self.args[i]._pretty(printer,\
                *args)))
            else:
                pform = prettyForm(*pform.right(printer._print(self.args[i],\
                *args)))
            if isinstance(self.args[i], (Add, Pow, QAdd)):
                pform = prettyForm(*pform.right(prettyForm(u'\u0029')))
            if i != length-1:
                pform = prettyForm(*pform.right(self.binopPretty))
        return pform


    @classmethod
    def _new_rawargs(cls, acts_like, hilbert_space, *args):
        """ Create new instance of own class with args exactly as provided by
            caller; Sets acts_like and hilbert_space attributes to those
            specified in inputs

            This is handy when we want to optimize things This is because we do
            not have to call the special '_qrules_*' to create a new object           
        """
        
        if len(args) == 1:
            return args[0]
        obj = Expr.__new__(cls, *args)  # NB no assumptions for Add/Mul
        obj.acts_like = acts_like
        obj.hilbert_space = hilbert_space
        return obj

