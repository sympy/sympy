from sympy.core.expr import Expr
from sympy.core.sympify import sympify
from sympy.printing.str import sstr
from sympy.core.add import Add
from sympy.core.power import Pow
from sympy.physics.qexpr import QExpr
from sympy.printing.pretty.stringpict import prettyForm

class QAssocOp(QExpr):
    """QAssocOp is Quantum version of Sympy's AssocOp.

    This is the base class of QAdd and QMul, the two associated binary
    operations for quantum expressions.
    """

    def __new__(cls, *args, **assumptions):
        if len(args) == 1:
            return sympify(args[0])
        if len(args) == 0:
            return sympify(1)
        return cls.eval(map(sympify, args))

    @classmethod
    def eval(cls, args):
        """Apply rules to each element in args.

        The rules check to see if the input is valid and then instantiate the
        objects if they can. The application of rules is where the
        ``acts_like`` and ``hilbert_space`` slots are computed. The
        application of rule also validates the quantum expression.
        """
        result = args[0]
        for i in range(len(args)-1):
            result = cls._apply_rules(result, args[i+1])
        return result

    @classmethod
    def _new_rawargs(cls, acts_like, hilbert_space, *args):
        """Create new instance of own class with args exactly as provided.

        This class method sets the ``acts_like`` and ``hilbert_space`` slots
        to those specified by the caller, rather than computing them
        dynamically.

        This is handy when we want to optimize things by not having to do
        a careful application of the rules and dynamically computing the
         ``acts_like`` and ``hilbert_space`` slots. If we
        know what it evaluates to and to what hilbert_space it will belong,
        one can hardcode that in. Thus, It's a quicker way of instantiating a
        QAssocOp If we know what type we will get. This was stolen from
        sympy.core.operations.

            >>> from sympy.physics.hilbert import HilbertSpace
            >>> from sympy.physics.quantum import Ket, Operator
            >>> from sympy.physics.qmul import QMul
            >>> QMul._new_rawargs(Ket, HilbertSpace(), Operator('psi'), Ket('a'))
            psi*|a>
            >>> a = _
            >>> a.acts_like
            <class 'sympy.physics.quantum.Ket'>
            >>> a.hilbert_space
            HilbertSpace()
        """

        if len(args) == 1:
            return args[0]
        obj = Expr.__new__(cls, *args)
        obj.acts_like = acts_like
        obj.hilbert_space = hilbert_space
        return obj

    def __getitem__(self, i):
        """This retrieves arguments from self.args"""
        return self.args[i]

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
