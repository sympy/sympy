from sympy.physics.qoperations import QAssocOp
from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.basic import S
from sympy.physics.qmul import QMul
from sympy.physics.quantum import QExpr
from sympy.physics.qexpr import QuantumError
from sympy.printing.pretty.stringpict import prettyForm
from sympy.physics.quantum import StateBase, Operator, Dagger, Commutator,\
KroneckerDelta, InnerProduct
from sympy.physics.qpow import QPow

class QAdd(QAssocOp):
    """
        Quantum Add operation
    """
    binop = ' + '
    binopPretty =  prettyForm(u' \u002B ')

    @classmethod
    def _rules_QAdd(cls, object1, object2):
        """
            This method is called by new to instantiate a QAdd class
            Applies rules of what can and can't be added together by checking
            types and hilbert spaces
            Returns an Add or QAdd objects on success
            Raises and exception if input violates quantum shape rules
        """
        if object2 is S.Zero:
            return object1
        elif object1 is S.Zero:
            return object2

        if (not isinstance(object1, QExpr) or\
        issubclass(object1.acts_like, InnerProduct))\
        and\
        (not isinstance(object2, (QExpr) or\
        issubclass(object2.acts_like, InnerProduct))):
                return Add(object1, object2)
        elif (not isinstance(object1, (StateBase, Operator, QAssocOp, QPow, \
        Dagger, Commutator, KroneckerDelta)))\
        or (not isinstance(object2, (StateBase, Operator, QAssocOp, QPow, \
        Dagger, Commutator, KroneckerDelta))):
            raise QuantumError("Can't add a %s and %s"\
            % (object1.__class__.__name__, object2.__class__.__name__))

        if object1.hilbert_space != object2.hilbert_space:
            raise QuantumError("Hilbert Spaces do not match")

        if issubclass(object1.acts_like, object2.acts_like)\
        or issubclass(object2.acts_like, object1.acts_like):
            retVal = cls.QAddflatten([object1, object2])
            retVal.hilbert_space = object1.hilbert_space
            retVal.acts_like = object1.acts_like
            return retVal
        else:
            raise QuantumError("Can't add (%s + %s)"\
            % (object1.acts_like.__name__, object2.acts_like.__name__))

    @classmethod
    def QAddflatten(cls, seq):
        """
            This simplifies expressions by combining like terms
            Assumes associativity and commutivity for addition
        """
        terms = {}      # term -> coeff
                        # e.g. x**2 -> 5   for ... + 5*x**2 + ...

        coeff = S.Zero  # standalone term
                        # e.g. 3 + ...
        order_factors = []

        for o in seq:
            # O(x)
            if o.is_Order:
                for o1 in order_factors:
                    if o1.contains(o):
                        o = None
                        break
                if o is None:
                    continue
                order_factors = [o]+[o1 for o1 in order_factors if not\
                o.contains(o1)]
                continue

            # 3
            elif not isinstance(o, QExpr):
                coeff += o
                continue

            # Add([...])
            elif isinstance(o, QAdd):
                # NB: here we assume Add is always commutative
                seq.extend(o.args)  # TODO zerocopy?
                continue

            # Mul([...])
            elif isinstance(o, QMul):
                c = o.args[0]

                # 3*...
                if not isinstance(c, QExpr):
                    if c is S.One:
                        s = o
                    else:
                        s = o.as_two_terms()[1]

                else:
                    c = S.One
                    s = o

            # everything else
            else:
                c = S.One
                s = o


            # now we have:
            # o = c*s, where
            #
            # c is a Number
            # s is an expression with number factor extracted

            # let's collect terms with the same s, so e.g.
            # 2*x**2 + 3*x**2  ->  5*x**2
            if s in terms:
                terms[s] += c
            else:
                terms[s] = c


        # now let's construct new args:
        # [2*x**2, x**3, 7*x**4, pi, ...]
        newseq = []
        noncommutative = False
        for s,c in terms.items():
            # 0*s
            if c is S.Zero:
                continue
            # 1*s
            elif c is S.One:
                newseq.append(s)
            # c*s
            else:
                if isinstance(s, QMul):
                    # Mul, already keeps its arguments in perfect order.
                    # so we can simply put c in slot0 and go the fast way.
                    cs = QMul._new_rawargs(s.acts_like, s.hilbert_space, *((c,)\
                    + s.args)) #figure out rawargs F
                    newseq.append(cs)

                else:
                    # alternatively we have to call all Mul's machinery (slow)
                    newseq.append(QMul(c,s))

            noncommutative = noncommutative or not s.is_commutative

        # nan
        if coeff is S.NaN:
            raise QuantumError("NaN for some reason")

        # oo, -oo
        elif (coeff is S.Infinity) or (coeff is S.NegativeInfinity):
            newseq = [f for f in newseq if not f.is_real]


        # process O(x)
        if order_factors:
            newseq2 = []
            for t in newseq:
                for o in order_factors:
                    # x + O(x) -> O(x)
                    if o.contains(t):
                        t = None
                        break
                # x + O(x**2) -> x + O(x**2)
                if t is not None:
                    newseq2.append(t)
            newseq = newseq2 + order_factors
            # 1 + O(1) -> O(1)
            for o in order_factors:
                if o.contains(coeff):
                    coeff = S.Zero
                    break


        # order args canonically
        # Currently we sort things using hashes, as it is quite fast. A better
        # solution is not to sort things at all - but this needs some more
        # fixing.
        newseq.sort(key=hash)

        # current code expects coeff to be always in slot-0
        if coeff is not S.Zero:
            newseq.insert(0, coeff)

        # we are done
        return Expr.__new__(cls, *newseq)

    @property
    def identity(self):
        return S.Zero

    def _eval_dagger(self):
        newargs = []
        for item in self.args:
            newargs.append(Dagger(item))
        return QAdd(*newargs)

    def _eval_expand_mul(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_mul'):
                newterm = term._eval_expand_mul(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

