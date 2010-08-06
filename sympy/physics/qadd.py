from sympy.physics.qoperations import QAssocOp
from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.basic import Basic, S, C
from sympy.physics.qmul import QMul
from sympy.physics.quantum import Dagger
from sympy.physics.quantum import QuantumBasic
from sympy.physics.quantumbasic import QuantumError

class QAdd(QAssocOp):
    binop = ' + '

    @classmethod    
    def _rules_QAdd(cls, Object1, Object2):
        from sympy.physics.quantum import StateBase, Operator, OuterProduct, KetBase, BraBase, OuterProduct, InnerProduct
        from sympy.core.add import Add    
        if (not isinstance(Object1, QuantumBasic)) and (not isinstance(Object2, QuantumBasic)):
            return Add(Object1, Object2)
        elif (not isinstance(Object1, QuantumBasic)) or (not isinstance(Object2, QuantumBasic)):
            raise Exception("Can't add a %s and %s" % (Object1.__class__.__name__, Object2.__class__.__name__))

        if Object1.hilbert_space != Object2.hilbert_space:
            raise Exception("Hilbert Spaces do not match")  

        if Object1.evaluates == Object2.evaluates:
            retVal = cls.QAddflatten([Object1, Object2])
            retVal.hilbert_space = Object1.hilbert_space
            retVal.evaluates = Object1.evaluates
            return retVal
        else:
            raise Exception("Can't add (%s + %s)" % (Object1.evaluates.__name__, Object2.evaluates.__name__))

    @classmethod
    def QAddflatten(cls, seq):
        # get Non-QUantum stuff out
        from sympy.physics.quantum import StateBase, Operator, OuterProduct, KetBase, BraBase, OuterProduct, InnerProduct
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
                order_factors = [o]+[o1 for o1 in order_factors if not o.contains(o1)]
                continue

            # 3
            elif not isinstance(o, QuantumBasic):
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
                if not isinstance(c, QuantumBasic):
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
                if s.is_Mul:
                    # Mul, already keeps its arguments in perfect order.
                    # so we can simply put c in slot0 and go the fast way.
                    cs = s._new_rawargs(*((c,) + s.args))
                    newseq.append(cs)

                else:
                    # alternatively we have to call all Mul's machinery (slow)
                    newseq.append(QMul(c,s))

            noncommutative = noncommutative or not s.is_commutative

        # nan
        if coeff is S.NaN:
            raise Exception("NaN for some reason")
            
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
        if len(newseq) == 1:
            return newseq[0]
        return Expr.__new__(cls, *newseq)

    @property
    def identity(self):
        return S.Zero

    def _eval_dagger(self):
        newargs = []
        for item in self.args:
            newargs.append(Dagger(item))
        return QAdd(*newargs)
