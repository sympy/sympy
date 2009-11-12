"""
A MathML printer.
"""

from sympy import Basic, sympify, C, S
from sympy.simplify import fraction
from printer import Printer
from conventions import split_super_sub


class MathMLPrinter(Printer):
    """Prints an expression to the MathML markup language

    Whenever possible tries to use Content markup and not Presentation markup.

    References: http://www.w3.org/TR/MathML2/
    """
    printmethod = "_mathml_"

    def __init__(self, *args, **kwargs):
        Printer.__init__(self, *args, **kwargs)
        from xml.dom.minidom import Document
        self.dom = Document()

    def doprint(self, e):
        mathML = Printer.doprint(self,expr)
        return mathML.toxml()

    def mathml_tag(self, e):
        """Returns the MathML tag for an expression."""
        translate = {
            'Add': 'plus',
            'Mul': 'times',
            'Derivative': 'diff',
            'Number': 'cn',
            'int': 'cn',
            'Pow': 'power',
            'Symbol': 'ci',
            'Integral': 'int',
            'sin': 'sin',
            'cos': 'cos',
            'tan': 'tan',
            'cot': 'cot',
            'asin': 'arcsin',
            'asinh': 'arcsinh',
            'acos': 'arccos',
            'acosh': 'arccosh',
            'atan': 'arctan',
            'atanh': 'arctanh',
            'acot': 'arccot',
            'atan2': 'arctan',
            'log': 'ln'
        }

        for cls in e.__class__.__mro__:
            n = cls.__name__
            if n in translate:
                return translate[n]
        # Not found in the MRO set
        n = e.__class__.__name__
        return n.lower()

    def _print_Mul(self, expr):
        coeff, terms  = expr.as_coeff_terms()

        if coeff.is_negative:
            x = self.dom.createElement('apply')
            x.appendChild(self.dom.createElement('minus'))
            x.appendChild(self._print_Mul(-expr))
            return x

        numer, denom = fraction(expr)

        if not denom is S.One:
            x = self.dom.createElement('apply')
            x.appendChild(self.dom.createElement('divide'))
            x.appendChild(self._print(numer))
            x.appendChild(self._print(denom))
            return x

        if coeff == 1 and len(terms) == 1:
            return self._print(terms[0])
        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement('times'))
        if(coeff != 1):
            x.appendChild(self._print(coeff))
        for term in terms:
            x.appendChild(self._print(term))
        return x

    # This is complicated because we attempt to order then results in order of Basic._compare_pretty
    # and use minus instead of negative
    def _print_Add(self, e):
        args = list(e.args)
        args.sort(Basic._compare_pretty)
        lastProcessed = self._print(args[0])
        args.pop(0)
        plusNodes = list()
        for i in range(0,len(args)):
            arg = args[i]
            coeff, terms = arg.as_coeff_terms()
            if(coeff.is_negative):
                #use minus
                x = self.dom.createElement('apply')
                x.appendChild(self.dom.createElement('minus'))
                x.appendChild(lastProcessed)
                x.appendChild(self._print(-arg))
                #invert expression  since this is now minused
                lastProcessed = x;
                if(arg == args[-1]):
                    plusNodes.append(lastProcessed)
            else:
                plusNodes.append(lastProcessed)
                lastProcessed = self._print(arg)
                if(arg == args[-1]):
                    plusNodes.append(self._print(arg))
        if len(plusNodes) == 1:
            return lastProcessed
        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement('plus'))
        while len(plusNodes) > 0:
            x.appendChild(plusNodes.pop(0))
        return x

    def _print_Matrix(self, m):
        x = self.dom.createElement('matrix')
        for i in range(m.lines):
            x_r = self.dom.createElement('matrixrow')
            for j in range(m.cols):
                x_r.appendChild(self._print(m[i,j]))
            x.appendChild(x_r)
        return x

    def _print_Rational(self, e):
        if e.q == 1:
            #don't divide
            x = self.dom.createElement('cn')
            x.appendChild(self.dom.createTextNode(str(e.p)))
            return x
        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement('divide'))
        #numerator
        xnum = self.dom.createElement('cn')
        xnum.appendChild(self.dom.createTextNode(str(e.p)))
        #denomenator
        xdenom = self.dom.createElement('cn')
        xdenom.appendChild(self.dom.createTextNode(str(e.q)))
        x.appendChild(xnum)
        x.appendChild(xdenom)
        return x

    def _print_Limit(self, e):
        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement(self.mathml_tag(e)))

        x_1 = self.dom.createElement('bvar')
        x_2 = self.dom.createElement('lowlimit')
        x_1.appendChild(self._print(e.args[1]))
        x_2.appendChild(self._print(e.args[2]))

        x.appendChild(x_1)
        x.appendChild(x_2)
        x.appendChild(self._print(e.args[0]))
        return x

    def _print_ImaginaryUnit(self,e):
        return self.dom.createElement('imaginaryi')

    def _print_EulerGamma(self,e):
        return self.dom.createElement('eulergamma')

    def _print_GoldenRatio(self,e):
        """We use unicode #x3c6 for greek letter phi as defined here
        http://www.w3.org/Math/characters/"""
        x = self.dom.createElement('cn')
        x.appendChild(self.dom.createTextNode(u"\u03c6"))
        return x

    def _print_Exp1(self,e):
        return self.dom.createElement('exponentiale')

    def _print_Pi(self, e):
        return self.dom.createElement('pi')

    def _print_Infinity(self, e):
        return self.dom.createElement('infinity')

    def _print_Negative_Infinity(self,e):
        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement('minus'))
        x.appendChild(self.dom.createElement('infinity'))
        return x

    def _print_Integral(self, e):
        def lime_recur(limits):
            x = self.dom.createElement('apply')
            x.appendChild(self.dom.createElement(self.mathml_tag(e)))
            bvar_elem = self.dom.createElement('bvar')
            bvar_elem.appendChild(self._print(limits[0][0]))
            x.appendChild(bvar_elem)

            if limits[0][1]:
                low_elem = self.dom.createElement('lowlimit')
                low_elem.appendChild(self._print(limits[0][1][0]))
                x.appendChild(low_elem)
                up_elem = self.dom.createElement('uplimit')
                up_elem.appendChild(self._print(limits[0][1][1]))
                x.appendChild(up_elem)
            if len(limits) == 1:
                x.appendChild(self._print(e.function))
            else:
                x.appendChild(lime_recur(limits[1:]))
            return x

        limits = list(e.limits)
        limits.reverse()
        return lime_recur(limits)

    def _print_Symbol(self, sym):
        ci = self.dom.createElement(self.mathml_tag(sym))

        def join(items):
            if len(items) > 1:
                mrow = self.dom.createElement('mml:mrow')
                for i, item in enumerate(items):
                    if i>0:
                        mo = self.dom.createElement('mml:mo')
                        mo.appendChild(self.dom.createTextNode(" "))
                        mrow.appendChild(mo)
                    mi = self.dom.createElement('mml:mi')
                    mi.appendChild(self.dom.createTextNode(item))
                    mrow.appendChild(mi)
                return mrow
            else:
                mi = self.dom.createElement('mml:mi')
                mi.appendChild(self.dom.createTextNode(items[0]))
                return mi

        name, supers, subs = split_super_sub(sym.name)
        mname = self.dom.createElement('mml:mi')
        mname.appendChild(self.dom.createTextNode(name))
        if len(supers) == 0:
            if len(subs) == 0:
                ci.appendChild(self.dom.createTextNode(name))
            else:
                msub = self.dom.createElement('mml:msub')
                msub.appendChild(mname)
                msub.appendChild(join(subs))
                ci.appendChild(msub)
        else:
            if len(subs) == 0:
                msup = self.dom.createElement('mml:msup')
                msup.appendChild(mname)
                msup.appendChild(join(supers))
                ci.appendChild(msup)
            else:
                msubsup = self.dom.createElement('mml:msubsup')
                msubsup.appendChild(mname)
                msubsup.appendChild(join(subs))
                msubsup.appendChild(join(supers))
                ci.appendChild(msubsup)
        return ci

    def _print_Pow(self, e):
        #Here we use root instead of power if the exponent is the reciprocal of an integer
        if e.exp.is_Rational and e.exp.p == 1:
            x = self.dom.createElement('apply')
            x.appendChild(self.dom.createElement('root'))
            if e.exp.q != 2:
                xmldeg = self.dom.createElement('degree')
                xmlci = self.dom.createElement('ci')
                xmlci.appendChild(self.dom.createTextNode(str(e.exp.q)))
                xmldeg.appendChild(xmlci)
                x.appendChild(xmldeg)
            x.appendChild(self._print(e.base))
            return x

        x = self.dom.createElement('apply')
        x_1 = self.dom.createElement(self.mathml_tag(e))
        x.appendChild(x_1)
        x.appendChild(self._print(e.base))
        x.appendChild(self._print(e.exp))
        return x

    def _print_Number(self, e):
        x = self.dom.createElement(self.mathml_tag(e))
        x.appendChild(self.dom.createTextNode(str(e)))
        return x

    def _print_Derivative(self, e):
        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement(self.mathml_tag(e)))

        x_1 = self.dom.createElement('bvar')
        for sym in e.symbols:
            x_1.appendChild(self._print(sym))

        x.appendChild(x_1)
        x.appendChild(self._print(e.expr))
        return x

    def _print_Function(self, e):
        x = self.dom.createElement("apply")
        x.appendChild(self.dom.createElement(self.mathml_tag(e)))
        for arg in e.args:
            x.appendChild(self._print(arg))
        return x

    def _print_Basic(self, e):
        x = self.dom.createElement(self.mathml_tag(e))
        for arg in e:
            x.appendChild(self._print(arg))
        return x

    def _print_AssocOp(self, e):
        x = self.dom.createElement('apply')
        x_1 = self.dom.createElement(self.mathml_tag(e))
        x.appendChild(x_1)
        for arg in e.args:
            x.appendChild(self._print(arg))
        return x

    def _print_list(self, seq):
        """MathML reference for the <list> element:
        http://www.w3.org/TR/MathML2/chapter4.html#contm.list"""
        dom_element = self.dom.createElement('list')
        for item in seq:
            dom_element.appendChild(self._print(item))
        return dom_element

    def _print_int(self, p):
        dom_element = self.dom.createElement(self.mathml_tag(p))
        dom_element.appendChild(self.dom.createTextNode(str(p)))
        return dom_element

def mathml(expr):
    """Returns the MathML representation of expr"""
    s = MathMLPrinter()
    return s._print(sympify(expr)).toxml(encoding="utf-8")

def print_mathml(expr):
    """
    Print's a pretty representation of the MathML code for expr

    >>> ##
    >>> from sympy.printing.mathml import print_mathml
    >>> from sympy.abc import x
    >>> print_mathml(x+1) #doctest: +NORMALIZE_WHITESPACE
    <apply>
        <plus/>
        <cn>
                1
        </cn>
        <ci>
                x
        </ci>
    </apply>

    """
    s = MathMLPrinter()
    print s._print(sympify(expr)).toprettyxml(encoding="utf-8")
