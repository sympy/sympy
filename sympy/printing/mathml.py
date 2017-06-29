"""
A MathML printer.
"""

from __future__ import print_function, division

from sympy import sympify, S, Mul
from sympy.core.function import _coeff_isneg
from sympy.core.alphabets import greeks
from sympy.core.compatibility import range
from .printer import Printer
from .pretty.pretty_symbology import greek_unicode
from .conventions import split_super_sub, requires_partial


class MathMLPrinter(Printer):
    """Prints an expression to the MathML markup language

    Whenever possible tries to use Content markup and not Presentation markup.

    References: https://www.w3.org/TR/MathML3/
    """
    printmethod = "_mathml"
    _default_settings = {
        "order": None,
        "encoding": "utf-8"
    }

    def __init__(self, settings=None):
        Printer.__init__(self, settings)
        from xml.dom.minidom import Document
        self.dom = Document()

    def doprint(self, expr):
        """
        Prints the expression as MathML.
        """
        mathML = Printer._print(self, expr)
        unistr = mathML.toxml()
        xmlbstr = unistr.encode('ascii', 'xmlcharrefreplace')
        res = xmlbstr.decode()
        return res

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
            'Sum': 'sum',
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
            'log': 'ln',
            'Equality': 'eq',
            'Unequality': 'neq',
            'GreaterThan': 'geq',
            'LessThan': 'leq',
            'StrictGreaterThan': 'gt',
            'StrictLessThan': 'lt',
        }

        for cls in e.__class__.__mro__:
            n = cls.__name__
            if n in translate:
                return translate[n]
        # Not found in the MRO set
        n = e.__class__.__name__
        return n.lower()

    def _print_Mul(self, expr):

        if _coeff_isneg(expr):
            x = self.dom.createElement('apply')
            x.appendChild(self.dom.createElement('minus'))
            x.appendChild(self._print_Mul(-expr))
            return x

        from sympy.simplify import fraction
        numer, denom = fraction(expr)

        if denom is not S.One:
            x = self.dom.createElement('apply')
            x.appendChild(self.dom.createElement('divide'))
            x.appendChild(self._print(numer))
            x.appendChild(self._print(denom))
            return x

        coeff, terms = expr.as_coeff_mul()
        if coeff is S.One and len(terms) == 1:
            # XXX since the negative coefficient has been handled, I don't
            # think a coeff of 1 can remain
            return self._print(terms[0])

        if self.order != 'old':
            terms = Mul._from_args(terms).as_ordered_factors()

        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement('times'))
        if(coeff != 1):
            x.appendChild(self._print(coeff))
        for term in terms:
            x.appendChild(self._print(term))
        return x

    def _print_Add(self, expr, order=None):
        args = self._as_ordered_terms(expr, order=order)
        lastProcessed = self._print(args[0])
        plusNodes = []
        for arg in args[1:]:
            if _coeff_isneg(arg):
                # use minus
                x = self.dom.createElement('apply')
                x.appendChild(self.dom.createElement('minus'))
                x.appendChild(lastProcessed)
                x.appendChild(self._print(-arg))
                # invert expression since this is now minused
                lastProcessed = x
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

    def _print_MatrixBase(self, m):
        x = self.dom.createElement('matrix')
        for i in range(m.rows):
            x_r = self.dom.createElement('matrixrow')
            for j in range(m.cols):
                x_r.appendChild(self._print(m[i, j]))
            x.appendChild(x_r)
        return x

    def _print_Rational(self, e):
        if e.q == 1:
            # don't divide
            x = self.dom.createElement('cn')
            x.appendChild(self.dom.createTextNode(str(e.p)))
            return x
        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement('divide'))
        # numerator
        xnum = self.dom.createElement('cn')
        xnum.appendChild(self.dom.createTextNode(str(e.p)))
        # denominator
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

    def _print_ImaginaryUnit(self, e):
        return self.dom.createElement('imaginaryi')

    def _print_EulerGamma(self, e):
        return self.dom.createElement('eulergamma')

    def _print_GoldenRatio(self, e):
        """We use unicode #x3c6 for Greek letter phi as defined here
        http://www.w3.org/2003/entities/2007doc/isogrk1.html"""
        x = self.dom.createElement('cn')
        x.appendChild(self.dom.createTextNode(u"\N{GREEK SMALL LETTER PHI}"))
        return x

    def _print_Exp1(self, e):
        return self.dom.createElement('exponentiale')

    def _print_Pi(self, e):
        return self.dom.createElement('pi')

    def _print_Infinity(self, e):
        return self.dom.createElement('infinity')

    def _print_Negative_Infinity(self, e):
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

            if len(limits[0]) == 3:
                low_elem = self.dom.createElement('lowlimit')
                low_elem.appendChild(self._print(limits[0][1]))
                x.appendChild(low_elem)
                up_elem = self.dom.createElement('uplimit')
                up_elem.appendChild(self._print(limits[0][2]))
                x.appendChild(up_elem)
            if len(limits[0]) == 2:
                up_elem = self.dom.createElement('uplimit')
                up_elem.appendChild(self._print(limits[0][1]))
                x.appendChild(up_elem)
            if len(limits) == 1:
                x.appendChild(self._print(e.function))
            else:
                x.appendChild(lime_recur(limits[1:]))
            return x

        limits = list(e.limits)
        limits.reverse()
        return lime_recur(limits)

    def _print_Sum(self, e):
        # Printer can be shared because Sum and Integral have the
        # same internal representation.
        return self._print_Integral(e)

    def _print_Symbol(self, sym):
        ci = self.dom.createElement(self.mathml_tag(sym))

        def join(items):
            if len(items) > 1:
                mrow = self.dom.createElement('mml:mrow')
                for i, item in enumerate(items):
                    if i > 0:
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

        # translate name, supers and subs to unicode characters
        greek_letters = set(greeks) # make a copy
        def translate(s):
            if s in greek_unicode:
                return greek_unicode.get(s)
            else:
                return s

        name, supers, subs = split_super_sub(sym.name)
        name = translate(name)
        supers = [translate(sup) for sup in supers]
        subs = [translate(sub) for sub in subs]

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
        # Here we use root instead of power if the exponent is the reciprocal of an integer
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
        diff_symbol = self.mathml_tag(e)
        if requires_partial(e):
            diff_symbol = 'partialdiff'
        x.appendChild(self.dom.createElement(diff_symbol))

        x_1 = self.dom.createElement('bvar')
        for sym in e.variables:
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

    def _print_Relational(self, e):
        x = self.dom.createElement('apply')
        x.appendChild(self.dom.createElement(self.mathml_tag(e)))
        x.appendChild(self._print(e.lhs))
        x.appendChild(self._print(e.rhs))
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

    def apply_patch(self):
        # Applying the patch of xml.dom.minidom bug
        # Date: 2011-11-18
        # Description: http://ronrothman.com/public/leftbraned/xml-dom-minidom-\
        #                   toprettyxml-and-silly-whitespace/#best-solution
        # Issue: http://bugs.python.org/issue4147
        # Patch: http://hg.python.org/cpython/rev/7262f8f276ff/

        from xml.dom.minidom import Element, Text, Node, _write_data

        def writexml(self, writer, indent="", addindent="", newl=""):
            # indent = current indentation
            # addindent = indentation to add to higher levels
            # newl = newline string
            writer.write(indent + "<" + self.tagName)

            attrs = self._get_attributes()
            a_names = list(attrs.keys())
            a_names.sort()

            for a_name in a_names:
                writer.write(" %s=\"" % a_name)
                _write_data(writer, attrs[a_name].value)
                writer.write("\"")
            if self.childNodes:
                writer.write(">")
                if (len(self.childNodes) == 1 and
                        self.childNodes[0].nodeType == Node.TEXT_NODE):
                    self.childNodes[0].writexml(writer, '', '', '')
                else:
                    writer.write(newl)
                    for node in self.childNodes:
                        node.writexml(
                            writer, indent + addindent, addindent, newl)
                    writer.write(indent)
                writer.write("</%s>%s" % (self.tagName, newl))
            else:
                writer.write("/>%s" % (newl))
        self._Element_writexml_old = Element.writexml
        Element.writexml = writexml

        def writexml(self, writer, indent="", addindent="", newl=""):
            _write_data(writer, "%s%s%s" % (indent, self.data, newl))
        self._Text_writexml_old = Text.writexml
        Text.writexml = writexml

    def restore_patch(self):
        from xml.dom.minidom import Element, Text
        Element.writexml = self._Element_writexml_old
        Text.writexml = self._Text_writexml_old


def mathml(expr, **settings):
    """Returns the MathML representation of expr"""
    return MathMLPrinter(settings).doprint(expr)


def print_mathml(expr, **settings):
    """
    Prints a pretty representation of the MathML code for expr

    Examples
    ========

    >>> ##
    >>> from sympy.printing.mathml import print_mathml
    >>> from sympy.abc import x
    >>> print_mathml(x+1) #doctest: +NORMALIZE_WHITESPACE
    <apply>
        <plus/>
        <ci>x</ci>
        <cn>1</cn>
    </apply>

    """
    s = MathMLPrinter(settings)
    xml = s._print(sympify(expr))
    s.apply_patch()
    pretty_xml = xml.toprettyxml()
    s.restore_patch()

    print(pretty_xml)
