# -*- coding: utf-8

from sympy.printing.mathml import DomPrinter
from sympy import sympify, S, Mul
from sympy.core.function import _coeff_isneg


class MathQuillPrinter(DomPrinter):

    def __init__(self, settings=None):
        DomPrinter.__init__(self, settings=settings)
        self.command_id_counter = 1
        self.block_id_counter = 1

    def get_mathquill_span(self, mqclass, op, counter="command_id"):
        span = self.dom.createElement("span")
        span.setAttribute("class", mqclass)  # "mq-binary-operator")
        if op is not None:
            span.appendChild(self.dom.createTextNode(op))
        if counter == "command_id":
            span.setAttribute("mathquill-command-id", "%i" % self.command_id_counter)
            self.command_id_counter += 1
        elif counter == "block_id":
            span.setAttribute("mathquill-block-id", "%i" % self.block_id_counter)
            self.block_id_counter += 1
        else:
            raise ValueError("unknown counter %s" % counter)
        return span

    def get_fragment(self):
        frag = self.dom.createDocumentFragment()

        def _frag_writexml(self, writer, indent="", addindent="", newl=""):
            for i in self.childNodes:
                i.writexml(writer, indent, addindent, newl)

        frag.__class__.writexml = _frag_writexml
        return frag


    def _print_Mul(self, expr):

        if _coeff_isneg(expr):
            x = self.get_fragment()
            x.appendChild(self.get_mathquill_span("", op="-"))
            x.appendChild(self._print_Mul(-expr))
            span = self.dom.createElement("span")
            span.setAttribute("style", "display:inline-block;width:0")
            # TODO: unescape
            span.appendChild(self.dom.createTextNode("​"))
            x.appendChild(span)
            return x

        from sympy.simplify import fraction
        numer, denom = fraction(expr)

        if denom is not S.One:
            x = self.get_mathquill_span("mq-fraction mq-non-leaf", None)
            mqNum = x.appendChild(self.get_mathquill_span("mq-numerator", op=None, counter="block_id"))
            mqNum.appendChild(self._print(numer))
            mqDen = x.appendChild(self.get_mathquill_span("mq-denominator", op=None, counter="block_id"))
            mqDen.appendChild(self._print(denom))
            return x

        coeff, terms = expr.as_coeff_mul()
        if coeff is S.One and len(terms) == 1:
            # XXX since the negative coefficient has been handled, I don't
            # think a coeff of 1 can remain
            return self._print(terms[0])

        if self.order != 'old':
            terms = Mul._from_args(terms).as_ordered_factors()

        frag = self.get_fragment()

        if(coeff != 1):
            frag.appendChild(self._print(coeff))
        for term in terms:
            frag.appendChild(self._print(term))
        return frag

    def _print_Add(self, expr, order=None):
        args = self._as_ordered_terms(expr, order=order)
        lastProcessed = self._print(args[0])
        x = self.get_fragment()
        if not _coeff_isneg(args[0]):
            x.appendChild(self._print(args[0]))
            args = args[1:]

        for arg in args:
            if _coeff_isneg(arg):
                # use minus
                x.appendChild(self.get_mathquill_span("mq-binary-operator", op='-', counter="command_id"))
            else:
                x.appendChild(self.get_mathquill_span("mq-binary-operator", op='+', counter="command_id"))
            x.appendChild(self._print(arg))
        return x

    def _print_Symbol(self, expr):
        v = self.dom.createElement("var")
        v.setAttribute("mathquill-command-id", str(self.command_id_counter))
        self.command_id_counter += 1
        v.appendChild(self.dom.createTextNode(expr.name))
        return v
