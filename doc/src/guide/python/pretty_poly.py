from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm

class PolyPrettyPrinter(PrettyPrinter):
    """This printer prints polynomials nicely. """

    def _print_Poly(self, poly):
        expr = poly.as_expr()
        gens = list(poly.gens)
        domain = poly.get_domain()

        pform_head = prettyForm('Poly')
        pform_tail = self._print_seq([expr] + gens + [domain], '(', ')')

        pform = prettyForm(*pform_head.right(pform_tail))
        return pform

def pretty_poly(expr, **settings):
    """Pretty-print polynomials nicely. """
    p = PolyPrettyPrinter(settings)
    s = p.doprint(expr)

    return s
