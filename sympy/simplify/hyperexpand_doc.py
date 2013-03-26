""" This module cooks up a docstring when imported. Its only purpose is to
    be displayed in the sphinx documentation. """

from sympy.simplify.hyperexpand import FormulaCollection
from sympy import latex, Eq, hyper

c = FormulaCollection()

doc = ""

for f in c.formulae:
    obj = Eq(hyper(f.indices.ap, f.indices.bq, f.z),
             f.closed_form.rewrite('nonrepsmall'))
    doc += ".. math::\n  %s\n" % latex(obj)

__doc__ = doc
