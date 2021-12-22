import os
from lxml import etree
from sympy.parsing.latex import parse_latex
from sympy import *
from sympy.printing.pycode import pycode

# this function will take mathml input 
# for eg. 
# mathml_input = """<math xmlns="http://www.w3.org/1998/Math/MathML">
#   <mrow>
#     <mfrac>
#       <mrow><mi>x</mi></mrow>
#       <mrow><mi>y</mi></mrow>
#     </mfrac>
#     <mrow><mi>+</mi></mrow>
#     <mrow><mi>z</mi></mrow>
#   </mrow>
# </math>"""
# and returns the python equation as output i.e (x/y)+z

def mathml2latex(equation):
    dom = etree.fromstring(equation)
    xslt_root = etree.parse("../../utilities/mathml/data/mmltex.xsl")
    transform = etree.XSLT(xslt_root)
    tex = str(transform(dom))
    tex = tex.replace('$', '')
    expr = parse_latex(tex)
    new_expr = pycode(expr)
    print(new_expr)






