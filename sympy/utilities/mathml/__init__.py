"""Module with some functions for MathML, like transforming MathML
content in MathML presentation.

To use this module, you will need lxml.
"""

from sympy.utilities.pkgdata import get_resource
import xml.dom.minidom


def add_mathml_headers(s):
    return """<math xmlns:mml="http://www.w3.org/1998/Math/MathML"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:schemaLocation="http://www.w3.org/1998/Math/MathML
        http://www.w3.org/Math/XMLSchema/mathml2/mathml2.xsd">""" + s + "</math>"


def apply_xsl(mml, xsl):
    """Apply a xsl to a MathML string
    @param mml: a string with MathML code
    @param xsl: a string representing a path to a xsl (xml stylesheet)
        file. This file name is relative to the PYTHONPATH
    >>> from sympy.utilities.mathml import apply_xsl
    >>> xsl = 'mathml/data/simple_mmlctop.xsl'
    >>> mml = '<apply> <inverse/> <sin/> </apply>'
    >>> res = apply_xsl(mml,xsl)
    >>> '-1' in res.splitlines()[6]
    True
    
    """
    from lxml import etree
    s = etree.XML(get_resource(xsl).read())
    transform = etree.XSLT(s)
    doc = etree.XML(mml)
    result = transform(doc)
    s = str(result)
    return s


def c2p(mml, simple=False):
    """Transforms a document in MathML content (like the one that sympy produces)
    in one document in MathML presentation, more suitable for printing, and more
    widely accepted
    """

    if not mml.startswith('<math'):
        mml = add_mathml_headers(mml)

    if simple:
        return apply_xsl(mml, 'mathml/data/simple_mmlctop.xsl')

    return apply_xsl(mml, 'mathml/data/mmlctop.xsl')
