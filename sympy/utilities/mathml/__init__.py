"""Module with some functions for MathML, like transforming MathML
content in MathML presentation.

To use this module, you will need lxml.
"""
from __future__ import annotations

from pathlib import Path

from sympy.utilities.decorator import deprecated


def _add_mathml_headers(s):
    return """<math xmlns:mml="http://www.w3.org/1998/Math/MathML"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:schemaLocation="http://www.w3.org/1998/Math/MathML
        http://www.w3.org/Math/XMLSchema/mathml2/mathml2.xsd">""" + s + "</math>"


@deprecated("""\
sympy.utilities.mathml.add_mathml_headers() is deprecated. Wrap the MathML
string in a <math> element directly instead, like

    '<math xmlns="http://www.w3.org/1998/Math/MathML">' + s + '</math>'""",
    deprecated_since_version="1.15",
    active_deprecations_target="deprecated-mathml-xsl")
def add_mathml_headers(s):
    """Wrap a MathML string in a ``<math>`` element.

    .. deprecated:: 1.15

       ``add_mathml_headers`` is deprecated. Wrap the MathML string in a
       ``<math>`` element directly instead. See
       :ref:`deprecated-mathml-xsl` for details.

    """
    return _add_mathml_headers(s)


def _read_binary(pkgname, filename):
    import sys

    if sys.version_info >= (3, 10):
        # files was added in Python 3.9 but only seems to work here in 3.10+
        from importlib.resources import files
        return files(pkgname).joinpath(filename).read_bytes()
    else:
        # read_binary was deprecated in Python 3.11
        from importlib.resources import read_binary
        return read_binary(pkgname, filename)


def _read_xsl(xsl):
    # Previously these values were allowed:
    if xsl == 'mathml/data/simple_mmlctop.xsl':
        xsl = 'simple_mmlctop.xsl'
    elif xsl == 'mathml/data/mmlctop.xsl':
        xsl = 'mmlctop.xsl'
    elif xsl == 'mathml/data/mmltex.xsl':
        xsl = 'mmltex.xsl'

    if xsl in ['simple_mmlctop.xsl', 'mmlctop.xsl', 'mmltex.xsl']:
        xslbytes = _read_binary('sympy.utilities.mathml.data', xsl)
    else:
        xslbytes = Path(xsl).read_bytes()

    return xslbytes


def _apply_xsl(mml, xsl):
    from lxml import etree

    parser = etree.XMLParser(resolve_entities=False)
    ac = etree.XSLTAccessControl.DENY_ALL

    s = etree.XML(_read_xsl(xsl), parser=parser)
    transform = etree.XSLT(s, access_control=ac)
    doc = etree.XML(mml, parser=parser)
    result = transform(doc)
    s = str(result)
    return s


@deprecated("""\
sympy.utilities.mathml.apply_xsl() is deprecated. Use lxml.etree.XSLT
directly to apply an XSLT stylesheet to a MathML string.""",
    deprecated_since_version="1.15",
    active_deprecations_target="deprecated-mathml-xsl")
def apply_xsl(mml, xsl):
    """Apply a xsl to a MathML string.

    .. deprecated:: 1.15

       ``apply_xsl`` is deprecated. Use ``lxml.etree.XSLT`` directly to apply
       an XSLT stylesheet to a MathML string. See
       :ref:`deprecated-mathml-xsl` for details.

    Parameters
    ==========

    mml
        A string with MathML code.
    xsl
        A string giving the name of an xsl (xml stylesheet) file which can be
        found in sympy/utilities/mathml/data. The following files are supplied
        with SymPy:

        - mmlctop.xsl
        - mmltex.xsl
        - simple_mmlctop.xsl

        Alternatively, a full path to an xsl file can be given.

    Examples
    ========

    >>> from sympy.utilities.mathml import apply_xsl
    >>> xsl = 'simple_mmlctop.xsl'
    >>> mml = '<apply> <plus/> <ci>a</ci> <ci>b</ci> </apply>'
    >>> res = apply_xsl(mml,xsl) # doctest: +SKIP
    >>> print(res) # doctest: +SKIP
    <?xml version="1.0"?>
    <mrow xmlns="http://www.w3.org/1998/Math/MathML">
      <mi>a</mi>
      <mo> + </mo>
      <mi>b</mi>
    </mrow>
    """
    return _apply_xsl(mml, xsl)


def _c2p(mml, simple=False):
    if not mml.startswith('<math'):
        mml = _add_mathml_headers(mml)

    if simple:
        return _apply_xsl(mml, 'mathml/data/simple_mmlctop.xsl')

    return _apply_xsl(mml, 'mathml/data/mmlctop.xsl')


@deprecated("""\
sympy.utilities.mathml.c2p() is deprecated. To get presentation MathML for
a SymPy expression, replace c2p(mathml(expr)) with

    mathml(expr, printer='presentation')

To convert other content MathML, apply an XSLT stylesheet with
lxml.etree.XSLT directly.""",
    deprecated_since_version="1.15",
    active_deprecations_target="deprecated-mathml-xsl")
def c2p(mml, simple=False):
    """Transforms a document in MathML content (like the one that sympy produces)
    in one document in MathML presentation, more suitable for printing, and more
    widely accepted

    .. deprecated:: 1.15

       ``c2p`` is deprecated. To get presentation MathML for a SymPy
       expression, use ``mathml(expr, printer='presentation')`` (see
       :func:`~sympy.printing.mathml.mathml`). See
       :ref:`deprecated-mathml-xsl` for details.

    Examples
    ========

    >>> from sympy.utilities.mathml import c2p
    >>> mml = '<apply> <exp/> <cn>2</cn> </apply>'
    >>> c2p(mml,simple=True) != c2p(mml,simple=False) # doctest: +SKIP
    True

    """
    return _c2p(mml, simple=simple)
