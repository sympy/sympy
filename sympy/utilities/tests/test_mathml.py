from __future__ import annotations
import os
from textwrap import dedent
from sympy.external import import_module
from sympy.testing.pytest import skip, warns_deprecated_sympy
from sympy.utilities.mathml import add_mathml_headers, apply_xsl, c2p



lxml = import_module('lxml')

path = os.path.abspath(os.path.join(os.path.dirname(__file__), "test_xxe.py"))


def test_xxe():
    assert os.path.isfile(path)
    if not lxml:
        skip("lxml not installed.")

    mml = dedent(
        rf"""
        <!--?xml version="1.0" ?-->
        <!DOCTYPE replace [<!ENTITY ent SYSTEM "file://{path}"> ]>
        <userInfo>
        <firstName>John</firstName>
        <lastName>&ent;</lastName>
        </userInfo>
        """
    )
    xsl = 'mathml/data/simple_mmlctop.xsl'

    with warns_deprecated_sympy():
        res = apply_xsl(mml, xsl)
    assert res == \
        '<?xml version="1.0"?>\n<userInfo>\n<firstName>John</firstName>\n<lastName/>\n</userInfo>\n'


def test_c2p():
    if not lxml:
        skip("lxml not installed.")

    mml = '<apply> <plus/> <ci>a</ci> <ci>b</ci> </apply>'
    with warns_deprecated_sympy():
        res = c2p(mml, simple=True)
    assert '<mrow>\n    <mi>a</mi>\n    <mo> + </mo>\n    <mi>b</mi>\n  </mrow>' in res
    with warns_deprecated_sympy():
        res = c2p(mml)
    assert '<mi>a</mi>' in res and '<mi>b</mi>' in res


def test_add_mathml_headers():
    with warns_deprecated_sympy():
        res = add_mathml_headers('<ci>a</ci>')
    assert res.startswith('<math xmlns:mml="http://www.w3.org/1998/Math/MathML"')
    assert res.endswith('><ci>a</ci></math>')
