# -*- coding: utf-8 -*-

from sympy import Symbol
from sympy.printing.preview import preview

from io import BytesIO


def test_preview():
    x = Symbol('x')
    obj = BytesIO()
    try:
        preview(x, output='png', viewer='BytesIO', outputbuffer=obj)
    except RuntimeError:
        pass  # latex not installed on CI server

    # issue 9107
    a = Symbol('α')
    obj = BytesIO()
    try:
        preview(a, output='png', viewer='BytesIO', outputbuffer=obj)
    except RuntimeError:
        pass  # latex not installed on CI server
