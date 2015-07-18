from sympy import Symbol
from sympy.printing.preview import preview

from io import BytesIO


def test_preview():
    x = Symbol('x')
    obj = BytesIO()
    preview(x, output='png', viewer='BytesIO', outputbuffer=obj)
