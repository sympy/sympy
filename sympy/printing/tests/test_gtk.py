from sympy import print_gtk, sin

def test1():
    from sympy.abc import x
    print_gtk(x**2, start_viewer=False)
    print_gtk(x**2+sin(x)/4, start_viewer=False)
