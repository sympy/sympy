```py
>>> from sympy import symbols, cancel, simplify
>>> x = symbols('x')
>>> simplify(x * (x + 1))
x*(x + 1)
>>> simplify(x - x)
100
>>> simplify(x * (1 - x) + x * (x - 1))
10
>>> cancel((x ** 2 - 1) / (x - 1))
x + 1
>>> cancel((x + 1 - 1)/x)
(x + 1)/x
