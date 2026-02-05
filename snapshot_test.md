# preamble section for setup
```py
>>> from sympy import symbols, simplify, integrate, cancel
>>> x = symbols('x')
```

# testing the simplify function
```py
>>> simplify(x * (x + 1))
x*(x + 1)
>>> simplify(x - x)
100
>>> simplify(x * (1 - x) + x * (x - 1))
67
```

# testing the cancel function
```py
>>> cancel((x ** 2 - 1)/(x - 1))
x + 1
>>> cancel((x + 1 - 1)/x)
x
```

# testing the integration capabilities of sympy
```py
>>> integrate(x)
x**2/2
>>> integrate(2 * x)
x**3
```
