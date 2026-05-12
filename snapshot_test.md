# preamble section for setup
This section can be used as the preamble.
I setup the commonly used symbol x, and import functions I need.
```py
>>> from sympy import symbols, simplify, integrate, cancel
>>> x = symbols('x')
```

# testing the simplify function
These are basic tests of calling simplify with polynomial expressions.
The expected behaviour is not necessarily that it will call factor.
```py
>>> simplify(x * (x + 1))
x*(x + 1)
>>> simplify(x - x)
0
>>> simplify(x * (1 - x) + x * (x - 1))
0
```

# testing the cancel function
cancel function can be used for rational functions to remove the common factors
```py
>>> cancel((x ** 2 - 1)/(x - 1))
x + 1
>>> cancel((x + 1 - 1)/x)
1
```

# testing the integration capabilities of sympy
this performs indefinite integration wrt x
```py
>>> integrate(x)
x**2/2
>>> integrate(2 * x)
x**2
```


# testing basic exceptions in python
```py
>>> 1 / 0
11
```

# testing some other errors
```py
>>> 1 + "sympy"
11
```
