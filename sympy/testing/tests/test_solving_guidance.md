# preamble
```
>>> from sympy import *
>>> from sympy.abc import x
```

# test_no_closed_form_solution
```
>>> from sympy import solve, cos
>>> from sympy.abc import x
>>> solve(cos(x) - x, x, dict=True)
Traceback (most recent call last):
    ...
NotImplementedError: multiple generators [x, cos(x)]
No algorithms are implemented to solve equation -x + cos(x)
```

# test_nsolve
```
>>> from sympy import nsolve, cos
>>> from sympy.abc import x
>>> nsolve(cos(x) - x, x, 2)
0.739085133215161
```

# test_polynomial_roots
```
>>> from sympy import solve
>>> from sympy.abc import x
>>> solutions = solve(x**5 - x - 1, x, dict=True)
>>> [solution[x].evalf(3) for solution in solutions]
[1.17, -0.765 - 0.352*I, -0.765 + 0.352*I, 0.181 - 1.08*I, 0.181 + 1.08*I]
```
