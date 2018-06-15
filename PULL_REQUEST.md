# Allow substitution of function symbols within an integral
<!-- Your title above should be a short description of what
was changed. Do not include the issue number in the title. -->

#### References to other Issues or PRs
<!-- If this pull request fixes an issue, write "Fixes #NNNN" in that exact
format, e.g. "Fixes #1234". See
https://github.com/blog/1506-closing-issues-via-pull-requests .-->
Fixes #14796

#### Brief description of what is fixed or changed
In the code for substitutions within an expression involving limits (such as an
integral), the code breaks if you try to substitute one function symbol for another.
The reason for this is found in the file `sympy/concrete/expr_with_limits.py` in
the implementation of the method `_eval_subs(self, old, new)` of the class `ExprWithLimits`. In
that method, the properties and methods `.free_symbols`, `.args`, `.atoms` are consulted.
However, those properties and methods do not return list (or iterables) in the case
where the `old` or `new` arguments are instances of `FunctionClass`.

A fix that works in the case that sparked the issue is to check if those properties
and methods return lists, and if not, proceed as if they returned empty lists. This
perhaps should be implemented in the `FunctionClass`.

#### Other comments
The `subs` function allows substitutions of subexpressions within an expression.
So for example:

```
from sympy import *
x, y = symbols("x y")
f, g = symbols("f g", cls=Function)
expr = f(x)
expr1 = expr.subs(x, y)
```

Afterwards, `expr1` will be the expression `f(y)`. You can even substitute function
symbols for other function symbols:

```
expr2 = expr.subs(f, cos)
```

Afterwards, `expr2` will be the expression `cos(x)`. However, if you try to do
the same substitution within an integral, `Integral(f(x), (x,0,1))`, (the integral
from x=0 to x=1 of f(x)) then instead of getting `Integral(cos(x), (x,0,1)`, the
code raises an error. This is because in substituting `cos` for `f`, it checks
the free variables of `f` and the args of `f` and the atoms of `cos`, and those
methods and properties are not defined on FunctionClass objects.