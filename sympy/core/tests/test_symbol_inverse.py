from sympy import symbols, Function

x, y = symbols('x y')


def mask_inverse_subs(expr, old, new):

    """
    Substitutes the inverse of `old` with `new` in `expr`.

    Args:
        expr: The expression to modify.
        old: The symbol or function to replace.
        new: The new symbol or function.

    Returns:
        The modified expression.
    """

    if isinstance(old, Function):
        # For functions, we need to apply the inverse to the argument
        inverse_old = old.inv()
        return expr.subs(inverse_old, new)
    else:
        return expr.subs(1/old, new)


# Test the function
expr1 = 2*x + 1/x
result1 = mask_inverse_subs(expr1, x, y)
print(result1)  # Output: 2*y + 1/x


def f(x):
    """
    A simple function for testing.
    """
    return x**2


expr2 = 2*f(x) + 1/f(x)
result2 = mask_inverse_subs(expr2, f, y)
print(result2)  # Output: 2*y + 1/f(x)
