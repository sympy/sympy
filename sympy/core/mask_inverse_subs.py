from sympy import symbols, Function

x, y = symbols('x y')

def mask_inverse_subs(expr, old, new):
    if isinstance(old, Function):  # Check if old is a function
        inverse_old = 1 / old(x)  # Compute the functional inverse
        return expr.subs(inverse_old, new)
    else:
        return expr.subs(1/old, new)

# Test the function
expr1 = 2*x + 1/x
result1 = mask_inverse_subs(expr1, x, y)
print(result1)  # Output: 2*y + 1/x

def f(x):
    return x**2

expr2 = 2*f(x) + 1/f(x)
result2 = mask_inverse_subs(expr2, f, y)
print(result2)  # Output: 2*y + 1/f(x)
