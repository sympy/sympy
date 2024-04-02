from sympy import symbols, Or, And, Not

# Define symbols
x, y, z = symbols('x y z')

# Custom XOR function with negations handling
def custom_xor(expr1, expr2):
    return Or(And(expr1, Not(expr2)), And(Not(expr1), expr2))

# Example usage
expr1 = x ^ y
expr2 = ~(x ^ z)
result = custom_xor(expr1, expr2)