from sympy import Boolean, Not, Or, And

class CustomBoolean(Boolean):
    def __xor__(self, other):
        # Check if other is a negated expression
        if isinstance(other, Not):
            return Or(And(self, other.args[0]), And(Not(self), other.args[0]))
        else:
            return super().__xor__(other)

# Define symbols
x, y, z = symbols('x y z')

# Use custom subclass for boolean expressions
x_expr = CustomBoolean(x)
y_expr = CustomBoolean(y)
z_expr = CustomBoolean(z)

# Example usage
expr1 = x_expr ^ y_expr
expr2 = ~(x_expr ^ z_expr)
result = expr1 ^ expr2

print(result)
