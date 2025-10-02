from sympy import *

x = symbols("x")
k = symbols("k")
y = Function("y")

# Test case 1: Your original equation
ode = Eq(x**k*y(x) + Derivative(y(x), (x, 2)), 0)
sol = dsolve(ode, func=y(x))
print("Solution for y'' + x^k*y = 0:")
print(sol)
# Test case 2: With a coefficient
ode2 = Eq(2*x**k*y(x) + Derivative(y(x), (x, 2)), 0)
sol2 = dsolve(ode2, func=y(x))
print("\nSolution for y'' + 2*x^k*y = 0:")
print(sol2)

# Test case 3: Numeric k (should still work with original method)
ode3 = Eq(x**2*y(x) + Derivative(y(x), (x, 2)), 0)
sol3 = dsolve(ode3, func=y(x))
print("\nSolution for y'' + x^2*y = 0:")
print(sol3)

