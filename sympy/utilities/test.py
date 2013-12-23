from sympy import symbols
import lambdify
x = symbols('x')
#lambdify._import("math") # This is working fine
lambdify._import("scipy")
lambdify._import("numpy")

