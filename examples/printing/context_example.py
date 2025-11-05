#!/usr/bin/env python
"""
Example demonstrating the ConTeXt printer in SymPy.

ConTeXt is a document preparation system based on TeX, similar to LaTeX.
This example shows how to convert SymPy expressions to ConTeXt format.
"""

from sympy import symbols, sqrt, sin, cos, exp, log, Integral, Sum, pi, I, Rational
from sympy import context, latex

def main():
    print("=" * 70)
    print("ConTeXt Printer Example")
    print("=" * 70)

    # Define symbols
    x, y, z, n = symbols('x y z n')

    print("\n1. Basic Expressions")
    print("-" * 70)
    expr1 = x**2 + y**2
    print(f"Expression: {expr1}")
    print(f"ConTeXt:    {context(expr1)}")
    print(f"LaTeX:      {latex(expr1)}")
    print(f"Same?       {context(expr1) == latex(expr1)}")

    print("\n2. Square Root and Fractions")
    print("-" * 70)
    expr2 = sqrt(x**2 + y**2)
    print(f"Expression: {expr2}")
    print(f"ConTeXt:    {context(expr2)}")

    expr3 = (x + 1) / (x - 1)
    print(f"\nExpression: {expr3}")
    print(f"ConTeXt:    {context(expr3)}")

    print("\n3. Trigonometric Functions")
    print("-" * 70)
    expr4 = sin(x)**2 + cos(x)**2
    print(f"Expression: {expr4}")
    print(f"ConTeXt:    {context(expr4)}")

    print("\n4. Exponential and Logarithm")
    print("-" * 70)
    expr5 = exp(I*pi) + 1
    print(f"Expression: {expr5}")
    print(f"ConTeXt:    {context(expr5)}")

    expr6 = log(x*y)
    print(f"\nExpression: {expr6}")
    print(f"ConTeXt:    {context(expr6)}")

    print("\n5. Integrals")
    print("-" * 70)
    expr7 = Integral(x**2, x)
    print(f"Expression: {expr7}")
    print(f"ConTeXt:    {context(expr7)}")

    expr8 = Integral(sin(x), (x, 0, pi))
    print(f"\nExpression: {expr8}")
    print(f"ConTeXt:    {context(expr8)}")

    print("\n6. Summations")
    print("-" * 70)
    expr9 = Sum(1/n**2, (n, 1, 10))
    print(f"Expression: {expr9}")
    print(f"ConTeXt:    {context(expr9)}")

    print("\n7. Different Modes")
    print("-" * 70)
    expr10 = x**2 + sqrt(y)

    print(f"Expression: {expr10}")
    print(f"\nPlain mode (default):")
    print(f"  {context(expr10, mode='plain')}")

    print(f"\nInline mode (with $ delimiters):")
    print(f"  {context(expr10, mode='inline')}")

    print(f"\nEquation mode (ConTeXt uses \\startformula...\\stopformula):")
    print(f"  {context(expr10, mode='equation')}")

    print(f"\nFor comparison, LaTeX equation mode:")
    print(f"  {latex(expr10, mode='equation')}")

    print("\n8. Key Differences: Equation Delimiters")
    print("-" * 70)
    expr11 = (2*x)**Rational(7, 2)

    print("LaTeX uses \\begin{equation}...\\end{equation}")
    print("ConTeXt uses \\startformula...\\stopformula")
    print()
    print(f"Expression: {expr11}")
    print(f"\nConTeXt equation mode:")
    print(f"  {context(expr11, mode='equation')}")
    print(f"\nLaTeX equation mode:")
    print(f"  {latex(expr11, mode='equation')}")

    print("\n9. Settings from LaTeX Printer")
    print("-" * 70)
    print("The ConTeXt printer inherits all settings from the LaTeX printer")
    print()

    expr12 = 3*x**2/y
    print(f"Expression: {expr12}")
    print(f"Default:           {context(expr12)}")
    print(f"fold_short_frac:   {context(expr12, fold_short_frac=True)}")

    expr13 = x**(Rational(3, 4))
    print(f"\nExpression: {expr13}")
    print(f"Default:           {context(expr13)}")
    print(f"fold_frac_powers:  {context(expr13, fold_frac_powers=True)}")

    print("\n" + "=" * 70)
    print("For more information, see the documentation:")
    print("  - sympy.printing.context.context()")
    print("  - sympy.printing.latex.latex()")
    print("=" * 70)


if __name__ == "__main__":
    main()
