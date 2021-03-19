#!/usr/bin/env python

"""Integration example
Demonstrates some integral operations.
"""

from sympy import symbols, integrate, pprint

def main():
    a,b = symbols("a,b")
    e = a ** 2 + a + 1

    # Evaluating Indefinite Integrals
    print("\nExpression : ")
    print()
    pprint(e)
    print("\n\nIntegrating w.r.t. a:")
    print()
    pprint(integrate(e,a))
    print("\n\nIntegrating w.r.t. b:")
    print()
    pprint(integrate(e,b))
    # Evaluating Definite Integrals
    print("\n\nIntegrating w.r.t. a with lower limit = 2, upper limit = 4:")
    print()
    pprint(integrate(e, (a,2,4)))

if __name__ == "__main__":
    main()
