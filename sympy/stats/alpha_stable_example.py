"""
Examples of using the alpha-stable distribution in SymPy.

The alpha-stable distribution is a generalization of the normal distribution
that allows for heavy tails and skewness.
"""

from sympy import Symbol, pprint
from sympy.stats import AlphaStable, density, characteristic_function

# Define symbolic variable
x = Symbol('x', real=True)
t = Symbol('t', real=True)

print("=" * 70)
print("Alpha-Stable Distribution Examples")
print("=" * 70)

# Example 1: Cauchy Distribution (alpha=1, beta=0)
print("\n1. Cauchy Distribution (alpha=1, beta=0)")
print("-" * 70)
X_cauchy = AlphaStable('X', 1, 0, 1, 0)
pdf_cauchy = density(X_cauchy)(x)
print("PDF:")
pprint(pdf_cauchy)
print("\nSimplified: {simplify(pdf_cauchy)}")

# Example 2: Normal Distribution (alpha=2)
print("\n2. Normal Distribution (alpha=2, beta=0)")
print("-" * 70)
Y_normal = AlphaStable('Y', 2, 0, 1, 0)
pdf_normal = density(Y_normal)(x)
print("PDF:")
pprint(pdf_normal)

# Example 3: Levy Distribution (alpha=0.5, beta=1)
print("\n3. Levy Distribution (alpha=0.5, beta=1)")
print("-" * 70)
from sympy import Rational
L_levy = AlphaStable('L', Rational(1, 2), 1, 1, 0)
print("Levy distribution created (PDF defined piecewise for x > 0)")

# Example 4: General Stable Distribution with Heavy Tails
print("\n4. Heavy-Tailed Stable Distribution (alpha=1.5, beta=0)")
print("-" * 70)
Z_heavy = AlphaStable('Z', 1.5, 0, 1, 0)
print("Distribution created with alpha=1.5 (heavier tails than normal)")
print("No closed-form PDF, but characteristic function available:")
cf = characteristic_function(Z_heavy)(t)
print("\nCharacteristic function phi(t):")
pprint(cf)

# Example 5: Skewed Stable Distribution
print("\n5. Skewed Stable Distribution (alpha=1.5, beta=0.5)")
print("-" * 70)
W_skewed = AlphaStable('W', 1.5, 0.5, 2, 1)
print("Parameters:")
print("  alpha (stability) = 1.5")
print("  alpha (skewness) = 0.5 (right-skewed)")
print("  scale = 2")
print("  location = 1")

# Example 6: Scaled and Shifted Cauchy
print("\n6. Scaled and Shifted Cauchy")
print("-" * 70)
C_scaled = AlphaStable('C', 1, 0, scale=3, location=5)
pdf_scaled = density(C_scaled)(x)
print("Cauchy with scale=3, location=5:")
print("PDF:")
pprint(pdf_scaled)

print("\n" + "=" * 70)
print("Applications:")
print("=" * 70)
print("""
Alpha-stable distributions are used in:
- Financial modeling (stock returns with fat tails)
- Signal processing (impulsive noise)
- Physics (anomalous diffusion)
- Network traffic modeling
- Extreme value analysis

Special cases:
- alpha=2: Normal distribution (finite variance)
- alpha=1, beta=0: Cauchy distribution (no moments)
- alpha=0.5, beta=1: Levy distribution (used in Levy flights""")
