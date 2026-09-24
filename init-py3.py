from __future__ import division
from sympy import *

from sympy.plotting import plot_parametric, plot3d

# Output of the last expression
_ = None

x, y, z, t = symbols('x,y,z,t')
k, m, n = symbols('k,m,n', integer=True)
f, g, h = symbols('f,g,h', cls=Function)

# lambda is a reserved word, pi is a built-in symbol

_greek = 'alpha beta gamma delta epsilon zeta eta ' \
         'theta iota kappa mu nu xi omicron rho '   \
         'sigma tau upsilon phi chi psi omega '     \
         'Alpha Beta Gamma Delta Epsilon Zeta Eta ' \
         'Theta Iota Kappa Mu Nu Xi Omicron Rho '   \
         'Sigma Tau Upsilon Phi Chi Psi Omega '     \
         'α β γ δ ε ζ η θ'     \
         'ι κ λ μ ν ξ ο π ρ '  \
         'σ τ υ φ χ ψ ω '      \
         'Α Β Γ Δ Ε Ζ Η Θ'     \
         'Ι Κ Λ Μ Ν Ξ Ο Π Ρ '  \
         'Σ Τ Υ Φ Χ Ψ Ω'

for _symbol in _greek.split(' '):
    exec("%s = Symbol('%s', real=True)" % (_symbol, _symbol))

del _symbol
