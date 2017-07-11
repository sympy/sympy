from __future__ import division, print_function

from sympy import Symbol, eye, solve_linear_system, zeros
from sympy.core.compatibility import range

N = 8
M = zeros(N, N + 1)
M[:, :N] = eye(N)
S = [Symbol('A%i' % i) for i in range(N)]


def timeit_linsolve_trivial():
    solve_linear_system(M, *S)
