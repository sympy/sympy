#!/usr/bin/env python

"""Vandermonde matrix example

Demonstrates matrix computations using the Vandermonde matrix.
  * http://en.wikipedia.org/wiki/Vandermonde_matrix
"""

from __future__ import division, print_function

from sympy import Matrix, pprint, Rational, sqrt, symbols, Symbol, zeros
from sympy.core.compatibility import range


def symbol_gen(sym_str):
    """Symbol generator

    Generates sym_str_n where n is the number of times the generator
    has been called.
    """
    n = 0
    while True:
        yield Symbol("%s_%d" % (sym_str, n))
        n += 1


def comb_w_rep(n, k):
    """Combinations with repetition

    Returns the list of k combinations with repetition from n objects.
    """
    if k == 0:
        return [[]]
    combs = [[i] for i in range(n)]
    for i in range(k - 1):
        curr = []
        for p in combs:
            for m in range(p[-1], n):
                curr.append(p + [m])
        combs = curr
    return combs


def vandermonde(order, dim=1, syms='a b c d'):
    """Computes a Vandermonde matrix of given order and dimension.

    Define syms to give beginning strings for temporary variables.

    Returns the Matrix, the temporary variables, and the terms for the
    polynomials.
    """
    syms = syms.split()
    n = len(syms)
    if n < dim:
        new_syms = []
        for i in range(dim - n):
            j, rem = divmod(i, n)
            new_syms.append(syms[rem] + str(j))
        syms.extend(new_syms)
    terms = []
    for i in range(order + 1):
        terms.extend(comb_w_rep(dim, i))
    rank = len(terms)
    V = zeros(rank)
    generators = [symbol_gen(syms[i]) for i in range(dim)]
    all_syms = []
    for i in range(rank):
        row_syms = [next(g) for g in generators]
        all_syms.append(row_syms)
        for j, term in enumerate(terms):
            v_entry = 1
            for k in term:
                v_entry *= row_syms[k]
            V[i*rank + j] = v_entry
    return V, all_syms, terms


def gen_poly(points, order, syms):
    """Generates a polynomial using a Vandermonde system"""
    num_pts = len(points)
    if num_pts == 0:
        raise ValueError("Must provide points")
    dim = len(points[0]) - 1
    if dim > len(syms):
        raise ValueError("Must provide at lease %d symbols for the polynomial" % dim)
    V, tmp_syms, terms = vandermonde(order, dim)
    if num_pts < V.shape[0]:
        raise ValueError(
            "Must provide %d points for order %d, dimension "
            "%d polynomial, given %d points" %
            (V.shape[0], order, dim, num_pts))
    elif num_pts > V.shape[0]:
        print("gen_poly given %d points but only requires %d, "\
            "continuing using the first %d points" % \
            (num_pts, V.shape[0], V.shape[0]))
        num_pts = V.shape[0]

    subs_dict = {}
    for j in range(dim):
        for i in range(num_pts):
            subs_dict[tmp_syms[i][j]] = points[i][j]
    V_pts = V.subs(subs_dict)
    V_inv = V_pts.inv()

    coeffs = V_inv.multiply(Matrix([points[i][-1] for i in range(num_pts)]))

    f = 0
    for j, term in enumerate(terms):
        t = 1
        for k in term:
            t *= syms[k]
        f += coeffs[j]*t
    return f


def main():
    order = 2
    V, tmp_syms, _ = vandermonde(order)
    print("Vandermonde matrix of order 2 in 1 dimension")
    pprint(V)

    print('-'*79)
    print("Computing the determinant and comparing to \sum_{0<i<j<=3}(a_j - a_i)")

    det_sum = 1
    for j in range(order + 1):
        for i in range(j):
            det_sum *= (tmp_syms[j][0] - tmp_syms[i][0])

    print("""
    det(V) = %(det)s
    \sum   = %(sum)s
           = %(sum_expand)s
    """ % {"det": V.det(),
            "sum": det_sum,
            "sum_expand": det_sum.expand(),
          })

    print('-'*79)
    print("Polynomial fitting with a Vandermonde Matrix:")
    x, y, z = symbols('x,y,z')

    points = [(0, 3), (1, 2), (2, 3)]
    print("""
    Quadratic function, represented by 3 points:
       points = %(pts)s
       f = %(f)s
    """ % {"pts": points,
            "f": gen_poly(points, 2, [x]),
          })

    points = [(0, 1, 1), (1, 0, 0), (1, 1, 0), (Rational(1, 2), 0, 0),
              (0, Rational(1, 2), 0), (Rational(1, 2), Rational(1, 2), 0)]
    print("""
    2D Quadratic function, represented by 6 points:
       points = %(pts)s
       f = %(f)s
    """ % {"pts": points,
            "f": gen_poly(points, 2, [x, y]),
          })

    points = [(0, 1, 1, 1), (1, 1, 0, 0), (1, 0, 1, 0), (1, 1, 1, 1)]
    print("""
    3D linear function, represented by 4 points:
       points = %(pts)s
       f = %(f)s
    """ % {"pts": points,
            "f": gen_poly(points, 1, [x, y, z]),
          })


if __name__ == "__main__":
    main()
