#!/usr/bin/env python

from __future__ import print_function

from sympy.galgebra import xdvi, Get_Program, Print_Function
from sympy.galgebra import MV, Format


def basic_multivector_operations_3D():
    Print_Function()

    (ex, ey, ez) = MV.setup('e*x|y|z')

    print('g_{ij} =', MV.metric)

    A = MV('A', 'mv')

    A.Fmt(1, 'A')
    A.Fmt(2, 'A')
    A.Fmt(3, 'A')

    A.even().Fmt(1, '%A_{+}')
    A.odd().Fmt(1, '%A_{-}')

    X = MV('X', 'vector')
    Y = MV('Y', 'vector')

    X.Fmt(1, 'X')
    Y.Fmt(1, 'Y')

    (X*Y).Fmt(2, 'X*Y')
    (X ^ Y).Fmt(2, 'X^Y')
    (X | Y).Fmt(2, 'X|Y')
    return


def basic_multivector_operations_2D():
    Print_Function()

    (ex, ey) = MV.setup('e*x|y')

    print('g_{ij} =', MV.metric)

    X = MV('X', 'vector')
    A = MV('A', 'spinor')

    X.Fmt(1, 'X')
    A.Fmt(1, 'A')

    (X | A).Fmt(2, 'X|A')
    (X < A).Fmt(2, 'X<A')
    (A > X).Fmt(2, 'A>X')
    return


def dummy():
    return


def main():
    Get_Program(True)
    Format()

    basic_multivector_operations_3D()
    basic_multivector_operations_2D()

    xdvi('simple_test_latex.tex')
    return

if __name__ == "__main__":
    main()
