#!/usr/bin/env python
"""
Calculates the Coupled-Cluster energy- and amplitude equations
See 'An Introduction to Coupled Cluster Theory' by
T. Daniel Crawford and Henry F. Schaefer III.
http://www.ccc.uga.edu/lec_top/cc/html/review.html
"""

from sympy.physics.secondquant import (AntiSymmetricTensor, wicks,
        F, Fd, NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator)
from sympy import (
    symbols, expand, pprint, Number, latex
)

pretty_dummies_dict={
        'above':'cdefgh',
        'below':'klmno',
        'general':'pqrstu'
        }


def get_CC_operators():
    """
    Returns a tuple (T1,T2) of unique operators.
    """
    i = symbols('i',below_fermi=True,dummy=True)
    a = symbols('a',above_fermi=True,dummy=True)
    t_ai = AntiSymmetricTensor('t',(a,),(i,))
    ai = NO(Fd(a)*F(i))
    i,j = symbols('ij',below_fermi=True,dummy=True)
    a,b = symbols('ab',above_fermi=True,dummy=True)
    t_abij = AntiSymmetricTensor('t',(a,b),(i,j))
    abji = NO(Fd(a)*Fd(b)*F(j)*F(i))

    T1 = t_ai*ai
    T2 = Number((1,4))*t_abij*abji
    return (T1,T2)

def main():
    print
    print "Calculates the Coupled-Cluster energy- and amplitude equations"
    print "See 'An Introduction to Coupled Cluster Theory' by"
    print "T. Daniel Crawford and Henry F. Schaefer III"
    print "http://www.ccc.uga.edu/lec_top/cc/html/review.html"
    print

    # setup hamiltonian
    p,q,r,s = symbols('pqrs',dummy=True)
    f = AntiSymmetricTensor('f',(p,),(q,))
    pr = NO((Fd(p)*F(q)))
    v = AntiSymmetricTensor('v',(p,q),(r,s))
    pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))

    H=f*pr

    # Uncomment the next line to use a 2-body hamiltonian:
    # H=f*pr + Number(1,4)*v*pqsr

    print "Using the hamiltonian:", latex(H)

    print "Calculating nested commutators"
    C = Commutator

    T1,T2 = get_CC_operators()
    T = T1+ T2
    print "comm1..."
    comm1 = wicks(C(H,T),simplify_dummies=True, simplify_kronecker_deltas=True)

    T1,T2 = get_CC_operators()
    T = T1+ T2
    print "comm2..."
    comm2 = wicks(C(comm1,T),simplify_dummies=True, simplify_kronecker_deltas=True)

    T1,T2 = get_CC_operators()
    T = T1+ T2
    print "comm3..."
    comm3 = wicks(C(comm2,T),simplify_dummies=True, simplify_kronecker_deltas=True)

    T1,T2 = get_CC_operators()
    T = T1+ T2
    print "comm4..."
    comm4 = wicks(C(comm3,T),simplify_dummies=True, simplify_kronecker_deltas=True)

    print "construct Hausdoff expansion..."
    eq = H + comm1+comm2/2  +comm3/6+comm4/24
    eq = eq.expand()
    eq = evaluate_deltas(eq)
    eq = substitute_dummies(eq, new_indices=True, reverse_order=False,
            pretty_indices=pretty_dummies_dict)
    print "*********************"
    print

    print "extracting CC equations from full Hbar"
    i,j,k,l = symbols('ijkl',below_fermi=True)
    a,b,c,d = symbols('abcd',above_fermi=True)
    print
    print "CC Energy:"
    print latex(wicks(eq, simplify_dummies=True,
        keep_only_fully_contracted=True))
    print
    print "CC T1:"
    eqT1 = wicks(NO(Fd(i)*F(a))*eq, simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
    eqT1 = substitute_dummies(eqT1,reverse_order=False)
    print latex(eqT1)
    print
    print "CC T2:"
    eqT2 = wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*eq,simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
    P = PermutationOperator
    eqT2 = simplify_index_permutations(eqT2,[P(a,b),P(i,j)])
    print latex(eqT2)

if __name__ == "__main__":
    main()
