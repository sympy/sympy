
"""
References:

  [1] xperm.c part of XPerm written by J. M. Martin-Garcia
        http://www.xact.es/index.html

  [2] test_xperm.c in cadabra by Kasper Peeters, http://cadabra.phi-sci.com/

Test that the tensor canonicalization in the case of tensor without
free indices gives the same result with xperm.c

In the end of this file choose the slot generators and `g` representing
the tensor.

From this directory run:
python test_with_xperm.py > test_xperm.cc

For instance with the choice put below one gets
(advanced)$ pyt test_with_xperm.py > test_xperm.cc
[1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 10, 48, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25, 28, 27, 30, 29, 32, 31, 34, 33, 36, 35, 38, 37, 40, 39, 42, 41, 44, 43, 46, 45, 49, 47, 50, 51, 52]

In stderr appears the result, in the form of the result
of running test_xperm, see example test_xperm.cc in [2];
compile test_xperm.c with xperm.c, see [2]

run ./test_xperm > test_output
if xperm.c gives different result from the one shown in the stderr,
an error message appears in stderr (except in the case in which
the result is the identity; so if there is an error message, check
in test_output if the result is the identity; if the stderr of the
python program is also the identity, the result is correct))
"""

import sys
sys.path.insert(0,'../..')
from time import time
from sympy.combinatorics.permutations import Permutation, cyclic, perm_af_invert
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.perm_algorithms import double_coset_can_rep


def str_head(gens, size):
    num_gens = len(gens)
    size = len(gens[0])
    s = \
"""
extern "C" {
#include "xperm.h"
}

#include <iostream>
#include <cassert>

void test()
{
   int gs[%d*%d]={
""" %(num_gens, size)
    for g in gens:
       for x in g:
           s += '%d,' % x
       s += '\n'
    s = s[:-2]
    s += '};\n'
    return s

def str_perm(num_gens, p, result):
    size = len(p)
    s1 = """
    int dummysetlabels[%d]={""" %(size)
    for x in range(size-2):
        s1 += '%d, ' % 1
    s1 = s1[:-2]
    s1 += '};\n'
    s1 += "int result[%d]={""" %(size)
    for i in range(size):
        s1 += '%d, ' % result[i]
    s1 = s1[:-2]
    s1 += '};\n'

    s1 += '    int dummies[%d]={' %(size-2)
    for x in p[:-2]:
        s1 += '%d, ' % x
    s1 = s1[:-2]
    s1 += '};\n'
    s1 += '    int perm[%d]={' %(size)
    for x in p[:-2]:
        s1 += '%d, ' % x
    s1 += '%d,%d};\n' %(size-1,size)
    s1 +='''
    int free_indices[2];
    int cperm[%d];
    canonical_perm(perm,0,0,0,gs,%d,%d,free_indices,0,dummies,%d,
                   dummysetlabels,1,1,cperm);
    ''' % (size, num_gens, size, size//2 - 1)
    s1 += '''
    std::cout << "Test" << std::endl;
    for(unsigned int i=0; i<%d; ++i) {
        std::cout << cperm[i] << " ";
        if(cperm[i] != result[i])
          std::cerr << "ERR i="<<i<<" "<<cperm[i] << " " << result[i]<< std::endl;
    }
    std::cout << std::endl;

}
    ''' % size
    s1 += '''
int main(int argc, char **argv)
{
  test();
}
'''

    return s1


def str_code(gens, g, result):
    size = len(g)
    # convert
    gens = [[y+1 for y in x] for x in gens]
    g = [y+1 for y in perm_af_invert(g)]
    s = str_head(gens, size) + str_perm(len(gens), g, result)
    return s

def run(gens, g):
    size = len(g)
    sgens = [Permutation(x) for x in gens]
    t0 = time()
    res = double_coset_can_rep(0, sgens, g)
    t1 = time()
    result = [x+1 for x in perm_af_invert(res)]
    sys.stderr.write('%s\n' % result)
    sys.stderr.write('%f\n' %(t1-t0))
    return str_code(gens, g, result)


if __name__ == '__main__':

    n = 90
    # The last two indices of the permutations represent the sign.
    # slot symmetry generators
    gens = [Permutation(cyclic(a, n)).array_form + [n,n+1] for a in
                        [ [(1,i)] for i in range(3, n - 2, 2)]]
    # choose a tensor with only dummy indices, represented as a
    # permutation
    i = 10
    g = Permutation.unrank_nonlex(n, i).array_form + [n, n+1]

    result = run(gens, g)
    print result
