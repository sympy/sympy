
"""
Canonicalization of Riemann invariants.
--------------------------------------

In [1] it is remarked that the products of tensors with the form
(- sign for covariant index)
   R(d1,d2,-d3,-d4)*R(d3,d4,-d5,-d6)*...*R(dn,dn,-d1,-d2)                 (1)

are the hardest Riemann invariants to canonicalize.

Here these invariants are computed; the input is equivalent to (1); it is
obtained by application of random slot and dummy symmetries.
The result is compared with xperm.c

Usage:
From this directory run:
python test_with_xperm.py nr random > test_xperm.cc

where nr is the number of Riemann tensors;
random=0     means that the canonical tensor has the form (1);
random != 0  random canonical tensor

In stderr appears the input and the output canonical tensor as products of
Riemann tensors, and the list in the form of the output of test_xperm

Compile the output test_xperm.cc with the version of xperm.c in [3]

run ./test_xperm > test_output
if xperm.c gives different result from the one shown in the stderr of the
Python program, an error message appears in stderr (except in the case in which
the result is the identity; so if there is an error message, check
in test_output if the result is the identity; if the stderr of the
python program is also the identity, the result is correct)


References:

  [1] J. M. Martin-Garcia, Comp. Phys. Commun. 179 (2008) 597-603 , arXiv: 0803.0862
  [2] xperm.c part of XPerm written by J. M. Martin-Garcia
      http://www.xact.es/index.html

  [3] cadabra by Kasper Peeters, http://cadabra.phi-sci.com/
      compile test_xperm.cc with the version of xperm.c dated 6 May 2006,
      which is included in cadabra-1.31/src/modules
"""

import sys
sys.path.insert(0, '../..')
from time import time
from sympy.combinatorics.permutations import Permutation, cyclic, perm_af_invert , perm_af_muln
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.perm_algorithms import double_coset_can_rep,dummy_sgs, PermutationGroup
from random import randint, seed, shuffle

seed(10)

def perm_af_direct_product(gens1, gens2, signed=False):
    """
    direct products of the generators gens1 and gens2
    """
    s = 2 if signed else 0
    n1 = len(gens1[0]) - s
    n2 = len(gens2[0]) - s
    start = range(n1)
    end = range(n1, n1 + n2)
    if signed:
        gens1 = [gen[:-2] + end + [gen[-2]+n2, gen[-1]+n2] for gen in gens1]
        gens2 = [start + [x + n1 for x in gen] for gen in gens2]
    else:
        gens1 = [gen + end for gen in gens1]
        gens2 = [start + [x + n1 for x in gen] for gen in gens2]

    return gens1 + gens2

riemann_sgs = [[1,0,2,3,5,4], [0,1,3,2,5,4], [2,3,0,1,4,5]]

class RiemannMonomial(object):
    """
    class to create the product of Riemann tensors with its slot symmetries
    due to the Riemann identities and commutativity
    """
    def __init__(self, indices, sgens, strong_base):
        self.indices = indices
        self.sgens = sgens
        self.strong_base = strong_base

    def __mul__(self, other):
        """
        other is required to be a Riemann tensor, not a product of them
        """
        if len(other.indices) > 4:
            raise NotImplementedError
        indices = self.indices + other.indices
        n = len(indices)
        sgens = perm_af_direct_product(self.sgens, other.sgens, 1)
        a = range(n-8)
        a += [n-4,n-3,n-2,n-1,n-8,n-7,n-6,n-5,n,n+1]
        sgens.append(a)
        n1 = len(self.indices)
        base = self.strong_base + [x + n1 for x in other.strong_base]
        return RiemannMonomial(indices, sgens, base)


########### code generation
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

############

def str_riemann(ind, g):
    """
    string representation of a product of Riemann tensors
    """
    size = len(g)
    if g[-1] == size - 1:
        s = 'R('
    else:
        s = '-R('
    for i in range(len(ind)):
        if i%4 != 3:
            s += '%s,' % ind[g[i]]
        else:
            if i != len(ind) - 1:
                s += '%s)*R(' % ind[g[i]]
            else:
                s += '%s)'% ind[g[i]]
    return s

def str_code(gens, g, result):
    size = len(g)
    # convert
    gens = [[y+1 for y in x] for x in gens]
    g = [y+1 for y in perm_af_invert(g)]
    s = str_head(gens, size) + str_perm(len(gens), g, result)
    return s

def run(gens, g, ind):
    """
    gens slot symmetry generators
    `g` permutation corresponding to the tensor, with the last two
    indices for the sign

    FIXME: it works only when the tensor appears with positive sign,
    because the sign is not passed correctly to xperm
    """
    size = len(g)
    #assert g[size-1] == size-1
    sgens = [Permutation(x) for x in gens]
    t0 = time()
    res = double_coset_can_rep(0, sgens, g)
    t1 = time()
    if res:
        sys.stderr.write('\noutput %s\n' % str_riemann(ind, res))
        result = [x+1 for x in perm_af_invert(res)]
        sys.stderr.write('%s\n' % result)
    else:
        sys.stderr.write('0\n')
        result = [0]*size
    sys.stderr.write('double_coset_can_rep: %f\n' %(t1-t0))
    return str_code(gens, g, result)

def riemann_products(nr, random_input):
    """
    construct a product of Riemann tensors
    """
    r = RiemannMonomial([0,2,5,7], riemann_sgs, [0,2])
    for i in range(1, nr):
        j = 4*i
        a = [j,j+2,j+5,j+7]
        if i == nr-1:
            a[2] = 1
            a[3] = 3
        r = r*RiemannMonomial(a, riemann_sgs, [0,2])
    gens = r.sgens
    n = len(r.indices)
    # slot symmetries
    sgens = [Permutation(x) for x in r.sgens]
    g = r.indices[:]

    if random_input:
        shuffle(g)
    g += [n, n+1]

    ind = []
    for i in range(n//2):
        ind.append('d%s' % (i+1))
        ind.append('-d%s' % (i+1))
    dummies = range(n)
    dsgs = dummy_sgs(dummies, 0, len(dummies))
    sgens = r.sgens
    sgensp = [Permutation(x) for x in sgens]
    ndsgs = len(dsgs)
    nsgens = len(sgens)
    # use random slot and dummy symmetries to bring the tensor to a random equivalent form
    dsgs = [Permutation(x) for x in dsgs]
    sgens = [Permutation(x) for x in sgens]
    t0 = time()
    S = PermutationGroup(sgens)
    S_ord = S.order()
    t1 = time()
    D = PermutationGroup(dsgs)
    D_ord = D.order()
    t2 = time()
    sys.stderr.write('setup Schreier-Sims S: %.2f D:%.2f\n' %(t1-t0,t2-t1))
    js = randint(0, S_ord-1)
    jd = randint(0, D_ord-1)
    s = S.coset_unrank(js).array_form
    d = D.coset_unrank(jd).array_form
    g1 = perm_af_muln(d, g, s)
    size = len(g1)
    if g1[-1] == size-2:
        g1[-2] = size-2
        g1[-1] = size-1
    sys.stderr.write('input %s\n' % str_riemann(ind, g1))
    result = run(gens, g1, ind)
    print result

if __name__ == '__main__':
    try:
        nr = int(sys.argv[1])
        random_input = int(sys.argv[2])
    except:
        sys.stderr.write(
'''enter nr random
         nr     =     number of Riemann tensors  > 1
         random = 0   use the lowest canonical tensor
                  1   use a random canonical tenso
         '''
)
        sys.exit()
    riemann_products(nr, random_input)
