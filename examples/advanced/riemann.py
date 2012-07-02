
"""
Canonicalization of Riemann invariants.
--------------------------------------

Canonicalize Riemann invariants using the monoterm symmetries of the
Riemann tensor, commutativity and dummy index symmetry to put
the Riemann invariant in the form in which the sequence of its indices
is the least possible according to the ordering of the indices.

For instance (- sign for covariant index)
R(-d4,-d2,-d6,-d3)*R(-d1,-d5,d2,d4)*R(d3,d6,d1,d5)
is equivalent to
R(d1,d2,d3,d4)*R(-d1,-d2,d5,d6)*R(-d3,-d4,-d5,-d6)
which has ordering of indices
(d1,d2,d3,d4,-d1,-d2,d5,d6,-d3,-d4,-d5,-d6
which is the least possible according to the ordering of indices
(d1,-d1,d2,-d2,d3,-d3,d4,-d4,d5,-d5,d6,-d6)

This is an example of the class of Riemann invariants equivalent to the form
   R(d1,d2,d3,d4)*R(-d3,-d4,d5,d6)*...*R(dn,dn,-d1,-d2)                 (1)

which in [1] were considered hard Riemann invariants to canonicalize.

Usage as a script
-----------------

From this directory run:
python riemann.py nr random_input random_regular > test_xperm.cc

where
  nr is the number of Riemann tensors;
  random_input=0     means that the canonical tensor has the form (1);
  random_input != 0  random canonical tensor, in general not of the form (1)

  random_regular != 0 use networkx to get a random regular graph to
                      determine the contractions of the Riemann tensors.

Example usage as a script:
python riemann.py 3 0 0 > test_xperm.cc
input R(-d4,-d2,-d6,-d3)*R(-d1,-d5,d2,d4)*R(d3,d6,d1,d5)

output R(d1,d2,d3,d4)*R(-d1,-d2,d5,d6)*R(-d3,-d4,-d5,-d6)
[1, 5, 2, 6, 3, 9, 4, 10, 7, 11, 8, 12, 13, 14]

In stderr appears the input and the output canonical tensor as products of
Riemann tensors, and the list in the form of the output of test_xperm

The stdout test_xperm.cc can be compiled with xperm.c [2] in the version
present in [3] to test the correctness of the result.

run ./test_xperm > test_output
if xperm.c gives different result from the one shown in the stderr of the
Python program, an error message appears in stderr (in the case in which
the result is the identity, there is an error message; in that case check
in test_output if the result is the identity; if the stderr of the
python program is also the identity, the result is correct)

Usage as a module
-----------------

>>> from riemann import RiemannMonomial
>>> s = 'R(-d12,d2,-d1,d6)*R(-d5,-d7,-d11,-d4)*R(d3,d8,-d2,d12)*R(d7,d5,d1,-d6)*R(d10,d9,d11,d4)*R(-d3,-d8,-d9,-d10)'
>>> r = RiemannMonomial.from_string(s)
>>> r.canonic()
-R(d1,d2,d3,d4)*R(-d1,-d2,d5,d6)*R(-d3,-d4,d7,d8)*R(-d5,-d6,d9,d10)*R(-d7,-d8,d11,d12)*R(-d9,-d10,-d11,-d12)


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
from sympy.combinatorics.permutations import Permutation, cyclic, perm_af_invert , perm_af_muln, perm_af_mul
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.perm_algorithms import double_coset_can_rep,dummy_sgs, PermutationGroup, orbit_transversal
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

riemann_gens = [[1,0,2,3,5,4], [0,1,3,2,5,4], [2,3,0,1,4,5]]

class RiemannMonomial(object):
    """
    class to create the product of Riemann tensors with its slot symmetries
    due to the Riemann identities and commutativity

    Use sgens = riemann_gens when creating a RiemannMonomial
    multiplying it by RiemannMonomial with only a Riemann tensor,
    the result has sgens equal to the SGS of the Riemann monomial

    Examples
    ========
    >>> from riemann import RiemannMonomial
    >>> r = RiemannMonomial([2,0,5,7])
    >>> r = r*RiemannMonomial([4,6,1,3])
    R(d1,d0,-d2,-d3)*R(d2,d3,-d0,-d1)
    >>> r.canonic()
    -R(d0,d1,d2,d3)*R(-d0,-d1,-d2,-d3)

    where the default index names 'd0, d1',..., are used
    indices are ordered according to
    ind = ['d0','-d0','d1','-d1','d2','-d2','d3','-d3']
    RiemannMonomial([2,0,5,7]) corresponds to R^{d1 d0}_{d2 d3}
    RiemannMonomial([4,6,1,3]) corresponds to R^{d2 d3}_{d0 d1}

    One can choose different index names:
    >>> r.set_index_names('a b c d'.split())
    >>> r
    R(b,a,-c,-d)*R(c,d,-a,-b)

    >>> s = 'R(-d4,-d2,-d6,-d3)*R(-d1,-d5,d2,d4)*R(d3,d6,d1,d5)'
    >>> r = RiemannMonomial.from_string(s)
    >>> r.canonic()
    R(d1,d2,d3,d4)*R(-d1,-d2,d5,d6)*R(-d3,-d4,-d5,-d6)

    """
    def __init__(self, indices, sgens=riemann_gens, sign=False):
        self.indices = indices
        self.sgens = sgens
        self.index_names = None
        self.sign = sign
        self.ind = None

    def set_index_names(self, index_names):
        self.index_names = index_names
        self.ind = get_indices(index_names)

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
        r = RiemannMonomial(indices, sgens)
        if self.ind:
            r.ind = self.ind
        return r

    def canonic(self):
        """
        canonicalize the RiemannMonomial

        Examples
        ========

        >>> from riemann import RiemannMonomial
        >>> s = 'R(-d12,d2,-d1,d6)*R(-d5,-d7,-d11,-d4)*R(d3,d8,-d2,d12)*R(d7,d5,d1,-d6)*R(d10,d9,d11,d4)*R(-d3,-d8,-d9,-d10)'
        >>> r = RiemannMonomial.from_string(s)
        >>> r.canonic()
        -R(d1,d2,d3,d4)*R(-d1,-d2,d5,d6)*R(-d3,-d4,d7,d8)*R(-d5,-d6,d9,d10)*R(-d7,-d8,d11,d12)*R(-d9,-d10,-d11,-d12)
        """
        g = self.indices[:]
        n = len(g)
        g = g + [n+1, n] if self.sign else g + [n, n+1]
        size = n + 2
        sgens = self.sgens
        sgens = [Permutation(x) for x in sgens]
        sgs = get_sgens(sgens, n)
        res = double_coset_can_rep(0, sgens, g, sgs)
        has_minus_sign = (res[-1] == size-2)
        nr = len(res)//4
        r = RiemannMonomial(res[:4])
        r.ind = self.ind
        for i in range(4, len(res)-2, 4):
            r = r*RiemannMonomial(res[i:i+4])
        r.sign = has_minus_sign
        return r

    def __str__(self):
        n = len(self.indices)
        if not self.ind:
            index_names = ['d%d' % i for i in range(n//2)]
            ind = get_indices(index_names)
        else:
            ind = self.ind
        if self.sign:
            return str_riemann(ind, self.indices+[n+1,n])
        else:
            return str_riemann(ind, self.indices+[n,n+1])

    __repr__ = __str__

    @classmethod
    def from_string(self, s):
        """
        return a RiemannMonomial from a string

        if the index names are of the form '\w*\d+', they are
        ordered lexicographically in the letter part, then by number
        orderind, e.g.
        index_names = ['d7','d3','d4','d6','d8','d10','d9','d1','d2','d5']
        is sorted to
        ['d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10']

        Examples
        ========

        >>> from riemann import RiemannMonomial
        >>> s = 'R(-d4,-d2,-d6,-d3)*R(-d1,-d5,d2,d4)*R(d3,d6,d1,d5)'
        >>> r = RiemannMonomial.from_string(s)
        >>> r.canonic()
        R(d1,d2,d3,d4)*R(-d1,-d2,d5,d6)*R(-d3,-d4,-d5,-d6)
        """
        if s[0] == '-':
            sign = True
            s = s[1:]
        else:
            sign = False
        s = s.replace('R','').replace(')','').replace('*','')[1:].replace('(',',')
        a = s.split(',')
        index_names = [x for x in a if x[0] != '-']
        # sort the index names
        index_names1 = []
        index_names1 = []
        for index_name in index_names:
            found = 0
            for i, c in enumerate(index_name):
                if c.isdigit():
                    index_names1.append((index_name[:i], int(index_name[i:])))
                    found = 1
                    break
            if not found:
                index_names1.append((index_name, -1))
        index_names1.sort()
        index_names = []
        for c, i in index_names1:
            if i >= 0:
                index_names.append('%s%d' %(c, i))
            else:
                index_names.append(c)
        ind = get_indices(index_names)
        indices = [ind.index(x) for x in a]
        r = RiemannMonomial(indices[:4])
        for i in range(4, len(indices), 4):
            r = r*RiemannMonomial(indices[i:i+4])
        r.set_index_names(index_names)
        r.sign = sign
        return r


########### code generation
def _str_head(gens, size):
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

def _str_perm(num_gens, p, result):
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
    if p[-1] == size:
        s1 += '%d,%d};\n' %(size-1,size)
    else:
        s1 += '%d,%d};\n' %(size,size-1)
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

def get_indices(ind_names):
    """
    ordered list of indices from index names

    Examples
    ========
    >>> from riemann import get_indices
    >>> get_indices('a b c d'.split())
    ['a', '-a', 'b', '-b', 'c', '-c', 'd', '-d']
    """
    ind = []
    for ind_name in ind_names:
        ind.append(ind_name)
        ind.append('-' + ind_name)
    return ind

def str_riemann(ind, g):
    """
    string representation of a product of Riemann tensors

    Examples
    ========

    >>> from riemann import get_indices, str_riemann
    >>> N = 2
    >>> n = 2*N
    >>> ind_names = ['d%d' % i for i in range(n)]
    >>> ind = get_indices(ind_names); ind
    ['d0', '-d0', 'd1', '-d1', 'd2', '-d2', 'd3', '-d3']
    >>> g = [0, 2, 5, 7, 4, 6, 1, 3]
    >>> str_riemann(ind, g+[2*n,2*n+1])
    'R(d0,d1,-d2,-d3)*R(d2,d3,-d0,-d1)'
    """
    size = len(g)
    if g[-1] == size - 1:
        s = 'R('
    else:
        s = '-R('
    n = size - 2
    for i in range(n):
        if i%4 != 3:
            s += '%s,' % ind[g[i]]
        else:
            if i != n - 1:
                s += '%s)*R(' % ind[g[i]]
            else:
                s += '%s)'% ind[g[i]]
    return s

def str_code(gens, g, result):
    """
    produce the code for test_xperm.cc
    """
    size = len(g)
    # convert
    gens = [[y+1 for y in x] for x in gens]
    g = [y+1 for y in perm_af_invert(g)]
    s = _str_head(gens, size) + _str_perm(len(gens), g, result)
    return s

def riemann_canon(sgens, g, ind, sgs):
    """
    sgens slot symmetry generators
    `g` permutation corresponding to the tensor, with the last two
    indices for the sign
    ind list of index names
    sgs = (S_cosets, b_S)
    """
    size = len(g)
    #assert g[size-1] == size-1
    t0 = time()
    res = double_coset_can_rep(0, sgens, g, sgs)
    t1 = time()
    if res:
        sys.stderr.write('\noutput %s\n' % str_riemann(ind, res))
        result = [x+1 for x in perm_af_invert(res)]
        sys.stderr.write('%s\n' % result)
    else:
        sys.stderr.write('0\n')
        result = [0]*size
    sys.stderr.write('double_coset_can_rep: %f\n' %(t1-t0))
    gens = [h.array_form for h in sgens]
    return gens, g, result

def random_dummy(n):
    """
    random permutation of the dummy indices

    ind = ['d1','-d1','d2','-d2',...]; n = len(ind)
    Return a random permutation `d`
    ind -> [ind[i] for i in d]
    """
    a = range(n)
    # random raising and lowering of dummy indices
    b = []
    for i in range(0, n, 2):
        j = randint(0,1)
        if j:
            b.append((a[i], a[i+1]))
        else:
            b.append((a[i+1], a[i]))
    # random permutation of indices
    p = range(n//2)
    shuffle(p)
    b = [b[i] for i in p]
    res = []
    for i in range(n//2):
        res.extend(b[i])
    return res

def get_sgens(sgens, n):
    """
    return (S.coset_repr(), S.strong_base())

    sgens is a strong generating set for a permutation group of
    signed permutations

    n size of permutations without sign
    """
    sgensx = [x.array_form for x in sgens]
    b_S = []
    S_cosets = []
    # generate S_cosets, b_S
    for ii in range(n-1):
        if not sgensx:
            S_cosets.append([range(n+2)])
            continue
        Sxtrav = orbit_transversal(sgensx, ii, af=True)
        if len(Sxtrav) > 1:
            b_S.append(ii)
        S_cosets.append(Sxtrav)
        nSxtrav = len(Sxtrav)
        sgensx = [h for h in sgensx if h[ii] == ii]
    return S_cosets, b_S

def get_random_S_element(S_cosets, n):
    random_s = range(n+2)
    for ii in range(n-1):
        Sxtrav = S_cosets[ii]
        nSxtrav = len(Sxtrav)
        random_Sxtrav_element = Sxtrav[randint(0, nSxtrav-1)]
        random_s = perm_af_mul(random_s, random_Sxtrav_element)
    return random_s

def riemann_products(nr, random_input, random_regular):
    """
    construct a product of Riemann tensors, canonicalize it

    return the C code to be run with xperm.c to check the result

    nr number of Riemann tensors

    random_input = 0 and random_regular = 0 : produce a Riemann invariant
    equivalent to
      R^{d1 d2}_{d3 d4} R^{d3 d4}_{d5 d6} ... R^{d_{n-1} dn}_{d1 d2}  (1)

    random_input  != 0  randomize to an inequivalent Riemann invariant

    random_regular != 0 produce a Riemann invariant with graph structure
                        which is a random regular graph without multilines
    """

    if not random_regular:
        # produce the Riemann invariant (1)
        r = RiemannMonomial([0,2,5,7], riemann_gens)
        for i in range(1, nr):
            j = 4*i
            a = [j,j+2,j+5,j+7]
            if i == nr-1:
                a[2] = 1
                a[3] = 3
            r = r*RiemannMonomial(a, riemann_gens)

    else:
        # produce a Riemann invariant with graph structure
        # which is a random regular graph without multilines
        import networkx as nx
        gx = nx.random_regular_graph(4, nr)

        en = 0
        vv = [[] for k in range(nr)]
        a = gx.adjacency_list()
        r = 0
        for v1, vs in enumerate(a):
            for v2 in vs:
                if v2 < v1:
                    continue
                vv[v1].append(en)
                vv[v2].append(en+1)
                en += 2
            if not r:
                r = RiemannMonomial(vv[v1], riemann_gens)
            else:
                r = r*RiemannMonomial(vv[v1], riemann_gens)

    gens = r.sgens
    n = len(r.indices)
    # slot symmetries
    sgens = [Permutation(x) for x in r.sgens]
    # permutation element without sign corrending to the Riemann invariant
    g = r.indices[:]

    if random_input == 2:
        shuffle(g)

    # add plus sign to it
    g += [n, n+1]

    # index names, contravariant and covariant, written in ordered list
    ind = []
    for i in range(n//2):
        ind.append('d%s' % (i+1))
        ind.append('-d%s' % (i+1))
    # SGS generators dsgs of the dummy symmetry group
    t0 = time()
    dummies = range(n)
    dsgs = dummy_sgs(dummies, 0, len(dummies))
    # SGS generators of the slot symmetry group
    sgens = r.sgens
    ndsgs = len(dsgs)
    # use random slot and dummy symmetries to bring the tensor to a random equivalent form
    dsgs = [Permutation(x) for x in dsgs]
    sgens = [Permutation(x) for x in sgens]
    sgs = get_sgens(sgens, n)
    s = get_random_S_element(sgs[0], n)
    t1 = time()
    sys.stderr.write('setup SGS: %.2f\n' %(t1-t0))
    d = random_dummy(n) + [n, n+1]
    g1 = perm_af_muln(d, g, s)
    return sgens, g1, ind, sgs

def ri2_code(ind, g1):
    """
    write file ri2.py
    """
    sys.stderr.write('input %s\n' % str_riemann(ind, g1))
    s11 = '''from riemann import RiemannMonomial
from time import time
t0 = time()
s = RiemannMonomial.from_string("%s").canonic()
t1 = time()
print (s)
print ('%%.3f' %%(t1-t0))
''' % str_riemann(ind, g1)
    open('ri2.py','w').write(s11)

if __name__ == '__main__':
    try:
        nr = int(sys.argv[1])
        random_input = int(sys.argv[2])
        random_regular = int(sys.argv[3])
    except:
        sys.stderr.write(
'''enter nr random_input random_regular
         nr           =     number of Riemann tensors  > 1
         random_input = 0   use the lowest canonical tensor
                        1   use a random canonical tensor
         random_regular  = 0  use a graph of the form equivalent to
                              R^{d1 d2 d3 d4}*R^(-d1 d2}^{d5,d6} ..
                              (sausage regular graph)
                           1  use a random regular graph;
                              networkx must be installed to run this
         '''
)
        sys.exit()
    if random_regular and nr < 5:
        sys.stderr.write("to get a random regular graph use nr >= 5\n")
        sys.exit()
    # generate random riemann products
    sgens, g1, ind, sgs = riemann_products(nr, random_input, random_regular)
    # write file ri2.py
    ri2_code(ind, g1)
    # canonize riemann product and write it to stderr
    sgens, g, result = riemann_canon(sgens, g1, ind, sgs)
    #print test_xperm.cc code
    print str_code(sgens, g, result)
