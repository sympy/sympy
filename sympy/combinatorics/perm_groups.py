from sympy.combinatorics import Permutation
from sympy.core import Basic
from sympy.combinatorics.permutations import perm_af_mul, \
 _new_from_array_form, perm_af_commutes_with, perm_af_invert, perm_af_muln
from random import randint, randrange, choice
from sympy.functions.combinatorial.factorials import factorial
from math import log
from sympy.ntheory import isprime, sieve

def _smallest_change(h, alpha):
    """
    find the smallest point not fixed by `h`
    """
    for i in range(alpha, len(h)):
        if h[i] != i:
            return i

class _Vertex(object):
    """
    vertex of JGraph

    neighbor list of neighbor vertices
    perm     list of permutations
    index_neighbor list of index position of vertices in neighbor
    if vertex j is a neighbor of vertex i, then
       vertex[i].index_neighbor[ vertex[i].neighbor[j] ] = j
    if vertex j is not a neighbor of vertex i,
    vertex[i].index_neighbor[j] = -1
    n size of permutation

    """
    def __init__(self, n):

      self.neighbor = []
      self.perm = []
      self.index_neighbor = [-1]*n

class _JGraph(object):
    """
    Jerrum graph

    vertex   vertices of the Jerrum graph for permutation group G < S(n);
    vertex[i] i-th vertex, with `i` in range(n)
    jg       array of Jerrums generators (edges of the graph)
    jgs number of occupied entries of jg
    freejg   stack of slots of jg. freejg[i] points to the
             i-th free position of jg
    To a directed edge (i, j) between vertices i, j
    it is associated the index of a permutation `g` satisfying
    g[i] = j
    g = jg[vertex[i].perm[vertex[i].index_neighbor[j]]]
      = jg[vertex[j].perm[vertex[j].index_neighbor[i]]]

    cycle  list of indices of vertices used to identify cycles
    G Permutation group
    n size of permutation
    r number of generators
    """
    def __init__(self, G):
        n = G._degree
        self.vertex = [_Vertex(n) for i in range(n)]
        self.gens = [p.array_form for p in G._generators]
        self.jg = [None]*n
        self.jgs = 0
        self.freejg = []
        self.cycle = []
        self.G = G
        self.n = G._degree
        self.r = G._r
        self.idn = range(n)

    def find_cycle(self, v, v1, v2, prev):
        """
        find if there is a cycle

        v      vertex from which start searching a cycle
        v1, v2 vertices through which the cycle must go
        prev   previous vertex
        """
        cycle = self.cycle
        neighbor = self.vertex[v].neighbor
        if v1 in neighbor:
            if v1 != prev:
                return True
        if v2 in neighbor:
            if v2 != prev:
                cycle.append(v2)
                if self.find_cycle(v2, v1, v2, v):
                    return True
                cycle.pop()
        for u in neighbor:
            if u != prev:
                cycle.append(u)
                if self.find_cycle(u, v1, v2, v):
                    return True
                cycle.pop()
        return False

    def insert_edge(self, g, i, ig):
        """
        insert edge (permutation g) moving i to ig (i < ig)
        """
        vertex = self.vertex
        self.jgs += 1
        jgslot = self.freejg.pop() # the last free generator place
        self.jg[jgslot] = g
        nn = len(vertex[i].neighbor)
        vertex[i].neighbor.append(ig)
        vertex[i].perm.append(jgslot)
        vertex[i].index_neighbor[ig] = nn
        nn = len(vertex[ig].neighbor)
        vertex[ig].neighbor.append(i)
        vertex[ig].perm.append(jgslot)
        vertex[ig].index_neighbor[i] = nn

    def remove_edge(self, i, ig):
        """
        remove edge from i to ig
        """
        vertex = self.vertex
        # remove the permutation labelling this edge
        self.jgs -= 1
        jgslot = vertex[i].perm[ vertex[i].index_neighbor[ig] ]
        self.jg[jgslot] = None
        self.freejg.append(jgslot) # now we gained a free place

        for i1, i2 in ((i, ig), (ig, i)):
            v = vertex[i1]
            j0 = v.index_neighbor[i2]
            v.neighbor.pop(j0)
            v.perm.pop(j0)
            # the index of the vertices >= j0  in vertex[i] has changed
            for j in range(j0, len(v.neighbor)):
                v.index_neighbor[ v.neighbor[j] ] = j
            v.index_neighbor[ig] = -1

    def insert(self, g, alpha):
        """
        insert permutation `g` in stabilizer chain at point alpha
        """
        n = len(g)
        if not g == self.idn:
            vertex = self.vertex
            jg = self.jg
            i = _smallest_change(g, alpha)
            ig = g[i]
            nn = vertex[i].index_neighbor[ig]
            if nn >= 0:  # if ig is already neighbor of i
                jginn = jg[vertex[i].perm[nn]]
                if g != jginn:
                    # cycle consisting of two edges;
                    # replace jginn by g and insert h = g**-1*jginn
                    g1 = perm_af_invert(g)
                    h = perm_af_mul(g1, jginn)
                    jg[ vertex[i].perm[nn] ] = g
                    self.insert(h, alpha)
            else:  # new edge
                self.insert_edge(g, i, ig)
                self.cycle = [i]
                if self.find_cycle(i, i, ig, -1):
                    cycle = self.cycle
                    cycle.append(cycle[0])
                    # find the smallest point (vertex) of the cycle
                    minn = min(cycle)
                    cmin = cycle.index(minn)

                    # now walk around the cycle starting from the smallest
                    # point, and multiply around the cycle to obtain h
                    # satisfying h[cmin] = cmin
                    ap = []
                    for c in range(cmin, len(cycle)-1) + range(cmin):
                        i = cycle[c]
                        j = cycle[c+1]
                        nn = vertex[i].index_neighbor[j]
                        p = jg[ vertex[i].perm[nn] ]

                        if i > j:
                            p = perm_af_invert(p)
                        ap.append(p)
                    ap.reverse()
                    h = perm_af_muln(*ap)
                    self.remove_edge(cycle[cmin], cycle[cmin + 1])
                    self.insert(h, alpha)

    def schreier_tree(self, alpha, gen):
        """
        traversal of the orbit of alpha

        Compute a traversal of the orbit of alpha, storing the values
        in G._coset_repr; G._coset_repr[i][alpha] = i if i belongs
        to the orbit of alpha.
        """
        G = self.G
        G._coset_repr[alpha] = gen
        G._coset_repr_n += 1
        genv = self.gens[:self.r]
        h = 0
        r = self.r
        stg = [gen]
        sta = [alpha]
        pos = [0]*self.n
        while 1:
            # backtrack when finished iterating over generators
            if pos[h] >= r:
                if h == 0:
                    return
                pos[h] = 0
                h -= 1
                sta.pop()
                stg.pop()
                continue
            g = genv[pos[h]]
            pos[h] += 1
            alpha = sta[-1]
            ag = g[alpha]

            if G._coset_repr[ag] == None:
                gen1 = perm_af_mul(g, stg[-1])
                G._coset_repr[ag] = gen1
                G._coset_repr_n += 1
                sta.append(ag)
                stg.append(gen1)
                h += 1

    def jerrum_filter(self, alpha, cri):
        """
        filter the generators of the stabilizer subgroup G_alpha

        alpha point for which the stabilizer is computed

        cri[i] inverse of G._coset_repr[i] if `i` is not None

        Schreier lemma: the stabilizer subgroup G_alpha of G
        is generated by the schreier generators
        h = cosrep[ p2[i] ]**-1 * g[j] * cosrep[i]
        where j=0,..,len(gens)-1 and i=0,..,n-1, where n is the degree.
        Proof that h belongs to G_alpha:
        cosrep[k][alpha] = k for all k; cosrep[k]**-1[k] = alpha
        p1 = cosrep[i]; p2 = g[j]
        p3 = cosrep[ p2[i] ]; p3[alpha] = p2[i]
        p3**-1[p2[i] = alpha
        p3**-1[p2[p1[alpha]] = alpha, so h[alpha] = alpha

        Using Jerrum's filter one can reduce the len(gens)*n generators
        of G_alpha produced by the Schreier lemma to at most n-1

        Jerrum's filter:
        (see Cameron 'Permutation groups', page 22)
        _JGraph has n-1 vertices; the edges (i, j) are labelled by
        group elements `g` with j = imin(g) = min(i | g[i] != i);
        define m(graph) = sum(imin(g) for g in graph)

        At the beginning the graph has no edges, so it is
        an acyclic graph.
        Insert a group element `g` produced by the Schreier lemma;
        introduce in _JGraph an edge (imin(g), g[imin(g));
        if the graph contains a cycle,
        let `i0` be the smallest point in the cycle, and `h` the
        product of the group elements labelling the edges in the cycle,
        starting from `i0`; h[j] = j for j <= i0;
        modify it eliminating the edge (i0, g0[i0])
        in the cycle; one obtains a new acyclic graph with
        m(graph_new) > m(graph). `g0` can be expressed as a product
        of `h` and the other elements in the cycle.
        Then insert `h` in the graph, and so on.
        Since m < n**2, this process ends after
        a finite number of times, so in the end one remains
        with an acyclic graph, with at most n-1 edges and
        the same number of generators.
        """
        n = self.n
        r = self.r
        G = self.G
        gens = self.gens
        cosrep = G._coset_repr
        self.jgs = 0
        for j in range(n):
            self.jg[j] = None
        self.freejg = range(n)
        for i in range(n):
            self.vertex[i].neighbor = []
            self.vertex[i].perm = []
        for i in range(n):
            for j in range(n):
                self.vertex[i].index_neighbor[j] = -1

        for i in range(n):
            if cosrep[i] != None:
                p1 = cosrep[i]
                for j in range(r):
                    p2 = gens[j]
                    p3 = cri[ p2[i] ]
                    h = [p3[p2[k]] for k in p1]
                    self.insert(h, alpha)
        r = 0
        for j in range(n):
            if self.jg[j] != None:
                gens[r] = self.jg[j]
                r += 1
        self.r = r

class PermutationGroup(Basic):
    """
    The class defining a Permutation group.

    Permutation(generator_list) returns the permutation group
    generated by permutation_list.

    >>> from sympy.combinatorics.permutations import Permutation
    >>> from sympy.combinatorics.perm_groups import PermutationGroup
    >>> a = Permutation([0, 2, 1])
    >>> b = Permutation([1, 0, 2])
    >>> G = PermutationGroup([a, b])
    >>> G
    PermutationGroup([Permutation([0, 2, 1]), Permutation([1, 0, 2])])
    """

    def __new__(cls, *args, **kw_args):
        """
        The default constructor.
        """
        obj = Basic.__new__(cls, *args, **kw_args)
        obj._generators = args[0]
        obj._order = None
        obj._center = []
        obj._is_abelian = None
        obj._is_transitive = None
        obj._is_sym = None
        obj._is_alt = None
        obj._is_primitive = None
        obj._max_div = None
        size = len(args[0][0].array_form)
        obj._r = len(obj._generators)
        if not all(len(args[0][i].array_form)==size for i in xrange(1, len(args[0]))):
                raise ValueError("Permutation group size is not correct")
        obj._degree = size

        # these attributes are assigned after running schreier_sims
        obj._base = []
        obj._coset_repr = []
        obj._coset_repr_n = []
        obj._stabilizers_gens = []

        # these attributes are assigned after running _pr_init
        obj._random_gens = []
        return obj

    def __mul__(self, other):
        """
        Return the direct product of two permutation groups as a permutation
        group.

        This implementation realizes the direct product by shifting
        the index set for the generators of the second group: so if we have
        G acting on n1 points and H acting on n2 points, G*H acts on n1+n2
        points.

        Examples
        ========

        >>> from sympy.combinatorics.perm_groups import (PermutationGroup,
        ... CyclicGroup)
        >>> G = CyclicGroup(5)
        >>> H = G*G
        >>> H
        PermutationGroup([Permutation([1, 2, 3, 4, 0, 5, 6, 7, 8, 9]),
        Permutation([0, 1, 2, 3, 4, 6, 7, 8, 9, 5])])
        >>> H.order()
        25

        """
        gens1 = [perm.array_form for perm in self.generators]
        gens2 = [perm.array_form for perm in other.generators]
        n1 = self.degree
        n2 = other.degree
        start = range(n1)
        end = range(n1, n1 + n2)
        for i in range(len(gens2)):
            gens2[i] = [x + n1 for x in gens2[i]]
        gens2 = [start + gen for gen in gens2]
        gens1 = [gen + end for gen in gens1]
        together = gens1 + gens2
        gens = [_new_from_array_form(x) for x in together]
        return PermutationGroup(gens)

    @property
    def degree(self):
        """
        Returns the size of the permutations in the group.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([1,0])
        >>> G = PermutationGroup([a])
        >>> G.degree
        2

        """
        return self._degree

    @property
    def generators(self):
        """
        Returns the generators of the group in array form.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.generators
        [Permutation([0, 2, 1]), Permutation([1, 0, 2])]

        """
        return self._generators

    @property
    def is_abelian(self):
        """
        Checks if the group is Abelian.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.is_abelian
        False
        >>> a = Permutation([0, 2, 1])
        >>> G = PermutationGroup([a])
        >>> G.is_abelian
        True

        """
        if self._is_abelian is not None:
            return self._is_abelian

        self._is_abelian = True
        gens = [p.array_form for p in self.generators]
        for x in gens:
            for y in gens:
                if y <= x:
                    continue
                if not perm_af_commutes_with(x, y):
                    self._is_abelian = False
                    return False
        return True

    def order(self):
        """
        return the order of the group

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.order()
        6

        """
        if self._order != None:
            return self._order
        if self._is_sym:
            n = self.degree
            self._order = factorial(n)
            return self._order
        if self._is_alt:
            n = self.degree
            self._order = factorial(n)/2
            return self._order
        self.schreier_sims()
        m = 1
        for x in self._coset_repr_n:
            m *= x
        return m

    def coset_repr(self):
        """
        Return the Schreier-Sims representation of the group.

        The Schreier-Sims representation is the list of the cosets of
        the chain of stabilizers, see schreier_sims.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.coset_repr()
        [[[0, 1, 2], [1, 0, 2], [2, 0, 1]], [[0, 1, 2], [0, 2, 1]]]

        """
        if not self._coset_repr:
            self.schreier_sims()
        return self._coset_repr

    def stabilizers_gens(self):
        """
        Schreier-Sims stabilizers generators

        Return the generators of the stabilizers chain in the
        Schreier-Sims representation.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.stabilizers_gens()
        [[0, 2, 1]]

        """

        if not self._coset_repr:
            self.schreier_sims()
        return self._stabilizers_gens

    def strong_base(self):
        """
        strong base in Schreier-Sims representation

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.strong_base()
        [0, 1]

        """
        if not self._coset_repr:
            self.schreier_sims()
        return self._base

    def schreier_sims(self):
        """
        Schreier-Sims algorithm.

        It computes the generators of the stabilizers chain
        G > G_{b_1} > .. > G_{b1,..,b_r} > 1
        in which G_{b_1,..,b_i} stabilizes b_1,..,b_i,
        and the corresponding `s` cosets.
        An element of the group can be written univoquely
        as the product h_1*..*h_s.

        We use Jerrum's filter in our implementation of the
        Schreier-Sims algorithm. It runs in polynomial time.

        This implementation is a translation of the C++ implementation in
        http://www.m8j.net

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.schreier_sims()
        >>> G.stabilizers_gens()
        [[0, 2, 1]]
        >>> G.coset_repr()
        [[[0, 1, 2], [1, 0, 2], [2, 0, 1]], [[0, 1, 2], [0, 2, 1]]]

        """
        if self._coset_repr:
            return
        JGr = _JGraph(self)
        alpha = 0
        n = JGr.n
        self._order = 1
        coset_repr = []
        num_generators = []
        generators = []
        gen = range(n)
        base = {}
        JGr.gens += [None]*(n - len(JGr.gens))
        while 1:
            self._coset_repr_n = 0
            self._coset_repr = [None]*n
            JGr.schreier_tree(alpha, gen)
            cri = []
            for p in self._coset_repr:
                if not p:
                    cri.append(p)
                else:
                    cri.append(perm_af_invert(p))
            JGr.jerrum_filter(alpha, cri)
            if self._coset_repr_n > 1:
                base[alpha] = self._coset_repr_n
            self._order *= self._coset_repr_n
            coset_repr.append([p for p in self._coset_repr if p])
            d = {}
            for p in self._coset_repr:
                if p:
                    d[p[alpha]] = p
            num_generators.append(JGr.r)
            if JGr.r:
                generators.extend(JGr.gens[:JGr.r])
            if JGr.r <= 0:
                break
            alpha += 1
        self._coset_repr = coset_repr
        a = []
        for p in generators:
            if p not in a:
                a.append(p)
        self._stabilizers_gens = a

        i = len(JGr.gens) - 1
        while not JGr.gens[i]:
            i -= 1
        JGr.gens = JGr.gens[:i+1]
        self._base = base.keys()
        self._coset_repr_n = base.values()

    def coset_decomposition(self, g):
        """
        Decompose `g` as h_0*...*h_{len(u)}

        The Schreier-Sims coset representation u of `G`
        gives a univoque decomposition of an element `g`
        as h_0*...*h_{len(u)}, where h_i belongs to u[i]

        Output: [h_0, .., h_{len(u)}] if `g` belongs to `G`
                False otherwise

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([[0, 1, 3, 7, 6, 4], [2, 5]])
        >>> b = Permutation([[0, 1, 3, 2], [4, 5, 7, 6]])
        >>> G = PermutationGroup([a, b])
        >>> c = Permutation([[0, 1, 2, 3, 4], [5, 6, 7]])
        >>> G.coset_decomposition(c)
        False
        >>> c = Permutation([[0, 6], [1, 7], [2, 4], [3, 5]])
        >>> G.coset_decomposition(c)
        [[6, 4, 2, 0, 7, 5, 3, 1], [0, 4, 1, 5, 2, 6, 3, 7], [0, 1, 2, 3, 4, 5, 6, 7]]
        >>> G.has_element(c)
        True

        """
        u = self.coset_repr()
        if isinstance(g, Permutation):
            g = g.array_form
        g1 = g
        n = len(u)
        a = []
        for i in range(n):
            x = g1[i]
            for h in u[i]:
                if h[i] == x:
                    a.append(h)
                    p2 = perm_af_invert(h)
                    g1 = perm_af_mul(p2, g1)
                    break
            else:
                return False
        if perm_af_muln(*a) == g:
            return a
        return False

    def coset_rank(self, g):
        """
        rank using Schreier-Sims representation

        The coset rank of `g` is the ordering number in which
        it appears in the lexicographic listing according to the
        coset decomposition, see coset_decomposition;
        the ordering is the same as in G.generate(method='coset').
        If `g` does not belong to the group it returns None

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([[0, 1, 3, 7, 6, 4], [2, 5]])
        >>> b = Permutation([[0, 1, 3, 2], [4, 5, 7, 6]])
        >>> G = PermutationGroup([a, b])
        >>> c = Permutation([[0, 1, 2, 3, 4], [5, 6, 7]])
        >>> G.coset_rank(c)
        >>> c = Permutation([[0, 6], [1, 7], [2, 4], [3, 5]])
        >>> G.coset_rank(c)
        40
        >>> G.coset_unrank(40, af=True)
        [6, 7, 4, 5, 2, 3, 0, 1]

        """
        u = self.coset_repr()
        if isinstance(g, Permutation):
            g = g.array_form
        g1 = g
        m = len(u)
        a = []

        un = self._coset_repr_n
        n = self.degree
        rank = 0
        base = [1]
        for i in un[m:0:-1]:
            base.append(base[-1]*i)
        base.reverse()

        a1 = [0]*m
        i1 = -1
        for i in self._base:
            i1 += 1
            x = g1[i]
            for j, h in enumerate(u[i]):
                if h[i] == x:
                    a.append(h)
                    a1[i] = j
                    rank += j*base[i1]
                    p2 = perm_af_invert(h)
                    g1 = perm_af_mul(p2, g1)
                    break
            else:
                return None
        if perm_af_muln(*a) == g:
            return rank
        return None

    def coset_unrank(self, rank, af=False):
        """
        unrank using Schreier-Sims representation

        coset_unrank is the inverse operation of coset_rank
        if 0 <= rank < order; otherwise it returns None.

        """
        u = self.coset_repr()
        if rank < 0 or rank >= self.order():
            return None
        un = self._coset_repr_n
        base = self._base
        m = len(u)
        nb = len(base)
        assert nb == len(un)
        v = [0]*m
        for i in range(nb-1, -1,-1):
            j = base[i]
            rank, c = divmod(rank, un[i])
            v[j] = c
        a = [u[i][v[i]] for i in range(m)]
        h = perm_af_muln(*a)
        if af:
            return h
        else:
            return _new_from_array_form(h)

    def random(self, af=False):
        """
        return a random group element

        """
        rank = randrange(0, self.order())
        return self.coset_unrank(rank, af)

    def has_element(self, g):
        """
        test if `g` belongs to G; see coset_decomposition

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1, 3])
        >>> b = Permutation([0, 2, 3, 1])
        >>> g = PermutationGroup([a, b])
        >>> g.has_element(Permutation([0, 1, 3, 2]))
        True
        >>> g.has_element(Permutation([1, 2, 3, 0]))
        False

        """
        return bool(self.coset_decomposition(g.array_form))

    def generate_dimino(self, af=False):
        """
        yield group elements using Dimino's algorithm

        If af == True it yields the array form of the permutations

        Reference:
        [1] The implementation of various algorithms for Permutation Groups in
        the Computer Algebra System: AXIOM, N.J. Doye, M.Sc. Thesis

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1, 3])
        >>> b = Permutation([0, 2, 3, 1])
        >>> g = PermutationGroup([a, b])
        >>> list(g.generate_dimino(af=True))
        [[0, 1, 2, 3], [0, 2, 1, 3], [0, 2, 3, 1], [0, 1, 3, 2], [0, 3, 2, 1], [0, 3, 1, 2]]

        """
        idn = range(self.degree)
        order = 0
        element_list = [idn]
        set_element_list = set([tuple(idn)])
        if af:
            yield idn
        else:
            yield _new_from_array_form(idn)
        gens = [p.array_form for p in self.generators]

        for i in xrange(len(gens)):
            # D elements of the subgroup G_i generated by gens[:i]
            D = element_list[:]
            N = [idn]
            while N:
                A = N
                N = []
                for a in A:
                    for g in gens[:i+1]:
                        ag = perm_af_mul(a, g)
                        if tuple(ag) not in set_element_list:
                            # produce G_i*g
                            for d in D:
                                order += 1
                                ap = perm_af_mul(d, ag)
                                if af:
                                    yield ap
                                else:
                                    p = _new_from_array_form(ap)
                                    yield p
                                element_list.append(ap)
                                set_element_list.add(tuple(ap))
                                N.append(ap)
        self._order = len(element_list)

    def generate_schreier_sims(self, af=False):
        """
        yield group elements using the Schreier-Sims representation

        If af = True it yields the array form of the permutations

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1, 3])
        >>> b = Permutation([0, 2, 3, 1])
        >>> g = PermutationGroup([a, b])
        >>> list(g.generate_schreier_sims(af=True))
        [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 3, 1], [0, 2, 1, 3], [0, 3, 2, 1], [0, 3, 1, 2]]

        """
        def get1(posmax):
            n = len(posmax) - 1
            for i in range(n,-1,-1):
                if posmax[i] != 1:
                    return i + 1
        n = self.degree
        u = self.coset_repr()
        # stg stack of group elements
        stg = [range(n)]
        # posmax[i] = len(u[i])
        posmax = [len(x) for x in u]
        n1 = get1(posmax)
        pos = [0]*n1
        posmax = posmax[:n1]
        h = 0
        while 1:
            # backtrack when finished iterating over coset
            if pos[h] >= posmax[h]:
                if h == 0:
                    raise StopIteration
                pos[h] = 0
                h -= 1
                stg.pop()
                continue
            p = perm_af_mul(stg[-1], u[h][pos[h]])
            pos[h] += 1
            stg.append(p)
            h += 1
            if h == n1:
                if af:
                    yield p
                else:
                    p1 = _new_from_array_form(p)
                    yield p1
                stg.pop()
                h -= 1

    def generate(self, method="coset", af=False):
        """
        return iterator to generate the elements of the group

        Iteration is done with one of these methods:
          method='coset'  using the Schreier-Sims coset representation
          method='dimino' using the Dimino method

        If af = True it yields the array form of the permutations

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1, 3])
        >>> b = Permutation([0, 2, 3, 1])
        >>> g = PermutationGroup([a, b])
        >>> list(g.generate(af=True))
        [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 3, 1], [0, 2, 1, 3], [0, 3, 2, 1], [0, 3, 1, 2]]

        """
        if method == "coset":
            return self.generate_schreier_sims(af)
        elif method == "dimino":
            return self.generate_dimino(af)
        else:
            raise ValueError('there is not this method')

    def orbit(self, alpha):
        r"""
        Compute the orbit of alpha `\{g(\alpha) | g \in G\}` as a set.

        The time complexity of the algorithm used here is `O(|Orb|*r)` where
        `|Orb|` is the size of the orbit and `r` is the number of generators
        of the group. For a proof of correctness, see [1], p.78.

        Examples
        ========

        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([1,2,0,4,5,6,3])
        >>> G = PermutationGroup([a])
        >>> G.orbit(0)
        set([0, 1, 2])
        >>> G.orbit(4)
        set([3, 4, 5, 6])

        See Also
        ========

        orbit_transversal

        References
        ==========

        [1] Holt, D., Eick, B., O'Brien, E.
        "Handbook of computational group theory"

        """
        orb = [alpha]
        used = [False]*self.degree
        used[alpha] = True
        gens = self.generators
        for b in orb:
            for gen in gens:
                temp = gen(b)
                if used[temp] == False:
                    orb.append(temp)
                    used[temp] = True
        return set(orb)

    def orbit_transversal(self, alpha, pairs=False):
        r"""
        Computes a transversal for the orbit of ``alpha`` as a set.

        For a permutation group `G`, a transversal for the orbit
        `Orb = \{g(\alpha) | g \in G\}` is a set
        `\{g_\beta | g_\beta(\alpha) = \beta\}` for `\beta \in Orb`.
        Note that there may be more than one possible transversal.
        If ``pairs`` is set to ``True``, it returns the list of pairs
        `(\beta, g_\beta)`. For a proof of correctness, see [1], p.79

        Examples
        ========

        >>> from sympy.combinatorics.perm_groups import (PermutationGroup,
        ... DihedralGroup)
        >>> G = DihedralGroup(6)
        >>> G.orbit_transversal(0)
        [Permutation([0, 1, 2, 3, 4, 5]), Permutation([1, 2, 3, 4, 5, 0]),
        Permutation([5, 4, 3, 2, 1, 0]), Permutation([2, 3, 4, 5, 0, 1]),
        Permutation([4, 3, 2, 1, 0, 5]), Permutation([3, 4, 5, 0, 1, 2])]

        See Also
        ========

        orbit

        References
        ==========

        [1] Holt, D., Eick, B., O'Brien, E.
        "Handbook of computational group theory"

        """
        n = self.degree
        tr = [(alpha, _new_from_array_form(range(n)))]
        used = [False]*n
        used[alpha] = True
        gens = self.generators
        for pair in tr:
            for gen in gens:
                temp = gen(pair[0])
                if used[temp] == False:
                    tr.append((temp, gen*pair[1]))
                    used[temp] = True
        if pairs:
            return tr
        return [pair[1] for pair in tr]

    def orbits(self, rep=False):
        """
        compute the orbits of G;
        if rep=False it returns a list of sets
        else it returns a list of representatives of the orbits

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.orbits()
        [set([0, 1, 2])]
        >>> G.orbits(rep=True)
        [0]

        """
        n = self._degree
        s1 = set(range(n))
        orbs = []
        while s1:
            i = s1.pop()
            si = self.orbit(i)
            if rep:
                orbs.append(i)
            else:
                orbs.append(si)
            s1 -= si
        return orbs

    @property
    def is_transitive(self):
        """
        test if the group is transitive

        A group is transitive if it has a single orbit.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1, 3])
        >>> b = Permutation([2, 0, 1, 3])
        >>> G1 = PermutationGroup([a, b])
        >>> G1.is_transitive
        False
        >>> c = Permutation([2, 3, 0, 1])
        >>> G2 = PermutationGroup([a, c])
        >>> G2.is_transitive
        True

        """
        if self._is_transitive is not None:
            return self._is_transitive

        ans = len(self.orbit(0)) == self.degree
        self._is_transitive = ans
        return ans

    def stabilizer(self, alpha):
        r"""
        Returns the stabilizer subgroup of ``alpha``.

        The stabilizer of `\alpha` is the group `G_\alpha =
        \{g \in G | g(\alpha) = \alpha\}`.
        For a proof of correctness, see [1], p.79.

        Examples
        ========

        >>> from sympy.combinatorics.perm_groups import (PermutationGroup,
        ... DihedralGroup)
        >>> G = DihedralGroup(6)
        >>> G.stabilizer(5)
        PermutationGroup([Permutation([4, 3, 2, 1, 0, 5]),
        Permutation([0, 1, 2, 3, 4, 5])])

        See Also
        ========

        orbit

        References
        ==========

        [1] Holt, D., Eick, B., O'Brien, E.
        "Handbook of computational group theory"

        """
        n = self.degree
        orb = [alpha]
        table = {alpha: _new_from_array_form(range(n))}
        used = [False]*n
        used[alpha] = True
        gens = self.generators
        stab_gens = []
        for b in orb:
            for gen in gens:
                temp = gen(b)
                if used[temp] == False:
                    gen_temp = gen*table[b]
                    orb.append(temp)
                    table[temp] = gen_temp
                    used[temp] = True
                else:
                    schreier_gen = (~table[temp])*gen*table[b]
                    if schreier_gen not in stab_gens:
                        stab_gens.append(schreier_gen)
        return PermutationGroup(list(stab_gens))

    def schreier_vector(self, alpha):
        """
        Computes the schreier vector for ``alpha``.

        The Schreier vector efficiently stores information
        about the orbit of ``alpha``. It can later be used to quickly obtain
        elements of the group that send ``alpha`` to a particular element
        in the orbit. Notice that the Schreier vector depends on the order
        in which the group generators are listed. For a definition, see [1].
        Since list indices start from zero, we adopt the convention to use
        "None" instead of 0 to signify that an element doesn't belong
        to the orbit.
        For the algorithm and its correctness, see [2], pp.78-80.

        Examples
        ========

        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([2,4,6,3,1,5,0])
        >>> b = Permutation([0,1,3,5,4,6,2])
        >>> G = PermutationGroup([a,b])
        >>> G.schreier_vector(0)
        [-1, None, 0, 1, None, 1, 0]

        See Also
        ========

        orbit

        References
        ==========

        [1] http://en.wikipedia.org/wiki/Schreier_vector

        [2] Holt, D., Eick, B., O'Brien, E.
        "Handbook of computational group theory"

        """
        n = self.degree
        v = [None]*n
        v[alpha] = -1
        orb = [alpha]
        used = [False]*n
        used[alpha] = True
        gens = self.generators
        r = len(gens)
        for b in orb:
            for i in range(r):
                temp = gens[i](b)
                if used[temp] == False:
                    orb.append(temp)
                    used[temp] = True
                    v[temp] = i
        return v

    def orbit_rep(self, alpha, beta):
        """
        Return a group element which sends ``alpha`` to ``beta``.

        If ``beta`` is not in the orbit of ``alpha``, the function returns
        ``False``. This implementation makes use of the schreier vector.
        For a proof of correctness, see [1], p.80

        Examples
        ========

        >>> from sympy.combinatorics.perm_groups import (PermutationGroup,
        ... AlternatingGroup)
        >>> G = AlternatingGroup(5)
        >>> G.orbit_rep(0,4)
        Permutation([4, 2, 3, 0, 1])

        See Also
        ========

        schreier_vector

        References
        ==========

        [1] Holt, D., Eick, B., O'Brien, E.
        "Handbook of computational group theory"

        """
        v = self.schreier_vector(alpha)
        if v[beta] == None:
            return False
        n = self.degree
        u = _new_from_array_form(range(n))
        k = v[beta]
        gens = self.generators
        while k != -1:
            u = u*gens[k]
            beta = (~gens[k])(beta)
            k = v[beta]
        return u

    def is_normal(self, gr):
        """
        test if G=self is a normal subgroup of gr

        G is normal in gr if
        for each g2 in G, g1 in gr, g = g1*g2*g1**-1 belongs to G
        It is sufficient to check this for each g1 in gr.generator and
        g2 g2 in G.generator

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([1, 2, 0])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G1 = PermutationGroup([a, Permutation([2, 0, 1])])
        >>> G1.is_normal(G)
        True

        """
        gens2 = [p.array_form for p in self.generators]
        gens1 = [p.array_form for p in gr.generators]
        for g1 in gens1:
            for g2 in gens2:
                p = perm_af_muln(g1, g2, perm_af_invert(g1))
                if not self.coset_decomposition(p):
                    return False
        return True

    def normal_closure(self, gens):
        """
        normal closure in self of a list gens2 of permutations

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([1, 2, 0])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.order()
        6
        >>> G1 = G.normal_closure([a])
        >>> list(G1.generate(af=True))
        [[0, 1, 2], [1, 2, 0], [2, 0, 1]]

        """
        G2 = PermutationGroup(gens)
        if G2.is_normal(self):
            return G2
        gens1 = [p.array_form for p in self.generators]
        b = 0
        while not b:
            gens2 = [p.array_form for p in G2.generators]
            b = 1
            for g1 in gens1:
                if not b:
                    break
                for g2 in gens2:
                    p = perm_af_muln(g1, g2, perm_af_invert(g1))
                    p = Permutation(p)
                    if not G2.has_element(p):
                        gens2 = G2.generators + [p]
                        G2 = PermutationGroup(gens2)
                        b = 0
                        break

        return G2

    def __eq__(self, gr):
        """
        test if two groups have the same elements

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = [[1,2,0,3,4,5], [1,0,2,3,4,5], [2,1,0,3,4,5], [1,2,0,3,4,5]]
        >>> a = [Permutation(p) for p in a]
        >>> g = Permutation([1,2,3,4,5,0])
        >>> G1,G2,G3 = [PermutationGroup(x) for x in [a[:2],a[2:4],[g,g**2]]]
        >>> assert G1.order() == G2.order() == G3.order() == 6
        >>> assert G1 == G2 and G1 != G3

        """
        if self.degree != gr.degree:
            return False
        if self.order() != gr.order():
            return False
        gens1 = self.generators
        for g in gens1:
            if not gr.has_element(g):
                return False
        return True

    def is_subgroup(self, gr):
        """
        test if self is a subgroup of gr

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([1,2,3,4,0])
        >>> b = Permutation([1,0,2,3,4])
        >>> G = PermutationGroup([a, b])
        >>> c = Permutation([1,0,3,2,4])
        >>> G1 = PermutationGroup([a, c])
        >>> G1.is_subgroup(G)
        True

        """
        if self.degree != gr.degree:
            return False
        if self.order() > gr.order():
            return False
        gens1 = self.generators
        for g in gens1:
            if not gr.has_element(g):
                return False
        return True


    def __ne__(self, gr):
        return not self == gr

    def commutator(self):
        """
        commutator subgroup

        The commutator subgroup is the subgroup generated by all
        commutators; it is equal to the normal closure of the set
        of commutators of the generators.

        see http://groupprops.subwiki.org/wiki/Derived_subgroup

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([1, 0, 2, 4, 3])
        >>> b = Permutation([0, 1, 3, 2, 4])
        >>> G = PermutationGroup([a, b])
        >>> C = G.commutator()
        >>> list(C.generate(af=True))
        [[0, 1, 2, 3, 4], [0, 1, 3, 4, 2], [0, 1, 4, 2, 3]]

        """
        r = self._r
        gens = [p.array_form for p in self.generators]
        gens_inv = [perm_af_invert(p) for p in gens]
        set_commutators = set()
        for i in range(r):
            for j in range(r):
                p1 = gens[i]
                p1inv = gens_inv[i]
                p2 = gens[j]
                p2inv = gens_inv[j]
                c = [p1[p2[p1inv[k]]] for k in p2inv]
                ct = tuple(c)
                if not ct in set_commutators:
                    set_commutators.add(ct)
        cms = [Permutation(p) for p in set_commutators]
        G2 = self.normal_closure(cms)
        return G2

    def is_solvable(self):
        """
        test if the group G is solvable

        G is solvable if the derived series
        G = G_0 < G_1 < ... < G_k = 1, with G_{i+1} = G.commutator()
        see http://en.wikipedia.org/wiki/Solvable_group

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([1,2,0])
        >>> b = Permutation([1,0,2])
        >>> G = PermutationGroup([a, b])
        >>> G.is_solvable()
        True

        """
        order = self.order()
        if order == 1:
            return True
        G = self
        while order > 1:

            G = G.commutator()
            order1 = G.order()
            if order1 == order:
                return False
            order = order1
        return True

    def _pr_init(self, r, n):
        deg = self.degree
        random_gens = self.generators[:]
        k = len(random_gens)
        if k < r:
            for i in range(k, r):
                random_gens.append(random_gens[i - k])
        acc = _new_from_array_form(range(deg))
        random_gens.append(acc)
        self._random_gens = random_gens
        for i in range(n):
            self.pr_random()

    def pr_random(self):
        if self._random_gens == []:
            self._pr_init(11, 50)
        random_gens = self._random_gens
        r = len(random_gens) - 1
        s = randrange(r)
        t = randrange(r - 1)
        if t == s:
            t = r - 1
        x = choice([1, 2])
        e = choice([-1, 1])
        if x == 1:
            random_gens[s] = random_gens[s]*(random_gens[t]**e)
            random_gens[r] = random_gens[r]*random_gens[s]
        else:
            random_gens[s] = (random_gens[t]**e)*random_gens[s]
            random_gens[r] = random_gens[s]*random_gens[r]
        return random_gens[r]

    def is_alt_sym(self, eps = 0.05):
        n = self.degree
        if n < 8:
            return False
        if not self.is_transitive:
            return False
        if n < 17:
            c_n = 0.34
        else:
            c_n = 0.57
        d_n = (c_n*log(2))/log(n)
        N_eps = int(-log(eps)/d_n)
        for i in range(N_eps):
            perm = self.pr_random()
            if _check_cycles_alt_sym(perm):
                return True
        return False

    def alt_or_sym(self):
        if self.is_alt_sym():
            gens = self.generators
            for perm in gens:
                if perm.is_odd:
                    self._is_sym = True
                    return 'S'
            self._is_alt = True
            return 'A'
        return False

    def minimal_block(self, points):
        if not self.is_transitive:
            return False
        n = self.degree
        gens = self.generators
        parents = range(n)
        ranks = [1]*n
        not_rep = []
        max_rank = self._max_div
        k = len(points)
        if k > max_rank:
            return False
        for i in xrange(k-1):
            parents[points[i+1]] = points[0]
            not_rep.append(points[i+1])
        ranks[points[0]] = k
        i = 0
        len_not_rep = k-1
        while i < len_not_rep:
            temp = not_rep[i]
            i += 1
            for gen in gens:
                delta = self._union_find_rep(temp, parents)
                len_not_rep += self._union_find_merge(gen(temp), gen(delta), ranks,\
                                                 parents, not_rep)
        for i in xrange(n):
            self._union_find_rep(i, parents)
        return parents

    def _union_find_rep(self, num, parents):
        rep, parent = num, parents[num]
        while parent != rep:
            rep = parent
            parent = parents[rep]
        #path compression
        temp, parent = num, parents[num]
        while parent != rep:
            parents[temp] = rep
            temp = parent
            parent = parents[temp]
        return rep

    def _union_find_merge(self, first, second, ranks, parents, not_rep):
        rep_first = self._union_find_rep(first, parents)
        rep_second = self._union_find_rep(second, parents)
        if rep_first != rep_second:
            #union by rank
            if ranks[rep_first] >= ranks[rep_second]:
                new1, new2 = rep_first, rep_second
            else:
                new1, new2 = rep_second, rep_first
            parents[new2] = new1
            ranks[new1] = ranks[new1] + ranks[new2]
            not_rep.append(new2)
            return True
        return False

    @property
    def max_div(self):
        if self._max_div != None:
            return self._max_div
        n = self.degree
        for x in sieve:
            if n % x == 0:
                d = n/x
                self._max_div = d
                return d

    @property
    def is_primitive(self):
        if self._is_primitive != None:
            return self._is_primitive
        stab = self.stabilizer(0)
        orbs = stab.orbits()


def SymmetricGroup(n):
    """
    Generates the symmetric group on ``n`` elements as a permutation group.

    The generators taken are the ``n``-cycle
    ``(0 1 2 ... n-1)`` and the transposition ``(0 1)`` (in cycle notation).
    (See [1]). After the group is generated, some of its basic properties
    are set.

    Examples
    ========

    >>> from sympy.combinatorics.perm_groups import SymmetricGroup
    >>> G = SymmetricGroup(4)
    >>> G.order()
    24
    >>> list(G.generate_schreier_sims(af=True))
    [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2],
    [0, 3, 2, 1], [1, 2, 3, 0], [1, 2, 0, 3], [1, 3, 2, 0],
    [1, 3, 0, 2], [1, 0, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [2, 3, 1, 0],
    [2, 0, 3, 1], [2, 0, 1, 3], [2, 1, 3, 0], [2, 1, 0, 3], [3, 0, 1, 2],
    [3, 0, 2, 1], [3, 1, 0, 2], [3, 1, 2, 0], [3, 2, 0, 1], [3, 2, 1, 0]]

    See Also
    ========

    CyclicGroup, DihedralGroup, AlternatingGroup

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Symmetric_group#Generators_and_relations

    """
    if n == 1:
        G = PermutationGroup([Permutation([0])])
    elif n == 2:
        G = PermutationGroup([Permutation([1, 0])])
    else:
        a = range(1,n)
        a.append(0)
        gen1 = _new_from_array_form(a)
        a = range(n)
        a[0], a[1] = a[1], a[0]
        gen2 = _new_from_array_form(a)
        G = PermutationGroup([gen1, gen2])

    if n<3:
        G._is_abelian = True
    else:
        G._is_abelian = False
    G._degree = n
    G._is_transitive = True
    G._is_sym = True
    return G

def CyclicGroup(n):
    """
    Generates the cyclic group of order ``n`` as a permutation group.

    The generator taken is the ``n``-cycle ``(0 1 2 ... n-1)``
    (in cycle notation). After the group is generated, some of its basic
    properties are set.

    Examples
    ========

    >>> from sympy.combinatorics.perm_groups import CyclicGroup
    >>> G = CyclicGroup(6)
    >>> G.order()
    6
    >>> list(G.generate_schreier_sims(af=True))
    [[0, 1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 0], [2, 3, 4, 5, 0, 1],
    [3, 4, 5, 0, 1, 2], [4, 5, 0, 1, 2, 3], [5, 0, 1, 2, 3, 4]]

    See Also
    ========

    SymmetricGroup, DihedralGroup, AlternatingGroup

    """
    a = range(1, n)
    a.append(0)
    gen = _new_from_array_form(a)
    G = PermutationGroup([gen])

    G._is_abelian = True
    G._degree = n
    G._is_transitive = True
    G._order = n
    return G

def DihedralGroup(n):
    r"""
    Generates the dihedral group `D_n` as a permutation group.

    The dihedral group `D_n` is the group of symmetries of the regular
    ``n``-gon. The generators taken are the ``n``-cycle ``a = (0 1 2 ... n-1)``
    (a rotation of the ``n``-gon) and ``b = (0 n-1)(1 n-2)...``
    (a reflection of the ``n``-gon) in cycle rotation. It is easy to see that
    these satisfy ``a**n = b**2 = 1`` and ``bab = ~a`` so they indeed generate
    `D_n` (See [1]). After the group is generated, some of its basic properties
    are set.

    Examples
    ========

    >>> from sympy.combinatorics.perm_groups import DihedralGroup
    >>> G = DihedralGroup(5)
    >>> a = list(G.generate_dimino())
    >>> [perm.cyclic_form for perm in a]
    [[[4], [3], [2], [1], [0]], [[0, 1, 2, 3, 4]], [[0, 2, 4, 1, 3]],
    [[0, 3, 1, 4, 2]], [[0, 4, 3, 2, 1]], [[2], [1, 3], [0, 4]],
    [[2, 3], [1, 4], [0]], [[3], [2, 4], [0, 1]], [[3, 4], [1], [0, 2]],
    [[4], [1, 2], [0, 3]]]

    See Also
    ========

    SymmetricGroup, CyclicGroup, AlternatingGroup

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Dihedral_group

    """
    # small cases are special
    if n == 1:
        return PermutationGroup([Permutation([1, 0])])
    if n == 2:
        return PermutationGroup([Permutation([1, 0, 3, 2]),
               Permutation([2, 3, 0, 1]), Permutation([3, 2, 1, 0])])

    a = range(1, n)
    a.append(0)
    gen1 = _new_from_array_form(a)
    a = range(n)
    a.reverse()
    gen2 = _new_from_array_form(a)
    G = PermutationGroup([gen1, gen2])

    G._is_abelian = False
    G._degree = n
    G._is_transitive = True
    G._order = 2*n
    return G

def AlternatingGroup(n):
    """
    Generates the alternating group on ``n`` elements as a permutation group.

    For ``n > 2``, the generators taken are ``(0 1 2), (0 1 2 ... n-1)`` for
    ``n`` odd
    and ``(0 1 2), (1 2 ... n-1)`` for ``n`` even (See [1], p.31, ex.6.9.).
    After the group is generated, some of its basic properties are set.
    The cases ``n = 1, 2`` are handled separately.

    Examples
    ========

    >>> from sympy.combinatorics.perm_groups import AlternatingGroup
    >>> G = AlternatingGroup(4)
    >>> a = list(G.generate_dimino())
    >>> len(a)
    12
    >>> [perm.is_even for perm in a]
    [True, True, True, True, True, True, True, True, True, True, True, True]

    See Also
    ========

    SymmetricGroup, CyclicGroup, DihedralGroup

    References
    ==========

    [1] Armstrong, M. "Groups and Symmetry"

    """
    # small cases are special
    if n == 1 or n == 2:
        return PermutationGroup([Permutation([0])])

    a = range(n)
    a[0], a[1], a[2] = a[1], a[2], a[0]
    gen1 = _new_from_array_form(a)
    if n % 2 == 1:
        a = range(1, n)
        a.append(0)
        gen2 = _new_from_array_form(a)
    else:
        a = range(2, n)
        a.append(1)
        gen2 = _new_from_array_form([0] + a)
    G = PermutationGroup([gen1, gen2])

    if n<4:
        G._is_abelian = True
    else:
        G._is_abelian = False
    G._degree = n
    G._is_transitive = True
    G._is_alt = True
    return G

def _check_cycles_alt_sym(perm):
    n = perm.size
    af = perm.array_form
    current_len = 0
    total_len = 0
    used = set()
    for i in xrange(n/2):
        if not i in used and i < n/2 - total_len:
            current_len = 1
            used.add(i)
            j = i
            while(af[j] != i):
                current_len += 1
                j = af[j]
                used.add(j)
            total_len += current_len
            if current_len > n/2 and current_len < n-2 and isprime(current_len):
                return True
    return False
