from sympy.combinatorics import Permutation
from sympy.core import Basic
from sympy.combinatorics.permutations import perm_af_mul, \
 _new_from_array_form, perm_af_commutes_with, perm_af_invert, perm_af_muln
from random import randint

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
        return obj

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
        rank = randint(0, self.order())
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
        """
        compute the orbit {g[i] for g in G}

        It returns the orbit as a set.

        Examples
        ========

        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([2, 0, 1])
        >>> b = Permutation([2, 1, 0])
        >>> g = PermutationGroup([a, b])
        >>> g.orbit(0)
        set([0, 1, 2])
        >>> g.orbits()
        [set([0, 1, 2])]
        >>> g.orbits(rep=True)
        [0]
        """
        genv = [p.array_form for p in self.generators]
        n = self._degree
        orb = set([alpha])
        h = 0
        r = len(genv)
        stg = [range(n)]
        sta = [alpha]
        pos = [0]*n
        while 1:
            # backtrack when finished iterating over generators
            if pos[h] >= r:
                if h == 0:
                    return orb
                pos[h] = 0
                h -= 1
                sta.pop()
                stg.pop()
                continue
            g = genv[pos[h]]
            pos[h] += 1
            alpha = sta[-1]
            ag = g[alpha]

            if ag not in orb:
                gen1 = perm_af_mul(g, stg[-1])
                orb.add(ag)
                sta.append(ag)
                stg.append(gen1)
                h += 1

    def orbit_traversal(self, alpha, af=False):
        """
        compute the orbit traversal

        Output: list of group elements; applying each to alpha
        one gets the orbit of alpha

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1])
        >>> b = Permutation([1, 0, 2])
        >>> G = PermutationGroup([a, b])
        >>> G.orbit_traversal(1, af=True)
        [[1, 0, 2], [0, 1, 2], [0, 2, 1]]
        """
        genv = [p.array_form for p in self.generators]
        n = self._degree
        coset_repr = [None]*n
        if af:
            coset_repr[alpha] = range(n)
        else:
            coset_repr[alpha] = Permutation(range(n))
        h = 0
        r = len(genv)
        stg = [range(n)]
        sta = [alpha]
        pos = [0]*n
        while 1:
            # backtrack when finished iterating over generators
            if pos[h] >= r:
                if h == 0:
                    return coset_repr
                pos[h] = 0
                h -= 1
                sta.pop()
                stg.pop()
                continue
            g = genv[pos[h]]
            pos[h] += 1
            alpha = sta[-1]
            ag = g[alpha]

            if coset_repr[ag] == None:
                gen1 = perm_af_mul(g, stg[-1])
                if af:
                    coset_repr[ag] = gen1
                else:
                    coset_repr[ag] = Permutation(gen1)
                sta.append(ag)
                stg.append(gen1)
                h += 1

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
        >>> G1.is_transitive()
        False
        >>> c = Permutation([2, 3, 0, 1])
        >>> G2 = PermutationGroup([a, c])
        >>> G2.is_transitive()
        True
        """
        return len(self.orbit(0)) == self.degree

    def stabilizer(self, alpha):
        """
        return the stabilizer subgroup leaving alpha fixed.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1, 3])
        >>> b = Permutation([2, 0, 1, 3])
        >>> g = PermutationGroup([a, b])
        >>> g2 = g.stabilizer(2)
        >>> g2
        PermutationGroup([Permutation([1, 0, 2, 3])])
        """
        if alpha == 0:
            self.schreier_sims()
            gens = self._stabilizers_gens
            gens = [_new_from_array_form(p) for p in gens]
            return PermutationGroup(gens)
        # h[alpha] = 0; h[0] = alpha
        n = self.degree
        h = range(n)
        h[0] = alpha
        h[alpha] = 0
        genv = [p.array_form for p in self.generators]
        # conjugate the group to breng alpha to 0
        gens = [_new_from_array_form(perm_af_muln(h, p, h)) for p in genv]
        G = PermutationGroup(gens)
        G.schreier_sims()
        # stabilizers for 0
        gens1 = G._stabilizers_gens
        # conjugate the group
        gens2 = [_new_from_array_form(perm_af_muln(h, p, h)) for p in gens1]
        return PermutationGroup(gens2)


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
