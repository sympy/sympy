from __future__ import annotations
from collections import defaultdict
import itertools
from sympy.combinatorics.coset_table import modified_coset_enumeration_r
from sympy.combinatorics.fp_groups import FpGroup, FpSubgroup, simplify_presentation
from sympy.combinatorics.free_groups import FreeGroup
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.core.intfunc import igcd
from sympy.functions.combinatorial.numbers import totient
from sympy.core.singleton import S

class GroupHomomorphism:
    '''
    A class representing group homomorphisms. Instantiate using `homomorphism()`.

    References
    ==========

    .. [1] Holt, D., Eick, B. and O'Brien, E. (2005). Handbook of computational group theory.

    '''

    def __init__(self, domain, codomain, images):
        self.domain = domain
        self.codomain = codomain
        self.images = images
        self._inverses = None
        self._kernel = None
        self._image = None
        self._surjective_cert = None
        self._surjective_cert_computed = False
        self._factorization = None

    def _invs(self):
        '''
        Return a dictionary with `{gen: inverse}` where `gen` is a rewriting
        generator of `codomain` (e.g. strong generator for permutation groups)
        and `inverse` is an element of its preimage

        '''
        image = self.image()
        inverses = {}
        for k in list(self.images.keys()):
            v = self.images[k]
            if not (v in inverses
                    or v.is_identity):
                inverses[v] = k
        if isinstance(self.codomain, PermutationGroup):
            gens = image.strong_gens
        else:
            gens = image.generators
        dtype = type(self.domain.identity)
        for g in gens:
            if g in inverses or g.is_identity:
                continue
            if isinstance(self.codomain, PermutationGroup):
                parts = image._strong_gens_slp[g][::-1]
            else:
                parts = g
            w = dtype.prod(
                inverses[s] if s in inverses else inverses[s**-1]**-1
                for s in parts
            )
            inverses[g] = w

        return inverses

    def invert(self, g):
        '''
        Return an element of the preimage of ``g`` or of each element
        of ``g`` if ``g`` is a list.

        Explanation
        ===========

        If the codomain is an FpGroup, the inverse for equal
        elements might not always be the same unless the FpGroup's
        rewriting system is confluent. However, making a system
        confluent can be time-consuming. If it's important, try
        `self.codomain.make_confluent()` first.

        '''
        from sympy.combinatorics import Permutation
        from sympy.combinatorics.free_groups import FreeGroupElement
        if isinstance(g, (Permutation, FreeGroupElement)):
            if isinstance(self.codomain, FpGroup):
                g = self.codomain.reduce(g)
            if self._inverses is None:
                self._inverses = self._invs()
            image = self.image()
            w = self.domain.identity
            if isinstance(self.codomain, PermutationGroup):
                if g not in image:
                    raise ValueError("Given element is not in the image of the homomorphism.")
                gens = image.generator_product(g)[::-1]
                dtype = type(w)
                return dtype.prod(
                    self._inverses[s] if s in self._inverses
                    else self._inverses[s**-1]**-1
                    for s in gens if not s.is_identity
                )

            else:
                current_coset = 0
                im_gens = list(self._inverses)
                C = modified_coset_enumeration_r(self.codomain, im_gens)
                dtype = type(w)
                factors = []
                preimages = list(self._inverses.values())
                temp_symbols = [gen.array_form[0][0] for gen in C._grp.generators]
                transl_map = dict(zip(temp_symbols, preimages))
                symbol_to_group_gen = {gen.array_form[0][0]: gen for gen in self.codomain.generators}
                for s,exp in g.array_form:
                    t=symbol_to_group_gen[s]
                    if exp<0:
                        j = C.A_dict[t**-1]
                    else:
                        j = C.A_dict[t]
                    for _ in range(abs(exp)):
                        word_intermediate = C.P[current_coset][j]
                        if word_intermediate is not None:
                            factors.extend(
                                transl_map[s]**exp
                                for s, exp in word_intermediate.array_form
                            )
                        current_coset = C.table[current_coset][j]
                if current_coset != 0:
                    raise ValueError("Given element is not in the image of the homomorphism.")
                return dtype.prod(factors)

        elif isinstance(g, list):
            return [self.invert(e) for e in g]

    def kernel(self):
        '''
        Compute the kernel of `self`.

        '''
        if self._kernel is None:
            self._kernel = self._compute_kernel()
        return self._kernel

    def _compute_kernel(self):
        G = self.domain
        H = self.codomain
        Ginf = G.order() is S.Infinity
        if Ginf and isinstance(H, (FpGroup, FreeGroup)):
            # TODO: Make this work for PermutationGroup codomain by replacing it
            # with FpGroup using PermutationGroup.presentation()
            if isinstance(H, FreeGroup):
                H = FpGroup(H, [])
            preimages = self._is_surjective_cert()
            if preimages is not None:
                symbol_to_preimage = {
                    g.ext_rep[0]: preimages[g] for g in H.generators
                }

                def _lift_to_domain(word):
                    dtype = type(G.identity)
                    return dtype.prod(
                        symbol_to_preimage[symbol]**power
                        for symbol, power in word.array_form
                    )

                kernel_gens = []
                for relator in H.relators:
                    word = _lift_to_domain(relator)
                    if not word.is_identity:
                        kernel_gens.append(word)

                for generator in G.generators:
                    lifted_image = _lift_to_domain(self.images[generator])
                    word = generator * lifted_image**-1
                    if not word.is_identity:
                        kernel_gens.append(word)

                return FpSubgroup(G, kernel_gens, normal=True)
            else:
                return self.factor()[0].kernel()

        if Ginf:
            raise NotImplementedError(
                "Kernel computation is not implemented for this group type")

        if isinstance(G, FpGroup):
            P, T = G._to_perm_group()
            perm_images = {g: self(T.invert(g)) for g in P.generators}
            perm_hom = GroupHomomorphism(P, H, perm_images)
            perm_kernel = perm_hom._compute_kernel()
            kernel_gens = T.invert(perm_kernel.generators)
            return FpSubgroup(G, kernel_gens, normal=True)

        if isinstance(G, PermutationGroup):
            return G.subgroup_search(lambda g: self(g).is_identity)

        raise NotImplementedError(
            "Kernel computation is not implemented for this group type")

    def image(self):
        '''
        Compute the image of `self`.

        '''
        if self._image is None:
            values = list(set(self.images.values()))
            if isinstance(self.codomain, PermutationGroup):
                self._image = self.codomain.subgroup(values)
            else:
                self._image = FpSubgroup(self.codomain, values)
        return self._image

    def _apply(self, elem):
        '''
        Apply `self` to `elem`.

        '''
        if elem not in self.domain:
            if isinstance(elem, (list, tuple)):
                return [self._apply(e) for e in elem]
            raise ValueError("The supplied element does not belong to the domain")
        if elem.is_identity:
            return self.codomain.identity
        else:
            images = self.images
            value = self.codomain.identity
            if isinstance(self.domain, PermutationGroup):
                gens = self.domain.generator_product(elem, original=True)
                for g in gens:
                    if g in self.images:
                        value = images[g]*value
                    else:
                        value = images[g**-1]**-1*value
            else:
                i = 0
                for _, p in elem.array_form:
                    if p < 0:
                        g = elem[i]**-1
                    else:
                        g = elem[i]
                    value = value*images[g]**p
                    i += abs(p)
        return value

    def __call__(self, elem):
        return self._apply(elem)

    def is_injective(self):
        '''
        Check if the homomorphism is injective

        '''
        return self.kernel().order() == 1

    def _is_surjective_cert(self):
        '''
        Return a certificate of surjectivity, i.e. a dictionary of
        preimages of the codomain generators, if the homomorphism
        is surjective. Otherwise, return None.

        '''
        if not self._surjective_cert_computed:
            preimages = None
            try:
                preimages = {g: self.invert(g) for g in self.codomain.generators}
            except ValueError:
                pass
            self._surjective_cert = preimages
            self._surjective_cert_computed = True
        return self._surjective_cert

    def is_surjective(self):
        '''
        Check if the homomorphism is surjective

        '''
        return self._is_surjective_cert() is not None

    def is_isomorphism(self):
        '''
        Check if `self` is an isomorphism.

        '''
        return self.is_injective() and self.is_surjective()

    def is_trivial(self):
        '''
        Check is `self` is a trivial homomorphism, i.e. all elements
        are mapped to the identity.

        '''
        return self.image().order() == 1

    def factor(self):
        '''
        Return a factorization of ``self`` as ``iota`` o ``pi`` where
        ``pi`` is surjective and ``iota`` is injective.

        Returns
        =======

        (pi, iota) : tuple of ``GroupHomomorphism``
            ``pi`` maps from the domain of ``self`` onto a presentation
            of the image of ``self`` and ``iota`` embeds that image into
            the codomain of ``self``.

        '''
        if self._factorization is not None:
            return self._factorization

        if isinstance(self.codomain, PermutationGroup):
            image = self.image()
            surj = homomorphism(
                self.domain,
                image,
                self.domain.generators,
                [self.images[g] for g in self.domain.generators],
                check=False,
            )
            inj = homomorphism(
                image,
                self.codomain,
                image.generators,
                image.generators,
            )
            self._factorization = (surj, inj)
            return self._factorization

        if not isinstance(self.codomain, (FpGroup, FreeGroup)):
            raise NotImplementedError(
                "Factorization is implemented only for permutation, "
                "finitely presented, and free codomains"
            )

        if isinstance(self.codomain, FpGroup):
            codomain = self.codomain
        else:
            codomain = FpGroup(self.codomain, [])

        image_gens = list(dict.fromkeys(self.images.values()))
        image_group, image_inclusion = codomain.subgroup(
            image_gens,
            homomorphism=True,
        )

        surj_images = [
            image_inclusion.invert(self.images[g]) for g in self.domain.generators
        ]
        surj = homomorphism(
            self.domain,
            image_group,
            self.domain.generators,
            surj_images,
            check=False,
        )

        inj = homomorphism(
            image_group,
            self.codomain,
            image_group.generators,
            image_inclusion(image_group.generators),
            check=False,
        )
        self._factorization = (surj, inj)
        return self._factorization

    def compose(self, other):
        '''
        Return the composition of `self` and `other`, i.e.
        the homomorphism phi such that for all g in the domain
        of `other`, phi(g) = self(other(g))

        '''
        if not other.image().is_subgroup(self.domain):
            raise ValueError("The image of `other` must be a subgroup of "
                    "the domain of `self`")
        images = {g: self(other(g)) for g in other.images}
        return GroupHomomorphism(other.domain, self.codomain, images)

    def restrict_to(self, H):
        '''
        Return the restriction of the homomorphism to the subgroup `H`
        of the domain.

        '''
        if not isinstance(H, PermutationGroup) or not H.is_subgroup(self.domain):
            raise ValueError("Given H is not a subgroup of the domain")
        domain = H
        images = {g: self(g) for g in H.generators}
        return GroupHomomorphism(domain, self.codomain, images)

    def invert_subgroup(self, H):
        '''
        Return the subgroup of the domain that is the inverse image
        of the subgroup ``H`` of the homomorphism image

        '''
        if not H.is_subgroup(self.image()):
            raise ValueError("Given H is not a subgroup of the image")
        gens = []
        P = PermutationGroup(self.image().identity)
        for h in H.generators:
            h_i = self.invert(h)
            if h_i not in P:
                gens.append(h_i)
                P = PermutationGroup(gens)
            for k in self.kernel().generators:
                if k*h_i not in P:
                    gens.append(k*h_i)
                    P = PermutationGroup(gens)
        return P

def homomorphism(domain, codomain, gens, images=(), check=True):
    '''
    Create (if possible) a group homomorphism from the group ``domain``
    to the group ``codomain`` defined by the images of the domain's
    generators ``gens``. ``gens`` and ``images`` can be either lists or tuples
    of equal sizes. If ``gens`` is a proper subset of the group's generators,
    the unspecified generators will be mapped to the identity. If the
    images are not specified, a trivial homomorphism will be created.

    If the given images of the generators do not define a homomorphism,
    an exception is raised.

    If ``check`` is ``False``, do not check whether the given images actually
    define a homomorphism.

    '''
    if not isinstance(domain, (PermutationGroup, FpGroup, FreeGroup)):
        raise TypeError("The domain must be a group")
    if not isinstance(codomain, (PermutationGroup, FpGroup, FreeGroup)):
        raise TypeError("The codomain must be a group")

    generators = domain.generators
    if not all(g in generators for g in gens):
        raise ValueError("The supplied generators must be a subset of the domain's generators")
    if not all(g in codomain for g in images):
        raise ValueError("The images must be elements of the codomain")

    if images and len(images) != len(gens):
        raise ValueError("The number of images must be equal to the number of generators")

    gens = list(gens)
    images = list(images)

    images.extend([codomain.identity]*(len(generators)-len(images)))
    gens.extend([g for g in generators if g not in gens])
    images = dict(zip(gens,images))

    if check and not _check_homomorphism(domain, codomain, images):
        raise ValueError("The given images do not define a homomorphism")
    return GroupHomomorphism(domain, codomain, images)

def _check_homomorphism(domain, codomain, images):
    """
    Check that a given mapping of generators to images defines a homomorphism.

    Parameters
    ==========
    domain : PermutationGroup, FpGroup, FreeGroup
    codomain : PermutationGroup, FpGroup, FreeGroup
    images : dict
        The set of keys must be equal to domain.generators.
        The values must be elements of the codomain.

    """
    pres = domain if hasattr(domain, 'relators') else domain.presentation()
    rels = pres.relators
    gens = pres.generators
    symbols = [g.ext_rep[0] for g in gens]
    symbols_to_domain_generators = dict(zip(symbols, domain.generators))
    identity = codomain.identity
    dtype = type(identity)

    def _image(r):
        return dtype.prod(
            images[symbols_to_domain_generators[symbol]]**power
            for symbol, power in r.array_form
        )

    for r in rels:
        if isinstance(codomain, FpGroup):
            s = codomain.equals(_image(r), identity)
            if s is None:
                # only try to make the rewriting system
                # confluent when it can't determine the
                # truth of equality otherwise
                success = codomain.make_confluent()
                s = codomain.equals(_image(r), identity)
                if s is None and not success:
                    raise RuntimeError("Can't determine if the images "
                        "define a homomorphism. Try increasing "
                        "the maximum number of rewriting rules "
                        "(group._rewriting_system.set_max(new_value); "
                        "the current value is stored in group._rewriting"
                        "_system.maxeqns)")
        else:
            s = _image(r).is_identity
        if not s:
            return False
    return True

def orbit_homomorphism(group, omega):
    '''
    Return the homomorphism induced by the action of the permutation
    group ``group`` on the set ``omega`` that is closed under the action.

    '''
    from sympy.combinatorics import Permutation
    from sympy.combinatorics.named_groups import SymmetricGroup
    codomain = SymmetricGroup(len(omega))
    identity = codomain.identity
    omega = list(omega)
    images = {g: identity*Permutation([omega.index(o^g) for o in omega]) for g in group.generators}
    group._schreier_sims(base=omega)
    H = GroupHomomorphism(group, codomain, images)
    if len(group.basic_stabilizers) > len(omega):
        H._kernel = group.basic_stabilizers[len(omega)]
    else:
        H._kernel = PermutationGroup([group.identity])
    return H

def block_homomorphism(group, blocks):
    '''
    Return the homomorphism induced by the action of the permutation
    group ``group`` on the block system ``blocks``. The latter should be
    of the same form as returned by the ``minimal_block`` method for
    permutation groups, namely a list of length ``group.degree`` where
    the i-th entry is a representative of the block i belongs to.

    '''
    from sympy.combinatorics import Permutation
    from sympy.combinatorics.named_groups import SymmetricGroup

    n = len(blocks)

    # number the blocks; m is the total number,
    # b is such that b[i] is the number of the block i belongs to,
    # p is the list of length m such that p[i] is the representative
    # of the i-th block
    m = 0
    p = []
    b = [None]*n
    for i in range(n):
        if blocks[i] == i:
            p.append(i)
            b[i] = m
            m += 1
    for i in range(n):
        b[i] = b[blocks[i]]

    codomain = SymmetricGroup(m)
    # the list corresponding to the identity permutation in codomain
    identity = range(m)
    images = {g: Permutation([b[p[i]^g] for i in identity]) for g in group.generators}
    H = GroupHomomorphism(group, codomain, images)
    return H

def group_isomorphism(G, H, isomorphism=True):
    '''
    Compute an isomorphism between 2 given groups.

    Parameters
    ==========

    G : A finite ``FpGroup`` or a ``PermutationGroup``.
        First group.

    H : A finite ``FpGroup`` or a ``PermutationGroup``
        Second group.

    isomorphism : bool
        This is used to avoid the computation of homomorphism
        when the user only wants to check if there exists
        an isomorphism between the groups.

    Returns
    =======

    If isomorphism = False -- Returns a boolean.
    If isomorphism = True  -- Returns a boolean and an isomorphism between `G` and `H`.

    Examples
    ========

    >>> from sympy.combinatorics import free_group, Permutation
    >>> from sympy.combinatorics.perm_groups import PermutationGroup
    >>> from sympy.combinatorics.fp_groups import FpGroup
    >>> from sympy.combinatorics.homomorphisms import group_isomorphism
    >>> from sympy.combinatorics.named_groups import DihedralGroup, AlternatingGroup

    >>> D = DihedralGroup(8)
    >>> p = Permutation(0, 1, 2, 3, 4, 5, 6, 7)
    >>> P = PermutationGroup(p)
    >>> group_isomorphism(D, P)
    (False, None)

    >>> F, a, b = free_group("a, b")
    >>> G = FpGroup(F, [a**3, b**3, (a*b)**2])
    >>> H = AlternatingGroup(4)
    >>> (check, T) = group_isomorphism(G, H)
    >>> check
    True
    >>> T(b*a*b**-1*a**-1*b**-1)
    (0 2 3)

    Notes
    =====

    Uses the approach suggested by Robert Tarjan to compute the isomorphism between two groups.
    First, the generators of ``G`` are mapped to the elements of ``H`` and
    we check if the mapping induces an isomorphism.
    If G is a permutation group only mappings between elements of the same order are checked.

    '''
    if not isinstance(G, (PermutationGroup, FpGroup)):
        raise TypeError("The group must be a PermutationGroup or an FpGroup")
    if not isinstance(H, (PermutationGroup, FpGroup)):
        raise TypeError("The group must be a PermutationGroup or an FpGroup")

    if isinstance(G, FpGroup) and isinstance(H, FpGroup):
        G = simplify_presentation(G)
        H = simplify_presentation(H)
        # Two infinite FpGroups with the same generators are isomorphic
        # when the relators are same but are ordered differently.
        if G.generators == H.generators and (G.relators).sort() == (H.relators).sort():
            if not isomorphism:
                return True
            return (True, homomorphism(G, H, G.generators, H.generators))

    #  `_H` is the permutation group isomorphic to `H`.
    _H = H
    g_order = G.order()
    h_order = H.order()

    if g_order is S.Infinity:
        raise NotImplementedError("Isomorphism methods are not implemented for infinite groups.")

    if isinstance(H, FpGroup):
        if h_order is S.Infinity:
            raise NotImplementedError("Isomorphism methods are not implemented for infinite groups.")
        _H, h_isomorphism = H._to_perm_group()

    if (g_order != h_order) or (G.is_abelian != H.is_abelian):
        if not isomorphism:
            return False
        return (False, None)

    if not isomorphism:
        # Two groups of the same cyclic numbered order
        # are isomorphic to each other.
        n = g_order
        if (igcd(n, totient(n))) == 1:
            return True

    # Match the generators of `G` with subsets of `_H`
    gens = list(G.generators)
    if isinstance(G, PermutationGroup):
        diz_gen = defaultdict(list)
        diz_imm = defaultdict(list)
        for g in gens:
            diz_gen[g.order()].append(g)
        for el in _H.generate():
            j = el.order()
            if j in diz_gen:
                diz_imm[j].append(el)

        # creating a list with inside lists that contain all possible permutations of elements of an order
        lista = []
        for key,g in diz_gen.items():
            sublist = list(itertools.permutations(diz_imm[key], len(g)))
            lista.append(sublist)

        # creating cartesian product of all permutations sorted by order
        all_matches = itertools.product(*lista)
        for m in all_matches:
            images = []
            # this is flattening out the product
            for i in range(len(diz_gen)):
                images.extend(m[i])
            _images = dict(zip(gens,images))
            if _check_homomorphism(G, _H, _images):
                if isinstance(H, FpGroup):
                    images = h_isomorphism.invert(images)
                T =  homomorphism(G, H, G.generators, images, check=False)
                if T.is_isomorphism():
                    # It is a valid isomorphism
                    if not isomorphism:
                        return True
                    return (True, T)

        if not isomorphism:
            return False
        return (False, None)

    else:
        for subset in itertools.permutations(_H.generate(), len(gens)):
            images = list(subset)
            images.extend([_H.identity]*(len(G.generators)-len(images)))
            _images = dict(zip(gens,images))
            if _check_homomorphism(G, _H, _images):
                if isinstance(H, FpGroup):
                    images = h_isomorphism.invert(images)
                T =  homomorphism(G, H, G.generators, images, check=False)
                if T.is_isomorphism():
                    # It is a valid isomorphism
                    if not isomorphism:
                        return True
                    return (True, T)

        if not isomorphism:
            return False
        return (False, None)

def is_isomorphic(G, H):
    '''
    Check if the groups are isomorphic to each other

    Parameters
    ==========

    G : A finite ``FpGroup`` or a ``PermutationGroup``
        First group.

    H : A finite ``FpGroup`` or a ``PermutationGroup``
        Second group.

    Returns
    =======

    boolean
    '''
    return group_isomorphism(G, H, isomorphism=False)
