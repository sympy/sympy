from __future__ import print_function, division

from sympy.combinatorics.fp_groups import FpGroup, FpSubgroup
from sympy.combinatorics.free_groups import FreeGroup, FreeGroupElement
from sympy.combinatorics.perm_groups import PermutationGroup

class GroupHomomorphism(object):
    '''
    A class representing group homomorphisms. Instantiate using `homomorphism()`.

    References
    ==========
    [1] Holt, D., Eick, B. and O'Brien, E. (2005). Handbook of computational group theory.

    '''

    def __init__(self, domain, codomain, images):
        self.domain = domain
        self.codomain = codomain
        self.images = images
        self._inverses = None
        self._kernel = None
        self._image = None

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
        gens = image.strong_gens
        for g in gens:
            if g in inverses or g.is_identity:
                continue
            w = self.domain.identity
            for s in image._strong_gens_slp[g]:
                if s in inverses:
                    w = inverses[s]*w
                else:
                    w = inverses[s**-1]**-1*w
            inverses[g] = w
        return inverses

    def invert(self, g):
        '''
        Return an element of the preimage of `g`

        '''
        if not isinstance(self.codomain, PermutationGroup):
            raise NotImplementedError(
                "Only elements of PermutationGroups can be inverted")
        if self._inverses is None:
            self._inverses = self._invs()
        image = self.image()
        w = self.domain.identity
        for g in image.generator_product(g):
            if g.is_identity:
                continue
            w = self._inverses[g]*w
        return w

    def kernel(self):
        '''
        Compute the kernel of `self`.

        '''
        if self._kernel is None:
            self._kernel = self._compute_kernel()
        return self._kernel

    def _compute_kernel(self):
        from sympy import S
        G = self.domain
        G_order = G.order()
        if G_order == S.Infinity:
            raise NotImplementedError(
                "Kernel computation is not implemented for infinite groups")
        gens = []
        K = FpSubgroup(G, gens)
        i = self.image().order()
        while K.order()*i != G_order:
            r = G.random_element()
            k = r*self.invert(self(r))
            if not k in K:
                gens.append(k)
                K = FpSubgroup(G, gens)
        return K

    def image(self):
        '''
        Compute the image of `self`.

        '''
        if self._image is None:
            values = list(set(self.images.values()))
            if isinstance(self.codomain, PermutationGroup):
                self._image = self.codomain.subgroup(values)
            elif isinstance(self.codomain, FpGroup):
                self._image = FpSubgroup(self.codomain, values)
            else:
                self._image = FreeSubgroup(self.codomain, values)
        return self._image

    def _apply(self, elem):
        '''
        Apply `self` to `elem`.

        '''
        if not elem in self.domain.free_group:
            raise ValueError("The supplied element doesn't belong to the domain")
        if elem.is_identity:
            return self.codomain.identity
        else:
            p = elem.array_form[0][1]
            if p < 0:
                g = elem[0]**-1
            else:
                g = elem[0]
            return self.images[g]**p*self._apply(elem.subword(abs(p), len(elem)))

    def __call__(self, elem):
        return self._apply(elem)

    def is_injective(self):
        '''
        Check if the homomorphism is injective

        '''
        return self.kernel().order() == 1

    def is_surjective(self):
        '''
        Check if the homomorphism is surjective

        '''
        from sympy import S
        im = self.image().order()
        oth = self.codomain.order()
        if im == S.Infinity and oth == S.Infinity:
            return None
        else:
            return im == oth

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

def homomorphism(domain, codomain, gens, images=[]):
    '''
    Create (if possible) a group homomorphism from the group `domain`
    to the group `codomain` defined by the images of the domain's
    generators `gens`. `gens` and `images` can be either lists or tuples
    of equal sizes. If `gens` is a proper subset of the group's generators,
    the unspecified generators will be mapped to the identity. If the
    images are not specified, a trivial homomorphism will be created.

    If the given images of the generators do not define a homomorphism,
    an exception is raised.

    '''
    if isinstance(domain, PermutationGroup):
        raise NotImplementedError("Homomorphisms from permutation groups are not currently implemented")
    elif not isinstance(domain, (FpGroup, FreeGroup)):
        raise TypeError("The domain must be a group")
    if not isinstance(codomain, (PermutationGroup, FpGroup, FreeGroup)):
        raise TypeError("The codomain must be a group")

    generators = domain.generators
    if any([g not in generators for g in gens]):
        raise ValueError("The supplied generators must be a subset of the domain's generators")
    if any([g not in codomain for g in images]):
        raise ValueError("The images must be elements of the codomain")

    if images and len(images) != len(gens):
        raise ValueError("The number of images must be equal to the number of generators")

    gens = list(gens)
    images = list(images)
    images.extend([codomain.identity]*(len(generators)-len(images)))
    gens.extend([g for g in generators if g not in gens])
    images = dict(zip(gens,images))

    if not _check_homomorphism(domain, images, codomain.identity):
        raise ValueError("The given images do not define a homomorphism")
    return GroupHomomorphism(domain, codomain, images)

def _check_homomorphism(domain, codomain, images):
    rels = domain.relators
    identity = codomain.identity
    def _image(r):
        if r.is_identity:
            return identity
        else:
            w = identity
            r_arr = r.array_form
            i = 0
            while i < len(r):
                power = r_arr[i][1]
                if r[i] in images:
                    w = w*images[r[i]]**power
                else:
                    w = w*images[r[i]**-1]**power
                i += abs(power)
            return w

    for r in rels:
        if isinstance(codomain, PermutationGroup):
            s = _image(r).is_identity
        else:
            codomain.make_confluent()
            s = codomain.equals(_image(r), identity)
        if not s:
            return False
    return True
