from __future__ import print_function, division

from sympy.combinatorics.fp_groups import FpGroup, FpSubgroup
from sympy.combinatorics.free_groups import FreeGroupElement
from sympy.combinatorics.perm_groups import PermutationGroup

class GroupHomomorphism(object):
    '''
    A class representing group homomorphisms. Instantiate using `homomorphism()`.

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
        generator of `codomain` and `inverse` is an element of its preimage

        '''
        vals = [v for v in self.images.values() if not v.is_identity]
        keys = self.images.keys()
        vs = list(set(vals))
        keys = [keys[vals.index(v)] for v in vs]
        inverses = dict(zip(vals,keys))
        gens = [g for g in self.image().strong_gens if g not in inverses]
        for g in gens:
            w = self.codomain.identity
            for s in self.codomain.coset_factor(g):
                if s.is_identity:
                    continue
                w = w*inverses[s]
            inverses[g] = w
        inverses[self.codomain.identity] = self.domain.free_group.identity
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
        w = self.domain.identity
        print(self.codomain.coset_factor(g))
        for s in self.codomain.coset_factor(g):
            w = w*self._inverses[s]
        return w

    def kernel(self):
        '''
        Compute the kernel of `self`.

        =======
        Example
        =======

        '''
        if self._kernel is None:
            self._kernel = self._compute_kernel()
        return self._kernel


    def _compute_kernel(self):
        from sympy import S
        G = self.domain
        if G.order() == S.Infinity:
            raise NotImplementedError(
                "Kernel computation is not implemented for infinite groups")
        gens = []
        K = FpSubgroup(G, gens)
        i = self.image().order()
        while K.order()*i != G.order():
            r = G.random_element()
            k = r*self.invert(self(r))
            if not k in K:
                gens.append(k)
                K = FpSubgroup(G, gens)
        return K
        
            
    def image(self):
        '''
        Compute the image of `self`.

        =======
        Example
        =======

        '''
        if self._image is None:
            values = list(set(self.images.values()))
            if isinstance(self.codomain, PermutationGroup):
                self._image = self.codomain.subgroup(values)
            else:
                self._image = FpSubgroup(self.codomain, values)
        return self._image

    def _apply(self, elem):
        '''Apply `self` to `elem`.'''
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

    def _compose(self, oth):
        '''
        Compose `self` with `oth`, so that the returned homomorphism
        is `self.oth(elem) == self(oth(elem))`

        '''
        if self.domain != oth.codomain:
            raise ValueError("The given homomorphism is incompatible with 'self'")
        domain = oth.domain
        codomain = self.codomain
        images = oth.images
        for g in images:
            images[g] = self(images[g])
        return GroupHomomorphism(domain, codomain, images)

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
            return self.image().order() == self.codomain.order()

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
    Create (if possible) a group homomorphism from the group `domain` defined by
    the images of its generators `gens`. The image elements must all belong to
    the same group. `gens` and `images` can be either lists or tuples of equal sizes.
    If `gens` is a proper subset of the group's generators, the unspecified generators
    will be mapped to the identity. Also, if the images are not specified, a trivial
    homomorphism will be created.

    If the given images of the generators do not define a homomorphism, an exception
    is raised.

    Examples
    ========

    '''
    if isinstance(domain, PermutationGroup):
        raise NotImplementedError("Homomorphisms from permutation groups are not currently implemented")
    elif not isinstance(domain, FpGroup):
        raise TypeError("The domain must be a group")
    if not(isinstance(codomain, PermutationGroup) or isinstance(codomain, FpGroup)):
        raise TypeError("The codomain must be a group")

    generators = domain.generators
    if any([g not in generators for g in gens]):
        raise ValueError("The supplied generators must be a subset of the domain's generators")
    if any([g not in codomain for g in images]):
        print([g for g in images if g not in codomain],g[0] in codomain, codomain)
        raise ValueError("The images must be elements of the codomain")

    if len(images) > len(gens):
        raise ValueError("Too many images were given")

    images.extend([codomain.identity]*(len(generators)-len(images)))
    gens.extend([g for g in generators if g not in gens])
    images = dict(zip(gens,images))
    
    if not _check_homomorphism(domain, images, codomain.identity):
        raise ValueError("The given images do not define a homomorphism")
    return GroupHomomorphism(domain, codomain, images)

def _check_homomorphism(domain, images, identity):
    rels = domain.relators
    def _image(r):
        if r.is_identity:
            return identity
        elif len(r) == 1:
            return images[r]
        else:
            power = r.array_form[0][1]
            return images[r[0]]**power*_image(r.subword(abs(power),len(r)))

    print([_image(r) for r in rels])
    if any([not _image(r).is_identity for r in rels]):
        print(
        return False
    else:
        return True
