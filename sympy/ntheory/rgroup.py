from sympy.core import Basic
from sympy.ntheory.residue_ntheory import Residue, relprimes

class Rgroup(Basic):
    """
    The class defining a residue group.
    """
    _modulus = None
    _elements = []
    _totient = None

    @property
    def modulus(self):
        """
        Returns the modulus.
        """
        return self._modulus

    @property
    def elements(self):
        """
        Returns the elements of the group.
        """
        return self._elements

    @property
    def totient(self):
        """
        Returns the totient.
        """
        return self._totient

    def __new__(cls, *args, **kw_args):
        """
        The default constructor for the residue group.
        """
        ret_obj = Basic.__new__(cls, *args, **kw_args)
        ret_obj.modulus = args[0]
        elements = kw_args["elements"]
        if elements == None:
            ret_obj._elements = [Residue(x, modulus) for x in relprimes(modulus)]
        else:
            ret_obj._elements = elements
            ret_obj._totient = len(ret_obj.elements)
        return ret_obj

    def exp(self,i):
        return [self[i]**x for x in xrange(1, self[i].ord() + 1)]

    def __getitem__(self,k):
        return self.elements[k]
