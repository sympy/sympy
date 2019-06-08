from sympy.core import Basic
from sympy import sieve
from sympy.combinatorics.perm_groups import PermutationGroup

class PcGroup(Basic):

    is_group = True
    is_solvable = True

    def __init__(self, _pcgs):
        self.perm_group = PermutationGroup(_pcgs)
        self.pc_series = self._pc_series()
        self.pcgs = self._compute_pcgs()

    def _pc_series(self):
        return self.perm_group.composition_series()

    def _compute_pcgs(self):
        # computes the generating sequence for polycyclic groups.
        series = self.pc_series
        pcgs = []
        for i in range(len(series)-1):
            for g in series[i].generators:
                if not g in series[i+1]:
                    pcgs.append(g)
        return pcgs

    def relative_orders(self):
        rel_orders = []
        for i in range(len(self.pc_series)-1):
            G = self.pc_series[i]
            H = self.pc_series[i+1]
            rel_orders.append(G.order()//H.order())
        return rel_orders

    def is_prime_order(self):
        for order in self.relative_orders():
            if order not in sieve:
                return False
        return True

    def length(self):
        return len(self.pcgs)

    def pc_element_exponent(self, element):
        series = self.pc_series
        pcgs = self.pcgs
        exponent = [0]*len(series)
        for i in range(len(series)):
            exp = 0
            if not element in series[i]:
                for j in range(len(pcgs)):
                    element = (pcgs[j]**-1)*element
                    exp = exp + 1
                    if element in series[i]:
                        exponent[i] = exp
                        break
        return exponent
