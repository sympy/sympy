from sympy.combinatorics.fp_groups import FpGroup, FpSubgroup, subgroup_quotient
from sympy import sieve

class PcGroup(FpGroup):
    '''
    Class implementing Polycyclic groups.
    To-do : Docstring completion.
            Add error and exception handling
    '''
    def __init__(fp_group, pc_series=None, pcgs=None):
        self.fp_group = fp_group
        self.pc_series = pc_series
        self.pcgs = pcgs
        self.rel_orders = self.relative_orders()

    def compute_pcgs(self):
        pcgs = []
        for i in range(0, len(self.pc_series) - 1):
            quotient = subgroup_quotient(self.fp_grp.free_group, self.fp_grp,self.pc_series[i], self.pc_series[i+1])
            pcgs.append(quotient.random())
        self.pcgs = pcgs

    def compute_pc_series(self):
        self.pc_series = self.fp_group.compute_polycyclic_series()

    def relative_orders(self):
        rel_orders = []
        for i in range(0, len(self.pc_series) - 1):
            first_grp = self.pc_series[i]
            sec_grp = self.pc_series[i+1]
            if isinstance(self.pc_series[i], list):
                first_grp = fp_group.subgroup(self.pc_series[i])
            if isinstance(self.pc_series[i+1], list):
                sec_grp = fp_group.subgroup(self.pc_series[i+1])
            rel_orders[i] = first_grp.order()/sec_grp.order()
        self.rel_orders = rel_orders

    def is_prime_order(self):
        for elem in self.relative_orders:
            if elem not in sieve:
                return False
        return True

    ##########
    #  To-Do #
    ##########

    # 1. Documentation
    # 2. Following methods.

    def compute_presentation(self):
        '''
        Compute the presentation of the PcGroup
        '''
        return

    def collection_algo(self):
        '''
        Implemenation of the collection algo
        '''
        return

    def is_consistent(self):
        '''
        Could possibly use the existing methods.
        '''
        return
