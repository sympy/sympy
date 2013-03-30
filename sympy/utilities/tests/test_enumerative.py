from itertools import izip_longest
from sympy.utilities.iterables import  multiset_partitions
from sympy.utilities.enumerative import (
    factoring_visitor,
    list_visitor,
    MultisetPartitionTraverser,
    multiset_partitions_taocp,
    PartComponent,
    part_key
    )

# first some functions only useful as test scaffolding - these provide
# straightforward, but slow reference implementations against which to
# compare the real versions, and also a comparision to verify that
# different versions are giving identical results.

def part_range_filter(partition_iterator, lb, ub):
    """
    Filters (on the number of parts) a multiset partition enumeration

    Arguments
    =========

    lb, and ub are a range (in the python slice sense) on the lpart
    variable returned from a multiset partition enumeration.  Recall
    that lpart is 0-based (it points to the topmost part on the part
    stack), so if you want to return parts of sizes 2,3,4,5 you would
    use lb=1 and ub=5.
    """
    for state in partition_iterator:
        f, lpart, pstack = state
        if lpart >= lb and lpart < ub:
            yield state

def multiset_partitions_slow(multiplicities):
    # todo cheesey version of smichr's code
    pass

def compare_multiset_states(s1, s2):
    """compare for equality two instances of multiset partition states

    This is useful for comparing different versions of the algorithm
    to verify correctness."""
    # Comparison is physical, the only use of semantics is to ignore
    # trash off the top of the stack.
    f1, lpart1, pstack1 = s1
    f2, lpart2, pstack2 = s2

    if (lpart1 == lpart2) and (f1[0:lpart1+1] == f2[0:lpart2+1]):
        if pstack1[0:f1[lpart1+1]] == pstack2[0:f2[lpart2+1]]:
            return True
    return False

def test_multiset_versions():
    """Compares Knuth-based versions of multiset_partitions"""
    multiplicities = [5,2,2,1]
    m = MultisetPartitionTraverser()
    for s1, s2 in izip_longest(m.enum_all(multiplicities),
                               multiset_partitions_taocp(multiplicities)):
        assert compare_multiset_states(s1, s2)

def test_subrange():
    """Compare filter-based and more optimized subrange implementations"""
    multiplicities = [4,4,2,1] # mississippi

    m = MultisetPartitionTraverser()
    assert m.count_partitions(multiplicities) == \
        m.count_partitions_slow(multiplicities)

    # Note - you can't have multiple traversals from the same
    # MultisetPartitionTraverser object going at the same time, hence
    # make several here.
    ma = MultisetPartitionTraverser()
    mc = MultisetPartitionTraverser()
    md = MultisetPartitionTraverser()

    #  Several paths to compute just the size two partitions
    a_it = ma.enum_range(multiplicities, 1, 2)
    b_it = part_range_filter(multiset_partitions_taocp(multiplicities),1,2)
    c_it = part_range_filter(mc.enum_small(multiplicities, 2), 1, 10)
    d_it = part_range_filter(md.enum_large(multiplicities, 1), 0, 2)

    for sa, sb, sc, sd in izip_longest(a_it, b_it, c_it, d_it):
        assert compare_multiset_states(sa, sb)
        assert compare_multiset_states(sa, sc)
        assert compare_multiset_states(sa, sd)
