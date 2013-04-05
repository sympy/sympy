#
#  Algorithms and classes to support enumerative combinatorics.
#  Currently just multiset partitions, but more will likely be added.
#

class PartComponent(object):
    """Internal class used in support of the multiset partitions
    enumerators and the associated visitor functions.

    Represents one component of one part of the current partition.

    A stack of these, plus an auxiliary frame array, f, represents a
    partition of the multiset.

    Knuth's psuedocode makes c, u, and v separate arrays.
    """

    __slots__ = ('c', 'u', 'v')
    def __init__(self):
        self.c = 0   # Component number
        self.u = 0   # The as yet unpartitioned amount in component c
                     # *before* it is allocated by this triple
        self.v = 0   # Amount of c component in the current part
                     # (v<=u).  An invariant of the representation is
                     # that the next higher triple for this component
                     # (if there is one) will have a value of u-v in
                     # its u attribute.
    def __repr__(self):
        "for debug/algorithm animation purposes"
        return  'c:%d u:%d v:%d' % (self.c,self.u,self.v)

    def __eq__(self,other):
        """Define  value oriented equality, which is useful for testers"""
        return (isinstance(other, self.__class__) and
                self.c == other.c and
                self.u == other.u and
                self.v == other.v)

    def __ne__(self, other):
        """Defined for consistency with __eq__"""
        return not self.__eq__(other)


# This function tries to be a faithful implementation of algorithm
# 7.1.2.5M in Volume 4A, Combinatoral Algorithms, Part 1, of The Art
# of Computer Programming, by Donald Knuth.  This includes using
# (mostly) the same variable names, etc.  This makes for rather
# low-level Python.

# Changes from Knuth's psuedocode include
# - use PartComponent struct/object instead of 3 arrays
# - make the function a generator
# - map (with some difficulty) the GOTOs to Python control structures.
# - Knuth uses 1-based numbering for components, this code is 0-based
# - renamed variable l to lpart.
# - flag variable x takes on values True/False instead of 1/0
#
def multiset_partitions_taocp(multiplicities):
    """Enumerates partions of a multiset. This is a low-level routine
    that most often will be called by multiset_partitions when appropriate.

    Input
    =====

    multiplicities: list of multiplicities of the components of the multiset.

    Yields
    ======

    [f, lpart, pstack]: an internal datastructure, which can be processed
    by a vistor function which takes input of the elements (not their
    multiplicities) and converts the output from this function into an actual
    partition of elements.  ### add a brief description of what f, lpart and pstack are.

    Usage
    =====

    XXX

    Examples
    =====

    XXX

    """

    # Important variables.
    # m is the number of components, i.e., number of distinct elements
    m = len(multiplicities)
    # n is the cardinality, total number of elements whether or not distinct
    n = sum(multiplicities)

    # The main data structure, f segments pstack into parts.  See
    # list_visitor() for example code indicating how this internal
    # state corresponds to a partition.

    # Note: allocation of space for stack is conservative.  Knuth's
    # exercise 7.2.1.5.68 gives some indication of how to tighten this
    # bound, but this is not implemented.
    pstack = [PartComponent() for i in xrange(n * m + 1)]
    f = [0] * (n + 1)

    # Step M1 in Knuth (Initialize)
    # Initial state - entire multiset in one part.
    for j in xrange(m):
        ps = pstack[j]
        ps.c = j
        ps.u = multiplicities[j]
        ps.v = multiplicities[j]

    # Other variables
    f[0] = 0
    a = 0
    lpart = 0
    f[1] = m
    b = m # in general, current stack frame is from a to b - 1

    while True:
        while True:
            # Step M2 (Subtract v from u)
            j = a
            k = b
            x = False
            while j < b:
                pstack[k].u =  pstack[j].u - pstack[j].v
                if  pstack[k].u == 0:
                    x = True
                elif not x:
                    pstack[k].c = pstack[j].c
                    pstack[k].v = min(pstack[j].v, pstack[k].u)
                    x = pstack[k].u < pstack[j].v
                    k = k + 1
                else:  # x is True
                    pstack[k].c = pstack[j].c
                    pstack[k].v = pstack[k].u
                    k = k + 1
                j = j + 1
                # Note: x is True iff v has changed

            # Step M3 (Push if nonzero.)
            if k > b:
                a = b
                b = k
                lpart = lpart + 1
                f[lpart + 1] = b
                # Return to M2
            else:
                break # Continue to M4

        # M4  Visit a partition
        state = [f, lpart, pstack]
        yield state

        # M5 (Decrease v)
        while True:
            j = b-1
            while (pstack[j].v == 0):
                j = j - 1
            if j == a and pstack[j].v == 1:
                # M6 (Backtrack)
                if lpart == 0 :
                    return
                lpart = lpart - 1
                b = a
                a = f[lpart]
                # Return to M5
            else:
                pstack[j].v = pstack[j].v - 1
                for k in xrange(j + 1, b):
                    pstack[k].v = pstack[k].u
                break # GOTO M2

# --------------- Visitor functions for multiset partitions ---------------
# A visitor takes a partition state generated by multiset_partitions_taocp
# and a description of the elements of the multiset and produces the actual
# partition.

def factoring_visitor(state, primes):
    """Use with multiset_partitions_taocp to show the divisors of a number.
    The exponents of the number are sent to partitioner while the
    corresponding primes are input here.

    Examples
    ========

    To enumerate the factorings of a number we can think of the elements to
    partition as being the prime factors and the multiplicities as being their
    exponents.

    >>> from sympy.utilities.enumerative import factoring_visitor
    >>> from sympy.utilities.enumerative import multiset_partitions_taocp
    >>> from sympy import factorint

    >>> primes, multiplicities = zip(*factorint(24).items())
    >>> primes
    (2, 3)
    >>> multiplicities
    (3, 1)

    >>> states = multiset_partitions_taocp(multiplicities)
    >>> list(factoring_visitor(state, primes) for state in states)
    [[24], [8, 3], [12, 2], [4, 6], [4, 2, 3], [6, 2, 2], [2, 2, 2, 3]]
    """
    f, lpart, pstack = state
    factoring = []
    for i in xrange(lpart + 1):
        factor = 1
        for ps in pstack[f[i]: f[i + 1]]:
            if ps.v > 0:
                factor *= primes[ps.c] ** ps.v
        factoring.append(factor)
    return factoring


def list_visitor(state, components):
    """Return a list of lists to represent the partition.

    Examples
    ========

    >>> from sympy.utilities.enumerative import list_visitor
    >>> from sympy.utilities.enumerative import multiset_partitions_taocp
    >>> states = multiset_partitions_taocp([1, 2, 1])
    >>> s = state.next()
    >>> list_visitor(s, 'abc')  # for multiset 'a b b c'
    [['a', 'b', 'b', 'c']]
    >>> s = states.next()
    >>> list_visitor(s, [1, 2, 3])  # for multiset '1 2 2 3
    [[1, 2, 2], [3]]
    """
    f, lpart, pstack = state

    partition = []
    for i in xrange(lpart+1):
        part = []
        for ps in pstack[f[i]:f[i+1]]:
            if ps.v > 0:
                part.extend([components[ps.c]] * ps.v)
        partition.append(part)

    return partition

class MultisetPartitionTraverser():
    """
    Has methods ``enumerate`` and ``count`` which ethe partitions of a multiset.

    This is a refactored and extended version of algorithm 7.1.2.5M
    in Knuth's "The Art of Computer Programming."

    See Also
    ========
    multiset_partitions_taocp


    Examples
    ========

    XXX todo

    References
    ==========

    ..[1] Algorithm 7.1.2.5M in Volume 4A, Combinatoral Algorithms,
          Part 1, of The Art of Computer Programming, by Donald Knuth.

    ..[2] On a Problem of Oppenheim concerning "Factorisatio Numerorum"
          E. R. CANFIELD, Paul Erdos, Carl Pomerance,
          JOURNAL OF NUMEER THEORY, Vol. 17, No. 1. August 1983
          (See section 7 for a description of an algorithm similar to
          Knuth's)

    ..[3] Generating Multiset Partitions, Brent Yorgey, The
          Monad.Reader, Issue 8, September 2007.

    """

    def __init__(self):
        self.debug = False
        # TRACING variables.  These are useful for gathering
        # statistics on the algorithm itself, but have no particular
        # benefit to a user of the code.
        self.k1 = 0
        self.k2 = 0
        self.p1 = 0


    def db_trace(self, msg):
        """Useful for usderstanding/debugging the algorithms.  Not
        generally activated in end-user code."""
        if self.debug:
            letters = 'abcdefghijklmnopqrstuvwxyz'
            state = [self.f, self.lpart, self.pstack]
            print "DBG:", msg, \
                ["".join(part) for part in list_visitor(state, letters)], \
                animation_visitor(state)

    #
    # Helper methods for enumeration
    #
    def initialize_enumeration(self, multiplicities):
        """Allocates and initializes the partition stack.

        This is called from the enumeration/counting routines, so
        there is no need to call it separately."""

        num_components = len(multiplicities)
        # n is the cardinality, total number of elements whether or not distinct
        cardinality = sum(multiplicities)

        # pstack is the partition stack, which is segmented by
        # f into parts.
        self.pstack = [PartComponent() for i in
                       xrange(num_components * cardinality + 1)]
        self.f = [0] * (cardinality + 1)

        # Initial state - entire multiset in one part.
        for j in xrange(num_components):
            ps = self.pstack[j]
            ps.c = j
            ps.u = multiplicities[j]
            ps.v = multiplicities[j]

        self.f[0] = 0
        self.f[1] = num_components
        self.lpart = 0

    # This is the method that gets changed if we want to return only
    # partitions with a restricted size range.
    # Corresponds to M5
    def decrement_part(self, part):
        """Decrements part (a subrange of pstack), if possible, returning
        True iff the part was successfully decremented.

        Examples
        ========

        XXX todo

        """
        plen = len(part)
        for j in xrange(plen - 1, -1, -1):
            if (j == 0 and part[j].v > 1) or (j > 0 and part[j].v > 0):
                # found val to decrement
                part[j].v -= 1
                # Reset trailing parts back to maximum
                for k in xrange(j + 1, plen):
                    part[k].v = part[k].u
                return True
        return False

    # Version to allow number of parts to be bounded from above.
    # Corresponds to (a modified) step M5.  This expands on (my
    # probably imperfect understanding of) the answer to problem 69.
    def decrement_part_small(self, part, ub):
        """Decrements part (a subrange of pstack), if possible, returning
        True iff the part was successfully decremented.

        Input
        =====

        XXX define other input


        Examples
        ========

        XXX todo

        Notes
        =====

        The goal of this modification of the ordinary decrement
        method is to fail when it can be proved that this part can
        only have child partitions which are larger than allowed by
        ``lb``.  <- XXX ub? If a decision is made to fail, it must be accurate,
        otherwise the enumeration will miss some partitions.  But, it
        is OK not to capture all the possible failures -- if a part is
        passed that shouldn't be, the resulting too-large partitions
        are filtered by the enumeration one level up.  However, as is
        usual in constrained enumerations, it is advantageous to fail
        as early as possible.

        The tests used by this method catch the most common cases,
        although this implementation is by no means the last word on
        this problem.  The tests include:

        1) ``lpart`` <- XXX part? must be less than ``ub`` by at least 2.  This is because
           once a a part which has been decremented, the partition
           will gain at least one child in the spread step.

        2) If the leading component of the part is about to be
           decremented, check for how many parts will be added in
           order to use up the unallocated multiplicity in that
           leading component, and fail if this number is greater than
           allowed by ``ub``.  (See code for the exact expression.)  This
           test is given in the answer to Knuth's problem 7.2.1.5.69.

        3) If there is *exactly* enough room to expand the leading
           component by the above test, check the next component (if
           it exists) once decrementing has finished.  If this has
           ``v = 0``, this next component will push the expansion over the
           limit by 1, so fail.
        """
        if self.lpart >= ub - 1:
            ### XXX here and below, should 'instrumentation' be 'increment*'?
            self.p1 += 1 # instrumentation to keep track of usefulness of tests
            return False
        plen = len(part)
        for j in xrange(plen - 1, -1, -1):
            # Knuth's mod, (answer to prob 69)
            if (j==0) and (part[0].v - 1)*(ub - self.lpart) < part[0].u:
                self.k1 += 1 # instrumentation
                return False

            if (j == 0 and part[j].v > 1) or (j > 0 and part[j].v > 0):
                # found val to decrement
                part[j].v -= 1
                # Reset trailing parts back to maximum
                for k in xrange(j + 1, plen):
                    part[k].v = part[k].u

                # Have now decremented part, but are we doomed to
                # failure when it is expanded?  Check one oddball case
                # that turns out to be surprisingly common - exactly
                # enough room to expand the leading component, but no
                # room for the second component, which has v=0.
                if (
                        plen > 1 and (part[1].v == 0) and
                        (part[0].u - part[0].v) ==
                        ((ub - self.lpart - 1) * part[0].v)) :
                    self.k2 += 1 # instrumentation
                    self.db_trace("Decrement fails test 3")
                    return False
                return True
        return False

    def decrement_part_large(self, part, amt, lb):
        """Decrements part, while respecting size constraint, returning
        True iff the part was successfully decremented.

        Input
        =====

        ``amt``: (0 or 1) is whether to always decrement.  Call with ``amt=0``
        when enforcing size constraint.

        XXX define other input


        Examples
        ========

        XXX todo

        Notes
        =====

        A part can have no children which are of sufficient size unless that
        part has sufficient unallocated multiplicity.  This method finds the
        next lower decrement of part (if possible) which has sufficient
        unallocated multiplicity.
        """

        if amt == 1:
            # TODO - might be advantageous to inline this, and combine
            # with other iterations.  For now, just do a regular
            # decrement.
            status = self.decrement_part(part)
            if not status:
                return False
        # now decrement as needed to maintain constraint on minimum
        # unallocated multiplicity
        min_unalloc = lb - self.lpart
        if min_unalloc <= 0:
            return True
        total_mult = sum(pc.u for pc in part)
        total_alloc = sum(pc.v for pc in part)
        if total_mult <= min_unalloc:
            return False

        deficit = min_unalloc - (total_mult - total_alloc)
        if deficit <=0:
            return True

        for i in xrange(len(part) - 1, -1, -1):
            if i == 0:
                if part[0].v > deficit:
                    part[0].v -= deficit
                    return True
                else:
                    return False # This shouldn't happen, due to above check
            else:
                 if part[i].v >= deficit:
                    part[i].v -= deficit
                    return True
                 else:
                    deficit -= part[i].v
                    part[i].v = 0

    def decrement_part_range(self, part, lb, ub):
        """Decrements part (a subrange of pstack), if possible, returning
        True iff the part was successfully decremented.

        Input
        =====

        XXX define other input


        Examples
        ========

        XXX todo

        Notes
        =====

        Combines the constraints of _small and _large decrement
        methods.  If returns success, part has been decremented at
        least once, but perhaps by quite a bit more if needed to meet
        the lb constraint.
        """

        status = self.decrement_part_small(part, ub)
        if status:
            status = self.decrement_part_large(part, 0, lb)

        # XXX Should we redo the small checks, in case the large
        # actually did some decrementing?
        return status

    def spread_part_multiplicity(self):
        """Returns True if a new part has been created, and
        adjusts pstack, f and lpart as needed.

        Examples
        ========

        XXX todo

        Notes
        =====

        Spreads unallocated multiplicity from the current top part
        into a new part created above the current on the stack.  This
        new part is constrained to be less than or equal to the old in
        terms of the part ordering.

        This call does nothing (and returns False) if the current top
        part has no unallocated multiplicity.

        """
        j = self.f[self.lpart] # base of current top part
        k = self.f[self.lpart + 1] # ub of current; potential base of next
        base = k  #  save for later comparison

        changed = False # Set to true when the new part (so far) is
                        # strictly less than (as opposed to less than
                        # or equal) to the old.
        for j in xrange(self.f[self.lpart], self.f[self.lpart + 1]):
            self.pstack[k].u =  self.pstack[j].u - self.pstack[j].v
            if  self.pstack[k].u == 0:
                changed = True
            else:
                self.pstack[k].c = self.pstack[j].c
                if changed: # Put all available multiplicity in this part
                    self.pstack[k].v = self.pstack[k].u
                else: # Still maintaining ordering constraint
                    if self.pstack[k].u < self.pstack[j].v:
                        self.pstack[k].v = self.pstack[k].u
                        changed = True
                    else:
                        self.pstack[k].v = self.pstack[j].v
                k = k + 1
        if k > base:
            # Adjust for the new part on stack
            self.lpart = self.lpart + 1
            self.f[self.lpart + 1] = k
            return True
        return False

    def top_part(self):
        """Return current top part on the stack, as a slice of pstack.

        Examples
        ========

        XXX todo

        """
        return self.pstack[self.f[self.lpart]:self.f[self.lpart + 1]]

    def enum_all(self, multiplicities):
        """Enumerate the partitions of a multiset.

        Examples
        ========

        XXX todo

        See also
        ========
        multiset_partitions_taocp()

        Notes
        =====

        The above function gives the same result as this method, but
        is about twice as fast.  Hence, this method is only useful for
        testing and as a base from which to develop the
        range-restricted enumerations.
        """

        self.initialize_enumeration(multiplicities)
        while True:
            while self.spread_part_multiplicity():
                pass

            # M4  Visit a partition
            state = [self.f, self.lpart, self.pstack]
            yield state

            # M5 (Decrease v)
            while not self.decrement_part(self.top_part()):
                # M6 (Backtrack)
                if self.lpart == 0:
                    return
                self.lpart -= 1

    def enum_small(self, multiplicities, ub):
        """Enumerate multiset partitions with no more than ub parts.

        Input
        =====

        XXX define other input


        Examples
        ========

        XXX todo

        References
        ==========
        Knuth, The Art of Computer Programming, Volume 4a,
        Section 7.2.1.5, exercise 69.

        """

        # Keep track of iterations which do not yield a partition.
        # Clearly, we would like to keep this number small.
        self.discarded = 0

        self.initialize_enumeration(multiplicities)
        while True:
            good_partition = True
            while self.spread_part_multiplicity():
                self.db_trace("spread 1")
                if self.lpart >= ub:
                    self.discarded += 1
                    good_partition = False
                    self.db_trace("  Discarding")
                    self.lpart = ub - 2
                    break

            # M4  Visit a partition
            if good_partition:
                state = [self.f, self.lpart, self.pstack]
                yield state

            # M5 (Decrease v)
            while not self.decrement_part_small(self.top_part(), ub):
                self.db_trace("Failed decrement, going to backtrack")
                # M6 (Backtrack)
                if self.lpart == 0 :
                    return
                self.lpart -= 1
                self.db_trace("Backtracked to")
            self.db_trace("decrement ok, about to expand")

    def enum_large(self, multiplicities, lb):
        """Enumerate the partitions of a multiset with lb < num(parts)

        Input
        =====

        XXX define other input


        Examples
        ========

        XXX todo

        """
        self.discarded = 0
        self.initialize_enumeration(multiplicities)
        self.decrement_part_large(self.top_part(), 0, lb)
        while True:
            good_partition = True
            while self.spread_part_multiplicity():
                if not self.decrement_part_large(self.top_part(), 0, lb):
                    # Failure here should be rare/impossible
                    self.discarded += 1
                    good_partition = False
                    break

            # M4  Visit a partition
            if good_partition:
                state = [self.f, self.lpart, self.pstack]
                yield state

            # M5 (Decrease v)
            while not self.decrement_part_large(self.top_part(), 1, lb):
                # M6 (Backtrack)
                if self.lpart == 0:
                    return
                self.lpart -= 1

    def enum_range(self, multiplicities, lb, ub):
        """Enumerate the partitions of a multiset with ``lb < num(parts) <= ub``.

        In particular, if partitions with exactly ``k`` parts are desired,
        call with ``(multiplicities, k - 1, k)``

        Input
        =====

        XXX define other input


        Examples
        ========

        XXX todo

        """
        # Code combines the constraints of the _large and _small enumerations.
        # This is starting as a mash-up of the two
        self.discarded = 0
        self.initialize_enumeration(multiplicities)
        self.decrement_part_large(self.top_part(), 0, lb)
        while True:
            good_partition = True
            while self.spread_part_multiplicity():
                self.db_trace("spread 1")
                if not self.decrement_part_large(self.top_part(), 0, lb):
                    # Failure here - possible in range case?
                    self.db_trace("  Discarding (large cons)")
                    self.discarded += 1
                    good_partition = False
                    break
                elif self.lpart >= ub:
                    self.discarded += 1
                    good_partition = False
                    self.db_trace("  Discarding small cons")
                    self.lpart = ub - 2
                    break

            # M4  Visit a partition
            if good_partition:
                state = [self.f, self.lpart, self.pstack]
                yield state

            # M5 (Decrease v)
            while not self.decrement_part_range(self.top_part(), lb, ub):
                self.db_trace("Failed decrement, going to backtrack")
                # M6 (Backtrack)
                if self.lpart == 0 :
                    return
                self.lpart -= 1
                self.db_trace("Backtracked to")
            self.db_trace("decrement ok, about to expand")

    def count_partitions_slow(self, multiplicities):
        """Returns the number of partitions of a multiset whose elements
        have the multiplicities given in ``multiplicities``.

        This is mostly for comparison purposes.  It follows the same path as
        enumerate, and counts, rather than generates, the partitions.

        Examples
        ========

        XXX todo

        """
        # number of partitions so far in the enumeration
        self.pcount = 0
        self.initialize_enumeration(multiplicities)
        while True:
            while self.spread_part_multiplicity():
                pass

            # M4  Visit (count) a partition
            self.pcount += 1

            # M5 (Decrease v)
            while not self.decrement_part(self.top_part()):
                # M6 (Backtrack)
                if self.lpart == 0:
                    return self.pcount
                self.lpart -= 1

    def counting_decrement(self, part):
        """Helper for count_partitions -- XXX say what it does

        Examples
        ========

        XXX todo

        """
        status = self.decrement_part(part)
        if status:
            # Check if decremented part is in the cache -- if so,
            # increment pcount and treat as failure - we have used up
            # this part.
            pkey = part_key(self.top_part())
            if pkey in self.dp_map:
                self.pcount += self.dp_map[pkey]
                status = False
        return status

    def count_partitions(self, multiplicities):
        """Returns the number of partitions of a multiset whose elements
        have the multiplicities given in ``multiplicities``.

        For larger counts, this method is much faster than calling one
        of the enumerators and counting the result.  Uses dynamic
        programming to cut down on the number of nodes actually
        explored.

        Examples
        ========

        XXX todo

        Notes
        =====

        One can think of an enumeration of multiset partitions as
        operating on a binary tree of parts.  A part has (up to) two
        children, the left child resulting from the spread operation,
        and the right child from the decrement operation.  The
        ordinary enumeration of multiset partitions is an orderer
        traversal of this tree, and with the partitions corresponding
        to paths from the root to the leaves. The mapping from paths
        to partitions is not quite direct, since the partition would
        contain only those parts which are leaves or the parents of a
        spread link, not those which are parents of a decrement link.

        For counting purposes, it is sufficient to count leaves, and
        this can be done with a recursive ordered traversal.  The
        number of leaves of a subtree rooted at a particular part is a
        function only of that part itself, so memoizing has the
        potential to speed up the counting dramatically.

        This method follows a computational approach which is similar
        to the hypothetical memoized recursive function, but with two
        differences:

        1) This method is iterative, borrowing its structure from the
           other enumerations and maintaining an explicit stack of
           parts which are in the process of being counted.  (It is at
           least possible that there are some multisets for which the
           partitions are countable by this routine, but which would
           overflow the default Python recursion limit with a
           recursive implementation.)

        2) Instead of using the part data structure directly, an
           explicit key is constructed.  This allows some states to
           coalesce, which would remain separate with a physical key.
        """
        # number of partitions so far in the enumeration
        self.pcount = 0
        # dp_stack is list of lists of (part_key, start_count) pairs
        self.dp_stack = []

        # dp_map is map part_key-> count, where count represents the
        # number of multiset which are descendants of a part with this
        # key, **or any of its decrements**

        # Thus, when we find a part in the map, we add its count
        # value to the running total, cut off the enumeration, and
        # backtrack

        if not hasattr(self, 'dp_map'):
            self.dp_map = {}

        self.initialize_enumeration(multiplicities)
        pkey = part_key(self.top_part())
        self.dp_stack.append([(pkey, 0),])
        while True:
            while self.spread_part_multiplicity():
                pkey = part_key(self.top_part())
                if pkey in self.dp_map:
                    # Already have a cached value for the count of the
                    # subtree rooted at this part.  Add it to the
                    # running counter, and break out of the spread
                    # loop.  The -1 below is to compensate for the
                    # leaf that this code path would otherwise find,
                    # and which gets incremented for below.

                    self.pcount += (self.dp_map[pkey] - 1)
                    self.lpart -= 1
                    break
                else:
                    self.dp_stack.append([(pkey, self.pcount),])

            # M4  count a leaf partition
            self.pcount += 1

            # M5 (Decrease v)
            while not self.counting_decrement(self.top_part()):
                # M6 (Backtrack)
                for key, oldcount in self.dp_stack.pop():
                    self.dp_map[key] = self.pcount - oldcount
                if self.lpart == 0 :
                    return self.pcount
                self.lpart -= 1

            # At this point have a successfully decremented the part on
            # the stack which does not appear in the cache.  It needs
            # to be added to the list at the top of dp_stack
            pkey = part_key(self.top_part())
            self.dp_stack[-1].append((pkey, self.pcount),)

def part_key(part):
    """Helper for MultisetPartitionTraverser.count_partitions that
    creates a key for ``part``, that only includes information which can
    affect the count for that part.  (Any irrelevant information just
    reduces the effectiveness of dynamic programming.)

    Examples
    ========

    XXX todo

    """
    # The component number is irrelevant for counting partitions, so
    # leave it out of the memo key.
    rval = []
    for ps in part:
        rval.append(ps.u)
        rval.append(ps.v)
    return tuple(rval)

#
#  End of multiset partitions functions/classes
#
