from sympy.core import Basic
from sympy.utilities.iterables import flatten
from sympy.ntheory.residue_ntheory import int_tested

from collections import defaultdict

class Prufer(Basic):
    """
    The Prufer correspondence is an algorithm that describes the
    bijection between labeled trees and the Prufer code. A Prufer
    code of a labeled tree is unique up to isomorphism and has
    a length of n - 2.

    Prufer sequences were first used by Heinz Prufer to give a
    proof of Cayley's formula.

    References
    ==========

    .. [1] http://mathworld.wolfram.com/LabeledTree.html

    """
    _prufer_repr = None
    _tree_repr = None
    _nodes = None
    _rank = None

    @property
    def prufer_repr(self):
        """Returns Prufer sequence for the Prufer object.

        This sequence is found by removing the highest numbered vertex,
        recording the node it was attached to, and continuuing until only
        two verices remain. The Prufer sequence is the list of recorded nodes.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> Prufer([[0, 3], [1, 3], [2, 3], [3, 4], [4, 5]]).prufer_repr
        [4, 3, 3, 3]
        >>> Prufer([1, 0, 0]).prufer_repr
        [1, 0, 0]

        See Also
        ========
        to_prufer

        """
        if self._prufer_repr is None:
            self._prufer_repr = self.to_prufer(self._tree_repr[:], self.nodes)
        return self._prufer_repr

    @property
    def tree_repr(self):
        """Returns the tree representation of the Prufer object.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> Prufer([[0, 3], [1, 3], [2, 3], [3, 4], [4, 5]]).tree_repr
        [[0, 3], [1, 3], [2, 3], [3, 4], [4, 5]]
        >>> Prufer([1, 0, 0]).tree_repr
        [[4, 1], [3, 0], [2, 0], [1, 0]]

        See Also
        ========
        to_tree

        """
        if self._tree_repr is None:
            self._tree_repr = self.to_tree(self._prufer_repr[:])
        return self._tree_repr

    @property
    def nodes(self):
        """Returns the number of nodes in the tree.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> Prufer([[0, 3], [1, 3], [2, 3], [3, 4], [4, 5]]).nodes
        6
        >>> Prufer([1, 0, 0]).nodes
        5

        """
        return self._nodes

    @property
    def rank(self):
        """Returns the rank of the Prufer sequence.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> p = Prufer([[0, 3], [1, 3], [2, 3], [3, 4], [4, 5]])
        >>> p.rank
        993
        >>> p.next().rank
        994
        >>> p.prev().rank
        992

        See Also
        ========
        prufer_rank, next, prev

        """
        if self._rank is None:
            self._rank = self.prufer_rank()
        return self._rank

    @staticmethod
    def to_prufer(tree, n):
        """Return the Prufer sequence for a tree given as a list of edges.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> a = Prufer([[0, 1], [0, 2], [0, 3]])
        >>> a.prufer_repr
        [0, 0]
        >>> Prufer.to_prufer([[0, 1], [0, 2], [0, 3]], 4)
        [0, 0]

        """
        d = defaultdict(int)
        L = defaultdict(int)
        for edge in tree:

            # Increment the value of the corresponding
            # node in the degree list as we encounter an
            # edge involving it.
            d[edge[0]] += 1
            d[edge[1]] += 1
        for i in xrange(n - 2):
            # find first 1 from the right
            x = n - 1
            while d[x] != 1:
                x -= 1
            y = n - 1
            while True:
                # Keep reducing the index till we hit
                # an edge that is in the tree.
                e = sorted([x, y])
                if e in tree:
                    break
                y -= 1
            L[i] = y
            d[x] -= 1
            d[y] -= 1
            # Get rid of the edge that we just took
            # into account while reducing the values in
            # the degree list.
            tree.remove(e)
        return [L[i] for i in xrange(len(L))]

    @staticmethod
    def to_tree(prufer):
        """Return the tree (as a list of edges) of the given Prufer sequence.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> a = Prufer([0, 2], 4)
        >>> a.tree_repr
        [[3, 0], [1, 2], [2, 0]]
        >>> Prufer.to_tree([0, 2])
        [[3, 0], [1, 2], [2, 0]]

        """
        tree = []
        prufer.append(0)
        n = len(prufer)
        d = defaultdict(lambda : 1)
        for i in xrange(n - 1):
            d[prufer[i]] += 1
        for i in xrange(n):
            x = n
            # Find the node whose degree is one and
            # get the edge associated with it.
            while d[x] != 1:
                x -= 1
            y = prufer[i]
            # Since the edge has been removed reduce
            # the degree of all the nodes by one.
            d[x] -= 1
            d[y] -= 1
            tree.append([x, y])
        return tree

    def prufer_rank(self):
        """Computes the rank of a Prufer sequence.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> a = Prufer([[0, 1], [0, 2], [0, 3]])
        >>> a.prufer_rank()
        0

        """
        r = 0
        p = 1
        for i in xrange(self.nodes - 3, -1, -1):
            r += p*self.prufer_repr[i]
            p *= self.nodes
        return r

    @classmethod
    def unrank(self, rank, n):
        """Finds the unranked Prufer sequence.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> Prufer.unrank(0, 4)
        Prufer([0, 0])

        """
        n = int_tested(n)
        rank = int_tested(rank)
        L = defaultdict(int)
        for i in xrange(n - 3, -1, -1):
            L[i] = rank % n
            rank = (rank - L[i])//n
        return Prufer([L[i] for i in xrange(len(L))])

    def __new__(cls, *args, **kw_args):
        """The constructor for the Prufer object.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer

        A Prufer object can be constructed from a list of edges:

        >>> a = Prufer([[0, 1], [0, 2], [0, 3]])
        >>> a.prufer_repr
        [0, 0]

        If the number of nodes is given, no checking of the nodes will
        be performed; it will be assumed that nodes 0 through n - 1 are
        present:

        >>> Prufer([[0, 1], [0, 2], [0, 3]], 4)
        Prufer([[0, 1], [0, 2], [0, 3]], 4)

        A Prufer object can be constructed from a Prufer sequence:

        >>> b = Prufer([1, 3])
        >>> b.tree_repr
        [[2, 1], [1, 3], [3, 0]]

        """
        ret_obj = Basic.__new__(cls, *args, **kw_args)
        if isinstance(args[0][0], list):
            if len(args) > 1:
                nnodes = args[1]
            else:
                nodes = set(flatten(args[0]))
                nnodes = max(nodes) + 1
                if nnodes != len(nodes):
                    missing = sorted(set(range(nnodes)) - nodes)
                    if len(missing) == 1:
                        msg = 'Node %s is missing' % missing[0]
                    else:
                        msg = 'Nodes %s are missing' % missing
                    raise ValueError(msg)
            ret_obj._tree_repr = args[0]
            ret_obj._nodes = nnodes
        else:
            ret_obj._prufer_repr = args[0]
            ret_obj._nodes = len(ret_obj._prufer_repr) + 2
        return ret_obj

    def next(self, delta=1):
        """Generates the Prufer sequence that is delta beyond the current one.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> a = Prufer([[0, 1], [0, 2], [0, 3]])
        >>> b = a.next(1)
        >>> b.tree_repr
        [[3, 0], [2, 1], [1, 0]]
        >>> b.rank
        1

        """
        return Prufer.unrank(self.rank + delta, self.nodes)

    def prev(self, delta=1):
        """Generates the Prufer sequence that is -delta before the current one.

        Examples
        ========
        >>> from sympy.combinatorics.prufer import Prufer
        >>> a = Prufer([[0, 1], [1, 2], [2, 3], [1, 4]])
        >>> a.rank
        36
        >>> b = a.prev()
        >>> b
        Prufer([1, 2, 0])
        >>> b.rank
        35
        """
        return Prufer.unrank(self.rank - delta, self.nodes)
