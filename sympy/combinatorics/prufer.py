from sympy.core import Basic

class Prufer(Basic):
    """
    The Prufer correspondence is an algorithm that describes the
    bijection between labelled trees and the Prufer code.
    """
    _prufer_repr = None

    @property
    def prufer_repr(self):
        return self._prufer_repr

    def __new__(cls, *args, **kw_args):
        """
        The constructor for the Prufer object.

        Examples:
        >>> from sympy.combinatorics.prufer import Prufer
        >>> a = Prufer([[0, 1], [0, 2], [0, 3]], 4)
        >>> a.prufer_repr
        [0, 0]
        """
        ret_obj = Basic.__new__(cls, *args, **kw_args)
        tree = args[0]
        n = args[1]
        d = [0] * n
        L = [0] * n
        for edge in tree:
            list.sort(edge)
            d[edge[0]] += 1
            d[edge[1]] += 1
        for i in xrange(0, n - 2):
            x = n - 1
            while d[x] != 1:
                x -= 1
            y = n - 1
            while True:
                e = [x, y]
                list.sort(e)
                if e in tree:
                    break
                y -= 1
            L[i] = y
            d[x] -= 1
            d[y] -= 1
            tree.remove(e)
        ret_obj._prufer_repr = L[:2]
        return ret_obj
