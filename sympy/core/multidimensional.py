"""Provides functionality for multidimensional usage of scalar-functions.

Read the vectorize docstring for more details.
"""

from __future__ import annotations

from functools import wraps
from typing import Any, Callable, Union


def apply_on_element(
    f: Callable[..., Any],
    args: list[Any],
    kwargs: dict[str, Any],
    n: Union[int, str]
) -> list[Any]:
    """
    Returns a structure with the same dimension as the specified argument,
    where each basic element is replaced by the function f applied on it. All
    other arguments stay the same.

    Parameters
    ----------
    f : Callable
        The function to apply to each element
    args : list
        Positional arguments
    kwargs : dict
        Keyword arguments
    n : int or str
        Index or key of the argument to treat multidimensionally

    Returns
    -------
    list
        Structure with function applied to all elements
    """
    # Get the specified argument.
    if isinstance(n, int):
        structure = args[n]
        is_arg = True
    elif isinstance(n, str):
        structure = kwargs[n]
        is_arg = False

    # Define reduced function that is only dependent on the specified argument.
    def f_reduced(x: Any) -> Any:
        if hasattr(x, "__iter__") and not isinstance(x, str):
            return list(map(f_reduced, x))
        else:
            if is_arg:
                args[n] = x
            else:
                kwargs[n] = x
            return f(*args, **kwargs)

    # f_reduced will call itself recursively so that in the end f is applied to
    # all basic elements.
    return list(map(f_reduced, structure))


def iter_copy(structure: list[Any]) -> list[Any]:
    """
    Returns a copy of an iterable object (also copying all embedded iterables).

    Parameters
    ----------
    structure : list
        An iterable structure to copy

    Returns
    -------
    list
        A deep copy of the iterable structure
    """
    return [iter_copy(i) if hasattr(i, "__iter__") and not isinstance(i, str) else i
            for i in structure]


def structure_copy(structure: Any) -> Any:
    """
    Returns a copy of the given structure (numpy-array, list, iterable, ..).

    Parameters
    ----------
    structure : Any
        The structure to copy

    Returns
    -------
    Any
        A copy of the structure
    """
    if hasattr(structure, "copy"):
        return structure.copy()
    return iter_copy(structure)


class vectorize:
    """
    Generalizes a function of a scalar to a function of an iterable.

    For example:

    >>> from sympy.core.multidimensional import vectorize
    >>> from sympy import sin
    >>> vsin = vectorize(0)(sin)
    >>> vsin([1, 2, 3])
    [sin(1), sin(2), sin(3)]
    >>> vsin([[1, 2], [3, 4]])
    [[sin(1), sin(2)], [sin(3), sin(4)]]

    Parameters
    ----------
    mdargs : int or list of int/str
        The indices or keys of the arguments that are to be treated
        multidimensionally.

    Examples
    --------

    >>> from sympy.core.multidimensional import vectorize
    >>> from sympy import diff, symbols, sin
    >>> x, y = symbols('x y')
    >>> f = vectorize(0, 1)(diff)
    >>> f(sin(x), x)
    cos(x)
    >>> f([sin(x), x], [x, y])
    [[cos(x), 0], [1, 0]]

    """
    def __init__(self, *mdargs: Union[int, str]) -> None:
        self.mdargs = mdargs

    def __call__(self, f: Callable[..., Any]) -> Callable[..., Any]:
        """
        Returns a wrapper for the function f.
        """
        @wraps(f)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            # We copy the arguments so that the original ones are not modified.
            # We only need to copy the arguments that are to be treated
            # multidimensionally.
            args_list = list(args)
            for n in self.mdargs:
                if isinstance(n, int):
                    args_list[n] = structure_copy(args_list[n])
                elif isinstance(n, str):
                    kwargs[n] = structure_copy(kwargs[n])

            # Apply the function.
            res = f(*args_list, **kwargs)
            for n in self.mdargs:
                res = apply_on_element(f, args_list, kwargs, n)
            return res
        return wrapper
