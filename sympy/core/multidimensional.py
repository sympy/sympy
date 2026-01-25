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
    Generalizes a function taking scalars to accept multidimensional arguments.
    
    This decorator allows a function that operates on scalar values to work with
    nested lists, arrays, or other iterable structures by applying the function
    element-wise.
    
    Examples
    ========
    >>> from sympy import vectorize, diff, sin, symbols, Function
    >>> x, y, z = symbols('x y z')
    >>> f, g, h = list(map(Function, 'fgh'))
    >>> @vectorize(0)
    ... def vsin(x):
    ...     return sin(x)
    >>> vsin([1, x, y])
    [sin(1), sin(x), sin(y)]
    
    >>> @vectorize(0, 1)
    ... def vdiff(f, y):
    ...     return diff(f, y)
    >>> vdiff([f(x, y, z), g(x, y, z), h(x, y, z)], [x, y, z])
    [[Derivative(f(x, y, z), x), Derivative(f(x, y, z), y), Derivative(f(x, y, z), z)], 
     [Derivative(g(x, y, z), x), Derivative(g(x, y, z), y), Derivative(g(x, y, z), z)], 
     [Derivative(h(x, y, z), x), Derivative(h(x, y, z), y), Derivative(h(x, y, z), z)]]
    """
    
    def __init__(self, *mdargs: Union[int, str]) -> None:
        """
        Initialize the vectorize decorator.
        
        Parameters
        ----------
        *mdargs : int or str
            The arguments that will be treated as data structures. Numbers indicate
            positional arguments, strings indicate keyword arguments.
            If no argument is given, everything is treated multidimensional.
            
        Raises
        ------
        TypeError
            If any argument is not an int or str
        """
        for a in mdargs:
            if not isinstance(a, (int, str)):
                raise TypeError("a is of invalid type")
        self.mdargs: tuple[Union[int, str], ...] = mdargs
    
    def __call__(self, f: Callable[..., Any]) -> Callable[..., Any]:
        """
        Returns a wrapper for the one-dimensional function that can handle
        multidimensional arguments.
        
        Parameters
        ----------
        f : Callable
            The function to decorate
            
        Returns
        -------
        Callable
            A wrapper function that handles multidimensional arguments
        """
        @wraps(f)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            # Get arguments that should be treated multidimensional
            if self.mdargs:
                mdargs = self.mdargs
            else:
                mdargs = tuple(range(len(args))) + tuple(kwargs.keys())
            
            arglength = len(args)
            
            for n in mdargs:
                if isinstance(n, int):
                    if n >= arglength:
                        continue
                    entry = args[n]
                    is_arg = True
                elif isinstance(n, str):
                    try:
                        entry = kwargs[n]
                    except KeyError:
                        continue
                    is_arg = False
                
                if hasattr(entry, "__iter__") and not isinstance(entry, str):
                    # Create now a copy of the given array and manipulate then
                    # the entries directly.
                    if is_arg:
                        args = list(args)
                        args[n] = structure_copy(entry)
                    else:
                        kwargs[n] = structure_copy(entry)
                    
                    result = apply_on_element(wrapper, args, kwargs, n)
                    return result
            
            return f(*args, **kwargs)
        
        return wrapper
