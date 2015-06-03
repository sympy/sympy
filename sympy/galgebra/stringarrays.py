# sympy/galgebra/stringarrays.py

"""
stringarrays.py are a group of helper functions to convert string
input to vector and multivector class function to arrays of SymPy
symbols.
"""

import operator

from sympy.core.compatibility import reduce
from itertools import combinations

from sympy import S, Symbol, Function
from sympy.core.compatibility import range


def str_array(base, n=None):
    """
    Generate one dimensional (list of strings) or two dimensional (list
    of list of strings) string array.

    For one dimensional arrays: -

        base is string of variable names separated by blanks such as
        base = 'a b c' which produces the string list ['a','b','c'] or
        it is a string with no blanks than in conjunction with the
        integer n generates -

            str_array('v',n=-3) = ['v_1','v_2','v_3']
            str_array('v',n=3) = ['v__1','v__2','v__3'].

        In the case of LaTeX printing the '_' would give a subscript and
        the '__' a super script.

    For two dimensional arrays: -

        base is string where elements are separated by spaces and rows by
        commas so that -

            str_array('a b,c d') = [['a','b'],['c','d']]

    """
    if n is None:
        if ',' in base:
            base_array = []
            base_split = base.split(',')
            for base_arg in base_split:
                base_array.append(list(filter(lambda x: x != '', base_arg.split(' '))))
            return base_array
        else:
            return base.split(' ')
    result = []
    if isinstance(n, str):
        if n[0] == '-':
            for index in n[1:].split(' '):
                result.append(base + '_' + index)
        if n[0] == '+':
            for index in n[1:].split(' '):
                result.append(base + '__' + index)
    if n > 0:
        for i in range(1, n + 1):
            result.append(base + '__' + str(i))
    if n < 0:
        for i in range(1, -n + 1):
            result.append(base + '_' + str(i))
    return result


def symbol_array(base, n=None):
    """
    Generates a string arrary with str_array and replaces each string in
    array with Symbol of same name.
    """
    symbol_str_lst = str_array(base, n)
    result = []
    for symbol_str in symbol_str_lst:
        result.append(S(symbol_str))
    return tuple(result)


def fct_sym_array(str_lst, coords=None):
    """
    Construct list of symbols or functions with names in 'str_lst'.  If
    'coords' are given (tuple of symbols) function list constructed,
    otherwise a symbol list is constructed.
    """
    if coords is None:
        fs_lst = []
        for sym_str in str_lst:
            fs_lst.append(Symbol(sym_str))
    else:
        fs_lst = []
        for fct_str in str_lst:
            fs_lst.append(Function(fct_str)(*coords))
    return fs_lst


def str_combinations(base, lst, rank=1, mode='_'):
    """
    Construct a list of strings of the form 'base+mode+indexes' where the
    indexes are formed by converting 'lst' to a list of strings and then
    forming the 'indexes' by concatenating combinations of elements from
    'lst' taken 'rank' at a time.
    """
    a1 = combinations([str(x) for x in lst], rank)
    a2 = [reduce(operator.add, x) for x in a1]
    str_lst = [base + mode + x for x in a2]
    return str_lst
