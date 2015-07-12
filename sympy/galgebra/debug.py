# sympy/galgebra/debug.py

from __future__ import division, print_function

from itertools import islice


def ostr(obj, dict_mode=False):
    """
    Recursively convert iterated object (list/tuple/dict/set) to string.
    """
    def ostr_rec(obj, dict_mode):
        global ostr_s
        if isinstance(obj, tuple):
            if len(obj) == 0:
                ostr_s += '(),'
            else:
                ostr_s += '('
                for obj_i in obj:
                    ostr_rec(obj_i, dict_mode)
                ostr_s = ostr_s[:-1] + '),'
        elif isinstance(obj, list):
            if len(obj) == 0:
                ostr_s += '[],'
            else:
                ostr_s += '['
                for obj_i in obj:
                    ostr_rec(obj_i, dict_mode)
                ostr_s = ostr_s[:-1] + '],'
        elif isinstance(obj, dict):
            if dict_mode:
                ostr_s += '\n'
                for key in obj.keys():
                    ostr_rec(key, dict_mode)
                    if ostr_s[-1] == ',':
                        ostr_s = ostr_s[:-1]
                    ostr_s += ' -> '
                    ostr_rec(obj[key], dict_mode)
                    if ostr_s[-1] == ',':
                        ostr_s = ostr_s[:-1]
                    ostr_s += '\n'
            else:
                ostr_s += '{'
                for key in obj.keys():
                    ostr_rec(key, dict_mode)
                    if ostr_s[-1] == ',':
                        ostr_s = ostr_s[:-1]
                    ostr_s += ':'
                    ostr_rec(obj[key], dict_mode)
                ostr_s = ostr_s[:-1] + '} '
        elif isinstance(obj, set):
            tmp_obj = list(obj)
            ostr_s += '{'
            for obj_i in tmp_obj:
                ostr_rec(obj_i, dict_mode)
            ostr_s = ostr_s[:-1] + '},'
        else:
            ostr_s += str(obj) + ','
        return
    global ostr_s
    ostr_s = ''
    if isinstance(obj, (tuple, list, dict, set)):
        ostr_rec(obj, dict_mode)
        return ostr_s[:-1]
    else:
        return str(obj)


def oprint(*args, **kwargs):
    """
    Debug printing for iterated (list/tuple/dict/set) objects. args is
    of form (title1,object1,title2,object2,...) and prints:

        title1 = object1
        title2 = object2
        ...

    If you only wish to print a title set object = None.
    """

    if 'dict_mode' in kwargs:
        dict_mode = kwargs['dict_mode']
    else:
        dict_mode = False

    if isinstance(args[0], str) or args[0] is None:
        titles = list(islice(args, None, None, 2))
        objs = tuple(islice(args, 1, None, 2))
        if len(args) > 2:
            if objs[0] is None:
                n = 0
            else:
                n = len(titles[0])
            for (title, obj) in zip(titles[1:], objs[1:]):
                if obj is not None:
                    if not (dict_mode and isinstance(obj, dict)):
                        n = max(n, len(title))
        else:
            n = len(titles[0])

        for (title, obj) in zip(titles, objs):
            if obj is None:
                print(title)
            else:
                npad = n - len(title)
                if isinstance(obj, dict):
                    print(title + ':' + ostr(obj, dict_mode))
                else:
                    print(title + npad * ' ' + ' = ' + ostr(obj, dict_mode))
    else:
        for arg in args:
            print(ostr(arg, dict_mode))
    return


def print_sub_table(title, keys, sdict, blade_rep=True):
    """
    Print substitution dictionary, sdict, according to order of keys in
    keys
    """
    if title is not None:
        print(title)
    for key in keys:
        print(str(key) + ' = ' + ostr(sdict[key]))
    return


def print_product_table(title, keys, pdict, op='*', blade_rep=True):
    """
    Print product dictionary, pdict, according to order of keys in keys
    """
    if title is not None:
        print(title)
    pop = ')' + op + '('
    for key1 in keys:
        for key2 in keys:
            print('(' + str(key1) + pop + str(key2) + ') = ' + ostr(pdict[(key1, key2)]))
    return
