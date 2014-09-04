"""
This module simplifies the saving and loading of sympy expressions.
"""

import re
import sympy


class SympyDict(dict):
    """
    This saves a dictonary of sympy expressions to a file
    in human readable form.
    >>> a, b = sympy.symbols('a, b')
    >>> d = SympyDict({'a':a, 'b':b})
    >>> d.save('name.sympy')
    >>> del d
    >>> d2 = SympyDict.load('name.sympy')
    """

    def __init__(self, *args, **kwargs):
        super(SympyDict, self).__init__(*args, **kwargs)

    def __repr__(self):
        d = dict(self)
        for key in d.keys():
            d[key] = sympy.srepr(d[key])
        # regex is just used here to insert a new line after
        # each dict key, value pair to make it more readable
        return re.sub('(: \"[^"]*\",)', r'\1\n',  d.__repr__())

    def save(self, file):
        with open(file, 'wb') as savefile:
            savefile.write(self.__repr__())

    @classmethod
    def load(cls, file_path):
        with open(file_path, 'r') as loadfile:
            exec('d =' + loadfile.read())
        d = locals()['d']
        for key in d.keys():
            d[key] = sympy.sympify(d[key])
        return cls(d)
