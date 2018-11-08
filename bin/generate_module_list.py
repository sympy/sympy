"""
Execute like this:

$ python bin/generate_module_list.py
modules = [
    'sympy.assumptions',
    'sympy.assumptions.handlers',
    'sympy.benchmarks',
    'sympy.calculus',
    'sympy.categories',
    'sympy.codegen',
    'sympy.combinatorics',
    'sympy.concrete',
    'sympy.core',
    'sympy.core.benchmarks',
    'sympy.crypto',
    'sympy.deprecated',
    'sympy.diffgeom',
    'sympy.external',
    'sympy.functions',
    'sympy.functions.combinatorial',
    'sympy.functions.elementary',
    'sympy.functions.elementary.benchmarks',
...
]

"""

from __future__ import print_function

from glob import glob


def get_paths(level=15):
    """
    Generates a set of paths for modules searching.

    Examples
    ========

    >>> get_paths(2)
    ['sympy/__init__.py', 'sympy/*/__init__.py', 'sympy/*/*/__init__.py']
    >>> get_paths(6)
    ['sympy/__init__.py', 'sympy/*/__init__.py', 'sympy/*/*/__init__.py',
    'sympy/*/*/*/__init__.py', 'sympy/*/*/*/*/__init__.py',
    'sympy/*/*/*/*/*/__init__.py', 'sympy/*/*/*/*/*/*/__init__.py']

    """
    wildcards = ["/"]
    for i in range(level):
        wildcards.append(wildcards[-1] + "*/")
    p = ["sympy" + x + "__init__.py" for x in wildcards]
    return p

def generate_module_list():
    g = []
    for x in get_paths():
        g.extend(glob(x))
    g = [".".join(x.split("/")[:-1]) for x in g]
    g = [i for i in g if not i.endswith('.tests')]
    g.remove('sympy')
    g = list(set(g))
    g.sort()
    return g

if __name__ == '__main__':
    g = generate_module_list()
    print("modules = [")
    for x in g:
        print("    '%s'," % x)
    print("]")
