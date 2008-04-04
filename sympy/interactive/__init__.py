
from sympy import *

x, y, z = symbols('xyz')
k, m, n = symbols('kmn', integer=True)

f = Function("f")

def init_ipython():
    import os, sys

    import IPython.ipapi
    ip = IPython.ipapi.get()

    ip.IP.compile("from __future__ \
        import division", "<input>", "single")

    def result_display(self, arg):
        """Pretty-printer display hook.

           Called for displaying pretty results to the user. Using this
           handler not only SymPy's  expression can be printed but also
           Python's lists, tuples and dictionaries.

           This function was adapted from:

             ipython/IPython/hooks.py:155

        """
        if self.rc.pprint:
            out = pretty(arg)

            if '\n' in out:
                print

            print out
        else:
            print repr(arg)

    ip.set_hook('result_display', result_display)

    from sympy import __version__ as sympy_version
    py_version = "%d.%d.%d" % sys.version_info[:3]

    welcome = "Python %s console for SymPy %s" \
        % (py_version, sympy_version)

    if os.getenv('SYMPY_USE_CACHE') == 'no':
        welcome += ' (cache: off)'

    def late_startup_hook(self):
        print welcome

    ip.set_hook('late_startup_hook', late_startup_hook)

    def shutdown_hook(self):
        print "Exiting ..."

    ip.set_hook('shutdown_hook', shutdown_hook)
