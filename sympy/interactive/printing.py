"""Tools for setting up printing in interactive sessions. """

def _init_python_printing(stringify_func):
    """Setup printing in Python interactive session. """
    import __builtin__, sys

    def displayhook(arg):
        """Python's pretty-printer display hook.

           This function was adapted from:

            http://www.python.org/dev/peps/pep-0217/

        """
        if arg is not None:
            __builtin__._ = None
            print stringify_func(arg)
            __builtin__._ = arg

    sys.displayhook = displayhook

def _init_ipython_printing(ip, stringify_func):
    """Setup printing in IPython interactive session. """

    def result_display(self, arg):
        """IPython's pretty-printer display hook.

           This function was adapted from:

            ipython/IPython/hooks.py:155

        """
        if self.rc.pprint:
            out = stringify_func(arg)

            if '\n' in out:
                print

            print out
        else:
            print repr(arg)

    ip.set_hook('result_display', result_display)

def init_printing(pretty_print=True, order=None, use_unicode=None, wrap_line=None):
    """Initializes pretty-printer depending on the environment. """
    from sympy.printing.printer import Printer

    if pretty_print in (True, False):
        if pretty_print:
            from sympy.printing import pretty as stringify_func
        else:
            from sympy.printing import sstrrepr as stringify_func
    else:
        stringify_func = pretty_print

    Printer.set_global_settings(order=order, use_unicode=use_unicode, wrap_line=wrap_line)

    try:
        import IPython
    except ImportError:
        _init_python_printing(stringify_func)
    else:
        ip = IPython.ipapi.get()

        if ip is not None:
            _init_ipython_printing(ip, stringify_func)
        else:
            _init_python_printing(stringify_func)
