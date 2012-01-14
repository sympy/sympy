"""Tools for setting up printing in interactive sessions. """

def _init_python_printing(stringify_func):
    """Setup printing in Python interactive session. """
    import __builtin__, sys

    def _displayhook(arg):
        """Python's pretty-printer display hook.

           This function was adapted from:

            http://www.python.org/dev/peps/pep-0217/

        """
        if arg is not None:
            __builtin__._ = None
            print stringify_func(arg)
            __builtin__._ = arg

    sys.displayhook = _displayhook

def _init_ipython_printing(ip, stringify_func):
    """Setup printing in IPython interactive session. """

    def _pretty_print(arg, p, cycle):
        """caller for pretty, for use in IPython 0.11"""
        p.text(stringify_func(arg))

    def _result_display(self, arg):
        """IPython's pretty-printer display hook, for use in IPython 0.10

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

    import IPython
    if IPython.__version__ >= '0.11':
        formatter = ip.display_formatter.formatters['text/plain']

        for cls in (object, tuple, list, set, frozenset, dict, str):
            formatter.for_type(cls, pretty_print)

        # this loads pretty printing for objects that inherit from Basic or Matrix:
        formatter.for_type_by_name(
            'sympy.core.basic', 'Basic', _pretty_print
        )
        formatter.for_type_by_name(
            'sympy.matrices.matrices', 'Matrix', _pretty_print
        )
    else:
        ip.set_hook('result_display', _result_display)

def init_printing(pretty_print=True, order=None, use_unicode=None, wrap_line=None, num_columns=None, no_global=False, ip=None):
    """
    Initializes pretty-printer depending on the environment.

    Parameters
    ==========

    pretty_print: boolean
        If true, use pretty print to stringify, if false, use sstrrepr to stringify.
    order: boolean or string
        Set to 'none' for long expressions if slow; default is None.
    use_unicode: boolean or none
        Use unicode characters instead of string characters.
    wrap_line: boolean
        Is line wrapping enabled or not?
    num_columns: boolean
        Number of columns before before line breaking;
        defaults to None which reads terminal width.
    no_global: boolean
        Whether to make this a global printer or not.
    ip: bool or None
        If true, printing is set up specifically for ipython;
        if false or None, printing is initialized for a normal python console.

    Examples
    ========
    >>> from sympy.interactive import init_printing
    >>> from sympy import Symbol
    >>> init_printing(pretty_print=True)
    >>> sqrt(5) # doctest: +SKIP
      ___
    \/ 5
    >>> init_printing(pretty_print=False)
    >>> sqrt(5) # doctest: +SKIP
    sqrt(5)
    >>> theta = Symbol('theta')
    >>> init_printing(use_unicode=True)
    >>> theta # doctest: +SKIP
    u'\u03b8'
    >>> init_printing(use_unicode=False)
    >>> theta # doctest: +SKIP
    theta
    """
    from sympy.printing.printer import Printer

    if pretty_print:
        from sympy.printing import pretty as stringify_func
    else:
        from sympy.printing import sstrrepr as stringify_func

    if not no_global:
        Printer.set_global_settings(order=order, use_unicode=use_unicode, wrap_line=wrap_line, num_columns=num_columns)
    else:
        _stringify_func = stringify_func

        if pretty_print:
            stringify_func = lambda expr: _stringify_func(expr, order=order, use_unicode=use_unicode, wrap_line=wrap_line, num_columns=num_columns)
        else:
            stringify_func = lambda expr: _stringify_func(expr, order=order)

    if ip is not None and ip.__module__.startswith('IPython'):
        _init_ipython_printing(ip, stringify_func)
    else:
        _init_python_printing(stringify_func)
