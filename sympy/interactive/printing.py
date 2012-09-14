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
            formatter.for_type(cls, _pretty_print)

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
        If True, use pretty_print to stringify;
        if False, use sstrrepr to stringify.
    order: string or None
        There are a few different settings for this parameter:
        lex (default), which is lexographic order;
        grlex, which is graded lexographic order;
        grevlex, which is reversed graded lexographic order;
        old, which is used for compatibility reasons and for long expressions;
        None, which sets it to lex.
    use_unicode: boolean or None
        If True, use unicode characters;
        if False, do not use unicode characters.
    wrap_line: boolean
        If True, lines will wrap at the end;
        if False, they will not wrap but continue as one line.
    num_columns: int or None
        If int, number of columns before wrapping is set to num_columns;
        if None, number of columns before wrapping is set to terminal width.
    no_global: boolean
        If True, the settings become system wide;
        if False, use just for this console/session.
    ip: An interactive console
        This can either be an instance of IPython,
        or a class that derives from code.InteractiveConsole.

    Examples
    ========
    >>> from sympy.interactive import init_printing
    >>> from sympy import Symbol, sqrt
    >>> from sympy.abc import x, y
    >>> sqrt(5)
    sqrt(5)
    >>> init_printing(pretty_print=True) # doctest: +SKIP
    >>> sqrt(5) # doctest: +SKIP
      ___
    \/ 5
    >>> theta = Symbol('theta') # doctest: +SKIP
    >>> init_printing(use_unicode=True) # doctest: +SKIP
    >>> theta # doctest: +SKIP
    \u03b8
    >>> init_printing(use_unicode=False) # doctest: +SKIP
    >>> theta # doctest: +SKIP
    theta
    >>> init_printing(order='lex') # doctest: +SKIP
    >>> str(y + x + y**2 + x**2) # doctest: +SKIP
    x**2 + x + y**2 + y
    >>> init_printing(order='grlex') # doctest: +SKIP
    >>> str(y + x + y**2 + x**2) # doctest: +SKIP
    x**2 + x + y**2 + y
    >>> init_printing(order='grevlex') # doctest: +SKIP
    >>> str(y * x**2 + x * y**2) # doctest: +SKIP
    x**2*y + x*y**2
    >>> init_printing(order='old') # doctest: +SKIP
    >>> str(x**2 + y**2 + x + y) # doctest: +SKIP
    x**2 + x + y**2 + y
    >>> init_printing(num_columns=10) # doctest: +SKIP
    >>> x**2 + x + y**2 + y # doctest: +SKIP
    x + y +
    x**2 + y**2
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
