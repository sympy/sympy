"""Tools for setting up interactive sessions. """

from sympy.interactive.printing import init_printing

preexec_source = """\
from __future__ import division
from sympy import *
x, y, z, t = symbols('x y z t')
k, m, n = symbols('k m n', integer=True)
f, g, h = symbols('f g h', cls=Function)
"""

verbose_message = """\
These commands were executed:
%(source)s
Documentation can be found at http://www.sympy.org
"""

no_ipython = """\
Couldn't locate IPython. Having IPython installed is greatly recommended.
See http://ipython.scipy.org for more details. If you use Debian/Ubuntu,
just install the 'ipython' package and start isympy again.
"""

def _make_message(ipython=True, quiet=False, source=None):
    """Create a banner for an interactive session. """
    from sympy import __version__ as sympy_version
    from sympy.polys.domains import GROUND_TYPES
    from sympy.utilities.misc import ARCH
    from sympy import SYMPY_DEBUG

    import sys
    import os

    python_version = "%d.%d.%d" % sys.version_info[:3]

    if ipython:
        shell_name = "IPython"
    else:
        shell_name = "Python"

    info = ['ground types: %s' % GROUND_TYPES]

    cache = os.getenv('SYMPY_USE_CACHE')

    if cache is not None and cache.lower() == 'no':
        info.append('cache: off')

    if SYMPY_DEBUG:
        info.append('debugging: on')

    args = shell_name, sympy_version, python_version, ARCH, ', '.join(info)
    message = "%s console for SymPy %s (Python %s-%s) (%s)\n" % args

    if not quiet:
        if source is None:
            source = preexec_source

        _source = ""

        for line in source.split('\n')[:-1]:
            if not line:
                _source += '\n'
            else:
                _source += '>>> ' + line + '\n'

        message += '\n' + verbose_message % {'source': _source}

    return message

def _init_ipython_session(argv=[], auto=False):
    """Construct new IPython session. """
    import IPython

    if IPython.__version__ >= '0.11':
        # use an app to parse the command line, and init config
        from IPython.frontend.terminal import ipapp
        app = ipapp.TerminalIPythonApp()

        # don't draw IPython banner during initialization:
        app.display_banner = False
        app.initialize(argv)

        import re
        re_nameerror = re.compile("name '(?P<symbol>[A-Za-z_][A-Za-z0-9_]*)' is not defined")

        def _handler(self, etype, value, tb, tb_offset=None):
            """Handle :exc:`NameError` exception and allow injection of missing symbols. """
            if etype is NameError and tb.tb_next and not tb.tb_next.tb_next:
                match = re_nameerror.match(str(value))

                if match is not None:
                    self.run_cell("%(symbol)s = Symbol('%(symbol)s')" %
                        {'symbol': match.group("symbol")}, store_history=False)

                    try:
                        code = self.user_ns_hidden['In'][-1]
                    except (KeyError, IndexError):
                        pass
                    else:
                        self.run_cell(code, store_history=False)
                        return None

            stb = self.InteractiveTB.structured_traceback(etype, value, tb, tb_offset=tb_offset)
            self._showtraceback(etype, value, stb)

        if auto:
            app.shell.set_custom_exc((NameError,), _handler)

        return app.shell
    else:
        from IPython.Shell import make_IPython
        return make_IPython(argv)

def _init_python_session():
    """Construct new Python session. """
    from code import InteractiveConsole

    class SymPyConsole(InteractiveConsole):
        """An interactive console with readline support. """

        def __init__(self):
            InteractiveConsole.__init__(self)

            try:
                import readline
            except ImportError:
                pass
            else:
                import os
                import atexit

                readline.parse_and_bind('tab: complete')

                if hasattr(readline, 'read_history_file'):
                    history = os.path.expanduser('~/.sympy-history')

                    try:
                        readline.read_history_file(history)
                    except IOError:
                        pass

                    atexit.register(readline.write_history_file, history)

    return SymPyConsole()

def init_session(ipython=None, pretty_print=True, order=None,
        use_unicode=None, quiet=False, auto=False, argv=[]):
    """
    Initialize an embedded IPython or Python session.

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
    quiet: boolean
        If True, init_session will not print messages regarding its status;
        if False, init_session will print messages regarding its status.
    auto: boolean
        If True, init_session will automatically create symbols for you;
        if False, it will not.
    ipython: boolean or None
        If True, printing will initialize for an IPython console;
        if False, printing will initialize for a normal console;
        The default is None, which does what False does.
    argv: list of arguments for IPython
        See sympy.bin.isympy for options that can be used to intialize IPython.

    See Also
    ========

    sympy.interactive.printing.init_printing: for examples and the rest of the parameters.


    Examples
    ========

    >>> from sympy import init_session, Symbol, sin, sqrt
    >>> sin(x) #doctest: +SKIP
    NameError: name 'x' is not defined
    >>> init_session() #doctest: +SKIP
    >>> sin(x) #doctest: +SKIP
    sin(x)
    >>> sqrt(5) #doctest: +SKIP
      ___
    \/ 5
    >>> init_session(pretty_print=False) #doctest: +SKIP
    >>> sqrt(5) #doctest: +SKIP
    sqrt(5)
    >>> y + x + y**2 + x**2 #doctest: +SKIP
    x**2 + x + y**2 + y
    >>> init_session(order='grlex') #doctest: +SKIP
    >>> y + x + y**2 + x**2 #doctest: +SKIP
    x**2 + y**2 + x + y
    >>> init_session(order='grevlex') #doctest: +SKIP
    >>> y * x**2 + x * y**2 #doctest: +SKIP
    x**2*y + x*y**2
    >>> init_session(order='old') #doctest: +SKIP
    >>> x**2 + y**2 + x + y #doctest: +SKIP
    x + y + x**2 + y**2
    >>> theta = Symbol('theta') #doctest: +SKIP
    >>> theta #doctest: +SKIP
    theta
    >>> init_session(use_unicode=True) #doctest: +SKIP
    >>> theta # doctest: +SKIP
    \u03b8
    """
    import sys

    in_ipython = False

    if ipython is False:
        ip = _init_python_session()
        mainloop = ip.interact
    else:
        try:
            import IPython
        except ImportError:
            if ipython is not True:
                if not quiet:
                    print no_ipython
                ip = _init_python_session()
                mainloop = ip.interact
            else:
                raise RuntimeError("IPython is not available on this system")
        else:
            ipython = True

            if IPython.__version__ >= '0.11':
                try:
                    ip = get_ipython()
                except NameError:
                    ip = None
            else:
                ip = IPython.ipapi.get()
                if ip:
                    ip = ip.IP

            if ip is not None:
                in_ipython = True
            else:
                ip = _init_ipython_session(argv=argv, auto=auto)

            if IPython.__version__ >= '0.11':
                # runsource is gone, use run_cell instead, which doesn't
                # take a symbol arg.  The second arg is `store_history`,
                # and False means don't add the line to IPython's history.
                ip.runsource = lambda src, symbol='exec': ip.run_cell(src, False)
                mainloop = ip.mainloop
            else:
                mainloop = ip.interact

    if auto and (not ipython or IPython.__version__ < '0.11'):
        raise RuntimeError("automatic construction of symbols is possible only in IPython 0.11 or above")

    _preexec_source = preexec_source

    ip.runsource(_preexec_source, symbol='exec')
    init_printing(pretty_print=pretty_print, order=order, use_unicode=use_unicode, ip=ip)

    message = _make_message(ipython, quiet, _preexec_source)

    if not in_ipython:
        mainloop(message)
        sys.exit('Exiting ...')
    else:
        ip.write(message)
        ip.set_hook('shutdown_hook', lambda ip: ip.write("Exiting ...\n"))
