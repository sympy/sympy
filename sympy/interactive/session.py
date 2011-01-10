"""Tools for setting up interactive sessions. """

from sympy.interactive.printing import init_printing

preexec_source = """\
from __future__ import division
from sympy import *
x, y, z = symbols('x,y,z')
k, m, n = symbols('k,m,n', integer=True)
f, g, h = symbols('f,g,h', cls=Function)
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

    args = shell_name, sympy_version, python_version, ', '.join(info)
    message = "%s console for SymPy %s (Python %s) (%s)\n" % args

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

def _init_ipython_session(argv=[]):
    """Construct new IPython session. """
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
        use_unicode=None, quiet=False, keep_sign=False, argv=[]):
    """Initialize an embedded IPython or Python session. """
    import sys

    in_ipython = False

    if ipython is False:
        ip = _init_python_session()
    else:
        try:
            import IPython
        except ImportError:
            if ipython is not True:
                print no_ipython
                ip = _init_python_session()
            else:
                raise RuntimeError("IPython is not available on this system")
        else:
            ip = IPython.ipapi.get()
            ipython = True

            if ip is not None:
                ip, in_ipython = ip.IP, True
            else:
                ip = _init_ipython_session(argv)

    _preexec_source = preexec_source

    if keep_sign:
        _preexec_source += "Basic.keep_sign = True\n"

    ip.runsource(_preexec_source, symbol='exec')
    init_printing(pretty_print=pretty_print, order=order, use_unicode=use_unicode)

    message = _make_message(ipython, quiet, _preexec_source)

    if not in_ipython:
        ip.interact(message)
        sys.exit('Exiting ...')
    else:
        ip.write(message)
        ip.set_hook('shutdown_hook', lambda ip: ip.write("Exiting ...\n"))
