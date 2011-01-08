"""Tools for setting up interactive sessions. """

from sympy.interactive.printing import init_printing

preexec_code = """\
from __future__ import division
from sympy import *
x, y, z = symbols('x,y,z')
k, m, n = symbols('k,m,n', integer=True)
f, g, h = symbols('f,g,h', cls=Function)
"""

verbose_message = """\
These commands were executed:
%(code)s
Documentation can be found at http://sympy.org/
"""

def _make_message(ipython=True, quiet=False):
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
        code = ""

        for line in preexec_code.split('\n')[:-1]:
            if not line:
                code += '\n'
            else:
                code += '>>> ' + line + '\n'

        message += '\n' + verbose_message % {'code': code}

    return message

def init_session(session="ipython", pretty_print=True, order=None, use_unicode=None, quiet=False, argv=[]):
    """Initialize embedded IPython or Python session. """
    import os, sys

    def init_IPython():
        return IPython.Shell.make_IPython(argv)

    def init_Python():
        import code

        class HistoryConsole(code.InteractiveConsole):
            def __init__(self):
                code.InteractiveConsole.__init__(self)

                history = os.path.expanduser('~/.sympy-history')

                try:
                    import readline, atexit

                    readline.parse_and_bind('tab: complete')

                    if hasattr(readline, 'read_history_file'):
                        try:
                            readline.read_history_file(history)
                        except IOError:
                            pass

                        atexit.register(readline.write_history_file, history)
                except ImportError:
                    pass

        return HistoryConsole()

    if session not in ['ipython', 'python']:
        raise ValueError("'%s' is not a valid session name" % session)

    in_ipyshell = False

    try:
        import IPython

        ip = IPython.ipapi.get()

        if ip is not None:
            if session == 'ipython':
                ip, in_ipyshell = ip.IP, True
            else:
                raise ValueError("Can't start Python shell from IPython")
        else:
            if session == 'ipython':
                ip = init_IPython()
            else:
                ip = init_Python()
    except ImportError:
        if session == 'ipython':
            raise
        else:
            ip = init_Python()

    ip.runsource(preexec_code, symbol='exec')
    init_printing(pretty_print=pretty_print, order=order, use_unicode=use_unicode)

    ipython = session == 'ipython'
    message = _make_message(ipython, quiet)

    if not in_ipyshell:
        ip.interact(message)
        sys.exit('Exiting ...')
    else:
        ip.write(message)
        ip.set_hook('shutdown_hook', lambda ip: ip.write("Exiting ...\n"))
