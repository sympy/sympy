"""Tools for setting up interactive sessions. """

verbose_message = """\
These commands were executed:
>>> from __future__ import division
>>> from sympy import *
>>> x, y, z = symbols('x,y,z')
>>> k, m, n = symbols('k,m,n', integer=True)
>>> f, g, h = symbols('f,g,h', cls=Function)

Documentation can be found at http://sympy.org/
"""

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

    ip.runcode(ip.compile("from __future__ import division"))
    ip.runcode(ip.compile("from sympy.interactive import *"))

    ip.runcode(ip.compile("init_printing(pretty_print=%s, order=%r, use_unicode=%s)" % (pretty_print, order, use_unicode)))

    from sympy import __version__ as sympy_version
    py_version = "%d.%d.%d" % sys.version_info[:3]

    if session == "ipython":
        py_name = "IPython"
    else:
        py_name = "Python"

    from sympy.polys.domains import GROUND_TYPES

    info = ['ground types: %s' % GROUND_TYPES]

    cache = os.getenv('SYMPY_USE_CACHE')

    if cache is not None and cache.lower() == 'no':
        info.append('cache: off')

    welcome = "%s console for SymPy %s (Python %s) (%s)" % \
        (py_name, sympy_version, py_version, ', '.join(info))

    if not quiet:
        message = welcome + '\n\n' + verbose_message
    else:
        message = welcome + '\n'

    if not in_ipyshell:
        ip.interact(message)
        sys.exit('Exiting ...')
    else:
        ip.write(message)
        ip.set_hook('shutdown_hook', lambda ip: ip.write("Exiting ...\n"))
