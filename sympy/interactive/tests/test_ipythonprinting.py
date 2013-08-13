"""Tests that the IPython printing module is properly loaded. """

from sympy.core.compatibility import u
from sympy.interactive.session import init_ipython_session
from sympy.external import import_module

# run_cell was added in IPython 0.11
ipython = import_module("IPython", min_module_version="0.11")

# disable tests if ipython is not present
if not ipython:
    disabled = True

def test_ipythonprinting():
    # Initialize and setup IPython session
    app = init_ipython_session()
    app.run_cell("ip = get_ipython()")
    app.run_cell("inst = ip.instance()")
    app.run_cell("format = inst.display_formatter.format")
    app.run_cell("from sympy import Symbol")

    # Printing without printing extension
    app.run_cell("a = format(Symbol('pi'))")
    app.run_cell("a2 = format(Symbol('pi')**2)")
    # Deal with API change starting at IPython 1.0
    if int(ipython.__version__.split(".")[0]) < 1:
        assert app.user_ns['a']['text/plain'] == "pi"
        assert app.user_ns['a2']['text/plain'] == "pi**2"
    else:
        assert app.user_ns['a'][0]['text/plain'] == "pi"
        assert app.user_ns['a2'][0]['text/plain'] == "pi**2"

    # Load printing extension
    app.run_cell("from sympy import init_printing")
    app.run_cell("init_printing()")
    # Printing with printing extension
    app.run_cell("a = format(Symbol('pi'))")
    app.run_cell("a2 = format(Symbol('pi')**2)")
    # Deal with API change starting at IPython 1.0
    if int(ipython.__version__.split(".")[0]) < 1:
        assert app.user_ns['a']['text/plain'] in (u('\u03c0'), 'pi')
        assert app.user_ns['a2']['text/plain'] in (u(' 2\n\u03c0 '), '  2\npi ')
    else:
        assert app.user_ns['a'][0]['text/plain'] in (u('\u03c0'), 'pi')
        assert app.user_ns['a2'][0]['text/plain'] in (u(' 2\n\u03c0 '), '  2\npi ')
