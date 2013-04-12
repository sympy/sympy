"""Tests that the IPython printing module is properly loaded. """

from sympy.interactive.session import init_ipython_session
from sympy.external import import_module

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
    assert app.user_ns['a']['text/plain'] == "pi"
    assert app.user_ns['a2']['text/plain'] == "pi**2"

    # Load printing extension
    app.run_cell("%load_ext sympy.interactive.ipythonprinting")
    # Printing with printing extension
    app.run_cell("a = format(Symbol('pi'))")
    app.run_cell("a2 = format(Symbol('pi')**2)")
    assert app.user_ns['a']['text/plain'] in (u'\u03c0', 'pi')
    assert app.user_ns['a2']['text/plain'] in (u' 2\n\u03c0 ', '  2\npi ')
