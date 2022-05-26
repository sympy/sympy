import pyodide_pkg from 'pyodide/pyodide.js'

const pyodide = await pyodide_pkg.loadPyodide()
await pyodide.loadPackage([
    'https://files.pythonhosted.org/packages/d4/cf/3965bddbb4f1a61c49aacae0e78fd1fe36b5dc36c797b31f30cf07dcbbb7/mpmath-1.2.1-py3-none-any.whl',  // latest mpmath on PyPI
    'https://cdn.jsdelivr.net/pyodide/v0.20.0/full/numpy-1.22.3-cp310-cp310-emscripten_wasm32.whl',  // latest numpy built by pyodide
    'http://localhost:8000/dist/sympy-1.11.dev0-py3-none-any.whl'  // git version sympy
])
pyodide.runPython(`
import sympy
if not sympy.test(subprocess=False):
    exit(1)
`)
