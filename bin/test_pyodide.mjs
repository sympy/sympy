import { argv } from 'process'
import { readdirSync } from 'fs'
import pyodide_pkg from '../pyodide/pyodide.js'

// Messing with fetch is needed to be able to use micropip in node.js. They
// prevent the error:
//
// File "/lib/python3.10/site-packages/pyodide/http.py", line 228, in pyfetch
//   from js import fetch as _jsfetch
// ImportError: cannot import name 'fetch' from 'js' (unknown location)
//
import fetch from 'node-fetch'
globalThis.fetch = fetch

let sympy
const fileNames = readdirSync('dist')
for (const fileName of fileNames) {
    if (fileName.startsWith('sympy') && fileName.endsWith('.whl')) {
        sympy = fileName
    }
}

const pyodide = await pyodide_pkg.loadPyodide()

await pyodide.loadPackage('micropip')

await pyodide.runPythonAsync(`
  import micropip
  await micropip.install('multipledispatch')
`)

await pyodide.loadPackage([
    'mpmath',  // provided by pyodide
    'numpy',  // built by pyodide
    `../dist/${sympy}`  // git version sympy
])

let split = 'None'
if (argv[2]) {
    if (argv[2].startsWith('--split=')) {
        split = `'${argv[2].slice(8)}'`
    }
}

pyodide.runPython(`
import sympy
if not sympy.test(split=${split}, subprocess=False):
    exit(1)
`)
