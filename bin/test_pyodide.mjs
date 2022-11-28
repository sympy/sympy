import { argv } from 'process'
import { readdirSync } from 'fs'
import pyodide_pkg from '../pyodide/pyodide.js'

let sympy
const fileNames = readdirSync('dist')
for (const fileName of fileNames) {
    if (fileName.startsWith('sympy') && fileName.endsWith('.whl')) {
        sympy = fileName
    }
}

const pyodide = await pyodide_pkg.loadPyodide()
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
