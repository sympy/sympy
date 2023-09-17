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
    'micropip',
    'mpmath',  // provided by Pyodide
    'numpy',  // built by Pyodide
    'pytest',  // provided by Pyodide
    `./dist/${sympy}`  // git version SymPy
])
console.log('Pyodide packages loaded successfully')

const micropip = pyodide.pyimport('micropip');
await micropip.install(['pytest-split', 'hypothesis']);
console.log('Micropip packages loaded successfully')

let group = 'None'
if (argv[2]) {
    if (argv[2].startsWith('--group=')) {
        group = argv[2].slice(8)
    }
}
let splits = 'None'
if (argv[3]) {
    if (argv[3].startsWith('--splits=')) {
        splits = argv[3].slice(9)
    }
}
let sympyInstallPath = '/lib/python3.10/site-packages/sympy'
console.log('Variables set successfully')

pyodide.runPython(`
import os
import sys

import pytest

os.chdir('${sympyInstallPath}')

args = [
    '--rootdir',
    '${sympyInstallPath}',
    '-m',
    'not slow',
    '--group',
    '${group}',
    '--splits',
    '${splits}',
    '--durations',
    '20',
    '--ignore',
    '${sympyInstallPath}/integrals/rubi/rubi_tests/tests',
    '--ignore',
    '${sympyInstallPath}/testing/tests/test_runtests_pytest.py',
    '-W',
    'ignore::pytest.PytestUnknownMarkWarning',
    '${sympyInstallPath}',
]
exit_code = pytest.main(args)

if exit_code != pytest.ExitCode.OK:
    print(f'Exiting with pytest exit code {exit_code}')
    sys.exit(exit_code.value)
else:
    print(f'pytest finished with exit code {exit_code}')
`)
console.log('`pyodide.runPython` completed successfully')
