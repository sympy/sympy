function StartTest($name) {
    Add-AppveyorTest -Name "$name Test" -Outcome Running
    Write-Output "Testing $name"
}

function EndTest($name) {
    If ($output | Select-String -Pattern "Exception: Tests failed") {
        Update-AppveyorTest -Name "$name Test" -Outcome Failed -StdOut $output
        # Uncomment this to fail the build when test fails
        # throw "$name test failed!"
    } Else {
        Update-AppveyorTest -Name "$name Test" -Outcome Passed -StdOut $output
    }
}

# lambdify with tensorflow and numexpr is tested here
If ($env:TEST_OPT_DEPENDENCY -match "numpy") {
    StartTest "NUMPY"
    @"
import sympy
if not (sympy.test('*numpy*', 'sympy/matrices/', 'sympy/physics/quantum/',
        'sympy/core/tests/test_numbers.py', 'sympy/core/tests/test_sympify.py',
        'sympy/utilities/tests/test_lambdify.py',
        blacklist=['sympy/physics/quantum/tests/test_circuitplot.py']) and sympy.
         doctest('sympy/matrices/', 'sympy/utilities/lambdify.py')):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
    EndTest "NUMPY"
}

If ($env:TEST_OPT_DEPENDENCY -match "scipy") {
    StartTest "SCIPY"
    @"
import sympy
# scipy matrices are tested in numpy testing
if not sympy.test('sympy/external/tests/test_scipy.py'):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
    EndTest "SCIPY"
}

If ($env:TEST_OPT_DEPENDENCY -match "llvmlite") {
    StartTest "LLVMJIT"
    @"
import sympy
if not (sympy.test('sympy/printing/tests/test_llvmjit.py')
        and sympy.doctest('sympy/printing/llvmjitcode.py')):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
    EndTest "LLVMJIT"
}

If ($env:TEST_OPT_DEPENDENCY -match "theano") {
    StartTest "THEANO"
    @"
import sympy
if not sympy.test('*theano*'):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
    EndTest "THEANO"
}

If ($env:TEST_OPT_DEPENDENCY -match "gmpy") {
    StartTest "GMPY"
    @"
import sympy
if not (sympy.test('sympy/polys/') and sympy.doctest('sympy/polys/')):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
    EndTest "GMPY"
}

If ($env:TEST_OPT_DEPENDENCY -match "matplotlib") {
    StartTest "MATPLOTLIB"
    @"
# Set matplotlib so that it works correctly in headless environment. We have to do
# this here because it doesn't work after the sympy plotting module is
# imported.
import matplotlib
matplotlib.use("Agg")
import sympy
# Unfortunately, we have to use subprocess=False so that the above will be
# applied, so no hash randomization here.
if not (sympy.test('sympy/plotting', 'sympy/physics/quantum/tests/test_circuitplot.py',
    subprocess=False) and sympy.doctest('sympy/plotting', subprocess=False)):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
    EndTest "MATPLOTLIB"
}

If ($env:TEST_OPT_DEPENDENCY -match "ipython") {
    StartTest "IPYTHON"
    @"
import sympy
if not sympy.test('*ipython*'):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
    EndTest "IPYTHON"
}

If ($env:TEST_OPT_DEPENDENCY -match "symengine") {
    StartTest "SYMENGINE"
    $env:USE_SYMENGINE=1
    @"
import sympy
if not sympy.test('sympy/physics/mechanics'):
    raise Exception('Tests failed')
if not sympy.test('sympy/liealgebras'):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
    Remove-Item Env:\USE_SYMENGINE
    EndTest "SYMENGINE"
}
