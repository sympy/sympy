Add-AppveyorTest -Name "Tests" -Outcome Running
If (!$env:TEST_SLOW) {
    @"
print('Testing SYMPY')
import sympy
if not sympy.test():
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
} Else {
    @"
print('Testing SLOW')
import sympy
if not sympy.test(slow=True, verbose=True):
    raise Exception('Tests failed')
"@ | python -We:invalid 2>&1 | Tee-Object -Variable output
}
$err = $output | ?{ $_ -is [System.Management.Automation.ErrorRecord] }
If (!$err) {
    Update-AppveyorTest -Name "Tests" -Outcome Passed -StdOut $output
} Else {
    Update-AppveyorTest -Name "Tests" -Outcome Failed -StdOut $output -StdErr $err
}
