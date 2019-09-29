import sys
from subprocess import Popen, PIPE

def test_debug_output():
    env = {'SYMPY_DEBUG':'True'}
    cmd = 'from sympy import *; x = Symbol("x"); print(integrate((1-cos(x))/x, x))'
    cmdline = [sys.executable, '-c', cmd]
    proc = Popen(cmdline, env=env, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    out = out.decode('ascii') # utf-8?
    err = err.decode('ascii')
    expected = 'substituted: -x*(cos(x) - 1), u: 1/x, u_var: _u'
    assert expected in err
