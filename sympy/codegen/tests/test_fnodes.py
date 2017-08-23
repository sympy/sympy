import os
from sympy import Symbol
from sympy.codegen.ast import (
    Assignment, Print, Declaration, FunctionDefinition, Return, real,
    FunctionCall, Variable, Element, integer
)
from sympy.codegen.fnodes import (
    isign, dsign, cmplx, kind, literal_dp, Program, Module, use, Subroutine,
    dimension, assumed_extent, intent_out, size, Do, SubroutineCall, sum_
)
from sympy.printing.fcode import fcode
from sympy.utilities._compilation import compile_run_strings, has_fortran
from sympy.utilities.pytest import skip


def test_size():
    x = Symbol('x', real=True)
    sx = size(x)
    assert fcode(sx, source_format='free') == 'size(x)'


def test_Program():
    x = Symbol('x', real=True)
    vx = Variable.deduced(x, 42)
    decl = Declaration(vx)
    prnt = Print([x, x+1])
    prog = Program('foo', [decl, prnt])
    if not has_fortran():
        skip("No fortran compiler found.")

    (stdout, stderr), info = compile_run_strings([('main.f90', fcode(prog, standard=90))], clean=True)
    assert '42' in stdout
    assert '43' in stdout
    assert stderr == ''
    assert info['exit_status'] == os.EX_OK


def test_Module():
    x = Symbol('x', real=True)
    v_x = Variable.deduced(x)
    sq = FunctionDefinition(real, 'sqr', [v_x], [Return(x**2)])
    mod_sq = Module('mod_sq', [], [sq])
    sq_call = FunctionCall('sqr', [42.])
    prg_sq = Program('foobar', [
        use('mod_sq', only=['sqr']),
        Print(['"Square of 42 = "', sq_call])
    ])
    if not has_fortran():
        skip("No fortran compiler found.")
    (stdout, stderr), info = compile_run_strings([
        ('mod_sq.f90', fcode(mod_sq, standard=90)),
        ('main.f90', fcode(prg_sq, standard=90))
    ])
    assert '42' in stdout
    assert str(42**2) in stdout
    assert stderr == ''


def test_Subroutine():
    # Code to generate the subroutine in the example from
    # http://www.fortran90.org/src/best-practices.html#arrays
    r = Symbol('r', real=True)
    i = Symbol('i', integer=True)
    v_r = Variable.deduced(r, attrs=(dimension(assumed_extent), intent_out))
    v_i = Variable.deduced(i)
    v_n = Variable('n', integer, value=size(r))
    do_loop = Do([
        Assignment(Element(r, [i]), literal_dp(1)/i**2)
    ], i, 1, v_n)
    sub = Subroutine("f", [r], [
        Declaration(v_n),
        Declaration(v_i),
        do_loop
    ])
    v_r3 = Variable.deduced(r, attrs=[dimension(3)])
    mod = Module('mymod', definitions=[sub])
    prog = Program('foo', [
        use(mod, only=[sub]),
        Declaration(v_r3),
        SubroutineCall(sub, [v_r3]),
        Print([sum_(v_r3), v_r3])
    ])
    (stdout, stderr), info = compile_run_strings([
        ('a.f90', fcode(mod, standard=90)),
        ('b.f90', fcode(prog, standard=90))
    ])
    ref = [1.0/i for i in range(1, 4)]
    assert str(sum(ref)) in stdout
    for _ in ref:
        assert str(_) in stdout
    assert stderr == ''


def test_isign():
    x = Symbol('x', integer=True)
    assert isign(1, x) == isign(1, x)
    assert fcode(isign(1, x), standard=95, source_format='free') == 'isign(1, x)'


def test_dsign():
    x = Symbol('x')
    assert dsign(1, x) == dsign(1, x)
    assert fcode(dsign(literal_dp(1), x), standard=95, source_format='free') == 'dsign(1d0, x)'


def test_cmplx():
    x = Symbol('x')
    assert cmplx(1, x) == cmplx(1, x)


def test_kind():
    x = Symbol('x')
    assert kind(x) == kind(x)


def test_literal_dp():
    assert fcode(literal_dp(0), source_format='free') == '0d0'
