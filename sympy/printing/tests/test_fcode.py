from sympy import sin, cos, atan2, gamma, conjugate, sqrt, Factorial, \
    Integral, Piecewise, Add, diff, symbols, S, Real, Dummy
from sympy import Catalan, EulerGamma, E, GoldenRatio, I, pi
from sympy import Function, Rational, Integer, Lambda

from sympy.printing.fcode import fcode, FCodePrinter
from sympy.tensor import IndexedBase, Idx
from sympy.utilities.lambdify import implemented_function
from sympy.utilities.pytest import raises


def test_printmethod():
    x = symbols('x')
    class nint(Function):
        def _fcode(self, printer):
            return "nint(%s)" % printer._print(self.args[0])
    assert fcode(nint(x)) == "      nint(x)"

def test_fcode_Pow():
    x, y = symbols('x,y')
    n = symbols('n', integer=True)

    assert fcode(x**3) == "      x**3"
    assert fcode(x**(y**3)) == "      x**(y**3)"
    assert fcode(1/(sin(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "      (3.5d0*sin(x))**(-x + y**x)/(y + x**2)"
    assert fcode(sqrt(x)) == '      sqrt(x)'
    assert fcode(sqrt(n)) == '      sqrt(dble(n))'
    assert fcode(x**0.5) == '      sqrt(x)'
    assert fcode(x**Rational(1,2)) == '      sqrt(x)'
    assert fcode(sqrt(10)) == '      sqrt(10.0d0)'

def test_fcode_Rational():
    assert fcode(Rational(3,7)) == "      3.0d0/7.0d0"
    assert fcode(Rational(18,9)) == "      2"
    assert fcode(Rational(3,-7)) == "      -3.0d0/7.0d0"
    assert fcode(Rational(-3,-7)) == "      3.0d0/7.0d0"

def test_fcode_Integer():
    assert fcode(Integer(67)) == "      67"
    assert fcode(Integer(-1)) == "      -1"

def test_fcode_Real():
    assert fcode(Real(42.0)) == "      42.0000000000000d0"
    assert fcode(Real(-1e20)) == "      -1.00000000000000d+20"

def test_fcode_functions():
    x, y = symbols('x,y')
    assert fcode(sin(x) ** cos(y)) == "      sin(x)**cos(y)"

def test_fcode_NumberSymbol():
    p = FCodePrinter()
    assert fcode(Catalan) == '      parameter (Catalan = 0.915965594177219d0)\n      Catalan'
    assert fcode(EulerGamma) == '      parameter (EulerGamma = 0.577215664901533d0)\n      EulerGamma'
    assert fcode(E) == '      parameter (E = 2.71828182845905d0)\n      E'
    assert fcode(GoldenRatio) == '      parameter (GoldenRatio = 1.61803398874989d0)\n      GoldenRatio'
    assert fcode(pi) == '      parameter (pi = 3.14159265358979d0)\n      pi'
    assert fcode(pi,precision=5) == '      parameter (pi = 3.1416d0)\n      pi'
    assert fcode(Catalan,human=False) == (set([(Catalan, p._print(Catalan.evalf(15)))]), set([]), '      Catalan')
    assert fcode(EulerGamma,human=False) == (set([(EulerGamma, p._print(EulerGamma.evalf(15)))]), set([]), '      EulerGamma')
    assert fcode(E,human=False) == (set([(E, p._print(E.evalf(15)))]), set([]), '      E')
    assert fcode(GoldenRatio,human=False) == (set([(GoldenRatio, p._print(GoldenRatio.evalf(15)))]), set([]), '      GoldenRatio')
    assert fcode(pi,human=False) == (set([(pi, p._print(pi.evalf(15)))]), set([]), '      pi')
    assert fcode(pi,precision=5,human=False) == (set([(pi, p._print(pi.evalf(5)))]), set([]), '      pi')

def test_fcode_complex():
    assert fcode(I) == "      cmplx(0,1)"
    x = symbols('x')
    assert fcode(4*I) == "      cmplx(0,4)"
    assert fcode(3+4*I) == "      cmplx(3,4)"
    assert fcode(3+4*I+x) == "      cmplx(3,4) + x"
    assert fcode(I*x) == "      cmplx(0,1)*x"
    assert fcode(3+4*I-x) == "      cmplx(3,4) - x"
    x = symbols('x', imaginary=True)
    assert fcode(5*x) == "      5*x"
    assert fcode(I*x) == "      cmplx(0,1)*x"
    assert fcode(3+x) == "      3 + x"

def test_implicit():
    x, y = symbols('x,y')
    assert fcode(sin(x)) == "      sin(x)"
    assert fcode(atan2(x,y)) == "      atan2(x, y)"
    assert fcode(conjugate(x)) == "      conjg(x)"

def test_not_fortran():
    x = symbols('x')
    g = Function('g')
    assert fcode(gamma(x)) == "C     Not Fortran:\nC     gamma(x)\n      gamma(x)"
    assert fcode(Integral(sin(x))) == "C     Not Fortran:\nC     Integral(sin(x), x)\n      Integral(sin(x), x)"
    assert fcode(g(x)) == "C     Not Fortran:\nC     g(x)\n      g(x)"

def test_user_functions():
    x = symbols('x')
    assert fcode(sin(x), user_functions={sin: "zsin"}) == "      zsin(x)"
    x = symbols('x')
    assert fcode(gamma(x), user_functions={gamma: "mygamma"}) == "      mygamma(x)"
    g = Function('g')
    assert fcode(g(x), user_functions={g: "great"}) == "      great(x)"
    n = symbols('n', integer=True)
    assert fcode(Factorial(n), user_functions={Factorial: "fct"}) == "      fct(n)"

def test_inline_function():
    x = symbols('x')
    g = implemented_function('g', Lambda(x, 2*x))
    assert fcode(g(x)) == "      2*x"
    g = implemented_function('g', Lambda(x, 2*pi/x))
    assert fcode(g(x)) == (
            "      parameter (pi = 3.14159265358979d0)\n"
            "      2*pi/x"
            )
    A = IndexedBase('A')
    i = Idx('i', symbols('n', integer=True))
    g = implemented_function('g', Lambda(x, x*(1 + x)*(2 + x)))
    assert fcode(g(A[i]), assign_to=A[i]) == (
            "      do i = 1, n\n"
            "         A(i) = (1 + A(i))*(2 + A(i))*A(i)\n"
            "      end do"
            )

def test_assign_to():
    x = symbols('x')
    assert fcode(sin(x), assign_to="s") == "      s = sin(x)"

def test_line_wrapping():
    x, y = symbols('x,y')
    assert fcode(((x+y)**10).expand(), assign_to="var") == (
        "      var = 45*x**8*y**2 + 120*x**7*y**3 + 210*x**6*y**4 + 252*x**5*y**5\n"
        "     @ + 210*x**4*y**6 + 120*x**3*y**7 + 45*x**2*y**8 + 10*x*y**9 + 10*y\n"
        "     @ *x**9 + x**10 + y**10"
    )
    e = [x**i for i in range(11)]
    assert fcode(Add(*e)) == (
        "      1 + x + x**2 + x**3 + x**4 + x**5 + x**6 + x**7 + x**8 + x**9 + x\n"
        "     @ **10"
    )

def test_fcode_Piecewise():
    x = symbols('x')
    code = fcode(Piecewise((x,x<1),(x**2,True)))
    expected = (
        "      if (x < 1) then\n"
        "         x\n"
        "      else\n"
        "         x**2\n"
        "      end if"
    )
    assert code == expected
    assert fcode(Piecewise((x,x<1),(x**2,True)), assign_to="var") == (
        "      if (x < 1) then\n"
        "         var = x\n"
        "      else\n"
        "         var = x**2\n"
        "      end if"
    )
    a = cos(x)/x
    b = sin(x)/x
    for i in xrange(10):
        a = diff(a, x)
        b = diff(b, x)
    expected = (
        "      if (x < 0) then\n"
        "         weird_name = -cos(x)/x - 1814400*cos(x)/x**9 - 604800*sin(x)/x\n"
        "     @ **8 - 5040*cos(x)/x**5 - 720*sin(x)/x**4 + 10*sin(x)/x**2 + 90*\n"
        "     @ cos(x)/x**3 + 30240*sin(x)/x**6 + 151200*cos(x)/x**7 + 3628800*\n"
        "     @ cos(x)/x**11 + 3628800*sin(x)/x**10\n"
        "      else\n"
        "         weird_name = -sin(x)/x - 3628800*cos(x)/x**10 - 1814400*sin(x)/\n"
        "     @ x**9 - 30240*cos(x)/x**6 - 5040*sin(x)/x**5 - 10*cos(x)/x**2 + 90\n"
        "     @ *sin(x)/x**3 + 720*cos(x)/x**4 + 151200*sin(x)/x**7 + 604800*cos(\n"
        "     @ x)/x**8 + 3628800*sin(x)/x**11\n"
        "      end if"
    )
    code = fcode(Piecewise((a,x<0),(b,True)), assign_to="weird_name")
    assert code == expected
    assert fcode(Piecewise((x,x<1),(x**2,x>1),(sin(x),True))) == (
        "      if (x < 1) then\n"
        "         x\n"
        "      else if (1 < x) then\n"
        "         x**2\n"
        "      else\n"
        "         sin(x)\n"
        "      end if"
    )
    assert fcode(Piecewise((x,x<1),(x**2,x>1),(sin(x),x>0))) == (
        "      if (x < 1) then\n"
        "         x\n"
        "      else if (1 < x) then\n"
        "         x**2\n"
        "      else if (0 < x) then\n"
        "         sin(x)\n"
        "      end if"
    )

def test_wrap_fortran():
    #   "########################################################################"
    printer = FCodePrinter()
    lines = [
        "C     This is a long comment on a single line that must be wrapped properly to produce nice output",
        "      this = is + a + long + and + nasty + fortran + statement + that * must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +  that * must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +   that * must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement + that*must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +   that*must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +    that*must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +     that*must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement + that**must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +  that**must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +   that**must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +    that**must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +     that**must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement(that)/must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran +     statement(that)/must + be + wrapped + properly",
    ]
    wrapped_lines = printer._wrap_fortran(lines)
    expected_lines = [
        "C     This is a long comment on a single line that must be wrapped",
        "C     properly to produce nice output",
        "      this = is + a + long + and + nasty + fortran + statement + that *",
        "     @ must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +  that *",
        "     @ must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +   that",
        "     @ * must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement + that*",
        "     @ must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +   that*",
        "     @ must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +    that",
        "     @ *must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +",
        "     @ that*must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement + that**",
        "     @ must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +  that**",
        "     @ must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +   that",
        "     @ **must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +    that",
        "     @ **must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement +",
        "     @ that**must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran + statement(that)/",
        "     @ must + be + wrapped + properly",
        "      this = is + a + long + and + nasty + fortran +     statement(that)",
        "     @ /must + be + wrapped + properly",
    ]
    for line in wrapped_lines:
        assert len(line) <= 72
    for w, e in zip(wrapped_lines, expected_lines):
        assert w == e
    assert len(wrapped_lines) == len(expected_lines)

def test_wrap_fortran_keep_d0():
    printer = FCodePrinter()
    lines = [
        '      this_variable_is_very_long_because_we_try_to_test_line_break=1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break =1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break  = 1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break   = 1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break    = 1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break = 10.0d0'
        ]
    expected = [
        '      this_variable_is_very_long_because_we_try_to_test_line_break=1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break =',
        '     @ 1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break  =',
        '     @ 1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break   =',
        '     @ 1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break    =',
        '     @ 1.0d0',
        '      this_variable_is_very_long_because_we_try_to_test_line_break =',
        '     @ 10.0d0'
        ]
    assert printer._wrap_fortran(lines) == expected

def test_settings():
    raises(TypeError, 'fcode(S(4), method="garbage")')

def test_free_form_code_line():
    x, y = symbols('x,y')
    assert fcode(cos(x) + sin(y), source_format='free') == "cos(x) + sin(y)"

def test_free_form_continuation_line():
    x, y = symbols('x,y')
    result = fcode(((cos(x) + sin(y))**(7)).expand(), source_format='free')
    expected = (
'7*cos(x)**6*sin(y) + 7*sin(y)**6*cos(x) + 21*cos(x)**5*sin(y)**2 + 35* &\n'
'      cos(x)**4*sin(y)**3 + 35*cos(x)**3*sin(y)**4 + 21*cos(x)**2*sin(y &\n'
'      )**5 + cos(x)**7 + sin(y)**7'
    )
    assert result == expected

def test_free_form_comment_line():
    printer = FCodePrinter({ 'source_format': 'free'})
    lines = [ "! This is a long comment on a single line that must be wrapped properly to produce nice output"]
    expected = [
        '! This is a long comment on a single line that must be wrapped properly',
        '! to produce nice output']
    assert printer._wrap_fortran(lines) == expected

def test_loops():
    n, m = symbols('n,m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)

    expected = (
            'do i = 1, m\n'
            '   y(i) = 0\n'
            'end do\n'
            'do i = 1, m\n'
            '   do j = 1, n\n'
            '      y(i) = %(rhs)s\n'
            '   end do\n'
            'end do'
            )

    code = fcode(A[i, j]*x[j], assign_to=y[i], source_format='free')
    assert (code == expected % {'rhs': 'A(i, j)*x(j) + y(i)'} or
            code == expected % {'rhs': 'x(j)*A(i, j) + y(i)'})

def test_dummy_loops():
    # the following line could also be
    # [Dummy(s, integer=True) for s in 'im']
    # or [Dummy(integer=True) for s in 'im']
    i, m = symbols('i m', integer=True, cls=Dummy)
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx(i, m)

    expected = (
            'do i_%(icount)i = 1, m_%(mcount)i\n'
            '   y(i_%(icount)i) = x(i_%(icount)i)\n'
            'end do'
            ) % {'icount': i.label.dummy_index, 'mcount': m.dummy_index}
    code = fcode(x[i], assign_to=y[i], source_format='free')
    assert code == expected

def test_derived_classes():
    class MyFancyFCodePrinter(FCodePrinter):
        _default_settings = FCodePrinter._default_settings.copy()
        _default_settings['assign_to'] = "bork"

    printer = MyFancyFCodePrinter()
    x = symbols('x')
    assert printer.doprint(sin(x)) == "      bork = sin(x)"

def test_indent():
    codelines = (
            'subroutine test(a)\n'
            'integer :: a, i, j\n'
            '\n'
            'do\n'
            'do \n'
            'do j = 1, 5\n'
            'if (a>b) then\n'
            'if(b>0) then\n'
            'a = 3\n'
            'donot_indent_me = 2\n'
            'do_not_indent_me_either = 2\n'
            'ifIam_indented_something_went_wrong = 2\n'
            'if_I_am_indented_something_went_wrong = 2\n'
            'end should not be unindented here\n'
            'end if\n'
            'endif\n'
            'end do\n'
            'end do\n'
            'enddo\n'
            'end subroutine\n'
            '\n'
            'subroutine test2(a)\n'
            'integer :: a\n'
            'do\n'
            'a = a + 1\n'
            'end do \n'
            'end subroutine\n'
            )
    expected = (
            'subroutine test(a)\n'
            'integer :: a, i, j\n'
            '\n'
            'do\n'
            '   do \n'
            '      do j = 1, 5\n'
            '         if (a>b) then\n'
            '            if(b>0) then\n'
            '               a = 3\n'
            '               donot_indent_me = 2\n'
            '               do_not_indent_me_either = 2\n'
            '               ifIam_indented_something_went_wrong = 2\n'
            '               if_I_am_indented_something_went_wrong = 2\n'
            '               end should not be unindented here\n'
            '            end if\n'
            '         endif\n'
            '      end do\n'
            '   end do\n'
            'enddo\n'
            'end subroutine\n'
            '\n'
            'subroutine test2(a)\n'
            'integer :: a\n'
            'do\n'
            '   a = a + 1\n'
            'end do \n'
            'end subroutine\n'
            )
    p = FCodePrinter({'source_format':'free'})
    result = p.indent_code(codelines)
    assert result == expected
