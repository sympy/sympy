from sympy.parsing.fortran.fortran_parser import src_to_sympy


imp_stmt = "from sympy import Symbol"

def test_declaration():
    src1 = "integer :: a,b"
    src2 = "real :: p,q"
    src3 = """\
    integer :: x
    real :: p
    integer :: y
    real :: q
    """
    res1 = src_to_sympy(src1)
    res2 = src_to_sympy(src2)
    res3 = src_to_sympy(src3)
    cmp1 = (
        imp_stmt + "\n" +
        "a = Symbol('a', integer=True)" + "\n" +
        "b = Symbol('b', integer=True)"
    )
    cmp2 = (
        imp_stmt + "\n" +
        "p = Symbol('p', real=True)" + "\n" +
        "q = Symbol('q', real=True)"
    )
    cmp3 = (
        imp_stmt + "\n" +
        "x = Symbol('x', integer=True)" + "\n" +
        "p = Symbol('p', real=True)" + "\n" +
        "y = Symbol('y', integer=True)" + "\n" +
        "q = Symbol('q', real=True)"
    )
    assert res1.strip() == cmp1
    assert res2.strip() == cmp2
    assert res3.strip() == cmp3


def test_assignment():
    #TODO: Arithmetic Assignments
    src = """\
    integer :: a, b
    a = b
    """
    res = src_to_sympy(src)
    cmp = (
        imp_stmt + "\n" +
        "a = Symbol('a', integer=True)" + "\n" +
        "b = Symbol('b', integer=True)" + "\n" +
        "a = b"
    )
    assert res.strip() == cmp

def test_assignment_add():
    src = """\
    real :: c, d, e
    e = c + d
    """
    res = src_to_sympy(src)
    cmp = (
        imp_stmt + "\n" +
        "c = Symbol('c', real=True)" + "\n" +
        "d = Symbol('d', real=True)" + "\n" +
        "e = Symbol('e', real=True)" + "\n" +
        "e = c + d"
    )
    assert res.strip() == cmp

def test_assignment_sub():
    src = """\
    real :: c, d, e
    e = c - d
    """
    res = src_to_sympy(src)
    cmp = (
        imp_stmt + "\n" +
        "c = Symbol('c', real=True)" + "\n" +
        "d = Symbol('d', real=True)" + "\n" +
        "e = Symbol('e', real=True)" + "\n" +
        "e = c - d"
    )
    assert res.strip() == cmp

def test_assignment_mul():
    src = """\
    real :: c, d, e
    e = c * d
    """
    res = src_to_sympy(src)
    cmp = (
        imp_stmt + "\n" +
        "c = Symbol('c', real=True)" + "\n" +
        "d = Symbol('d', real=True)" + "\n" +
        "e = Symbol('e', real=True)" + "\n" +
        "e = c * d"
    )
    assert res.strip() == cmp

def test_assignment_div():
    src = """\
    real :: c, d, e
    e = c / d
    """
    res = src_to_sympy(src)
    cmp = (
        imp_stmt + "\n" +
        "c = Symbol('c', real=True)" + "\n" +
        "d = Symbol('d', real=True)" + "\n" +
        "e = Symbol('e', real=True)" + "\n" +
        "e = c / d"
    )
    assert res.strip() == cmp

def test_binop():
    src1 = """\
    real :: c, d, e, f, g
    g = c + d - e * f
    """
    src2 = """\
    real :: c, d, e, f, g
    g = c - d + (e * f)
    """
    src3 = """\
    real :: c, d, e, f, g
    g = (c + e * (f / d)) - (f + (e * c))
    """
    res1 = src_to_sympy(src1)
    res2 = src_to_sympy(src2)
    res3 = src_to_sympy(src3)
    cmp1 = (
        imp_stmt + "\n" +
        "c = Symbol('c', real=True)" + "\n" +
        "d = Symbol('d', real=True)" + "\n" +
        "e = Symbol('e', real=True)" + "\n" +
        "f = Symbol('f', real=True)" + "\n" +
        "g = Symbol('g', real=True)" + "\n" +
        "g = c + d - e * f"
    )
    cmp2 = (
        imp_stmt + "\n" +
        "c = Symbol('c', real=True)" + "\n" +
        "d = Symbol('d', real=True)" + "\n" +
        "e = Symbol('e', real=True)" + "\n" +
        "f = Symbol('f', real=True)" + "\n" +
        "g = Symbol('g', real=True)" + "\n" +
        "g = c - d + e * f"
    )
    cmp3 = (
        imp_stmt + "\n" +
        "c = Symbol('c', real=True)" + "\n" +
        "d = Symbol('d', real=True)" + "\n" +
        "e = Symbol('e', real=True)" + "\n" +
        "f = Symbol('f', real=True)" + "\n" +
        "g = Symbol('g', real=True)" + "\n" +
        "g = c + e * (f / d) - (f + e * c)"
    )
    assert res1.strip() == cmp1
    assert res2.strip() == cmp2
    assert res3.strip() == cmp3




def test_function():
    src = """\
    function abc(a ,b) result(r)
    integer, intent(in) :: a,b
    r = a * b
    end function
    """
    res = src_to_sympy(src)
    cmp = (
        imp_stmt + "\n" +
        "a = Symbol('a', integer=True)" + "\n" +
        "b = Symbol('b', integer=True)" + "\n" +
        "r = Symbol('r', integer=True)" + "\n" +
        "\n" + "\n" +
        "def abc(a, b):" + "\n" +
        "    " + "r = a * b" + "\n" +
        "    " + "return r"
    )
    assert res.strip() == cmp
