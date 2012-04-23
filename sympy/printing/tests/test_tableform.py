
from sympy import TableForm, S
from sympy.printing.latex import latex
from sympy.abc import x
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import sin
from sympy.utilities.pytest import raises

from textwrap import dedent

def test_TableForm():
    s = str(TableForm([["a", "b"], ["c", "d"], ["e", 0]],
        headings="automatic"))
    assert s == """\
  | 1 2\n\
------\n\
1 | a b\n\
2 | c d\n\
3 | e  \n\
"""
    s = str(TableForm([["a", "b"], ["c", "d"], ["e", 0]],
        headings="automatic", wipe_zeros=False))
    assert s == """\
  | 1 2\n\
------\n\
1 | a b\n\
2 | c d\n\
3 | e 0\n\
"""
    s = str(TableForm([[x**2, "b"], ["c", x**2], ["e", "f"]],
            headings=("automatic", None)))
    assert s == """\
1 | x**2 b   \n\
2 | c    x**2\n\
3 | e    f   \n\
"""
    s = str(TableForm([["a", "b"], ["c", "d"], ["e", "f"]],
            headings=(None, "automatic")))
    assert s == """\
1 2\n\
--\n\
a b\n\
c d\n\
e f\n\
"""
    s = str(TableForm([[5, 7], [4, 2], [10, 3]],
            headings=[["Group A", "Group B", "Group C"], ["y1", "y2"]]))
    assert s == """\
        | y1 y2\n\
--------------\n\
Group A | 5  7 \n\
Group B | 4  2 \n\
Group C | 10 3 \n\
"""
    raises(ValueError,
            dedent('''
            TableForm([[5, 7], [4, 2], [10, 3]],
            headings=[["Group A", "Group B", "Group C"], ["y1", "y2"]],
            alignment="middle")''')
            )
    s = str(TableForm([[5, 7], [4, 2], [10, 3]],
            headings=[["Group A", "Group B", "Group C"], ["y1", "y2"]],
            alignment="right"))
    assert s == """\
        | y1 y2\n\
--------------\n\
Group A |  5  7\n\
Group B |  4  2\n\
Group C | 10  3\n\
"""

    s = latex(TableForm([[0, x**3], ["c", S(1)/4], [sqrt(x), sin(x**2)]],
            wipe_zeros=True, headings=("automatic", "automatic")))
    assert s == """\
\\begin{tabular}{l l l}\n\
 & 1 & 2 \\\\\n\
\\hline\n\
1 &   & $x^{3}$ \\\\\n\
2 & $c$ & $\\frac{1}{4}$ \\\\\n\
3 & $\\sqrt{x}$ & $\\sin{\\left (x^{2} \\right )}$ \\\\\n\
\\end{tabular}\
"""
    s = latex(TableForm([["a", x**3], ["c", S(1)/4], [sqrt(x), sin(x**2)]],
            headings=("automatic", "automatic")))
    assert s == """\
\\begin{tabular}{l l l}\n\
 & 1 & 2 \\\\\n\
\\hline\n\
1 & $a$ & $x^{3}$ \\\\\n\
2 & $c$ & $\\frac{1}{4}$ \\\\\n\
3 & $\\sqrt{x}$ & $\\sin{\\left (x^{2} \\right )}$ \\\\\n\
\\end{tabular}\
"""
    s = latex(TableForm([["a", x**3], ["c", S(1)/4], [sqrt(x), sin(x**2)]],
            column_formats=['(%s)', None], headings=("automatic", "automatic")))
    assert s == r"""\begin{tabular}{l l l}
 & 1 & 2 \\
\hline
1 & (a) & $x^{3}$ \\
2 & (c) & $\frac{1}{4}$ \\
3 & (sqrt(x)) & $\sin{\left (x^{2} \right )}$ \\
\end{tabular}"""
    s = latex(TableForm([["a", x**3], ["c", S(1)/4], [sqrt(x), sin(x**2)]]))
    assert s == """\
\\begin{tabular}{l l}\n\
$a$ & $x^{3}$ \\\\\n\
$c$ & $\\frac{1}{4}$ \\\\\n\
$\\sqrt{x}$ & $\\sin{\\left (x^{2} \\right )}$ \\\\\n\
\\end{tabular}\
"""
