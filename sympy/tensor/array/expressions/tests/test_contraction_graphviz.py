from textwrap import dedent

from sympy import tensorcontraction, tensorproduct, MatrixSymbol, Trace
from sympy.tensor.array.expressions import ArraySymbol, ArrayDiagonal, ArrayContraction, ArrayTensorProduct, \
    PermuteDims, convert_matrix_to_array
from sympy.tensor.array.expressions.contraction_graphviz import ArrayDotGraphPrinter


def _get_graphviz_dot_content(expr):
    adgp = ArrayDotGraphPrinter()
    output = adgp.doprint(expr)
    return output


def test_output_contraction_graphviz():
    A = ArraySymbol('A', (3, 4, 5))
    B = ArraySymbol('B', (5, 4))

    expr: ArrayContraction = tensorcontraction(tensorproduct(A, B), (1, 4), (2, 3))
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            B1 [label="{ <p0> [5] | <p1> [4] }", xlabel="B"];
            A2 [label="{ <p0> 0 [3] | <p1> [4] | <p2> [5] }", xlabel="A"];

            B1:p0 -- A2:p2;
            B1:p1 -- A2:p1;
        }
    """).strip()

    expr: ArrayContraction = ArrayContraction(ArrayTensorProduct(A, B), (1, 4), (2, 3))
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 0 [3] | <p1> [4] | <p2> [5] }", xlabel="A"];
            B2 [label="{ <p0> [5] | <p1> [4] }", xlabel="B"];

            A1:p1 -- B2:p1;
            A1:p2 -- B2:p0;
        }
    """).strip()

    A = ArraySymbol("A", (3, 3))

    expr = ArrayDiagonal(A, (0, 1))
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 0 [3] | <p1> 0 [3] }", xlabel="A"];

        }
    """).strip()

    expr = ArrayDiagonal(ArrayTensorProduct(A, A), (0, 1))
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 2 [3] | <p1> 2 [3] }", xlabel="A"];
            A2 [label="{ <p0> 2 [3] | <p1> 3 [3] }", xlabel="A"];

        }
    """).strip()

    expr = ArrayDiagonal(ArrayTensorProduct(A, A), (0, 1), (2, 3))
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 0 [3] | <p1> 0 [3] }", xlabel="A"];
            A2 [label="{ <p0> 1 [3] | <p1> 1 [3] }", xlabel="A"];

        }
    """).strip()

    expr = ArrayDiagonal(ArrayTensorProduct(A, A), (0, 2), (1, 3))
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 0 [3] | <p1> 1 [3] }", xlabel="A"];
            A2 [label="{ <p0> 0 [3] | <p1> 1 [3] }", xlabel="A"];

        }
    """).strip()

    expr = ArrayTensorProduct(A, A)
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 0 [3] | <p1> 1 [3] }", xlabel="A"];
            A2 [label="{ <p0> 2 [3] | <p1> 3 [3] }", xlabel="A"];

        }
    """).strip()

    expr = ArrayContraction(ArrayTensorProduct(A, A), (1, 2))
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 0 [3] | <p1> [3] }", xlabel="A"];
            A2 [label="{ <p0> [3] | <p1> 1 [3] }", xlabel="A"];

            A1:p1 -- A2:p0;
        }
    """).strip()

    expr = ArrayContraction(ArrayTensorProduct(A, A), (1, 2), (0, 3))
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> [3] | <p1> [3] }", xlabel="A"];
            A2 [label="{ <p0> [3] | <p1> [3] }", xlabel="A"];

            A1:p0 -- A2:p1;
            A1:p1 -- A2:p0;
        }
    """).strip()

    expr = PermuteDims(A, [1, 0])
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 1 [3] | <p1> 0 [3] }", xlabel="A"];

        }
    """).strip()

    A = ArraySymbol("A", (1, 2, 3))
    expr = PermuteDims(A, index_order_old="ijk", index_order_new="kij")
    assert _get_graphviz_dot_content(expr) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 1 [1] | <p1> 2 [2] | <p2> 0 [3] }", xlabel="A"];

        }
    """).strip()

    M = MatrixSymbol("M", 3, 3)
    N = MatrixSymbol("N", 3, 3)
    P = MatrixSymbol("P", 3, 3)

    A = ArraySymbol("A", (3, 3, 3))

    assert _get_graphviz_dot_content(ArrayDiagonal(A, (0, 1, 2))) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            A1 [label="{ <p0> 0 [3] | <p1> 0 [3] | <p2> 0 [3] }", xlabel="A"];

        }
    """).strip()

    assert _get_graphviz_dot_content(convert_matrix_to_array(M * N)) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            M1 [label="{ <p0> 0 [3] | <p1> [3] }", xlabel="M"];
            N2 [label="{ <p0> [3] | <p1> 1 [3] }", xlabel="N"];

            M1:p1 -- N2:p0;
        }""").strip()
    assert _get_graphviz_dot_content(convert_matrix_to_array(M * N * P)) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            N1 [label="{ <p0> [3] | <p1> [3] }", xlabel="N"];
            M2 [label="{ <p0> 0 [3] | <p1> [3] }", xlabel="M"];
            P3 [label="{ <p0> [3] | <p1> 1 [3] }", xlabel="P"];

            N1:p0 -- M2:p1;
            N1:p1 -- P3:p0;
        }""").strip()
    assert _get_graphviz_dot_content(convert_matrix_to_array(M * N.T * P)) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            N1 [label="{ <p0> [3] | <p1> [3] }", xlabel="N"];
            M2 [label="{ <p0> 0 [3] | <p1> [3] }", xlabel="M"];
            P3 [label="{ <p0> [3] | <p1> 1 [3] }", xlabel="P"];

            N1:p0 -- P3:p0;
            N1:p1 -- M2:p1;
        }""").strip()
    assert _get_graphviz_dot_content(convert_matrix_to_array(Trace(M.T * N.T * P))) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            M1 [label="{ <p0> [3] | <p1> [3] }", xlabel="M"];
            N2 [label="{ <p0> [3] | <p1> [3] }", xlabel="N"];
            P3 [label="{ <p0> [3] | <p1> [3] }", xlabel="P"];

            M1:p0 -- N2:p1;
            M1:p1 -- P3:p1;
            N2:p0 -- P3:p0;
        }""").strip()
    assert _get_graphviz_dot_content(convert_matrix_to_array(Trace(M.T * N.T * P))) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            M1 [label="{ <p0> [3] | <p1> [3] }", xlabel="M"];
            N2 [label="{ <p0> [3] | <p1> [3] }", xlabel="N"];
            P3 [label="{ <p0> [3] | <p1> [3] }", xlabel="P"];

            M1:p0 -- N2:p1;
            M1:p1 -- P3:p1;
            N2:p0 -- P3:p0;
        }""").strip()
    assert _get_graphviz_dot_content(convert_matrix_to_array(Trace(M))) == dedent("""
        graph ArrayExpr {
            rankdir = LR;
            splines = true;

            node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];

            M1 [label="{ <p0> [3] | <p1> [3] }", xlabel="M"];

            M1:p0 -- M1:p1;
        }""").strip()
