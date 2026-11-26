# -*- coding: utf-8 -*-
from __future__ import annotations

from io import BytesIO
from pathlib import Path
from subprocess import CalledProcessError
import tempfile
from unittest.mock import patch

from sympy.core.relational import Eq
from sympy.core.symbol import Symbol
from sympy.functions.elementary.piecewise import Piecewise
from sympy.printing.preview import _get_latex_main, preview
from sympy.testing.pytest import raises


def test_preview():
    x = Symbol('x')
    obj = BytesIO()
    try:
        preview(x, output='png', viewer='BytesIO', outputbuffer=obj)
    except RuntimeError:
        pass  # latex not installed on CI server


def test_preview_unicode_symbol():
    # issue 9107
    a = Symbol('α')
    obj = BytesIO()
    try:
        preview(a, output='png', viewer='BytesIO', outputbuffer=obj)
    except RuntimeError:
        pass  # latex not installed on CI server


def test_preview_latex_construct_in_expr():
    # see PR 9801
    x = Symbol('x')
    pw = Piecewise((1, Eq(x, 0)), (0, True))
    obj = BytesIO()
    try:
        preview(pw, output='png', viewer='BytesIO', outputbuffer=obj)
    except RuntimeError:
        pass  # latex not installed on CI server


def test_get_latex_main_default():
    x = Symbol('x')
    tex = _get_latex_main(x + 1)
    assert r"\documentclass[varwidth,12pt]{standalone}" in tex
    assert r"\usepackage{amsmath}" in tex
    assert r"\usepackage{amsfonts}" in tex
    assert r"\usepackage{euler}" in tex
    assert r"\begin{document}" in tex
    assert r"$\displaystyle x + 1$" in tex
    assert tex.endswith(r"\end{document}")


def test_get_latex_main_context():
    x = Symbol('x')
    tex = _get_latex_main(x + 1, latex_cmd='context')
    assert r"\setupbodyfont[12pt]" in tex
    assert r"\starttext" in tex
    assert r"\startTEXpage" in tex
    assert r"$\displaystyle x + 1$" in tex
    assert r"\stopTEXpage" in tex
    assert tex.endswith(r"\stoptext")
    assert r"\documentclass" not in tex
    assert r"\usepackage" not in tex
    assert r"\begin{document}" not in tex
    assert r"\end{document}" not in tex


def test_get_latex_main_context_options():
    x = Symbol('x')
    tex = _get_latex_main(
        x,
        packages=('tikz', 'units'),
        extra_preamble='\n\\setupcolors[state=start]',
        fontsize=10,
        latex_cmd='/usr/local/bin/context',
    )
    assert r"\setupbodyfont[10pt]" in tex
    assert r"\usemodule[tikz]" in tex
    assert r"\usemodule[units]" in tex
    assert r"\setupcolors[state=start]" in tex
    assert r"$\displaystyle x$" in tex

    custom_tex = _get_latex_main(
        r"$\displaystyle x$",
        preamble="\\setupbodyfont[10pt]\n\\starttext",
        latex_cmd='context',
    )
    assert custom_tex.startswith("\\setupbodyfont[10pt]\n\\starttext")
    assert r"\stopTEXpage" not in custom_tex
    assert custom_tex.endswith(r"\stoptext")

    raises(
        ValueError,
        lambda: _get_latex_main(
            x,
            preamble=r"\starttext",
            packages=('tikz',),
            latex_cmd='context',
        ),
    )


def test_preview_context_pdf_mock():
    x = Symbol('x')
    obj = BytesIO()
    calls = []

    def fake_which(cmd):
        return "/usr/bin/context" if cmd == "context" else None

    def fake_check_output(cmd, cwd=None, **kwargs):
        calls.append(cmd)
        Path(cwd, "texput.pdf").write_bytes(b"%PDF-1.4\n")
        return b""

    with tempfile.TemporaryDirectory() as tmpdir:
        tex_out = Path(tmpdir) / "out.tex"
        with (
            patch("sympy.printing.preview.shutil.which", side_effect=fake_which),
            patch(
                "sympy.printing.preview._check_output_no_window",
                side_effect=fake_check_output,
            ),
        ):
            preview(
                x + 1,
                output='pdf',
                viewer='BytesIO',
                outputbuffer=obj,
                outputTexFile=str(tex_out),
                latex_cmd='context',
            )
        tex_content = tex_out.read_text(encoding='utf-8')

    assert calls == [["context", "--nonstopmode", "texput.tex"]]
    assert obj.getvalue() == b"%PDF-1.4\n"
    assert r"\starttext" in tex_content
    assert r"\stoptext" in tex_content


def test_preview_context_png_mock():
    x = Symbol('x')
    obj = BytesIO()
    calls = []

    def fake_which(cmd):
        if cmd == "context":
            return "/usr/bin/context"
        if cmd == "pdftoppm":
            return "/usr/bin/pdftoppm"
        return None

    def fake_check_output(cmd, cwd=None, **kwargs):
        calls.append(cmd)
        if cmd[0] == "context":
            Path(cwd, "texput.pdf").write_bytes(b"%PDF-1.4\n")
        else:
            Path(cwd, "texput.png").write_bytes(b"\x89PNG\r\n")
        return b""

    with (
        patch("sympy.printing.preview.shutil.which", side_effect=fake_which),
        patch(
            "sympy.printing.preview._check_output_no_window",
            side_effect=fake_check_output,
        ),
    ):
        preview(
            x,
            output='png',
            viewer='BytesIO',
            outputbuffer=obj,
            latex_cmd='context',
        )

    assert calls == [
        ["context", "--nonstopmode", "texput.tex"],
        ["/usr/bin/pdftoppm", "-r", "150", "-png", "-singlefile", "texput.pdf", "texput"],
    ]
    assert obj.getvalue() == b"\x89PNG\r\n"


def test_preview_context_errors_mock():
    x = Symbol('x')
    obj = BytesIO()

    with patch("sympy.printing.preview.shutil.which", return_value=None):
        raises(
            RuntimeError,
            lambda: preview(
                x, output='pdf', viewer='BytesIO', outputbuffer=obj,
                latex_cmd='context',
            ),
        )

    def fake_which(cmd):
        return "/usr/bin/context" if cmd == "context" else None

    with (
        patch("sympy.printing.preview.shutil.which", side_effect=fake_which),
        patch("sympy.printing.preview._check_output_no_window", return_value=b""),
    ):
        raises(
            ValueError,
            lambda: preview(
                x, output='dvi', viewer='BytesIO', outputbuffer=obj,
                latex_cmd='context',
            ),
        )

    with (
        patch("sympy.printing.preview.shutil.which", side_effect=fake_which),
        patch(
            "sympy.printing.preview._check_output_no_window",
            side_effect=CalledProcessError(1, "context", output=b"TeX error"),
        ),
    ):
        raises(
            RuntimeError,
            lambda: preview(
                x, output='pdf', viewer='BytesIO', outputbuffer=obj,
                latex_cmd='context',
            ),
        )


def test_preview_pdflatex_mock():
    x = Symbol('x')
    obj = BytesIO()
    calls = []

    def fake_which(cmd):
        return "/usr/bin/pdflatex" if cmd == "pdflatex" else None

    def fake_check_output(cmd, cwd=None, **kwargs):
        calls.append(cmd)
        Path(cwd, "texput.pdf").write_bytes(b"%PDF-1.4\n")
        return b""

    with tempfile.TemporaryDirectory() as tmpdir:
        tex_out = Path(tmpdir) / "out.tex"
        with (
            patch("sympy.printing.preview.shutil.which", side_effect=fake_which),
            patch(
                "sympy.printing.preview._check_output_no_window",
                side_effect=fake_check_output,
            ),
        ):
            preview(
                x,
                output='pdf',
                viewer='BytesIO',
                outputbuffer=obj,
                outputTexFile=str(tex_out),
                latex_cmd='pdflatex',
            )
        tex_content = tex_out.read_text(encoding='utf-8')

    assert calls == [
        ["pdflatex", "-halt-on-error", "-interaction=nonstopmode", "texput.tex"]
    ]
    assert obj.getvalue() == b"%PDF-1.4\n"
    assert r"\begin{document}" in tex_content
    assert r"\end{document}" in tex_content
