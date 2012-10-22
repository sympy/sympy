from __future__ import with_statement

import os
import time
import tempfile

from latex import latex

def preview(expr, output='png', viewer=None, euler=True, packages=(), **latex_settings):
    r"""
    View expression or LaTeX markup in PNG, DVI, PostScript or PDF form.

    If the expr argument is an expression, it will be exported to LaTeX and
    then compiled using available the TeX distribution.  The first argument,
    'expr', may also be a LaTeX string.  The function will then run the
    appropriate viewer for the given output format or use the user defined
    one. By default png output is generated.

    By default pretty Euler fonts are used for typesetting (they were used to
    typeset the well known "Concrete Mathematics" book). For that to work, you
    need the 'eulervm.sty' LaTeX style (in Debian/Ubuntu, install the
    texlive-fonts-extra package). If you prefer default AMS fonts or your
    system lacks 'eulervm' LaTeX package then unset the 'euler' keyword
    argument.

    To use viewer auto-detection, lets say for 'png' output, issue

    >>> from sympy import symbols, preview, Symbol
    >>> x, y = symbols("x,y")

    >>> preview(x + y, output='png') # doctest: +SKIP

    This will choose 'pyglet' by default. To select a different one, do

    >>> preview(x + y, output='png', viewer='gimp') # doctest: +SKIP

    The 'png' format is considered special. For all other formats the rules
    are slightly different. As an example we will take 'dvi' output format. If
    you would run

    >>> preview(x + y, output='dvi') # doctest: +SKIP

    then 'view' will look for available 'dvi' viewers on your system
    (predefined in the function, so it will try evince, first, then kdvi and
    xdvi). If nothing is found you will need to set the viewer explicitly.

    >>> preview(x + y, output='dvi', viewer='superior-dvi-viewer') # doctest: +SKIP

    This will skip auto-detection and will run user specified
    'superior-dvi-viewer'. If 'view' fails to find it on your system it will
    gracefully raise an exception. You may also enter 'file' for the viewer
    argument. Doing so will cause this function to return a file object in
    read-only mode.

    Currently this depends on pexpect, which is not available for windows.

    Additional keyword args will be passed to the latex call, e.g., the
    symbol_names flag.

    >>> phidd = Symbol('phidd')
    >>> preview(phidd, symbol_names={phidd:r'\ddot{\varphi}'}) # doctest: +SKIP

    """

    # we don't want to depend on anything not in the
    # standard library with SymPy by default
    import pexpect

    special = [ 'pyglet' ]

    if viewer is None:
        if output == "png":
            viewer = "pyglet"
        else:
            # sorted in order from most pretty to most ugly
            # very discussable, but indeed 'gv' looks awful :)
            candidates = {
                "dvi" : [ "evince", "okular", "kdvi", "xdvi" ],
                "ps"  : [ "evince", "okular", "gsview", "gv" ],
                "pdf" : [ "evince", "okular", "kpdf", "acroread", "xpdf", "gv" ],
            }

            try:
                for candidate in candidates[output]:
                    if pexpect.which(candidate):
                        viewer = candidate
                        break
                else:
                    raise SystemError("No viewers found for '%s' output format." % output)
            except KeyError:
                raise SystemError("Invalid output format: %s" % output)
    else:
        if viewer not in special and not pexpect.which(viewer):
            raise SystemError("Unrecognized viewer: %s" % viewer)

    actual_packages = packages + ("amsmath", "amsfonts")
    if euler:
        actual_packages += ("euler",)
    package_includes = "\n".join(["\\usepackage{%s}" % p
                                  for p in actual_packages])

    format = r"""\documentclass[12pt]{article}
                 %s
                 \begin{document}
                 \pagestyle{empty}
                 %s
                 \vfill
                 \end{document}
              """ % (package_includes, "%s")

    if isinstance(expr, str):
        latex_string = expr
    else:
        latex_string = latex(expr, mode='inline', **latex_settings)


    tmp = tempfile.mktemp()

    with open(tmp + ".tex", "w") as tex:
        tex.write(format % latex_string)

    cwd = os.getcwd()
    os.chdir(tempfile.gettempdir())

    if os.system("latex -halt-on-error %s.tex" % tmp) != 0:
        raise SystemError("Failed to generate DVI output.")

    os.remove(tmp + ".tex")
    os.remove(tmp + ".aux")
    os.remove(tmp + ".log")

    if output != "dvi":
        command = {
            "ps"  : "dvips -o %s.ps %s.dvi",
            "pdf" : "dvipdf %s.dvi %s.pdf",
            "png" : "dvipng -T tight -z 9 " + \
                    "--truecolor -o %s.png %s.dvi",
        }

        try:
            if os.system(command[output] % (tmp, tmp)) != 0:
                raise SystemError("Failed to generate '%s' output." % output)
            else:
                os.remove(tmp + ".dvi")
        except KeyError:
            raise SystemError("Invalid output format: %s" % output)

    src = "%s.%s" % (tmp, output)
    src_file = None

    if viewer == "file":
        src_file = open(src, 'rb')
    elif viewer == "pyglet":
        try:
            from pyglet import window, image, gl
            from pyglet.window import key
        except ImportError:
            raise ImportError("pyglet is required for plotting.\n visit http://www.pyglet.org/")

        if output == "png":
            from pyglet.image.codecs.png import PNGImageDecoder
            img = image.load(src, decoder=PNGImageDecoder())
        else:
            raise SystemError("pyglet preview works only for 'png' files.")

        offset = 25

        win = window.Window(
            width = img.width + 2*offset,
            height = img.height + 2*offset,
            caption = "sympy",
            resizable = False
        )

        win.set_vsync(False)

        try:
            def on_close():
                win.has_exit = True

            win.on_close = on_close

            def on_key_press(symbol, modifiers):
                if symbol in [key.Q, key.ESCAPE]:
                    on_close()

            win.on_key_press = on_key_press

            def on_expose():
                gl.glClearColor(1.0, 1.0, 1.0, 1.0)
                gl.glClear(gl.GL_COLOR_BUFFER_BIT)

                img.blit(
                    (win.width - img.width) / 2,
                    (win.height - img.height) / 2
                )

            win.on_expose = on_expose

            while not win.has_exit:
                win.dispatch_events()
                win.flip()
        except KeyboardInterrupt:
            pass

        win.close()
    else:
        os.system("%s %s &> /dev/null &" % (viewer, src))
        time.sleep(2) # wait for the viewer to read data

    os.remove(src)
    os.chdir(cwd)

    if src_file is not None:
        return src_file
