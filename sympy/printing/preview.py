
import os
import time
import tempfile

from latex import latex

def preview(expr, output='ps', viewer=None, euler=True):
    """View expression in PNG, DVI, PostScript or PDF form.

       This will generate LaTeX representation of the given expression
       and compile it using available TeX distribiution. Then it will
       run appropriate viewer for the given output format or use the
       user defined one. If you prefer not to use external viewer
       then you can use combination of 'png' output and 'pyglet'
       viewer. By default PostScript output is generated.

       By default pretty Euler fonts are used for typesetting (they
       were used to typeset the well known "Concrete Mathematics"
       book). If you prefer default AMS fonts or your system lacks
       'eulervm' LaTeX package then unset 'euler' keyword argument.

       To use viewer auto-detection, lets say for 'png' output, issue::

           >> from sympy import *
           >> x, y = symbols("xy")

           >> preview(x + y, output='png')

       This will choose 'pyglet by default. To select different one::

           >> preview(x + y, output='png', viewer='gimp')

       The 'png' format is considered special. For all other formats
       the rules are sligtly different. As an example we will take
       'dvi' output format. If you would run::

           >> preview(x + y, output='dvi')

       then 'view' will look for available 'dvi' viewers on your
       system (predefined in the function, so it will try evince,
       first, then kdvi and xdvi). If nothing is found you will
       need to set the viewer explicitely::

           >> preview(x + y, output='dvi', viewer='superior-dvi-viwer')

       This will skip auto-dection and will run user specified
       'superior-dvi-viwer'. If 'view' fails to find it on
       your system it will gracefully raise an exception.

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
                "dvi" : [ "evince", "kdvi", "xdvi" ],
                "ps"  : [ "evince", "gsview", "gv" ],
                "pdf" : [ "evince", "kpdf", "acroread", "xpdf", "gv" ],
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

    if not euler:
        format = r"""\documentclass[12pt]{article}
                     \usepackage{amsmath}
                     \begin{document}
                     \pagestyle{empty}
                     %s
                     \vfill
                     \end{document}
                 """
    else:
        format = r"""\documentclass[12pt]{article}
                     \usepackage{amsmath}
                     \usepackage{eulervm}
                     \begin{document}
                     \pagestyle{empty}
                     %s
                     \vfill
                     \end{document}
                 """

    if viewer == "pyglet":
        # import pyglet before we change the current dir, because after that it
        # would fail:
        from sympy.thirdparty import import_thirdparty
        pyglet = import_thirdparty("pyglet")
    tmp = tempfile.mktemp()

    tex = open(tmp + ".tex", "w")
    tex.write(format % latex(expr, inline=False))
    tex.close()

    cwd = os.getcwd()
    os.chdir(tempfile.gettempdir())

    def execute(cmd):
        return pexpect.run(cmd, withexitstatus=True)[1]

    if execute("latex -halt-on-error %s.tex" % tmp) != 0:
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
            if execute(command[output] % (tmp, tmp)) != 0:
                raise SystemError("Failed to generate '%s' output." % output)
            else:
                os.remove(tmp + ".dvi")
        except KeyError:
            raise SystemError("Invalid output format: %s" % output)

    src = "%s.%s" % (tmp, output)

    if viewer == "pyglet":
        from pyglet import window, image, gl
        from pyglet.window import key

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

# some convenience wrappers for 'preview'

view = preview

def pngview(expr, viewer=None, euler=True):
    """Calls preview(expr, output="png", ...)

    This is a convenience function. See the docstring of preview for more
    information.
    """
    return preview(expr, output="png", viewer=viewer, euler=euler)

def pdfview(expr, viewer=None, euler=True):
    """Calls preview(expr, output="pdf", ...)

    This is a convenience function. See the docstring of preview for more
    information.
    """
    return preview(expr, output="pdf", viewer=viewer, euler=euler)

def dviview(expr, viewer=None, euler=True):
    """Calls preview(expr, output="dvi", ...)

    This is a convenience function. See the docstring of preview for more
    information.
    """
    return preview(expr, output="dvi", viewer=viewer, euler=euler)
