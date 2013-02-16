from __future__ import with_statement

import os
from os.path import isabs, join
import time
from subprocess import Popen, check_call, PIPE, STDOUT
import tempfile
import shutil
from cStringIO import StringIO

from sympy.utilities.exceptions import SymPyDeprecationWarning
from latex import latex


def preview(expr, output='png', viewer=None, euler=True, packages=(),
            filename=None, outputbuffer=None, formatstr=None, dvioptions=None,
            **latex_settings):
    r"""
    View expression or LaTeX markup in PNG, DVI, PostScript or PDF form.

    If the expr argument is an expression, it will be exported to LaTeX and
    then compiled using the available TeX distribution.  The first argument,
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
    gracefully raise an exception. Currently this depends on pexpect, which
    is not available for windows.

    You may also enter 'file' for the viewer argument. Doing so will cause
    this function to return a file object in read-only mode, if 'filename'
    is unset. However, if it was set, then 'preview' writes the genereted
    file to this filename instead.

    There is also support for writing to a StringIO like object, which needs
    to be passed to the 'outputbuffer' argument.

    >>> from StringIO import StringIO
    >>> obj = StringIO()
    >>> preview(x + y, output='png', viewer='StringIO',
    ...         outputbuffer=obj) # doctest: +SKIP

    The template for the LaTeX code, which is processed by the LaTeX
    interpreter can customized with the 'formatstr' argument. This can be
    used, e.g., to set a differnt font size or use a differnt documentclass.

    >>> fmt =r"\documentclass[10pt]{extarticle}%s\begin{document}%s\end{document}"
    >>> preview(x + y, output='png', formatstr=fmt) # doctest: +SKIP

    If the value of 'output' is different from 'dvi' then command line
    options can be set ('dvioptions' argument) for the execution of the
    'dvi'+output conversion tool. These options have to be in the form of a
    list of strings (see subprocess.Popen).

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
                "dvi": [ "evince", "okular", "kdvi", "xdvi" ],
                "ps": [ "evince", "okular", "gsview", "gv" ],
                "pdf": [ "evince", "okular", "kpdf", "acroread", "xpdf", "gv" ],
            }

            try:
                for candidate in candidates[output]:
                    if pexpect.which(candidate):
                        viewer = candidate
                        break
                else:
                    raise SystemError(
                        "No viewers found for '%s' output format." % output)
            except KeyError:
                raise SystemError("Invalid output format: %s" % output)
    else:
        if viewer == "file":
            if filename is None:
                SymPyDeprecationWarning(feature="Using viewer=\"file\" without a "
                    "specified filename ", last_supported_version="0.7.3",
                    use_instead="viewer=\"file\" and filename=\"desiredname\"")
        elif viewer == "StringIO":
            if outputbuffer is None:
                raise ValueError("outputbuffer has to be a StringIO "
                                 "compatible object if viewer=\"StringIO\"")
        elif viewer not in special and not pexpect.which(viewer):
            raise SystemError("Unrecognized viewer: %s" % viewer)

    actual_packages = packages + ("amsmath", "amsfonts")
    if euler:
        actual_packages += ("euler",)
    package_includes = "\n" + "\n".join(["\\usepackage{%s}" % p
                                         for p in actual_packages])

    if formatstr is None:
        formatstr = r"""\documentclass[12pt]{article}
                        %s
                        \begin{document}
                        \pagestyle{empty}
                        %s
                        \vfill
                        \end{document}
                     """
    formatstr = formatstr % (package_includes, "%s")

    if isinstance(expr, str):
        latex_string = expr
    else:
        latex_string = latex(expr, mode='inline', **latex_settings)

    try:
        workdir = tempfile.mkdtemp()

        p = Popen(['latex', '-halt-on-error'], cwd=workdir, stdin=PIPE,
                  stdout=PIPE)
        _, perr = p.communicate(formatstr % latex_string)
        if p.returncode != 0:
            raise SystemError("Failed to generate DVI output: %s", perr)

        if output != "dvi":
            defaultoptions = {
                "ps": [],
                "pdf": [],
                "png": ["-T", "tight", "-z", "9", "--truecolor"]
            }

            command_end = {
                "ps": ["-o", "texput.ps", "texput.dvi"],
                "pdf": ["texput.dvi", "texput.pdf"],
                "png": ["-o", "texput.png", "texput.dvi"]
            }

            cmd = ["dvi" + output]
            try:
                if dvioptions is not None:
                    cmd.extend(dvioptions)
                else:
                    cmd.extend(defaultoptions[output])
                cmd.extend(command_end[output])
            except KeyError:
                raise SystemError("Invalid output format: %s" % output)

            with open(os.devnull, 'w') as devnull:
                check_call(cmd, cwd=workdir, stdout=devnull, stderr=STDOUT)

        src = "texput.%s" % (output)

        if viewer == "file":
            if filename is None:
                buffer = StringIO()
                with open(join(workdir, src), 'rb') as fh:
                    buffer.write(fh.read())
                return buffer
            else:
                if not isabs(filename):
                    raise ValueError("Provided filename has to be an absolute "
                                     "path")
                shutil.move(join(workdir,src), filename)
        elif viewer == "StringIO":
            with open(join(workdir, src), 'rb') as fh:
                outputbuffer.write(fh.read())
        elif viewer == "pyglet":
            try:
                from pyglet import window, image, gl
                from pyglet.window import key
            except ImportError:
                raise ImportError("pyglet is required for plotting.\n visit http://www.pyglet.org/")

            if output == "png":
                from pyglet.image.codecs.png import PNGImageDecoder
                img = image.load(join(workdir, src), decoder=PNGImageDecoder())
            else:
                raise SystemError("pyglet preview works only for 'png' files.")

            offset = 25

            win = window.Window(
                width=img.width + 2*offset,
                height=img.height + 2*offset,
                caption="sympy",
                resizable=False
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
            with open(os.devnull, 'w') as devnull:
                check_call([viewer, src], cwd=workdir, stdout=devnull,
                           stderr=STDOUT)
    finally:
        try:
            shutil.rmtree(workdir) # delete directory
        except OSError, e:
            if e.errno != 2: # code 2 - no such file or directory
                raise
