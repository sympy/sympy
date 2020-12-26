import os
from os.path import join
import shutil
import tempfile

try:
    from subprocess import STDOUT, CalledProcessError, check_output
except ImportError:
    pass

from sympy.utilities.decorator import doctest_depends_on
from .latex import latex

__doctest_requires__ = {('preview',): ['pyglet']}


def _check_output_no_window(*args, **kwargs):
    # Avoid showing a cmd.exe window when running this
    # on Windows
    if os.name == 'nt':
        creation_flag = 0x08000000 # CREATE_NO_WINDOW
    else:
        creation_flag = 0 # Default value
    return check_output(*args, creationflags=creation_flag, **kwargs)


def _run_pyglet(fname, fmt):
    from pyglet import window, image, gl
    from pyglet.window import key

    if fmt == "png":
        from pyglet.image.codecs.png import PNGImageDecoder
        img = image.load(fname, decoder=PNGImageDecoder())
    else:
        raise ValueError("pyglet preview works only for 'png' files.")

    offset = 25

    config = gl.Config(double_buffer=False)
    win = window.Window(
        width=img.width + 2*offset,
        height=img.height + 2*offset,
        caption="sympy",
        resizable=False,
        config=config
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


@doctest_depends_on(exe=('latex', 'dvipng'), modules=('pyglet',),
            disable_viewers=('evince', 'gimp', 'superior-dvi-viewer'))
def preview(expr, output='png', viewer=None, euler=True, packages=(),
            filename=None, outputbuffer=None, preamble=None, dvioptions=None,
            outputTexFile=None, **latex_settings):
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

    >>> preview(x + y, output='png')

    This will choose 'pyglet' by default. To select a different one, do

    >>> preview(x + y, output='png', viewer='gimp')

    The 'png' format is considered special. For all other formats the rules
    are slightly different. As an example we will take 'dvi' output format. If
    you would run

    >>> preview(x + y, output='dvi')

    then 'view' will look for available 'dvi' viewers on your system
    (predefined in the function, so it will try evince, first, then kdvi and
    xdvi). If nothing is found you will need to set the viewer explicitly.

    >>> preview(x + y, output='dvi', viewer='superior-dvi-viewer')

    This will skip auto-detection and will run user specified
    'superior-dvi-viewer'. If 'view' fails to find it on your system it will
    gracefully raise an exception.

    You may also enter 'file' for the viewer argument. Doing so will cause
    this function to return a file object in read-only mode, if 'filename'
    is unset. However, if it was set, then 'preview' writes the genereted
    file to this filename instead.

    There is also support for writing to a BytesIO like object, which needs
    to be passed to the 'outputbuffer' argument.

    >>> from io import BytesIO
    >>> obj = BytesIO()
    >>> preview(x + y, output='png', viewer='BytesIO',
    ...         outputbuffer=obj)

    The LaTeX preamble can be customized by setting the 'preamble' keyword
    argument. This can be used, e.g., to set a different font size, use a
    custom documentclass or import certain set of LaTeX packages.

    >>> preamble = "\\documentclass[10pt]{article}\n" \
    ...            "\\usepackage{amsmath,amsfonts}\\begin{document}"
    >>> preview(x + y, output='png', preamble=preamble)

    If the value of 'output' is different from 'dvi' then command line
    options can be set ('dvioptions' argument) for the execution of the
    'dvi'+output conversion tool. These options have to be in the form of a
    list of strings (see subprocess.Popen).

    Additional keyword args will be passed to the latex call, e.g., the
    symbol_names flag.

    >>> phidd = Symbol('phidd')
    >>> preview(phidd, symbol_names={phidd:r'\ddot{\varphi}'})

    For post-processing the generated TeX File can be written to a file by
    passing the desired filename to the 'outputTexFile' keyword
    argument. To write the TeX code to a file named
    "sample.tex" and run the default png viewer to display the resulting
    bitmap, do

    >>> preview(x + y, outputTexFile="sample.tex")


    """
    special = [ 'pyglet' ]

    if viewer is None:
        if output == "png":
            viewer = "pyglet"
        else:
            # sorted in order from most pretty to most ugly
            # very discussable, but indeed 'gv' looks awful :)
            # TODO add candidates for windows to list
            candidates = {
                "dvi": [ "evince", "okular", "kdvi", "xdvi" ],
                "ps": [ "evince", "okular", "gsview", "gv" ],
                "pdf": [ "evince", "okular", "kpdf", "acroread", "xpdf", "gv" ],
            }

            try:
                candidate_viewers = candidates[output]
            except KeyError:
                raise ValueError("Invalid output format: %s" % output) from None

            for candidate in candidate_viewers:
                path = shutil.which(candidate)
                if path is not None:
                    viewer = path
                    break
            else:
                raise OSError(
                    "No viewers found for '%s' output format." % output)
    else:
        if viewer == "file":
            if filename is None:
                raise ValueError("filename has to be specified if viewer=\"file\"")
        elif viewer == "BytesIO":
            if outputbuffer is None:
                raise ValueError("outputbuffer has to be a BytesIO "
                                 "compatible object if viewer=\"BytesIO\"")
        elif viewer not in special and not shutil.which(viewer):
            raise OSError("Unrecognized viewer: %s" % viewer)


    if preamble is None:
        actual_packages = packages + ("amsmath", "amsfonts")
        if euler:
            actual_packages += ("euler",)
        package_includes = "\n" + "\n".join(["\\usepackage{%s}" % p
                                             for p in actual_packages])

        preamble = r"""\documentclass[varwidth,12pt]{standalone}
%s

\begin{document}
""" % (package_includes)
    else:
        if packages:
            raise ValueError("The \"packages\" keyword must not be set if a "
                             "custom LaTeX preamble was specified")

    if isinstance(expr, str):
        latex_string = expr
    else:
        latex_string = ('$\\displaystyle ' +
                        latex(expr, mode='plain', **latex_settings) +
                        '$')

    latex_main = preamble + '\n' + latex_string + '\n\n' + r"\end{document}"

    with tempfile.TemporaryDirectory() as workdir:
        with open(join(workdir, 'texput.tex'), 'w', encoding='utf-8') as fh:
            fh.write(latex_main)

        if outputTexFile is not None:
            shutil.copyfile(join(workdir, 'texput.tex'), outputTexFile)

        if not shutil.which('latex'):
            raise RuntimeError("latex program is not installed")

        try:
            _check_output_no_window(
                ['latex', '-halt-on-error', '-interaction=nonstopmode',
                 'texput.tex'],
                cwd=workdir,
                stderr=STDOUT)
        except CalledProcessError as e:
            raise RuntimeError(
                "'latex' exited abnormally with the following output:\n%s" %
                e.output)

        src = "texput.%s" % (output)

        if output != "dvi":
            # in order of preference
            commandnames = {
                "ps": ["dvips"],
                "pdf": ["dvipdfmx", "dvipdfm", "dvipdf"],
                "png": ["dvipng"],
                "svg": ["dvisvgm"],
            }
            try:
                cmd_variants = commandnames[output]
            except KeyError:
                raise ValueError("Invalid output format: %s" % output) from None

            # find an appropriate command
            for cmd_variant in cmd_variants:
                cmd_path = shutil.which(cmd_variant)
                if cmd_path:
                    cmd = [cmd_path]
                    break
            else:
                if len(cmd_variants) > 1:
                    raise RuntimeError("None of %s are installed" % ", ".join(cmd_variants))
                else:
                    raise RuntimeError("%s is not installed" % cmd_variants[0])

            defaultoptions = {
                "dvipng": ["-T", "tight", "-z", "9", "--truecolor"],
                "dvisvgm": ["--no-fonts"],
            }

            commandend = {
                "dvips": ["-o", src, "texput.dvi"],
                "dvipdf": ["texput.dvi", src],
                "dvipdfm": ["-o", src, "texput.dvi"],
                "dvipdfmx": ["-o", src, "texput.dvi"],
                "dvipng": ["-o", src, "texput.dvi"],
                "dvisvgm": ["-o", src, "texput.dvi"],
            }

            if dvioptions is not None:
                cmd.extend(dvioptions)
            else:
                cmd.extend(defaultoptions.get(cmd_variant, []))
            cmd.extend(commandend[cmd_variant])

            try:
                _check_output_no_window(cmd, cwd=workdir, stderr=STDOUT)
            except CalledProcessError as e:
                raise RuntimeError(
                    "'%s' exited abnormally with the following output:\n%s" %
                    (' '.join(cmd), e.output))


        if viewer == "file":
            shutil.move(join(workdir, src), filename)
        elif viewer == "BytesIO":
            with open(join(workdir, src), 'rb') as fh:
                outputbuffer.write(fh.read())
        elif viewer == "pyglet":
            try:
                import pyglet  # noqa: F401
            except ImportError:
                raise ImportError("pyglet is required for preview.\n visit http://www.pyglet.org/")

            return _run_pyglet(join(workdir, src), fmt=output)
        else:
            try:
                _check_output_no_window(
                    [viewer, src], cwd=workdir, stderr=STDOUT)
            except CalledProcessError as e:
                raise RuntimeError(
                    "'%s %s' exited abnormally with the following output:\n%s" %
                    (viewer, src, e.output))
