from __future__ import print_function, division

from os.path import join
import tempfile
import shutil
from subprocess import STDOUT, CalledProcessError

from sympy.core.compatibility import cStringIO as StringIO
from sympy.core.compatibility import check_output
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.utilities.misc import find_executable
from sympy.utilities.decorator import doctest_depends_on
from .latex import latex


@doctest_depends_on(exe=('latex', 'dvipng'), modules=('pyglet', 'matplotlib'),
            disable_viewers=('evince', 'gimp', 'superior-dvi-viewer'))
def preview(expr, output='png', viewer=None, euler=True, packages=(),
            filename=None, outputbuffer=None, preamble=None, dvioptions=None,
            outputtexfile=None, **latex_settings):
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

    This will choose 'matplotlib' by default. To select a different one, do

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

    There is also support for writing to a StringIO like object, which needs
    to be passed to the 'outputbuffer' argument.

    >>> from StringIO import StringIO
    >>> obj = StringIO()
    >>> preview(x + y, output='png', viewer='StringIO',
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
    passing the desired filename to the 'outputtexfile' keyword
    argument. To write the TeX code to a file named
    "sample.tex" and run the default png viewer to display the resulting
    bitmap, do

    >>> preview(x + y, outputtexfile="sample.tex")

    """

    if viewer == "file" and filename is None:
        SymPyDeprecationWarning(feature="Using viewer=\"file\" without a "
            "specified filename", deprecated_since_version="0.7.3",
            useinstead="viewer=\"file\" and filename=\"desiredname\"",
            issue=3919).warn()
    if viewer == "StringIO" and outputbuffer is None:
        raise ValueError("outputbuffer has to be a StringIO "
                         "compatible object if viewer=\"StringIO\"")
    use_matplotlib = False
    if viewer is None or viewer in ('file', 'StringIO', 'matplotlib'):
        if output == "png":
            try:
                import matplotlib.pyplot as plt
                use_matplotlib = True
                if viewer is None:
                    viewer = 'matplotlib'
            except ImportError:
                if viewer == 'matplotlib':
                    raise ImportError("matplotlib is not installed.\n visit http://www.matplotlib.org/")
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
                    path = find_executable(candidate)
                    if path is not None:
                        viewer = path
                        break
                else:
                    raise SystemError(
                        "No viewers found for '%s' output format." % output)
            except KeyError:
                raise SystemError("Invalid output format: %s" % output)
    elif viewer != 'pyglet'  and not find_executable(viewer):
        raise SystemError("Unrecognized viewer: %s" % viewer)

    if isinstance(expr, str):
        latex_string = expr
    else:
        latex_string = latex(expr, mode='inline', **latex_settings)

    workdir = tempfile.mkdtemp()
    try:
        src = "texput.%s" % (output)
        if use_matplotlib:
            plt.figure(figsize=(1, 1), frameon=False, dpi=50)
            plt.axes(frameon=0)
            plt.text(0.01, 0.8, latex_string, fontsize=50)
            plt.xticks(())
            plt.yticks(())
            if viewer == 'matplotlib':
                plt.show()
                return
            plt.savefig(join(workdir, src), bbox_inches='tight')
            plt.close()
        else:
            _render_with_latex(latex_string, output, workdir, preamble,
                               packages, euler, outputtexfile, dvioptions)
        fullpath = join(workdir, src)
        if viewer == "file":
            if filename is None:
                buffer = StringIO()
                with open(fullpath, 'rb') as fh:
                    buffer.write(fh.read())
                return buffer
            else:
                shutil.move(fullpath, filename)
        elif viewer == "StringIO":
            with open(fullpath, 'rb') as fh:
                outputbuffer.write(fh.read())
        elif viewer == "pyglet":
            try:
                from pyglet import window, image, gl
                from pyglet.window import key
            except ImportError:
                raise ImportError("pyglet is required for preview.\n visit http://www.pyglet.org/")

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
            try:
                check_output([viewer, src], cwd=workdir, stderr=STDOUT)
            except CalledProcessError as e:
                raise RuntimeError(
                    "'%s %s' exited abnormally with the following output:\n%s" %
                    (viewer, src, e.output))
    finally:
        try:
            shutil.rmtree(workdir) # delete directory
        except OSError as e:
            if e.errno != 2: # code 2 - no such file or directory
                raise

def _render_with_latex(latex_string, output, workdir, preamble, packages, euler,
                       outputtexfile, dvioptions):
    if preamble is None:
        actual_packages = packages + ("amsmath", "amsfonts")
        if euler:
            actual_packages += ("euler",)
        package_includes = "\n" + "\n".join(["\\usepackage{%s}" % p
                                             for p in actual_packages])

        preamble = r"""\documentclass[12pt]{article}
\pagestyle{empty}
%s

\begin{document}
""" % (package_includes)
    else:
        if len(packages) > 0:
            raise ValueError("The \"packages\" keyword must not be set if a "
                             "custom LaTeX preamble was specified")
    latex_main = preamble + '\n%s\n\n' + r"\end{document}"

    with open(join(workdir, 'texput.tex'), 'w') as fh:
        fh.write(latex_main % latex_string)

    if outputtexfile is not None:
        shutil.copyfile(join(workdir, 'texput.tex'), outputtexfile)
    if not find_executable('latex'):
        raise RuntimeError("latex program is not installed")

    try:
        check_output(['latex', '-halt-on-error', '-interaction=nonstopmode',
                      'texput.tex'], cwd=workdir, stderr=STDOUT)
    except CalledProcessError as e:
        raise RuntimeError(
            "'latex' exited abnormally with the following output:\n%s" %
            e.output)

    if output != "dvi":
        defaultoptions = {
            "ps": [],
            "pdf": [],
            "png": ["-T", "tight", "-z", "9", "--truecolor"]
        }

        commandend = {
            "ps": ["-o", "texput.ps", "texput.dvi"],
            "pdf": ["texput.dvi", "texput.pdf"],
            "png": ["-o", "texput.png", "texput.dvi"]
        }

        cmd = ["dvi" + output]
        if not find_executable(cmd[0]):
            raise RuntimeError("%s is not installed" % cmd[0])
        try:
            if dvioptions is not None:
                cmd.extend(dvioptions)
            else:
                cmd.extend(defaultoptions[output])
            cmd.extend(commandend[output])
        except KeyError:
            raise SystemError("Invalid output format: %s" % output)

        try:
            check_output(cmd, cwd=workdir, stderr=STDOUT)
        except CalledProcessError as e:
            raise RuntimeError(
                "'%s' exited abnormally with the following output:\n%s" %
                (' '.join(cmd), e.output))
