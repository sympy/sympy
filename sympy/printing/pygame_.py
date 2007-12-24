"""
install TeX and these Debian packages: python-pygame, python-pexpect, dvipng

To view the equation in the evince:

>>> from sympy import *
>>> import sympy.printing as printing
>>> x = Symbol('x')

>>> printing.view(1/log(x))
>>>

You can use any other viewer:

>>> printing.view(1/log(x), psviewer="kpdf")
>>>

Finally, you can view the equation in the pygame window:

>>> printing.print_pygame(1/log(x))
>>>
"""

import tempfile
from sympy.printing import latex

def print_pygame(st):
    try:
        import pygame
    except ImportError:
        print "Pygame is not installed. In Debian, install the " \
            "python-pygame package."
        return

    from pygame import QUIT, KEYDOWN, K_ESCAPE, K_q

    st = latex(st)
    pygame.font.init()
    size = 640, 240
    screen = pygame.display.set_mode(size)
    screen.fill((255, 255, 255))
    font = pygame.font.Font(None, 24)
    text = font.render(st, True, (0, 0, 0))
    textpos = text.get_rect(centerx=screen.get_width()/2)
    screen.blit(text, textpos)
    pygame.display.flip()

    image = tex2png(st,pygame)
    imagepos = image.get_rect(centerx=screen.get_width()/2).move((0,30))
    screen.blit(image, imagepos)
    pygame.display.flip()

    while 1:
        for event in pygame.event.get():
            if event.type == QUIT:
                return
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                return
            elif event.type == KEYDOWN and event.key == K_q:
                return

tex_str = r"""\documentclass{article}
\begin{document}
%s
\vfill
\end{document}"""

def tex2png(eq, pygame):
    """
    Accepts a latex equation in "eq" and returns an image with this equation.
    """
    #http://www.fauskes.net/nb/htmleqII/
    import os
    import pexpect

    x = tempfile.mktemp()
    tmp1 = '%s.tex'%x

    # create a LaTeX document and insert equations
    f = open(tmp1,'w')
    f.write(tex_str % eq)
    f.close()

    # compile LaTeX document. A DVI file is created
    cwd = os.getcwd()
    os.chdir(tempfile.gettempdir())
    pexpect.run('latex %s' % tmp1)

    # Run dvipng on the generated DVI file. Use tight bounding box.
    # Magnification is set to 1200
    # currently, the dvipng is broken on debian.....
    cmd = "dvipng -T tight -x 1728 -z 9 -bg transparent " \
    + "-o %s.png %s.dvi" % (x,x)
    pexpect.run(cmd)
    image = pygame.image.load("%s.png" % x)

    #remove temporary files
    os.remove("%s.tex" % x)
    os.remove("%s.dvi" % x)
    os.remove("%s.log" % x)
    os.remove("%s.png" % x)
    os.chdir(cwd)

    return image

def preview(expr, output='ps', viewer='gv', euler=True):
    """Preview expression in DVI, PostScript or PDF form.

       This will generate .dvi, .ps or .pdf file respetively in system
       temporary directory using available latex distribution. It will
       make the document using Euler fonts (they were used to typeset
       the well known "Concrete Mathematics" book).

       You can also select the tool, which will be used to preview
       the generated file. By default 'gv' is the tool of choice
       as it is the most lighweight software for this purpose.

    """

    if not euler:
        format = r"""\documentclass[12pt]{article}
                     \begin{document}
                     \begin{equation*}
                     %s
                     \end{equation*}
                     \vfill
                     \end{document}
                 """
    else:
        format = r"""\documentclass[12pt]{article}
                     \usepackage{amsmath}
                     \usepackage{eulervm}
                     \begin{document}
                     \begin{equation*}
                     %s
                     \end{equation*}
                     \vfill
                     \end{document}
                 """

    import os
    import pexpect

    tmp = tempfile.mktemp()

    tex = open(tmp + ".tex", "w")
    tex.write(format % latex(expr)[1:-1])
    tex.close()

    cwd = os.getcwd()
    os.chdir(tempfile.gettempdir())

    pexpect.run("latex -halt-on-error %s.tex" % tmp)

    os.remove(tmp + ".tex")
    os.remove(tmp + ".aux")
    os.remove(tmp + ".log")

    if output != "dvi":
        command = {
            "ps"  : "dvips -o %s.ps %s.dvi",
            "pdf" : "dvipdf %s.dvi %s.pdf",
        }

        try:
            pexpect.run(command[output] % (tmp, tmp))
            os.remove(tmp + ".dvi")
        except KeyError:
            raise "Invalid output format: %s" % output

    os.system("%s %s.%s &" % (viewer, tmp, output))

    os.chdir(cwd)
