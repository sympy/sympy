#!/usr/bin/python2.7
#printer.py

import os
import sys

preamble = \
"""
\\documentclass{article}
\\pagestyle{empty}
\\usepackage[english]{babel}
\\usepackage[T1]{fontenc}
\\usepackage{listings}
\\usepackage{lmodern}
\\renewcommand{\\familydefault}{\\sfdefault}
\\usepackage{blindtext}
\\usepackage[pdftex]{geometry,graphicx,color}
\\usepackage{bm}
\\usepackage{amsmath}
\\usepackage{mathtools}
\\usepackage{tensor}
\\usepackage{bbding}
\\usepackage{tensor}
\\usepackage{amsfonts}
\\usepackage{srcltx}
\\usepackage{amssymb}
\\usepackage{pdfpages}
\\usepackage{rotating}
\\usepackage{setspace}
\\usepackage{sectsty}
\\usepackage{eufrak}
\\usepackage{makecell}
\\usepackage{longtable}
\\newcommand{\\tbi}[1]{\\textbf{\\textit{#1}}}
\\newcommand{\\tb}[1]{\\textbf{#1}}
\\newcommand{\\ti}[1]{\\textit{#1}}
\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}
\\newcommand{\\lp}{\\left (}
\\newcommand{\\rp}{\\right )}
\\newcommand{\\half}{\\frac{1}{2}}
\\newcommand{\\llt}{\\left <}
\\newcommand{\\rgt}{\\right >}
\\newcommand{\\abs}[1]{\\left |{#1}\\right |}
\\newcommand{\\pdiff}[2]{\\bfrac{\\partial {#1}}{\\partial {#2}}}
\\newcommand{\\pdifftwo}[3]{\\bfrac{\\partial^{2} {#1}}{\\partial {#2}\\partial {#3}}}
\\newcommand{\\lbrc}{\\left \\{}
\\newcommand{\\rbrc}{\\right \\}}
\\newcommand{\\set}[1]{\\lbrc {#1} \\rbrc}
\\newcommand{\\W}{\\wedge}
\\newcommand{\\R}{\\dagger}
\\newcommand{\\lbrk}{\\left [}
\\newcommand{\\rbrk}{\\right ]}
\\newcommand{\\com}[1]{\\lbrk {#1} \\rbrk}
\\newcommand{\\proj}[2]{\\llt {#1} \\rgt_{#2}}
%\\newcommand{\\bm}{\\boldsymbol}
\\newcommand{\\braces}[1]{\\left \\{ {#1} \\right \\}}
\\newcommand{\\grade}[1]{\\left < {#1} \\right >}
\\newcommand{\\f}[2]{{#1}\\lp {#2} \\rp}
\\newcommand{\\paren}[1]{\\lp {#1} \\rp}
\\newcommand{\\eval}[2]{\\left . {#1} \\right |_{#2}}
\\newcommand{\\prm}[1]{{#1}'}
\\newcommand{\\ddt}[1]{\\bfrac{d{#1}}{dt}}
\\newcommand{\\deriv}[3]{\\bfrac{d^{#3}#1}{d{#2}^{#3}}}
\\newcommand{\\be}{\\begin{equation}}
\\newcommand{\\ee}{\\end{equation}}
\\newcommand{\\eb}{\\bm{e}}
\\newcommand{\\ehb}{\\bm{\\hat{e}}}
\\newcommand{\\Tn}[2]{\\f{\\mathcal{T}_{#2}}{#1}}
\\newcommand{\\tr}{\\mbox{tr}}
\\newcommand{\\T}[1]{\\texttt{#1}}
\\newcommand{\\grd}{\\bm{\\nabla}}
\\newcommand{\\ubh}{\\bm{\\hat{u}}}
\\newcommand{\\ebh}{\\bm{\\hat{e}}}
\\newcommand{\\ebf}{\\bm{e}}
\\newcommand{\\mat}[1]{\\left [ {#1} \\right ]}
\\begin{document}
"""

postscript = '\\end{document}\n'

paper_size = '\\documentclass{article}'


texfile = open(sys.argv[1],'r')

latex_str = texfile.read()

latex_str = preamble + latex_str + postscript

print 'latex file =', latex_str

latex_file = open('tmp_latex.tex', 'w')
latex_file.write(latex_str)
latex_file.close()
os.system('pdflatex tmp_latex.tex')
os.system('pdfcrop tmp_latex.pdf')
#os.system('evince tmp_latex-crop.pdf')
os.system('rm tmp_latex.*')
os.system('Pdf2Png tmp_latex-crop')
os.system('mv tmp_latex-crop.png '+sys.argv[1][:-4]+'.png')
os.system('rm tmp_latex-crop.pdf')
os.system('eog '+sys.argv[1][:-4]+'.png')
