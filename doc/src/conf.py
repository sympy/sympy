# -*- coding: utf-8 -*-
#
# SymPy documentation build configuration file, created by
# sphinx-quickstart.py on Sat Mar 22 19:34:32 2008.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default value; values that are commented out
# serve to show the default value.

import sys

# If your extensions are in another directory, add it here.
sys.path.extend(['../sympy', 'ext'])

# General configuration
# ---------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.addons.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.pngmath', 'math_dollar']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['.templates']

# The suffix of source filenames.
source_suffix = '.txt'

# The master toctree document.
master_doc = 'index'

# General substitutions.
project = 'SymPy'
copyright = '2008, 2009, 2010 SymPy Development Team'

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
#
# The short X.Y version.
version = '0.6.7'
# The full version, including alpha/beta/rc tags.
release = '0.6.7'

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# List of documents that shouldn't be included in the build.
#unused_docs = []

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# Options for HTML output
# -----------------------

# The style sheet to use for HTML and HTML Help pages. A file of that name
# must exist either in Sphinx' static/ path, or in one of the custom paths
# given in html_static_path.
html_style = 'default.css'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

html_logo = '_static/sympylogo.png'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Content template for the index page.
#html_index = ''

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_use_modindex = True

# If true, the reST sources are included in the HTML build as _sources/<name>.
#html_copy_source = True

# Output file base name for HTML help builder.
htmlhelp_basename = 'SymPydoc'


# Options for LaTeX output
# ------------------------

# The paper size ('letter' or 'a4').
#latex_paper_size = 'letter'

# The font size ('10pt', '11pt' or '12pt').
#latex_font_size = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [('index', 'sympy.tex', 'SymPy Documentation',
                        'SymPy Development Team', 'manual')]

# Additional stuff for the LaTeX preamble.
#latex_preamble = ''

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_use_modindex = True

default_role = 'math'
pngmath_divpng_args = ['-gamma 1.5','-D 110']
pngmath_latex_preamble =  '\\usepackage{amsmath}\n'+\
              '\\usepackage{bm}\n'+\
              '\\usepackage{amsfonts}\n'+\
              '\\usepackage{amssymb}\n'+\
              '\\setlength{\\parindent}{0pt}\n'+\
              '\\newcommand{\\bfrac}[2]{\\displaystyle\\frac{#1}{#2}}\n'+\
              '\\newcommand{\\lp}{\\left (}\n'+\
              '\\newcommand{\\rp}{\\right )}\n'+\
              '\\newcommand{\\half}{\\frac{1}{2}}\n'+\
              '\\newcommand{\\llt}{\\left <}\n'+\
              '\\newcommand{\\rgt}{\\right >}\n'+\
              '\\newcommand{\\abs}[1]{\\left |{#1}\\right | }\n'+\
              '\\newcommand{\\pdiff}[2]{\\bfrac{\\partial {#1}}{\\partial {#2}}}\n'+\
              '\\newcommand{\\lbrc}{\\left \\{}\n'+\
              '\\newcommand{\\rbrc}{\\right \\}}\n'+\
              '\\newcommand{\\W}{\\wedge}\n'+\
              '\\newcommand{\\R}{\\dagger}\n'+\
              '\\newcommand{\\lbrk}{\\left [}\n'+\
              '\\newcommand{\\rbrk}{\\right ]}\n'+\
              '\\newcommand{\\proj}[2]{\\llt {#1} \\rgt_{#2}}\n'+\
              '\\newcommand{\\bs}{$\\backslash$}\n'+\
              '\\newcommand{\\sinf}[1]{\\sin\\lp{#1}\\rp}\n'+\
              '\\newcommand{\\cosf}[1]{\\cos\\lp{#1}\\rp}\n'+\
              '\\newcommand{\\ebh}{\\hat{\\bm{e}}}\n'

