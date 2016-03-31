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
import sympy

# If your extensions are in another directory, add it here.
sys.path = ['../sympy', 'ext'] + sys.path

# General configuration
# ---------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.addons.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.viewcode', 'sphinx.ext.mathjax',
              'numpydoc', 'sympylive', 'sphinx.ext.graphviz', ]

# Use this to use pngmath instead
#extensions = ['sphinx.ext.autodoc', 'sphinx.ext.viewcode', 'sphinx.ext.pngmath', ]

# MathJax file, which is free to use.  See http://www.mathjax.org/docs/2.0/start.html
mathjax_path = 'http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML-full'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['.templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General substitutions.
project = 'SymPy'
copyright = '2015 SymPy Development Team'

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
#
# The short X.Y version.
version = sympy.__version__
# The full version, including alpha/beta/rc tags.
release = version

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
html_favicon = '../_build/logo/sympy-notailtext-favicon.ico'
# See http://sphinx-doc.org/theming.html#builtin-themes.


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
html_domain_indices = ['py-modindex']

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
# (source start file, target name, title, author, document class [howto/manual], toctree_only).
# toctree_only is set to True so that the start file document itself is not included in the
# output, only the documents referenced by it via TOC trees.  The extra stuff in the master
# document is intended to show up in the HTML, but doesn't really belong in the LaTeX output.
latex_documents = [('index', 'sympy-%s.tex' % release, 'SymPy Documentation',
                    'SymPy Development Team', 'manual', True)]

# Additional stuff for the LaTeX preamble.
# Tweaked to work with XeTeX.
latex_elements = {
    'babel':     '',
    'fontenc': r'''
\usepackage{bm}
\usepackage{amssymb}
\usepackage{fontspec}
\usepackage[english]{babel}
\defaultfontfeatures{Mapping=tex-text}
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}
''',
    'fontpkg':   '',
    'inputenc':  '',
    'utf8extra': '',
    'preamble':  r'''
% redefine \LaTeX to be usable in math mode
\expandafter\def\expandafter\LaTeX\expandafter{\expandafter\text\expandafter{\LaTeX}}
'''
}

# SymPy logo on title page
html_logo = '_static/sympylogo.png'
latex_logo = '_static/sympylogo_big.png'

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# Show page numbers next to internal references
latex_show_pagerefs = True

# We use False otherwise the module index gets generated twice.
latex_use_modindex = False

default_role = 'math'
pngmath_divpng_args = ['-gamma 1.5', '-D 110']
# Note, this is ignored by the mathjax extension
# Any \newcommand should be defined in the file
pngmath_latex_preamble = '\\usepackage{amsmath}\n' \
    '\\usepackage{bm}\n' \
    '\\usepackage{amsfonts}\n' \
    '\\usepackage{amssymb}\n' \
    '\\setlength{\\parindent}{0pt}\n'

texinfo_documents = [
    (master_doc, 'sympy', 'SymPy Documentation', 'SymPy Development Team',
   'SymPy', 'Computer algebra system (CAS) in Python', 'Programming', 1),
]

# Use svg for graphviz

graphviz_output_format = 'svg'
