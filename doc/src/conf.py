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
import inspect
import os
import subprocess
from datetime import datetime

import sympy

# If your extensions are in another directory, add it here.
sys.path = ['ext'] + sys.path

# General configuration
# ---------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.addons.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.linkcode', 'sphinx_math_dollar',
              'sphinx.ext.mathjax', 'numpydoc', 'sympylive',
              'sphinx.ext.graphviz', 'matplotlib.sphinxext.plot_directive']

# Use this to use pngmath instead
#extensions = ['sphinx.ext.autodoc', 'sphinx.ext.viewcode', 'sphinx.ext.pngmath', ]

# Enable warnings for all bad cross references. These are turned into errors
# with the -W flag in the Makefile.
nitpicky = True

# To stop docstrings inheritance.
autodoc_inherit_docstrings = False

# MathJax file, which is free to use.  See https://www.mathjax.org/#gettingstarted
# As explained in the link using latest.js will get the latest version even
# though it says 2.7.5.
mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/latest.js?config=TeX-AMS_HTML-full'

# See https://www.sympy.org/sphinx-math-dollar/
mathjax_config = {
    'tex2jax': {
        'inlineMath': [ ["\\(","\\)"] ],
        'displayMath': [["\\[","\\]"] ],
    },
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

suppress_warnings = ['ref.citation', 'ref.footnote']

# General substitutions.
project = 'SymPy'
copyright = '{} SymPy Development Team'.format(datetime.utcnow().year)

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

# Don't show the source code hyperlinks when using matplotlib plot directive.
plot_html_show_source_link = False

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

html_theme = 'classic'

html_logo = '_static/sympylogo.png'
html_favicon = '../_build/logo/sympy-notailtext-favicon.ico'
# See http://www.sphinx-doc.org/en/master/theming.html#builtin-themes


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
% Define version of \LaTeX that is usable in math mode
\let\OldLaTeX\LaTeX
\renewcommand{\LaTeX}{\text{\OldLaTeX}}

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


# Requried for linkcode extension.
# Get commit hash from the external file.
commit_hash_filepath = '../commit_hash.txt'
commit_hash = None
if os.path.isfile(commit_hash_filepath):
    with open(commit_hash_filepath) as f:
        commit_hash = f.readline()

# Get commit hash from the external file.
if not commit_hash:
    try:
        commit_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'])
        commit_hash = commit_hash.decode('ascii')
        commit_hash = commit_hash.rstrip()
    except:
        import warnings
        warnings.warn(
            "Failed to get the git commit hash as the command " \
            "'git rev-parse HEAD' is not working. The commit hash will be " \
            "assumed as the SymPy master, but the lines may be misleading " \
            "or nonexistent as it is not the correct branch the doc is " \
            "built with. Check your installation of 'git' if you want to " \
            "resolve this warning.")
        commit_hash = 'master'

fork = 'sympy'
blobpath = \
    "https://github.com/{}/sympy/blob/{}/sympy/".format(fork, commit_hash)


def linkcode_resolve(domain, info):
    """Determine the URL corresponding to Python object."""
    if domain != 'py':
        return

    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        return

    obj = submod
    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except Exception:
            return

    # strip decorators, which would resolve to the source of the decorator
    # possibly an upstream bug in getsourcefile, bpo-1764286
    try:
        unwrap = inspect.unwrap
    except AttributeError:
        pass
    else:
        obj = unwrap(obj)

    try:
        fn = inspect.getsourcefile(obj)
    except Exception:
        fn = None
    if not fn:
        return

    try:
        source, lineno = inspect.getsourcelines(obj)
    except Exception:
        lineno = None

    if lineno:
        linespec = "#L%d-L%d" % (lineno, lineno + len(source) - 1)
    else:
        linespec = ""

    fn = os.path.relpath(fn, start=os.path.dirname(sympy.__file__))
    return blobpath + fn + linespec
