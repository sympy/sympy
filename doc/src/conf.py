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

# Make sure we import sympy from git
sys.path.insert(0, os.path.abspath('../..'))

import sympy

# If your extensions are in another directory, add it here.
sys.path = ['ext'] + sys.path

# General configuration
# ---------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.addons.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.linkcode',
              'sphinx_math_dollar', 'sphinx.ext.mathjax', 'numpydoc',
              'sphinx_reredirects', 'sphinx_copybutton',
              'sphinx.ext.graphviz', 'matplotlib.sphinxext.plot_directive',
              'myst_parser', 'convert-svg-to-pdf', 'sphinx.ext.intersphinx',
              ]

# Add redirects here. This should be done whenever a page that is in the
# existing release docs is moved somewhere else so that the URLs don't break.
# The format is

# "page/path/without/extension": "../relative_path_with.html"

# Note that the html path is relative to the redirected page. Always test the
# redirect manually (they aren't tested automatically). See
# https://documatt.gitlab.io/sphinx-reredirects/usage.html

redirects = {
    "guides/getting_started/install": "../../install.html",
    "documentation-style-guide": "contributing/documentation-style-guide.html",
    "gotchas": "explanation/gotchas.html",
    "special_topics/classification": "../explanation/classification.html",
    "special_topics/finite_diff_derivatives": "../explanation/finite_diff_derivatives.html",
    "special_topics/intro": "../explanation/index.html",
    "special_topics/index": "../explanation/index.html",
    "modules/index": "../reference/index.html",
    "modules/physics/index": "../../reference/public/physics/index.html",

    "guides/contributing/index": "../../contributing/index.html",
    "guides/contributing/dev-setup": "../../contributing/dev-setup.html",
    "guides/contributing/dependencies": "../../contributing/dependencies.html",
    "guides/contributing/build-docs": "../../contributing/build-docs.html",
    "guides/contributing/debug": "../../contributing/debug.html",
    "guides/contributing/docstring": "../../contributing/docstring.html",
    "guides/documentation-style-guide": "../../contributing/contributing/documentation-style-guide.html",
    "guides/make-a-contribution": "../../contributing/make-a-contribution.html",
    "guides/contributing/deprecations": "../../contributing/deprecations.html",

    "tutorial/preliminaries": "../tutorials/intro-tutorial/preliminaries.html",
    "tutorial/intro": "../tutorials/intro-tutorial/intro.html",
    "tutorial/index": "../tutorials/intro-tutorial/index.html",
    "tutorial/gotchas": "../tutorials/intro-tutorial/gotchas.html",
    "tutorial/features": "../tutorials/intro-tutorial/features.html",
    "tutorial/next": "../tutorials/intro-tutorial/next.html",
    "tutorial/basic_operations": "../tutorials/intro-tutorial/basic_operations.html",
    "tutorial/printing": "../tutorials/intro-tutorial/printing.html",
    "tutorial/simplification": "../tutorials/intro-tutorial/simplification.html",
    "tutorial/calculus": "../tutorials/intro-tutorial/calculus.html",
    "tutorial/solvers": "../tutorials/intro-tutorial/solvers.html",
    "tutorial/matrices": "../tutorials/intro-tutorial/matrices.html",
    "tutorial/manipulation": "../tutorials/intro-tutorial/manipulation.html",

}

html_baseurl = "https://docs.sympy.org/latest/"

# Configure Sphinx copybutton (see https://sphinx-copybutton.readthedocs.io/en/latest/use.html)
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

# Use this to use pngmath instead
#extensions = ['sphinx.ext.autodoc', 'sphinx.ext.viewcode', 'sphinx.ext.pngmath', ]

# Enable warnings for all bad cross references. These are turned into errors
# with the -W flag in the Makefile.
nitpicky = True

nitpick_ignore = [
    ('py:class', 'sympy.logic.boolalg.Boolean')
]

# To stop docstrings inheritance.
autodoc_inherit_docstrings = False

# See https://www.sympy.org/sphinx-math-dollar/
mathjax3_config = {
  "tex": {
    "inlineMath": [['\\(', '\\)']],
    "displayMath": [["\\[", "\\]"]],
  }
}

# Myst configuration (for .md files). See
# https://myst-parser.readthedocs.io/en/latest/syntax/optional.html
myst_enable_extensions = ["dollarmath", "linkify"]
myst_heading_anchors = 6
# myst_update_mathjax = False

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
sys.path.append(os.path.abspath("./_pygments"))
pygments_style = 'styles.SphinxHighContrastStyle'
pygments_dark_style = 'styles.NativeHighContrastStyle'

# Don't show the source code hyperlinks when using matplotlib plot directive.
plot_html_show_source_link = False

# Options for HTML output
# -----------------------

# The style sheet to use for HTML and HTML Help pages. A file of that name
# must exist either in Sphinx' static/ path, or in one of the custom paths
# given in html_static_path.
# html_style = 'default.css'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# was classic
# html_theme = "classic"
html_theme = "furo"

# Adjust the sidebar so that the entire sidebar is scrollable
html_sidebars = {
    "**": [
        "sidebar/scroll-start.html",
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/navigation.html",
        "sidebar/versions.html",
        "sidebar/scroll-end.html",
    ],
}

common_theme_variables = {
    # Main "SymPy green" colors. Many things uses these colors.
    "color-brand-primary": "#52833A",
    "color-brand-content": "#307748",

    # The left sidebar.
    "color-sidebar-background": "#3B5526",
    "color-sidebar-background-border": "var(--color-background-primary)",
    "color-sidebar-link-text": "#FFFFFF",
    "color-sidebar-brand-text": "var(--color-sidebar-link-text--top-level)",
    "color-sidebar-link-text--top-level": "#FFFFFF",
    "color-sidebar-item-background--hover": "var(--color-brand-primary)",
    "color-sidebar-item-expander-background--hover": "var(--color-brand-primary)",

    "color-link-underline--hover": "var(--color-link)",
    "color-api-keyword": "#000000bd",
    "color-api-name": "var(--color-brand-content)",
    "color-api-pre-name": "var(--color-brand-content)",
    "api-font-size": "var(--font-size--normal)",
    "color-foreground-secondary": "#53555B",

    # TODO: Add the other types of admonitions here if anyone uses them.
    "color-admonition-title-background--seealso": "#CCCCCC",
    "color-admonition-title--seealso": "black",
    "color-admonition-title-background--note": "#CCCCCC",
    "color-admonition-title--note": "black",
    "color-admonition-title-background--warning": "var(--color-problematic)",
    "color-admonition-title--warning": "white",
    "admonition-font-size": "var(--font-size--normal)",
    "admonition-title-font-size": "var(--font-size--normal)",

    # Note: this doesn't work. If we want to change this, we have to set
    # it as the .highlight background in custom.css.
    "color-code-background": "hsl(80deg 100% 95%)",

    "code-font-size": "var(--font-size--small)",
    "font-stack--monospace": 'DejaVu Sans Mono,"SFMono-Regular",Menlo,Consolas,Monaco,Liberation Mono,Lucida Console,monospace;'
    }

html_theme_options = {
    "light_css_variables": common_theme_variables,
    # The dark variables automatically inherit values from the light variables
    "dark_css_variables": {
        **common_theme_variables,
        "color-brand-primary": "#33CB33",
        "color-brand-content": "#1DBD1D",

        "color-api-keyword": "#FFFFFFbd",
        "color-api-overall": "#FFFFFF90",
        "color-api-paren": "#FFFFFF90",

        "color-sidebar-item-background--hover": "#52833A",
        "color-sidebar-item-expander-background--hover": "#52833A",
        # This is the color of the text in the right sidebar
        "color-foreground-secondary": "#9DA1AC",

        "color-admonition-title-background--seealso": "#555555",
        "color-admonition-title-background--note": "#555555",
        "color-problematic": "#B30000",
    },
    # See https://pradyunsg.me/furo/customisation/footer/
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/sympy/sympy",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}

# Add a header for PR preview builds. See the Circle CI configuration.
if os.environ.get("CIRCLECI") == "true":
    PR_NUMBER = os.environ.get('CIRCLE_PR_NUMBER')
    SHA1 = os.environ.get('CIRCLE_SHA1')
    html_theme_options['announcement'] = f"""This is a preview build from
SymPy pull request <a href="https://github.com/sympy/sympy/pull/{PR_NUMBER}">
#{PR_NUMBER}</a>. It was built against <a
href="https://github.com/sympy/sympy/pull/{PR_NUMBER}/commits/{SHA1}">{SHA1[:7]}</a>.
If you aren't looking for a PR preview, go to <a
href="https://docs.sympy.org/">the main SymPy documentation</a>. """

# custom.css contains changes that aren't possible with the above because they
# aren't specified in the Furo theme as CSS variables
html_css_files = ['custom.css']

# html_js_files = []

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
# html_copy_source = True

# Output file base name for HTML help builder.
htmlhelp_basename = 'SymPydoc'

language = 'en'

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
html_favicon = '../_build/logo/sympy-notailtext-favicon.ico'

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

# Enable links to other packages
intersphinx_mapping = {
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'mpmath': ('https://mpmath.org/doc/current/', None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}
# Require :external: to reference intersphinx. Prevents accidentally linking
# to something from matplotlib.
intersphinx_disabled_reftypes = ['*']

# Required for linkcode extension.
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
