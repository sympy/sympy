# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import sys
import inspect
import os
import subprocess
from datetime import datetime

import sympy

# If your extensions are in another directory, add it here.
sys.path = ['ext'] + sys.path
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = '2'
copyright = '2021, 2'
author = '2'

# The short X.Y version
version = '2'

# The full version, including alpha/beta/rc tags
release = '2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary',
              'sphinx.ext.napoleon', 'sphinx.ext.intersphinx']

napoleon_custom_sections = ["Explanation"]

autodoc_mock_imports = ['sympy.plotting', 'sympy.testing.benchmarking',
                         'sympy.utilities.benchmarking', 'sympy.conftest', 'sympy.galgebra', 
                         'sympy.solvers', 'sympy.integrals.rubi.rubi_test', 
                         'sympy.integrals.rubi.rubi_tests.tests.test_trinomials',
                         'sympy.integrals.rubi.rubi_tests.tests.test_tangent',
                         'sympy.integrals.rubi.rubi_tests.tests.test_1_2',
                         'sympy.integrals.rubi.rubi_tests.tests.test_1_3',
                         'sympy.integrals.rubi.rubi_tests.tests.test_1_4',
                         'sympy.integrals.rubi.rubi_tests.tests.test_exponential',
                         'sympy.integrals.rubi.rubi_tests.tests.test_hyperbolic_sine',
                         'sympy.integrals.rubi.rubi_tests.tests.test_inverse_sine',
                         'sympy.integrals.rubi.rubi_tests.tests.test_logarithms',
                         'sympy.integrals.rubi.rubi_tests.tests.test_miscellaneous_algebra',
                         'sympy.integrals.rubi.rubi_tests.tests.test_secant',
                         'sympy.integrals.rubi.rubi_tests.tests.test_sine',
                         'sympy.integrals.rubi.rubi_tests.tests.test_special_functions',
                         'sympy.integrals.rubi.rubi_tests.tests.test_inverse_hyperbolic_sine']

# To stop docstrings inheritance.
autodoc_inherit_docstrings = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

exclude_patterns = ['_build', '_templates', 'reference/*']

autosummary_generate = True

#The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

suppress_warnings = ['ref.citation', 'ref.footnote']

autosummary_imported_members = True

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

html_theme = "classic"

html_logo = '_static/sympylogo.png'
#html_favicon = '../_build/logo/sympy-notailtext-favicon.ico'
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
