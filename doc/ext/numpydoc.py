"""
========
numpydoc
========

Sphinx extension that handles docstrings in the Numpy standard format. [1]

It will:

- Convert Parameters etc. sections to field lists.
- Convert See Also section to a See also entry.
- Renumber references.
- Extract the signature from the docstring, if it can't be determined
  otherwise.

.. [1] https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

"""

import re
import pydoc
import sphinx
import inspect
from collections.abc import Callable

if sphinx.__version__ < '1.0.1':
    raise RuntimeError("Sphinx 1.0.1 or newer is required")

from docscrape_sphinx import get_doc_object, SphinxDocString


def mangle_docstrings(app, what, name, obj, options, lines,
                      reference_offset=[0]):

    cfg = {'use_plots': app.config.numpydoc_use_plots,
           'show_class_members': app.config.numpydoc_show_class_members,
           'show_inherited_class_members':
           app.config.numpydoc_show_inherited_class_members,
           'class_members_toctree': app.config.numpydoc_class_members_toctree}

    u_NL = '\n'
    if what == 'module':
        # Strip top title
        pattern = '^\\s*[#*=]{4,}\\n[a-z0-9 -]+\\n[#*=]{4,}\\s*'
        title_re = re.compile(pattern, re.IGNORECASE | re.DOTALL)
        lines[:] = title_re.sub('', u_NL.join(lines)).split(u_NL)
    else:
        # jupyterlite-sphinx uses the "!! processed by numpydoc !!" tag to identify
        # section boundaries when parsing docstrings. However, if this tag appears
        # inconsistently or in the middle of content, it can cause index errors
        # during section parsing. Here, we ensure the tag appears exactly once at
        # the end of the docstring with consistent formatting (i.e., '...' on its
        # own line followed by the indented tag). This helps jupyterlite-sphinx
        # correctly identify where the Examples section ends per docstring.
        try:
            doc = get_doc_object(obj, what, u_NL.join(lines), config=cfg)
            if doc is None:
                return

            # Split into lines and strip the processed tag if it
            # is present.
            doc_str = str(doc)
            new_lines = [line for line in doc_str.split(u_NL)
                        if not line.strip().endswith('!! processed by numpydoc !!')]

            # Ensure the tag appears exactly once at the end
            new_lines.extend(['..', '    !! processed by numpydoc !!'])

            lines[:] = new_lines
        except Exception as e:
            import warnings
            warning_msg = f"Failed to process docstring for {name}: {str(e)}"
            warnings.warn(warning_msg)
            return

    # replace reference numbers so that there are no duplicates
    references = []
    for line in lines:
        line = line.strip()
        m = re.match('^.. \\[([a-z0-9_.-])\\]', line, re.IGNORECASE)
        if m:
            references.append(m.group(1))

    # start renaming from the longest string, to avoid overwriting parts
    if references:
        references.sort(key=lambda x: -len(x))
        for i, line in enumerate(lines):
            for r in references:
                if re.match('^\\d+$', r):
                    new_r = "R%d" % (reference_offset[0] + int(r))
                else:
                    new_r = "%s%d" % (r, reference_offset[0])
                lines[i] = lines[i].replace('[%s]_' % r,
                                            '[%s]_' % new_r)
                lines[i] = lines[i].replace('.. [%s]' % r,
                                            '.. [%s]' % new_r)

        reference_offset[0] += len(references)


def mangle_signature(app, what, name, obj, options, sig, retann):
    # Do not try to inspect classes that don't define `__init__`
    if (inspect.isclass(obj) and
        (not hasattr(obj, '__init__') or
            'initializes x; see ' in pydoc.getdoc(obj.__init__))):
        return '', ''

    if not (isinstance(obj, Callable) or
            hasattr(obj, '__argspec_is_invalid_')):
        return

    if not hasattr(obj, '__doc__'):
        return

    doc = SphinxDocString(pydoc.getdoc(obj))
    if doc['Signature']:
        sig = re.sub("^[^(]*", "", doc['Signature'])
        return sig, ''


def setup(app, get_doc_object_=get_doc_object):
    if not hasattr(app, 'add_config_value'):
        return  # probably called by nose, better bail out

    global get_doc_object
    get_doc_object = get_doc_object_

    app.connect('autodoc-process-docstring', mangle_docstrings)
    app.connect('autodoc-process-signature', mangle_signature)
    app.add_config_value('numpydoc_edit_link', None, False)
    app.add_config_value('numpydoc_use_plots', None, False)
    app.add_config_value('numpydoc_show_class_members', True, True)
    app.add_config_value('numpydoc_show_inherited_class_members', True, True)
    app.add_config_value('numpydoc_class_members_toctree', True, True)

    # Extra mangling domains
    app.add_domain(NumpyPythonDomain)
    app.add_domain(NumpyCDomain)

# ------------------------------------------------------------------------------
# Docstring-mangling domains
# ------------------------------------------------------------------------------

from docutils.statemachine import ViewList
from sphinx.domains.c import CDomain
from sphinx.domains.python import PythonDomain


class ManglingDomainBase:
    directive_mangling_map = {}

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.wrap_mangling_directives()

    def wrap_mangling_directives(self):
        for name, objtype in list(self.directive_mangling_map.items()):
            self.directives[name] = wrap_mangling_directive(
                self.directives[name], objtype)


class NumpyPythonDomain(ManglingDomainBase, PythonDomain):
    name = 'np'
    directive_mangling_map = {
        'function': 'function',
        'class': 'class',
        'exception': 'class',
        'method': 'function',
        'classmethod': 'function',
        'staticmethod': 'function',
        'attribute': 'attribute',
    }
    indices = []


class NumpyCDomain(ManglingDomainBase, CDomain):
    name = 'np-c'
    directive_mangling_map = {
        'function': 'function',
        'member': 'attribute',
        'macro': 'function',
        'type': 'class',
        'var': 'object',
    }


def wrap_mangling_directive(base_directive, objtype):
    class directive(base_directive):
        def run(self):
            env = self.state.document.settings.env

            name = None
            if self.arguments:
                m = re.match(r'^(.*\s+)?(.*?)(\(.*)?', self.arguments[0])
                name = m.group(2).strip()

            if not name:
                name = self.arguments[0]

            lines = list(self.content)
            mangle_docstrings(env.app, objtype, name, None, None, lines)
            self.content = ViewList(lines, self.content.parent)

            return base_directive.run(self)

    return directive
