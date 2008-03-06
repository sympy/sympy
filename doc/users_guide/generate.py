#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Generate SymPy Documentation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Generates a bunch of html files containing the documentation.

    :copyright: 2006-2007 by Armin Ronacher, Georg Brandl.
    :license: BSD, see LICENSE for more details.

    Taken from the Jinja project and modified to produce nice SymPy docs.

    Prerequisites - on Debian install: python-pygments python-jinja
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import re
import inspect
from datetime import datetime
from cgi import escape

from docutils import nodes
from docutils.parsers.rst import directives
from docutils.core import publish_parts
from docutils.writers import html4css1

from jinja import Environment

from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import HtmlFormatter

def generate_list_of_filters():
    from jinja.filters import FILTERS
    result = []

    filters = {}
    for name, f in FILTERS.iteritems():
        if not f in filters:
            filters[f] = ([name], inspect.getdoc(f))
        else:
            filters[f][0].append(name)
    for names, _ in filters.itervalues():
        names.sort(key=lambda x: -len(x))

    for names, doc in sorted(filters.values(), key=lambda x: x[0][0].lower()):
        name = names[0]
        if len(names) > 1:
            aliases = '\n\n    :Aliases: %s\n' % ', '.join(names[1:])
        else:
            aliases = ''

        doclines = []
        for line in doc.splitlines():
            doclines.append('    ' + line)
        doc = '\n'.join(doclines)
        result.append('`%s`\n%s%s' % (name, doc, aliases))

    return '\n'.join(result)

def generate_list_of_tests():
    from jinja.tests import TESTS
    result = []

    tests = {}
    for name, f in TESTS.iteritems():
        if not f in tests:
            tests[f] = ([name], inspect.getdoc(f))
        else:
            tests[f][0].append(name)
    for names, _ in tests.itervalues():
        names.sort(key=lambda x: -len(x))

    for names, doc in sorted(tests.values(), key=lambda x: x[0][0].lower()):
        name = names[0]
        if len(names) > 1:
            aliases = '\n\n    :Aliases: %s\n' % ', '.join(names[1:])
        else:
            aliases = ''

        doclines = []
        for line in doc.splitlines():
            doclines.append('    ' + line)
        doc = '\n'.join(doclines)
        result.append('`%s`\n%s%s' % (name, doc, aliases))

    return '\n'.join(result)

def generate_list_of_loaders():
    from jinja import loaders as loader_module

    result = []
    loaders = []
    for item in loader_module.__all__:
        loaders.append(getattr(loader_module, item))
    loaders.sort(key=lambda x: x.__name__.lower())

    for loader in loaders:
        doclines = []
        for line in inspect.getdoc(loader).splitlines():
            doclines.append('    ' + line)
        result.append('`%s`\n%s' % (loader.__name__, '\n'.join(doclines)))

    return '\n\n'.join(result)

def generate_list_of_baseloaders():
    from jinja import loaders as loader_module

    result = []
    loaders = []
    for item in dir(loader_module):
        obj = getattr(loader_module, item)
        try:
            if issubclass(obj, loader_module.BaseLoader) and \
               obj.__name__ != 'BaseLoader' and \
               obj.__name__ not in loader_module.__all__:
                loaders.append(obj)
        except TypeError:
            pass
    loaders.sort(key=lambda x: x.__name__.lower())

    for loader in loaders:
        doclines = []
        for line in inspect.getdoc(loader).splitlines():
            doclines.append('    ' + line)
        result.append('`%s`\n%s' % (loader.__name__, '\n'.join(doclines)))

    return '\n\n'.join(result)

def generate_environment_doc():
    from jinja.environment import Environment
    return '%s\n\n%s' % (
        inspect.getdoc(Environment),
        inspect.getdoc(Environment.__init__)
    )

e = Environment()

PYGMENTS_FORMATTER = HtmlFormatter(style='pastie', cssclass='syntax')

LIST_OF_FILTERS = generate_list_of_filters()
LIST_OF_TESTS = generate_list_of_tests()
LIST_OF_LOADERS = generate_list_of_loaders()
LIST_OF_BASELOADERS = generate_list_of_baseloaders()
ENVIRONMENT_DOC = generate_environment_doc()
#CHANGELOG = file(os.path.join(os.path.dirname(__file__), os.pardir, 'CHANGES'))\
CHANGELOG = ""

FULL_TEMPLATE = e.from_string('''\
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN"
   "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
  <title>{{ title }} &mdash; SymPy Documentation</title>
  <meta http-equiv="content-type" content="text/html; charset=utf-8">
  <link rel="stylesheet" href="style.css" type="text/css">
  <style type="text/css">
    {{ style|e }}
  </style>
</head>
<body>
  <div id="content">
    {% if file_id == 'index' %}
      <div id="sympylogo"></div>
      <h2 class="subheading plain">{{ title }}</h2>
    {% else %}
      <h1 class="heading"><span>SymPy</span></h1>
      <h2 class="subheading">{{ title }}</h2>
    {% endif %}
    {% if file_id != 'index' or toc %}
    <div id="toc">
      <h2>Navigation</h2>
      <ul>
        <li><a href="index.html">back to index</a></li>
      </ul>
      {% if toc %}
        <h2>Contents</h2>
        <ul class="contents">
        {% for key, value in toc %}
          <li><a href="{{ key }}">{{ value }}</a></li>
        {% endfor %}
        </ul>
      {% endif %}
    </div>
    {% endif %}
    <div id="contentwrapper">
      {{ body }}
    </div>
  </div>
</body>
<!-- generated on: {{ generation_date }}
     file id: {{ file_id }} -->
</html>\
''')

PREPROC_TEMPLATE = e.from_string('''\
<!-- TITLE -->{{ title }}<!-- ENDTITLE -->
<!-- TOC -->{% for key, value in toc %}<li><a href="{{
    key }}">{{ value }}</a></li>{% endfor %}<!-- ENDTOC -->
<!-- BODY -->{{ body }}<!-- ENDBODY -->\
''')

def pygments_directive(name, arguments, options, content, lineno,
                      content_offset, block_text, state, state_machine):
    try:
        lexer = get_lexer_by_name(arguments[0])
    except ValueError:
        # no lexer found
        lexer = get_lexer_by_name('text')
    parsed = highlight(u'\n'.join(content), lexer, PYGMENTS_FORMATTER)
    return [nodes.raw('', parsed, format="html")]
pygments_directive.arguments = (1, 0, 1)
pygments_directive.content = 1
directives.register_directive('sourcecode', pygments_directive)


def create_translator(link_style):
    class Translator(html4css1.HTMLTranslator):
        def visit_reference(self, node):
            refuri = node.get('refuri')
            if refuri is not None and '/' not in refuri and refuri.endswith('.txt'):
                node['refuri'] = link_style(refuri[:-4])
            html4css1.HTMLTranslator.visit_reference(self, node)
    return Translator


class DocumentationWriter(html4css1.Writer):

    def __init__(self, link_style):
        html4css1.Writer.__init__(self)
        self.translator_class = create_translator(link_style)

    def translate(self):
        html4css1.Writer.translate(self)
        # generate table of contents
        contents = self.build_contents(self.document)
        contents_doc = self.document.copy()
        contents_doc.children = contents
        contents_visitor = self.translator_class(contents_doc)
        contents_doc.walkabout(contents_visitor)
        self.parts['toc'] = self._generated_toc

    def build_contents(self, node, level=0):
        sections = []
        i = len(node) - 1
        while i >= 0 and isinstance(node[i], nodes.section):
            sections.append(node[i])
            i -= 1
        sections.reverse()
        toc = []
        for section in sections:
            try:
                reference = nodes.reference('', '', refid=section['ids'][0], *section[0])
            except IndexError:
                continue
            ref_id = reference['refid']
            text = escape(reference.astext().encode('utf-8'))
            toc.append((ref_id, text))

        self._generated_toc = [('#%s' % href, caption) for href, caption in toc]
        # no further processing
        return []


def generate_documentation(data, link_style):
    writer = DocumentationWriter(link_style)
    data = data.replace('[[list_of_filters]]', LIST_OF_FILTERS)\
               .replace('[[list_of_tests]]', LIST_OF_TESTS)\
               .replace('[[list_of_loaders]]', LIST_OF_LOADERS)\
               .replace('[[list_of_baseloaders]]', LIST_OF_BASELOADERS)\
               .replace('[[environment_doc]]', ENVIRONMENT_DOC)\
               .replace('[[changelog]]', CHANGELOG)
    parts = publish_parts(
        data,
        writer=writer,
        settings_overrides={
            'initial_header_level': 2,
            'field_name_limit': 50,
        }
    )
    return {
        'title':        parts['title'].encode('utf-8'),
        'body':         parts['body'].encode('utf-8'),
        'toc':          parts['toc']
    }


def handle_file(filename, fp, dst, preproc):
    now = datetime.now()
    title = os.path.basename(filename)[:-4]
    content = fp.read().decode('utf-8')
    suffix = not preproc and '.html' or ''
    parts = generate_documentation(content, (lambda x: './%s%s' % (x, suffix)))
    result = file(os.path.join(dst, title + '.html'), 'w')
    c = dict(parts)
    c['style'] = PYGMENTS_FORMATTER.get_style_defs('.syntax')
    c['generation_date'] = now
    c['file_id'] = title
    if preproc:
        tmpl = PREPROC_TEMPLATE
    else:
        tmpl = FULL_TEMPLATE
    result.write(tmpl.render(c).encode('utf-8'))
    result.close()


def run(dst, preproc, sources=(), handle_file=handle_file):
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'src'))
    if not sources:
        sources = [os.path.join(path, fn) for fn in os.listdir(path)]
    for fn in sources:
        if not os.path.isfile(fn):
            continue
        print 'Processing %s' % fn
        f = open(fn)
        try:
            handle_file(fn, f, dst, preproc)
        finally:
            f.close()


def main(dst='html/', preproc=False, *sources):
    run(os.path.realpath(dst), str(preproc).lower() == 'true', sources)


if __name__ == '__main__':
    main(*sys.argv[1:])
