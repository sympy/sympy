# -*- coding: utf-8 -*-
"""
    sphinx.ext.mathjax
    ~~~~~~~~~~~~~~~~~~

    Allow `MathJax <http://mathjax.org/>`_ to be used to display math in
    Sphinx's HTML writer -- requires the MathJax JavaScript library on your
    webserver/computer.

    :copyright: Copyright 2007-2011 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from docutils import nodes

from sphinx.application import ExtensionError
from sphinx.ext.mathbase import setup_math as mathbase_setup


def html_visit_math(self, node):
    self.body.append(self.starttag(node, 'span', '', CLASS='math'))
    self.body.append(self.builder.config.mathjax_inline[0] +
                     self.encode(node['latex']) +
                     self.builder.config.mathjax_inline[1] + '</span>')
    raise nodes.SkipNode

def html_visit_displaymath(self, node):
    self.body.append(self.starttag(node, 'div', CLASS='math'))
    if node['nowrap']:
        self.body.append(self.builder.config.mathjax_display[0] +
                         node['latex'] +
                         self.builder.config.mathjax_display[1])
        self.body.append('</div>')
        raise nodes.SkipNode

    parts = [prt for prt in node['latex'].split('\n\n') if prt.strip()]
    for i, part in enumerate(parts):
        part = self.encode(part)
        if i == 0:
            # necessary to e.g. set the id property correctly
            if node['number']:
                self.body.append('<span class="eqno">(%s)</span>' %
                                 node['number'])
        if '&' in part or '\\\\' in part:
            self.body.append(self.builder.config.mathjax_display[0] +
                             '\\begin{split}' + part + '\\end{split}' +
                             self.builder.config.mathjax_display[1])
        else:
            self.body.append(self.builder.config.mathjax_display[0] + part +
                             self.builder.config.mathjax_display[1])
    self.body.append('</div>\n')
    raise nodes.SkipNode

def builder_inited(app):
    if not app.config.mathjax_path:
        raise ExtensionError('mathjax_path config value must be set for the '
                             'mathjax extension to work')
    app.add_javascript(app.config.mathjax_path)


def setup(app):
    mathbase_setup(app, (html_visit_math, None), (html_visit_displaymath, None))
    app.add_config_value('mathjax_path',
                         'http://cdn.mathjax.org/mathjax/latest/MathJax.js?'
                         'config=TeX-AMS-MML_HTMLorMML', False)
    app.add_config_value('mathjax_inline', [r'\(', r'\)'], 'html')
    app.add_config_value('mathjax_display', [r'\[', r'\]'], 'html')
    app.connect('builder-inited', builder_inited)
