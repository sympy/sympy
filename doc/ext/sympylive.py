"""
    sympylive
    ~~~~~~~~~

    Allow `SymPy Live <http://live.sympy.org/>`_ to be used for interactive
    evaluation of SymPy's code examples.

    :copyright: Copyright 2014 by the SymPy Development Team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""


def builder_inited(app):
    if not app.config.sympylive_url:
        raise ExtensionError('sympylive_url config value must be set'
                             ' for the sympylive extension to work')

    app.add_javascript(app.config.sympylive_url + '/static/utilities.js')
    app.add_javascript(app.config.sympylive_url + '/static/external/classy.js')

    app.add_stylesheet(app.config.sympylive_url + '/static/live-core.css')
    app.add_stylesheet(app.config.sympylive_url +
                       '/static/live-autocomplete.css')
    app.add_stylesheet(app.config.sympylive_url + '/static/live-sphinx.css')

    app.add_javascript(app.config.sympylive_url + '/static/live-core.js')
    app.add_javascript(app.config.sympylive_url +
                       '/static/live-autocomplete.js')
    app.add_javascript(app.config.sympylive_url + '/static/live-sphinx.js')


def setup(app):
    app.add_config_value('sympylive_url', 'http://live.sympy.org', False)
    app.connect('builder-inited', builder_inited)
