from os import path

from docutils import nodes
from sphinx.environment.collectors import EnvironmentCollector
from sphinx import addnodes
from sphinx.util.osutil import relative_uri

__version__ = (0, 8, 0)


class SimpleTocTreeCollector(EnvironmentCollector):
    """A TocTree collector that saves toctrees in a simple dict.

    sphinx.environment.collectors.toctree.TocTreeCollector saves
    TocTree as docutils.nodes which are hard to work with...

    Executed once per document/page, at sphinx's "read" phase.

    Saved data example:
    >>> {
    >>>  'sections': [{'title': 'Demo', 'href': '#demo'}],
    >>>  'toctrees': [<toctree: >]
    >>> }
    """
    def enable(self, app):
        super().enable(app)
        # env is populated from cache, if not cache create/initalize attibute
        if not hasattr(app.env, 'toc_dict'):
            app.env.toc_dict = {}

    def clear_doc(self, app, env, docname):
        env.toc_dict.pop(docname, None)

    def merge_other(self, app, env, docnames, other):
        for docname in docnames:
            env.toc_dict[docname] = other.toc_dict[docname]


    def process_doc(self, app, doctree):
        docname = app.env.docname # sphinx mutates this, ouch!!!

        # print(f"================ Collector\n{docname}\n============\n")
        # get 1 level document toc (sections)
        section_nodes = [s for s in doctree if isinstance(s, nodes.section)]
        # if first level is a single section,
        # ignore it and use second level of sections
        if len(section_nodes) == 1:
            section2_nodes = [s for s in section_nodes[0]
                              if isinstance(s, nodes.section)]
            if section2_nodes: # do not replace with level-2 sections if None
                section_nodes = section2_nodes

        sections = []
        for node in section_nodes:
            sections.append({
                'title': node[0].astext(),
                'href': '#{}'.format(node['ids'][0]),
            })

        app.env.toc_dict[docname] = {
            'sections': sections,
            'toctrees': doctree.traverse(addnodes.toctree)
        }



def add_toctree_data(app, pagename, templatename, context, doctree):
    """Create toctree_data, used to build sidebar navigation

    :param pagename: The name of the page
    :type pagename: str
    :param templatename: The name of the templatename
    :type templatename: str
    :param context: The context
    :type context: dict
    :param doctree: A doctree
    :type doctree: docutils.nodes.document

    Add to `toctree_data` to `context` that will be available on templates.
    Although data is "global", it is called once per page because current
    page is "highlighted", and some part of TOC might be collapsed.

    :return: None
    """
    # print(f"---------- Context\n{pagename}\n-------------\n")

    # start from master_doc
    master = app.env.get_doctree(app.env.config.master_doc)

    # each toctree will create navigation section
    res = [] # list of top level toctrees in master_doc
    for tree in master.traverse(addnodes.toctree):

        # special case for toctree that includes a single item
        # that contains a nested toctree.
        # In this case, just use the referenced toctree directly
        if len(tree['entries']) == 1:
            entry_docname = tree['entries'][0][1]
            toctrees = list(app.env.toc_dict[entry_docname]['toctrees'])

            if toctrees:
                # FIXME
                assert len(toctrees) == 1, "Press: Not supported more then one toctree on nested toctree"
                tree = toctrees[0]

        current0 = False # same page might have multiple tocs

        # add toc tree items, expand one more level if toctree is current page
        entries = []
        # entries contain a pair (title, name|external url)
        # if name of document is given and not title, title is taken from document
        for title, name in tree['entries']:
            if not title:
                title = app.env.titles[name].astext()

            # check if entry is an external resource
            ext_resource = '://' in name

            current1 = (pagename == name)
            children = []
            if current1:
                current0 = True
                # if current, add another level
                children = app.env.toc_dict[name]['sections']
            # add page_toc for current page
            entries.append({
                'name': name,
                'title': title,
                'current': current1,
                'children': children,
                'ext_resource': ext_resource,
            })

        toc_docname = tree['parent'] # docname where this toc appears
        title = tree['caption']

        # Anchor element is the section containing the toc,
        # as the toc itself does not contain ID.
        anchor_id = ''

        # tree.parent is the parent docutils node.
        # First parent is "compound" node toctree-wrapper,
        # second parent is the section containing the toctree
        toc_section = tree.parent.parent
        if toc_section['ids']: # no id means toc actually not in a section
            # TODO: should we be strict about toc being inside a section
            anchor_id = toc_section['ids'][0]
            if not title:
                title = toc_section['names'][0]

        # sphinx `pathto` does not play nice with anchors when
        # `allow_sharp_as_current_path` is True
        baseuri = app.builder.get_target_uri(pagename).rsplit('#', 1)[0]
        toc_uri = app.builder.get_target_uri(toc_docname).rsplit('#', 1)[0]
        toc_href = '{}#{}'.format(relative_uri(baseuri, toc_uri), anchor_id)
        res.append({
            'docname': toc_docname,
            'href': toc_href,
            'title': title,
            'current': current0,
            'entries': entries,
        })
    context['toctree_data'] = res



def setup(app):
    app.add_env_collector(SimpleTocTreeCollector)
    app.connect('html-page-context', add_toctree_data)
    app.add_html_theme('press', path.abspath(path.dirname(__file__)))
