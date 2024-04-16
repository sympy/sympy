import sys
import re
import inspect
import textwrap
import pydoc
import sphinx
import collections

from docscrape import NumpyDocString, FunctionDoc, ClassDoc


class SphinxDocString(NumpyDocString):
    def __init__(self, docstring, config={}):
        NumpyDocString.__init__(self, docstring, config=config)
        self.load_config(config)

    def load_config(self, config):
        self.use_plots = config.get('use_plots', False)
        self.class_members_toctree = config.get('class_members_toctree', True)

    # string conversion routines
    def _str_header(self, name, symbol='`'):
        return ['.. rubric:: ' + name, '']

    def _str_field_list(self, name):
        return [':' + name + ':']

    def _str_indent(self, doc, indent=4):
        out = []
        for line in doc:
            out += [' '*indent + line]
        return out

    def _str_signature(self):
        return ['']
        if self['Signature']:
            return ['``%s``' % self['Signature']] + ['']
        else:
            return ['']

    def _str_summary(self):
        return self['Summary'] + ['']

    def _str_extended_summary(self):
        return self['Extended Summary'] + ['']

    def _str_returns(self, name='Returns'):
        out = []
        if self[name]:
            out += self._str_field_list(name)
            out += ['']
            for param, param_type, desc in self[name]:
                if param_type:
                    out += self._str_indent(['**{}** : {}'.format(param.strip(),
                                                              param_type)])
                else:
                    out += self._str_indent([param.strip()])
                if desc:
                    out += ['']
                    out += self._str_indent(desc, 8)
                out += ['']
        return out

    def _str_param_list(self, name):
        out = []
        if self[name]:
            out += self._str_field_list(name)
            out += ['']
            for param, param_type, desc in self[name]:
                if param_type:
                    out += self._str_indent(['**{}** : {}'.format(param.strip(),
                                                              param_type)])
                else:
                    out += self._str_indent(['**%s**' % param.strip()])
                if desc:
                    out += ['']
                    out += self._str_indent(desc, 8)
                out += ['']
        return out

    @property
    def _obj(self):
        if hasattr(self, '_cls'):
            return self._cls
        elif hasattr(self, '_f'):
            return self._f
        return None

    def _str_member_list(self, name):
        """
        Generate a member listing, autosummary:: table where possible,
        and a table where not.

        """
        out = []
        if self[name]:
            out += ['.. rubric:: %s' % name, '']
            prefix = getattr(self, '_name', '')

            if prefix:
                prefix = '~%s.' % prefix

            # Lines that are commented out are used to make the
            # autosummary:: table. Since SymPy does not use the
            # autosummary:: functionality, it is easiest to just comment it
            # out.
            # autosum = []
            others = []
            for param, param_type, desc in self[name]:
                param = param.strip()

                # Check if the referenced member can have a docstring or not
                param_obj = getattr(self._obj, param, None)
                if not (callable(param_obj)
                        or isinstance(param_obj, property)
                        or inspect.isgetsetdescriptor(param_obj)):
                    param_obj = None

                # if param_obj and (pydoc.getdoc(param_obj) or not desc):
                #     # Referenced object has a docstring
                #     autosum += ["   %s%s" % (prefix, param)]
                # else:
                others.append((param, param_type, desc))

            # if autosum:
            #     out += ['.. autosummary::']
            #     if self.class_members_toctree:
            #         out += ['   :toctree:']
            #     out += [''] + autosum

            if others:
                maxlen_0 = max(3, max(len(x[0]) for x in others))
                hdr = "="*maxlen_0 + "  " + "="*10
                fmt = '%%%ds  %%s  ' % (maxlen_0,)
                out += ['', '', hdr]
                for param, param_type, desc in others:
                    desc = " ".join(x.strip() for x in desc).strip()
                    if param_type:
                        desc = "({}) {}".format(param_type, desc)
                    out += [fmt % (param.strip(), desc)]
                out += [hdr]
            out += ['']
        return out

    def _str_section(self, name):
        out = []
        if self[name]:
            out += self._str_header(name)
            out += ['']
            content = textwrap.dedent("\n".join(self[name])).split("\n")
            out += content
            out += ['']
        return out

    def _str_see_also(self, func_role):
        out = []
        if self['See Also']:
            see_also = super()._str_see_also(func_role)
            out = ['.. seealso::', '']
            out += self._str_indent(see_also[2:])
        return out

    def _str_warnings(self):
        out = []
        if self['Warnings']:
            out = ['.. warning::', '']
            out += self._str_indent(self['Warnings'])
        return out

    def _str_index(self):
        idx = self['index']
        out = []
        if len(idx) == 0:
            return out

        out += ['.. index:: %s' % idx.get('default', '')]
        for section, references in idx.items():
            if section == 'default':
                continue
            elif section == 'refguide':
                out += ['   single: %s' % (', '.join(references))]
            else:
                out += ['   {}: {}'.format(section, ','.join(references))]
        return out

    def _str_references(self):
        out = []
        if self['References']:
            out += self._str_header('References')
            if isinstance(self['References'], str):
                self['References'] = [self['References']]
            out.extend(self['References'])
            out += ['']
            # Latex collects all references to a separate bibliography,
            # so we need to insert links to it
            if sphinx.__version__ >= "0.6":
                out += ['.. only:: latex', '']
            else:
                out += ['.. latexonly::', '']
            items = []
            for line in self['References']:
                m = re.match(r'.. \[([a-z0-9._-]+)\]', line, re.I)
                if m:
                    items.append(m.group(1))
            out += ['   ' + ", ".join(["[%s]_" % item for item in items]), '']
        return out

    def _str_examples(self):
        examples_str = "\n".join(self['Examples'])

        if (self.use_plots and 'import matplotlib' in examples_str
                and 'plot::' not in examples_str):
            out = []
            out += self._str_header('Examples')
            out += ['.. plot::', '']
            out += self._str_indent(self['Examples'])
            out += ['']
            return out
        else:
            return self._str_section('Examples')

    def __str__(self, indent=0, func_role="obj"):
        out = []
        out += self._str_signature()
        out += self._str_index() + ['']
        out += self._str_summary()
        out += self._str_extended_summary()
        out += self._str_param_list('Parameters')
        out += self._str_returns('Returns')
        out += self._str_returns('Yields')
        for param_list in ('Other Parameters', 'Raises', 'Warns'):
            out += self._str_param_list(param_list)
        out += self._str_warnings()
        for s in self._other_keys:
            out += self._str_section(s)
        out += self._str_see_also(func_role)
        out += self._str_references()
        out += self._str_member_list('Attributes')
        out = self._str_indent(out, indent)
        return '\n'.join(out)


class SphinxFunctionDoc(SphinxDocString, FunctionDoc):
    def __init__(self, obj, doc=None, config={}):
        self.load_config(config)
        FunctionDoc.__init__(self, obj, doc=doc, config=config)


class SphinxClassDoc(SphinxDocString, ClassDoc):
    def __init__(self, obj, doc=None, func_doc=None, config={}):
        self.load_config(config)
        ClassDoc.__init__(self, obj, doc=doc, func_doc=None, config=config)


class SphinxObjDoc(SphinxDocString):
    def __init__(self, obj, doc=None, config={}):
        self._f = obj
        self.load_config(config)
        SphinxDocString.__init__(self, doc, config=config)


def get_doc_object(obj, what=None, doc=None, config={}):
    if inspect.isclass(obj):
        what = 'class'
    elif inspect.ismodule(obj):
        what = 'module'
    elif callable(obj):
        what = 'function'
    else:
        what = 'object'

    if what == 'class':
        return SphinxClassDoc(obj, func_doc=SphinxFunctionDoc, doc=doc,
                              config=config)
    elif what in ('function', 'method'):
        return SphinxFunctionDoc(obj, doc=doc, config=config)
    else:
        if doc is None:
            doc = pydoc.getdoc(obj)
        return SphinxObjDoc(obj, doc, config=config)
