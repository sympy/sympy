# -*- coding: utf-8 -*-

# Copyright (c) 2012-2013 by Christoph Reller. All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.

#    2. Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY CHRISTOPH RELLER ''AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
# EVENT SHALL CHRISTOPH RELLER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are
# those of the authors and should not be interpreted as representing official
# policies, either expressed or implied, of Christoph Reller.

"""
    sphinxcontrib.tikz
    ~~~~~~~~~~~~~~~~~~

    Draw pictures with the `TikZ/PGF LaTeX package.

    See README.rst file for details

    Author: Christoph Reller <christoph.reller@gmail.com>
    Version: 0.4.1
"""    

import tempfile
import posixpath
import shutil
import sys
from os import path, getcwd, chdir, mkdir, system
from subprocess import Popen, PIPE, call
try:
    from hashlib import sha1 as sha
except ImportError:
    from sha import sha

from docutils import nodes, utils
from docutils.parsers.rst import directives

from sphinx.errors import SphinxError
try:
    from sphinx.util.osutil import ensuredir, ENOENT, EPIPE
except:
    from sphinx.util import ensuredir, ENOENT, EPIPE
    
from sphinx.util.compat import Directive

class TikzExtError(SphinxError):
    category = 'Tikz extension error'

class tikzinline(nodes.Inline, nodes.Element):
    pass

def tikz_role(role, rawtext, text, lineno, inliner, option={}, content=[]):
    tikz = utils.unescape(text, restore_backslashes=True)
    return [tikzinline(tikz=tikz)], []

class tikz(nodes.Part, nodes.Element):
    pass

class TikzDirective(Directive):
    has_content = True
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec = {'libs':directives.unchanged,'stringsubst':directives.flag}

    def run(self):
        node = tikz()
        if not self.content:
            node['caption'] = ''
            node['tikz'] = '\n'.join(self.arguments)
        else:
            node['tikz'] = '\n'.join(self.content)
            node['caption'] = '\n'.join(self.arguments)
        node['libs'] = self.options.get('libs', '')
        if 'stringsubst' in self.options:
            node['stringsubst'] = True
        else:
            node['stringsubst'] = False
        return [node]

DOC_HEAD = r'''
\documentclass[12pt]{article}
\usepackage[utf8]{inputenc}
\usepackage{tikz}
\usetikzlibrary{%s}
\pagestyle{empty}
'''

DOC_BODY = r'''
\begin{document}
\begin{tikzpicture}
%s
\end{tikzpicture}
\end{document}
'''

def render_tikz(self,tikz,libs='',stringsubst=False):
    hashkey = tikz.encode('utf-8')
    fname = 'tikz-%s.png' % (sha(hashkey).hexdigest())
    relfn = posixpath.join(self.builder.imgpath, fname)
    outfn = path.join(self.builder.outdir, '_images', fname)

    if path.isfile(outfn):
        return relfn

    if hasattr(self.builder, '_tikz_warned'):
        return None
    
    ensuredir(path.dirname(outfn))
    curdir = getcwd()

    latex = DOC_HEAD % libs
    latex += self.builder.config.tikz_latex_preamble
    if stringsubst:
        tikz = tikz % {'wd': curdir}
    latex += DOC_BODY % tikz
    if isinstance(latex, unicode):
        latex = latex.encode('utf-8')

    if not hasattr(self.builder, '_tikz_tempdir'):
        tempdir = self.builder._tikz_tempdir = tempfile.mkdtemp()
    else:
        tempdir = self.builder._tikz_tempdir

    chdir(tempdir)

    tf = open('tikz.tex', 'wb')
    tf.write(latex)
    tf.close()

    try:
        try:
            p = Popen(['pdflatex', '--interaction=nonstopmode', 'tikz.tex'],
                      stdout=PIPE, stderr=PIPE)
        except OSError, err:
            if err.errno != ENOENT:   # No such file or directory
                raise
            self.builder.warn('LaTeX command cannot be run')
            self.builder._tikz_warned = True
            return None
    finally:
        chdir(curdir)

    stdout, stderr = p.communicate()
    if p.returncode != 0:
        raise TikzExtError('Error (tikz extension): latex exited with error:\n'
                           '[stderr]\n%s\n[stdout]\n%s' % (stderr, stdout))

    chdir(tempdir)

    # the following does not work for pdf patterns
    # p1 = Popen(['convert', '-density', '120', '-colorspace', 'rgb',
    #             '-trim', 'tikz.pdf', outfn], stdout=PIPE, stderr=PIPE)
    # stdout, stderr = p1.communicate()

    try:
        p = Popen(['pdftoppm', '-r', '120', 'tikz.pdf', 'tikz'],
                  stdout=PIPE, stderr=PIPE)
    except OSError, e:
        if e.errno != ENOENT:   # No such file or directory
            raise
        self.builder.warn('pdftoppm command cannot be run')
        self.builder.warn(err)
        self.builder._tikz_warned = True
        chdir(curdir)
        return None
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        self.builder._tikz_warned = True
        raise TikzExtError('Error (tikz extension): pdftoppm exited with error:'
                           '\n[stderr]\n%s\n[stdout]\n%s' % (stderr, stdout))

    if self.builder.config.tikz_proc_suite == 'ImageMagick':
        convert_args = []
        if self.builder.config.tikz_transparent:
            convert_args = ['-fuzz', '2%', '-transparent', 'white']

        try:
            p1 = Popen(['convert', '-trim'] + convert_args +
                       ['tikz-1.ppm', outfn],
                       stdout=PIPE, stderr=PIPE)
        except OSError, e:
            if e.errno != ENOENT:   # No such file or directory
                raise
            self.builder.warn('convert command cannot be run')
            self.builder.warn(err)
            self.builder._tikz_warned = True
            chdir(curdir)
            return None
        stdout, stderr = p1.communicate()
        if p1.returncode != 0:
            self.builder._tikz_warned = True
            chdir(curdir)
            raise TikzExtError('Error (tikz extension): convert exited with '
                               'error:\n[stderr]\n%s\n[stdout]\n%s'
                               % (stderr, stdout))

    elif self.builder.config.tikz_proc_suite == 'Netpbm':
        try:
            p1 = Popen(['pnmcrop', 'tikz-1.ppm'], stdout=PIPE, stderr=PIPE)
        except OSError, err:
            if err.errno != ENOENT:   # No such file or directory
                raise
            self.builder.warn('pnmcrop command cannot be run:')
            self.builder.warn(err)
            self.builder._tikz_warned = True
            chdir(curdir)
            return None

        pnm_args = []
        if self.builder.config.tikz_transparent:
            pnm_args = ['-transparent', 'white']
    
        try:
            p2 = Popen(['pnmtopng'] + pnm_args, stdin=p1.stdout,
                       stdout=PIPE, stderr=PIPE)
        except OSError, err:
            if err.errno != ENOENT:   # No such file or directory
                raise
            self.builder.warn('pnmtopng command cannot be run:')
            self.builder.warn(err)
            self.builder._tikz_warned = True
            chdir(curdir)
            return None
    
        pngdata, stderr2 = p2.communicate()
        dummy, stderr1 = p1.communicate()
        if p1.returncode != 0:
            self.builder._tikz_warned = True
            raise TikzExtError('Error (tikz extension): pnmcrop exited with '
                               'error:\n[stderr]\n%s' % (stderr1))
        if p2.returncode != 0:
            self.builder._tikz_warned = True
            raise TikzExtError('Error (tikz extension): pnmtopng exited with '
                               'error:\n[stderr]\n%s' % (stderr2))
        f = open(outfn,'wb')
        f.write(pngdata)
        f.close()

    else:
        self.builder._tikz_warned = True
        chdir(curdir)
        raise TikzExtError('Error (tikz extension): Invalid configuration '
                           'value for tikz_proc_suite')

    chdir(curdir)
    return relfn

def html_visit_tikzinline(self,node):
    libs = self.builder.config.tikz_tikzlibraries
    libs = libs.replace(' ', '').replace('\t', '').strip(', ')
    try:
        fname = render_tikz(self,node['tikz'],libs);
    except TikzExtError, exc:
        info = str(exc)[str(exc).find('!'):-1]
        sm = nodes.system_message(info, type='WARNING', level=2,
                                  backrefs=[], source=node['tikz'])
        sm.walkabout(self)
        self.builder.warn('display latex %r: \n' % node['tikz'] + str(exc))
        raise nodes.SkipNode
    if fname is None:
        # something failed -- use text-only as a bad substitute
        self.body.append('<span class="math">%s</span>' %
                         self.encode(node['tikz']).strip())
    else:
        self.body.append('<img class="math" src="%s" alt="%s"/>' %
                         (fname, self.encode(node['tikz']).strip()))
        raise nodes.SkipNode

def html_visit_tikz(self,node):
    libs = self.builder.config.tikz_tikzlibraries + ',' + node['libs']
    libs = libs.replace(' ', '').replace('\t', '').strip(', ')

    try:
        fname = render_tikz(self,node['tikz'],libs,node['stringsubst'])
    except TikzExtError, exc:
        info = str(exc)[str(exc).find('!'):-1]
        sm = nodes.system_message(info, type='WARNING', level=2,
                                  backrefs=[], source=node['tikz'])
        sm.walkabout(self)
        self.builder.warn('display latex %r: \n' % node['tikz'] + str(exc))
        raise nodes.SkipNode
    if fname is None:
        # something failed -- use text-only as a bad substitute
        self.body.append('<span class="math">%s</span>' %
                         self.encode(node['tikz']).strip())
    else:
        self.body.append(self.starttag(node, 'div', CLASS='figure'))
        self.body.append('<p>')
        self.body.append('<img src="%s" alt="%s" /></p>\n' %
                         (fname, self.encode(node['tikz']).strip()))
        if node['caption']:
            self.body.append('<p class="caption">%s</p>' %
                             self.encode(node['caption']).strip())
        self.body.append('</div>')
        raise nodes.SkipNode

def latex_visit_tikzinline(self, node):
    tikz = node['tikz']
    if tikz[0] == '[':
        cnt,pos = 1,1
        while cnt > 0 and cnt < len(tikz):
            if tikz[pos] == '[':
                cnt = cnt + 1
            if tikz[pos] == ']':
                cnt = cnt - 1
            pos = pos + 1
        tikz = tikz[:pos] + '{' + tikz[pos:]
    else:
        tikz = '{' + tikz
    self.body.append('\\tikz' + tikz + '}')
    raise nodes.SkipNode

def latex_visit_tikz(self, node):
    if node['caption']:
        latex = '\\begin{figure}[htp]\\centering\\begin{tikzpicture}' + \
                node['tikz'] + '\\end{tikzpicture}' + '\\caption{' + \
                self.encode(node['caption']).strip() + '}\\end{figure}'
    else:
        latex = '\\begin{center}\\begin{tikzpicture}' + node['tikz'] + \
            '\\end{tikzpicture}\\end{center}'
    self.body.append(latex)

def depart_tikz(self,node):
    pass

def cleanup_tempdir(app, exc):
    if exc:
        return
    if not hasattr(app.builder, '_tikz_tempdir'):
        return
    try:
        shutil.rmtree(app.builder._tikz_tempdir)
    except Exception:
        pass

def setup(app):
    app.add_node(tikz,
                 html=(html_visit_tikz, depart_tikz),
                 latex=(latex_visit_tikz, depart_tikz))
    app.add_node(tikzinline,
                 html=(html_visit_tikzinline, depart_tikz),
                 latex=(latex_visit_tikzinline, depart_tikz))
    app.add_role('tikz', tikz_role)
    app.add_directive('tikz', TikzDirective)
    app.add_config_value('tikz_latex_preamble', '', 'html')
    app.add_config_value('tikz_tikzlibraries', '', 'html')
    app.add_config_value('tikz_transparent', True, 'html')
    app.add_config_value('tikz_proc_suite', 'Netpbm', 'html')
    app.connect('build-finished', cleanup_tempdir)
