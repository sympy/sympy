# ----------------------------------------------------------------------------
# pyglet
# Copyright (c) 2006-2008 Alex Holkner
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions 
# are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright 
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#  * Neither the name of pyglet nor the names of its
#    contributors may be used to endorse or promote products
#    derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# ----------------------------------------------------------------------------

'''Load fonts and render text.

This is a fairly-low level interface to text rendering.  Obtain a font using
`load`::

    from pyglet import font
    arial = font.load('Arial', 14, bold=True, italic=False)

pyglet will load any system-installed fonts.  You can add additional fonts
(for example, from your program resources) using `add_file` or
`add_directory`.

Obtain a list of `Glyph` objects for a string of text using the `Font`
object::

    text = 'Hello, world!'
    glyphs = arial.get_glyphs(text)

The most efficient way to render these glyphs is with a `GlyphString`::

    glyph_string = GlyphString(text, glyphs)
    glyph_string.draw()

There are also a variety of methods in both `Font` and
`GlyphString` to facilitate word-wrapping.

A convenient way to render a string of text is with a `Text`::

    text = Text(font, text)
    text.draw()

See the `pyglet.font.base` module for documentation on the base classes used
by this package.
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 2136 2008-07-24 00:50:15Z Alex.Holkner $'

import sys
import os
import math
import weakref

import pyglet
from pyglet.gl import *
from pyglet import gl
from pyglet import image
from pyglet import window

class GlyphString(object):
    '''An immutable string of glyphs that can be rendered quickly.

    This class is ideal for quickly rendering single or multi-line strings
    of text that use the same font.  To wrap text using a glyph string,
    call `get_break_index` to find the optimal breakpoint for each line,
    the repeatedly call `draw` for each breakpoint.

    :deprecated: Use `pyglet.text.layout` classes.
    '''

    def __init__(self, text, glyphs, x=0, y=0):
        '''Create a glyph string.
        
        The `text` string is used to determine valid breakpoints; all glyphs
        must have already been determined using
        `pyglet.font.base.Font.get_glyphs`.  The string
        will be positioned with the baseline of the left-most glyph at the
        given coordinates.
        
        :Parameters:
            `text` : str or unicode
                String to represent.
            `glyphs` : list of `pyglet.font.base.Glyph`
                Glyphs representing `text`.
            `x` : float
                X coordinate of the left-side bearing of the left-most glyph.
            `y` : float
                Y coordinate of the baseline.

        '''
        # Create an interleaved array in GL_T2F_V3F format and determine
        # state changes required.
        
        lst = []
        texture = None
        self.text = text
        self.states = []
        self.cumulative_advance = [] # for fast post-string breaking
        state_from = 0
        state_length = 0
        for i, glyph in enumerate(glyphs):
            if glyph.owner != texture:
                if state_length:
                    self.states.append((state_from, state_length, texture))
                texture = glyph.owner
                state_from = i
                state_length = 0
            state_length += 1
            t = glyph.tex_coords
            lst += [t[0], t[1], t[2], 1.,
                    x + glyph.vertices[0], y + glyph.vertices[1], 0., 1.,
                    t[3], t[4], t[5], 1.,
                    x + glyph.vertices[2], y + glyph.vertices[1], 0., 1.,
                    t[6], t[7], t[8], 1.,
                    x + glyph.vertices[2], y + glyph.vertices[3], 0., 1.,
                    t[9], t[10], t[11], 1.,
                    x + glyph.vertices[0], y + glyph.vertices[3], 0., 1.]
            x += glyph.advance
            self.cumulative_advance.append(x)
        self.states.append((state_from, state_length, texture))

        self.array = (c_float * len(lst))(*lst)
        self.width = x

    def get_break_index(self, from_index, width):
        '''Find a breakpoint within the text for a given width.
        
        Returns a valid breakpoint after `from_index` so that the text
        between `from_index` and the breakpoint fits within `width` pixels.

        This method uses precomputed cumulative glyph widths to give quick
        answer, and so is much faster than 
        `pyglet.font.base.Font.get_glyphs_for_width`.  

        :Parameters:
            `from_index` : int
                Index of text to begin at, or 0 for the beginning of the
                string. 
            `width` : float
                Maximum width to use.

        :rtype: int
        :return: the index of text which will be used as the breakpoint, or
            `from_index` if there is no valid breakpoint.
        '''
        to_index = from_index
        if from_index >= len(self.text):
            return from_index
        if from_index:
            width += self.cumulative_advance[from_index-1]
        for i, (c, w) in enumerate(
                zip(self.text[from_index:], 
                    self.cumulative_advance[from_index:])):
            if c in u'\u0020\u200b':
                to_index = i + from_index + 1
            if c == '\n':
                return i + from_index + 1
            if w > width:
                return to_index 
        return to_index

    def get_subwidth(self, from_index, to_index):
        '''Return the width of a slice of this string.

        :Parameters:
            `from_index` : int
                The start index of the string to measure.
            `to_index` : int
                The end index (exclusive) of the string to measure.

        :rtype: float
        '''
        if to_index <= from_index:
            return 0
        width = self.cumulative_advance[to_index-1] 
        if from_index:
            width -= self.cumulative_advance[from_index-1]
        return width

    def draw(self, from_index=0, to_index=None):
        '''Draw a region of the glyph string.
        
        Assumes texture state is enabled.  To enable the texture state::

            from pyglet.gl import *
            glEnable(GL_TEXTURE_2D)

        :Parameters:
            `from_index` : int
                Start index of text to render.
            `to_index` : int
                End index (exclusive) of text to render.

        '''
        if from_index >= len(self.text) or \
           from_index == to_index or \
           not self.text:
            return

        # XXX Safe to assume all required textures will use same blend state I
        # think.  (otherwise move this into loop)
        self.states[0][2].apply_blend_state()

        if from_index:
            glPushMatrix()
            glTranslatef(-self.cumulative_advance[from_index-1], 0, 0)
        if to_index is None:
            to_index = len(self.text)

        glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT)
        glInterleavedArrays(GL_T4F_V4F, 0, self.array)
        for state_from, state_length, texture in self.states:
            if state_from + state_length < from_index:
                continue
            state_from = max(state_from, from_index)
            state_length = min(state_length, to_index - state_from)
            if state_length <= 0:
                break
            glBindTexture(GL_TEXTURE_2D, texture.id)
            glDrawArrays(GL_QUADS, state_from * 4, state_length * 4)
        glPopClientAttrib()

        if from_index:
            glPopMatrix()

class _TextZGroup(pyglet.graphics.Group):
    z = 0

    def set_state(self):
        glTranslatef(0, 0, self.z)

    def unset_state(self):
        glTranslatef(0, 0, -self.z)

class Text(object):
    '''Simple displayable text.

    This is a convenience class for rendering strings of text.  It takes
    care of caching the vertices so the text can be rendered every frame with
    little performance penalty.

    Text can be word-wrapped by specifying a `width` to wrap into.  If the
    width is not specified, it gives the width of the text as laid out.

    :Ivariables:
        `x` : int
            X coordinate of the text
        `y` : int
            Y coordinate of the text

    :deprecated: Use `pyglet.text.Label`.
    '''

    # Alignment constants

    #: Align the left edge of the text to the given X coordinate.
    LEFT = 'left'
    #: Align the horizontal center of the text to the given X coordinate.
    CENTER = 'center'
    #: Align the right edge of the text to the given X coordinate.
    RIGHT = 'right'
    #: Align the bottom of the descender of the final line of text with the
    #: given Y coordinate.
    BOTTOM = 'bottom'
    #: Align the baseline of the first line of text with the given Y
    #: coordinate.
    BASELINE = 'baseline'
    #: Align the top of the ascender of the first line of text with the given
    #: Y coordinate.
    TOP = 'top'

    # None: no multiline
    # 'width': multiline, wrapped to width
    # 'multiline': multiline, no wrap
    _wrap = None

    # Internal bookkeeping for wrap only.
    _width = None

    def __init__(self, font, text='', x=0, y=0, z=0, color=(1,1,1,1),
            width=None, halign=LEFT, valign=BASELINE):
        '''Create displayable text.

        :Parameters:
            `font` : `Font`
                Font to render the text in.
            `text` : str
                Initial string to render.
            `x` : float
                X coordinate of the left edge of the text.
            `y` : float
                Y coordinate of the baseline of the text.  If the text is
                word-wrapped, this refers to the first line of text.
            `z` : float
                Z coordinate of the text plane.
            `color` : 4-tuple of float
                Color to render the text in.  Alpha values can be specified
                in the fourth component.
            `width` : float
                Width to limit the rendering to. Text will be word-wrapped
                if necessary.
            `halign` : str
                Alignment of the text.  See `Text.halign` for details.
            `valign` : str
                Controls positioning of the text based off the y coordinate.
                One of BASELINE, BOTTOM, CENTER or TOP. Defaults to BASELINE.
        '''
        multiline = False
        if width is not None:
            self._width = width
            self._wrap = 'width'
            multiline = True
        elif '\n' in text:
            self._wrap = 'multiline'
            multiline = True

        self._group = _TextZGroup()
        self._document = pyglet.text.decode_text(text)
        self._layout = pyglet.text.layout.TextLayout(self._document, 
                                              width=width,
                                              multiline=multiline,
                                              dpi=font.dpi,
                                              group=self._group)

        self._layout.begin_update()
        if self._wrap == 'multiline':
            self._document.set_style(0, len(text), dict(wrap=False))
        self.font = font
        self.color = color
        self._x = x
        self.y = y
        self.z = z
        self.width = width
        self.halign = halign
        self.valign = valign
        self._update_layout_halign()
        self._layout.end_update()

    def _get_font(self):
        return self._font

    def _set_font(self, font):
        self._font = font
        self._layout.begin_update()
        self._document.set_style(0, len(self._document.text), {
            'font_name': font.name,
            'font_size': font.size,
            'bold': font.bold,
            'italic': font.italic,
        })
        self._layout._dpi = font.dpi
        self._layout.end_update()

    font = property(_get_font, _set_font)

    def _get_color(self):
        color = self._document.get_style('color')
        if color is None:
            return (1., 1., 1., 1.)
        return tuple([c/255. for c in color])

    def _set_color(self, color):
        color = [int(c * 255) for c in color]
        self._document.set_style(0, len(self._document.text), {
            'color': color,
        })

    color = property(_get_color, _set_color)

    def _update_layout_halign(self):
        if self._layout.multiline:
            # TextLayout has a different interpretation of halign that doesn't
            # consider the width to be a special factor; here we emulate the
            # old behaviour by fudging the layout x value.
            if self._layout.anchor_x == 'left':
                self._layout.x = self.x
            elif self._layout.anchor_x == 'center':
                self._layout.x = self.x + self._layout.width - \
                    self._layout.content_width // 2
            elif self._layout.anchor_x == 'right': 
                self._layout.x = self.x + 2 * self._layout.width - \
                    self._layout.content_width
        else:
            self._layout.x = self.x

    def _get_x(self):
        return self._x

    def _set_x(self, x):
        self._x = x
        self._update_layout_halign()

    x = property(_get_x, _set_x)

    def _get_y(self):
        return self._layout.y

    def _set_y(self, y):
        self._layout.y = y

    y = property(_get_y, _set_y)

    def _get_z(self):
        return self._group.z

    def _set_z(self, z):
        self._group.z = z

    z = property(_get_z, _set_z)

    def _update_wrap(self):
        if self._width is not None:
            self._wrap = 'width'
        elif '\n' in self.text:
            self._wrap = 'multiline'

        self._layout.begin_update()
        if self._wrap == None:
            self._layout.multiline = False
        elif self._wrap == 'width':
            self._layout.multiline = True
            self._layout.width = self._width
            self._document.set_style(0, len(self.text), dict(wrap=True))
        elif self._wrap == 'multiline':
            self._layout.multiline = True
            self._document.set_style(0, len(self.text), dict(wrap=False))
        self._update_layout_halign()
        self._layout.end_update()

    def _get_width(self):
        if self._wrap == 'width':
            return self._layout.width
        else:
            return self._layout.content_width

    def _set_width(self, width):
        self._width = width
        self._update_wrap()

    width = property(_get_width, _set_width, 
        doc='''Width of the text.

        When set, this enables word-wrapping to the specified width.
        Otherwise, the width of the text as it will be rendered can be
        determined.
        
        :type: float
        ''')

    def _get_height(self):
        return self._layout.content_height

    height = property(_get_height,
        doc='''Height of the text.
        
        This property is the ascent minus the descent of the font, unless
        there is more than one line of word-wrapped text, in which case
        the height takes into account the line leading.  Read-only.

        :type: float
        ''')

    def _get_text(self):
        return self._document.text

    def _set_text(self, text):
        self._document.text = text
        self._update_wrap()

    text = property(_get_text, _set_text,
        doc='''Text to render.

        The glyph vertices are only recalculated as needed, so multiple
        changes to the text can be performed with no performance penalty.
        
        :type: str
        ''')

    def _get_halign(self):
        return self._layout.anchor_x

    def _set_halign(self, halign):
        self._layout.anchor_x = halign
        self._update_layout_halign()

    halign = property(_get_halign, _set_halign,
        doc='''Horizontal alignment of the text.

        The text is positioned relative to `x` and `width` according to this
        property, which must be one of the alignment constants `LEFT`,
        `CENTER` or `RIGHT`.

        :type: str
        ''')

    def _get_valign(self):
        return self._layout.anchor_y

    def _set_valign(self, valign):
        self._layout.anchor_y = valign

    valign = property(_get_valign, _set_valign,
        doc='''Vertical alignment of the text.

        The text is positioned relative to `y` according to this property,
        which must be one of the alignment constants `BOTTOM`, `BASELINE`,
        `CENTER` or `TOP`.

        :type: str
        ''')

    def _get_leading(self):
        return self._document.get_style('leading') or 0

    def _set_leading(self, leading):
        self._document.set_style(0, len(self._document.text), {
            'leading': leading,
        })

    leading = property(_get_leading, _set_leading,
        doc='''Vertical space between adjacent lines, in pixels.

        :type: int
        ''')

    def _get_line_height(self):
        return self._font.ascent - self._font.descent + self.leading

    def _set_line_height(self, line_height):
        self.leading = line_height - (self._font.ascent - self._font.descent)

    line_height = property(_get_line_height, _set_line_height,
        doc='''Vertical distance between adjacent baselines, in pixels.

        :type: int
        ''')

    def draw(self):
        self._layout.draw()

if not getattr(sys, 'is_epydoc', False):
    if sys.platform == 'darwin':
        from pyglet.font.carbon import CarbonFont
        _font_class = CarbonFont
    elif sys.platform in ('win32', 'cygwin'):
        if pyglet.options['font'][0] == 'win32':
            from pyglet.font.win32 import Win32Font
            _font_class = Win32Font
        elif pyglet.options['font'][0] == 'gdiplus':
            from pyglet.font.win32 import GDIPlusFont
            _font_class = GDIPlusFont
        else:
            assert False, 'Unknown font driver'
    else:
        from pyglet.font.freetype import FreeTypeFont
        _font_class = FreeTypeFont

def load(name=None, size=None, bold=False, italic=False, dpi=None):
    '''Load a font for rendering.

    :Parameters:
        `name` : str, or list of str
            Font family, for example, "Times New Roman".  If a list of names
            is provided, the first one matching a known font is used.  If no
            font can be matched to the name(s), a default font is used.  In
            pyglet 1.1, the name may be omitted.
        `size` : float
            Size of the font, in points.  The returned font may be an exact
            match or the closest available.  In pyglet 1.1, the size may be
            omitted, and defaults to 12pt.
        `bold` : bool
            If True, a bold variant is returned, if one exists for the given
            family and size.
        `italic` : bool
            If True, an italic variant is returned, if one exists for the given
            family and size.
        `dpi` : float
            The assumed resolution of the display device, for the purposes of
            determining the pixel size of the font.  Defaults to 96. 

    :rtype: `Font`
    '''
    # Arbitrary default size
    if size is None:
        size = 12

    if dpi is None:
        dpi = 96

    # Find first matching name
    if type(name) in (tuple, list):
        for n in name:
            if _font_class.have_font(n):
                name = n
                break
        else:
            name = None

    # Locate or create font cache   
    shared_object_space = gl.current_context.object_space
    if not hasattr(shared_object_space, 'pyglet_font_font_cache'):
        shared_object_space.pyglet_font_font_cache = \
            weakref.WeakValueDictionary()
        shared_object_space.pyglet_font_font_hold = []
    font_cache = shared_object_space.pyglet_font_font_cache
    font_hold = shared_object_space.pyglet_font_font_hold

    # Look for font name in font cache
    descriptor = (name, size, bold, italic, dpi)
    if descriptor in font_cache:
        return font_cache[descriptor]

    # Not in cache, create from scratch
    font = _font_class(name, size, bold=bold, italic=italic, dpi=dpi)

    # Save parameters for new-style layout classes to recover
    font.name = name
    font.size = size
    font.bold = bold
    font.italic = italic
    font.dpi = dpi

    # Cache font in weak-ref dictionary to avoid reloading while still in use
    font_cache[descriptor] = font

    # Hold onto refs of last three loaded fonts to prevent them being
    # collected if momentarily dropped.
    del font_hold[3:]
    font_hold.insert(0, font)
    
    return font

def add_file(font):
    '''Add a font to pyglet's search path.

    In order to load a font that is not installed on the system, you must
    call this method to tell pyglet that it exists.  You can supply
    either a filename or any file-like object.

    The font format is platform-dependent, but is typically a TrueType font
    file containing a single font face.  Note that to load this file after
    adding it you must specify the face name to `load`, not the filename.

    :Parameters:
        `font` : str or file
            Filename or file-like object to load fonts from.

    '''
    if type(font) in (str, unicode):
        font = open(font, 'rb')
    if hasattr(font, 'read'):
        font = font.read()
    _font_class.add_font_data(font)


def add_directory(dir):
    '''Add a directory of fonts to pyglet's search path.

    This function simply calls `add_file` for each file with a ``.ttf``
    extension in the given directory.  Subdirectories are not searched.

    :Parameters:
        `dir` : str
            Directory that contains font files.

    '''
    import os
    for file in os.listdir(dir):
        if file[-4:].lower() == '.ttf':
            add_file(os.path.join(dir, file))

