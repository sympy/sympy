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

'''
'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: freetype.py 2084 2008-05-27 12:42:19Z Alex.Holkner $'

import ctypes
from ctypes import *
from warnings import warn

import pyglet.lib
from pyglet.font import base
from pyglet import image
from pyglet.font.freetype_lib import *

# fontconfig library definitions
fontconfig = pyglet.lib.load_library('fontconfig')

FcResult = c_int

fontconfig.FcPatternBuild.restype = c_void_p
fontconfig.FcFontMatch.restype = c_void_p
fontconfig.FcFreeTypeCharIndex.restype = c_uint

FC_FAMILY = 'family'
FC_SIZE = 'size'
FC_SLANT = 'slant'
FC_WEIGHT = 'weight'
FC_FT_FACE = 'ftface'
FC_FILE = 'file'

FC_WEIGHT_REGULAR = 80
FC_WEIGHT_BOLD = 200

FC_SLANT_ROMAN = 0
FC_SLANT_ITALIC = 100

FT_STYLE_FLAG_ITALIC = 1
FT_STYLE_FLAG_BOLD = 2

(FT_RENDER_MODE_NORMAL,
 FT_RENDER_MODE_LIGHT,
 FT_RENDER_MODE_MONO,
 FT_RENDER_MODE_LCD,
 FT_RENDER_MODE_LCD_V) = range(5)

def FT_LOAD_TARGET_(x):
    return (x & 15) << 16

FT_LOAD_TARGET_NORMAL = FT_LOAD_TARGET_(FT_RENDER_MODE_NORMAL)
FT_LOAD_TARGET_LIGHT = FT_LOAD_TARGET_(FT_RENDER_MODE_LIGHT)
FT_LOAD_TARGET_MONO = FT_LOAD_TARGET_(FT_RENDER_MODE_MONO)
FT_LOAD_TARGET_LCD = FT_LOAD_TARGET_(FT_RENDER_MODE_LCD)
FT_LOAD_TARGET_LCD_V = FT_LOAD_TARGET_(FT_RENDER_MODE_LCD_V)

(FT_PIXEL_MODE_NONE,
 FT_PIXEL_MODE_MONO,
 FT_PIXEL_MODE_GRAY,
 FT_PIXEL_MODE_GRAY2,
 FT_PIXEL_MODE_GRAY4,
 FT_PIXEL_MODE_LCD,
 FT_PIXEL_MODE_LCD_V) = range(7)

(FcTypeVoid,
 FcTypeInteger,
 FcTypeDouble, 
 FcTypeString, 
 FcTypeBool,
 FcTypeMatrix,
 FcTypeCharSet,
 FcTypeFTFace,
 FcTypeLangSet) = range(9)
FcType = c_int

(FcMatchPattern,
 FcMatchFont) = range(2)
FcMatchKind = c_int

class _FcValueUnion(Union):
    _fields_ = [
        ('s', c_char_p),
        ('i', c_int),
        ('b', c_int),
        ('d', c_double),
        ('m', c_void_p),
        ('c', c_void_p),
        ('f', c_void_p),
        ('p', c_void_p),
        ('l', c_void_p),
    ]

class FcValue(Structure):
    _fields_ = [
        ('type', FcType),
        ('u', _FcValueUnion)
    ]

# End of library definitions

def f16p16_to_float(value):
    return float(value) / (1 << 16)

def float_to_f16p16(value):
    return int(value * (1 << 16))

def f26p6_to_float(value):
    return float(value) / (1 << 6)

def float_to_f26p6(value):
    return int(value * (1 << 6))

class FreeTypeGlyphRenderer(base.GlyphRenderer):
    def __init__(self, font):
        super(FreeTypeGlyphRenderer, self).__init__(font)
        self.font = font

    def render(self, text):
        face = self.font.face
        FT_Set_Char_Size(face, 0, self.font._face_size, 
                         self.font._dpi, self.font._dpi)
        glyph_index = fontconfig.FcFreeTypeCharIndex(byref(face), ord(text[0]))
        error = FT_Load_Glyph(face, glyph_index, FT_LOAD_RENDER)
        if error != 0:
            raise base.FontException(
                'Could not load glyph for "%c"' % text[0], error) 
        glyph_slot = face.glyph.contents
        width = glyph_slot.bitmap.width
        height = glyph_slot.bitmap.rows
        baseline = height - glyph_slot.bitmap_top
        lsb = glyph_slot.bitmap_left
        advance = int(f26p6_to_float(glyph_slot.advance.x))
        mode = glyph_slot.bitmap.pixel_mode
        pitch = glyph_slot.bitmap.pitch
        
        if mode == FT_PIXEL_MODE_MONO:
            # BCF fonts always render to 1 bit mono, regardless of render
            # flags. (freetype 2.3.5)
            bitmap_data = cast(glyph_slot.bitmap.buffer, 
                               POINTER(c_ubyte * (pitch * height))).contents
            data = (c_ubyte * (pitch * 8 * height))()
            data_i = 0
            for byte in bitmap_data:
                # Data is MSB; left-most pixel in a byte has value 128.
                data[data_i + 0] = (byte & 0x80) and 255 or 0
                data[data_i + 1] = (byte & 0x40) and 255 or 0
                data[data_i + 2] = (byte & 0x20) and 255 or 0
                data[data_i + 3] = (byte & 0x10) and 255 or 0
                data[data_i + 4] = (byte & 0x08) and 255 or 0
                data[data_i + 5] = (byte & 0x04) and 255 or 0
                data[data_i + 6] = (byte & 0x02) and 255 or 0
                data[data_i + 7] = (byte & 0x01) and 255 or 0
                data_i += 8
            pitch <<= 3
        elif mode == FT_PIXEL_MODE_GRAY:
            # Usual case
            data = glyph_slot.bitmap.buffer
        else:
            raise base.FontException('Unsupported render mode for this glyph')

        # pitch should be negative, but much faster to just swap tex_coords
        img = image.ImageData(width, height, 'A', data, pitch)
        glyph = self.font.create_glyph(img) 
        glyph.set_bearings(baseline, lsb, advance)
        t = list(glyph.tex_coords)
        glyph.tex_coords = t[9:12] + t[6:9] + t[3:6] + t[:3]

        return glyph

class FreeTypeMemoryFont(object):
    def __init__(self, data):
        self.buffer = (ctypes.c_byte * len(data))()
        ctypes.memmove(self.buffer, data, len(data))

        ft_library = ft_get_library()
        self.face = FT_Face()
        r = FT_New_Memory_Face(ft_library, 
            self.buffer, len(self.buffer), 0, self.face)
        if r != 0:
            raise base.FontException('Could not load font data')

        self.name = self.face.contents.family_name
        self.bold = self.face.contents.style_flags & FT_STYLE_FLAG_BOLD != 0
        self.italic = self.face.contents.style_flags & FT_STYLE_FLAG_ITALIC != 0

        # Replace Freetype's generic family name with TTF/OpenType specific
        # name if we can find one; there are some instances where Freetype
        # gets it wrong.
        if self.face.contents.face_flags & FT_FACE_FLAG_SFNT:
            name = FT_SfntName()
            for i in range(FT_Get_Sfnt_Name_Count(self.face)):
                result = FT_Get_Sfnt_Name(self.face, i, name)
                if result != 0:
                    continue
                if not (name.platform_id == TT_PLATFORM_MICROSOFT and
                        name.encoding_id == TT_MS_ID_UNICODE_CS):
                    continue
                if name.name_id == TT_NAME_ID_FONT_FAMILY:
                    string = string_at(name.string, name.string_len)
                    self.name = string.decode('utf-16be', 'ignore')

    def __del__(self):
        try:
            FT_Done_Face(self.face)
        except:
            pass

class FreeTypeFont(base.Font):
    glyph_renderer_class = FreeTypeGlyphRenderer

    # Map font (name, bold, italic) to FreeTypeMemoryFont
    _memory_fonts = {}

    def __init__(self, name, size, bold=False, italic=False, dpi=None):
        super(FreeTypeFont, self).__init__()

        if dpi is None:
            dpi = 96  # as of pyglet 1.1; pyglet 1.0 had 72.

        # Check if font name/style matches a font loaded into memory by user
        lname = name and name.lower() or ''
        if (lname, bold, italic) in self._memory_fonts:
            font = self._memory_fonts[lname, bold, italic]
            self._set_face(font.face, size, dpi)
            return

        # Use fontconfig to match the font (or substitute a default).
        ft_library = ft_get_library()

        match = self.get_fontconfig_match(name, size, bold, italic)
        if not match:
            raise base.FontException('Could not match font "%s"' % name)

        f = FT_Face()
        if fontconfig.FcPatternGetFTFace(match, FC_FT_FACE, 0, byref(f)) != 0:
            value = FcValue()
            result = fontconfig.FcPatternGet(match, FC_FILE, 0, byref(value))
            if result != 0:
                raise base.FontException('No filename or FT face for "%s"' % \
                                         name)
            result = FT_New_Face(ft_library, value.u.s, 0, byref(f))
            if result:
                raise base.FontException('Could not load "%s": %d' % \
                                         (name, result))

        fontconfig.FcPatternDestroy(match)

        self._set_face(f, size, dpi)

    def _set_face(self, face, size, dpi):
        self.face = face.contents
        self._face_size = float_to_f26p6(size)
        self._dpi = dpi

        FT_Set_Char_Size(self.face, 0, float_to_f26p6(size), dpi, dpi)
        metrics = self.face.size.contents.metrics
        if metrics.ascender == 0 and metrics.descender == 0:
            # Workaround broken fonts with no metrics.  Has been observed with
            # courR12-ISO8859-1.pcf.gz: "Courier" "Regular"
            #
            # None of the metrics fields are filled in, so render a glyph and
            # grab its height as the ascent, and make up an arbitrary
            # descent.
            i = fontconfig.FcFreeTypeCharIndex(byref(self.face), ord('X'))
            FT_Load_Glyph(self.face, i, FT_LOAD_RENDER)
            self.ascent = self.face.available_sizes.contents.height
            self.descent = -self.ascent // 4  # arbitrary.
        else:
            self.ascent = int(f26p6_to_float(metrics.ascender))
            self.descent = int(f26p6_to_float(metrics.descender))

    @staticmethod
    def get_fontconfig_match(name, size, bold, italic):
        if bold:
            bold = FC_WEIGHT_BOLD
        else:
            bold = FC_WEIGHT_REGULAR

        if italic:
            italic = FC_SLANT_ITALIC
        else:
            italic = FC_SLANT_ROMAN

        fontconfig.FcInit()

        if isinstance(name, unicode):
            name = name.encode('utf8')

        pattern = fontconfig.FcPatternCreate()
        fontconfig.FcPatternAddDouble(pattern, FC_SIZE, c_double(size))
        fontconfig.FcPatternAddInteger(pattern, FC_WEIGHT, bold)
        fontconfig.FcPatternAddInteger(pattern, FC_SLANT, italic)
        fontconfig.FcPatternAddString(pattern, FC_FAMILY, name)
        fontconfig.FcConfigSubstitute(0, pattern, FcMatchPattern)
        fontconfig.FcDefaultSubstitute(pattern)

        # Look for a font that matches pattern
        result = FcResult()
        match = fontconfig.FcFontMatch(0, pattern, byref(result))
        fontconfig.FcPatternDestroy(pattern)
        
        return match

    @classmethod
    def have_font(cls, name):
        value = FcValue()
        match = cls.get_fontconfig_match(name, 12, False, False)
        result = fontconfig.FcPatternGet(match, FC_FAMILY, 0, byref(value))
        if value.u.s == name:
            return True
        else:
            name = name.lower()
            for font in cls._memory_fonts.values():
                if font.name.lower() == name:
                    return True
        return False
    
    @classmethod
    def add_font_data(cls, data):
        font = FreeTypeMemoryFont(data)
        cls._memory_fonts[font.name.lower(), font.bold, font.italic] = font
