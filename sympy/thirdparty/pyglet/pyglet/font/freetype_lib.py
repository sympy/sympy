# ----------------------------------------------------------------------------
# pyglet
# Copyright (c) 2006-2007 Alex Holkner
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
#  * Neither the name of the pyglet nor the names of its
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
__version__ = '$Id: freetype_lib.py 1322 2007-10-23 12:58:03Z Alex.Holkner $'

from ctypes import *

import pyglet.lib

_libfreetype = pyglet.lib.load_library('freetype')

_font_data = {}

def _get_function(name, argtypes, rtype):
    try:
        func = getattr(_libfreetype, name)
        func.argtypes = argtypes
        func.restype = rtype
        return func
    except AttributeError, e:
            raise ImportError(e)

FT_Done_FreeType = _get_function('FT_Done_FreeType', [c_void_p], None)
FT_Done_Face = _get_function('FT_Done_Face', [c_void_p], None)


class FT_LibraryRec(Structure):
    _fields_ = [
        ('dummy', c_int),
    ]

    def __del__(self, byref=byref, FT_Done_FreeType=FT_Done_FreeType):
        global _library
        FT_Done_FreeType(byref(self))
        _library = None
FT_Library = POINTER(FT_LibraryRec)

class FT_Glyph_Metrics(Structure):
    _fields_ = [
        ('width', c_long),
        ('height', c_long),
        ('horiBearingX', c_long),
        ('horiBearingY', c_long),
        ('horiAdvance', c_long),
        ('vertBearingX', c_long),
        ('vertBearingY', c_long),
        ('vertAdvance', c_long),
    ]

    def dump(self):
        for (name, type) in self._fields_:
            print 'FT_Glyph_Metrics', name, `getattr(self, name)`

class FT_Generic(Structure):
    _fields_ = [('data', c_void_p), ('finalizer', c_void_p)]

class FT_BBox(Structure):
    _fields_ = [('xMin', c_long), ('yMin', c_long), ('xMax', c_long),
        ('yMax', c_long)]

class FT_Vector(Structure):
    _fields_ = [('x', c_long), ('y', c_long)]

class FT_Bitmap(Structure):
    _fields_ = [
        ('rows', c_int),
        ('width', c_int),
        ('pitch', c_int),
        # declaring buffer as c_char_p confuses ctypes, poor dear
        ('buffer', POINTER(c_ubyte)),
        ('num_grays', c_short),
        ('pixel_mode', c_ubyte),
        ('palette_mode', c_char),
        ('palette', c_void_p),
    ]

class FT_Outline(Structure):
    _fields_ = [
        ('n_contours', c_short),      # number of contours in glyph
        ('n_points', c_short),        # number of points in the glyph
        ('points', POINTER(FT_Vector)),  # the outline's points
        ('tags', c_char_p),            # the points flags
        ('contours', POINTER(c_short)),  # the contour end points
        ('flags', c_int),             # outline masks
    ]

class FT_GlyphSlotRec(Structure):
    _fields_ = [
        ('library', FT_Library),
        ('face', c_void_p),
        ('next', c_void_p),
        ('reserved', c_uint),
        ('generic', FT_Generic),

        ('metrics', FT_Glyph_Metrics),
        ('linearHoriAdvance', c_long),
        ('linearVertAdvance', c_long),
        ('advance', FT_Vector),

        ('format', c_int),

        ('bitmap', FT_Bitmap),
        ('bitmap_left', c_int),
        ('bitmap_top', c_int),

        ('outline', FT_Outline),
        ('num_subglyphs', c_uint),
        ('subglyphs', c_void_p),

        ('control_data', c_void_p),
        ('control_len', c_long),

        ('lsb_delta', c_long),
        ('rsb_delta', c_long),
        ('other', c_void_p),
        ('internal', c_void_p),
    ]
FT_GlyphSlot = POINTER(FT_GlyphSlotRec)

class FT_Size_Metrics(Structure):
    _fields_ = [
        ('x_ppem', c_ushort),    # horizontal pixels per EM
        ('y_ppem', c_ushort),    # vertical pixels per EM

        ('x_scale', c_long),     # two scales used to convert font units
        ('y_scale', c_long),     # to 26.6 frac. pixel coordinates

        ('ascender', c_long),    # ascender in 26.6 frac. pixels
        ('descender', c_long),   # descender in 26.6 frac. pixels
        ('height', c_long),      # text height in 26.6 frac. pixels
        ('max_advance', c_long), # max horizontal advance, in 26.6 pixels
    ]

class FT_SizeRec(Structure):
    _fields_ = [
        ('face', c_void_p),
        ('generic', FT_Generic),
        ('metrics', FT_Size_Metrics),
        ('internal', c_void_p),
    ]
FT_Size = POINTER(FT_SizeRec)

class FT_Bitmap_Size(Structure):
    _fields_ = [
        ('height', c_ushort),
        ('width', c_ushort),
        ('size', c_long),
        ('x_ppem', c_long),
        ('y_ppem', c_long),
    ]

# face_flags values
FT_FACE_FLAG_SCALABLE          = 1 <<  0
FT_FACE_FLAG_FIXED_SIZES       = 1 <<  1
FT_FACE_FLAG_FIXED_WIDTH       = 1 <<  2
FT_FACE_FLAG_SFNT              = 1 <<  3
FT_FACE_FLAG_HORIZONTAL        = 1 <<  4
FT_FACE_FLAG_VERTICAL          = 1 <<  5
FT_FACE_FLAG_KERNING           = 1 <<  6
FT_FACE_FLAG_FAST_GLYPHS       = 1 <<  7
FT_FACE_FLAG_MULTIPLE_MASTERS  = 1 <<  8
FT_FACE_FLAG_GLYPH_NAMES       = 1 <<  9
FT_FACE_FLAG_EXTERNAL_STREAM   = 1 << 10
FT_FACE_FLAG_HINTER            = 1 << 11

class FT_FaceRec(Structure):
    _fields_ = [
          ('num_faces', c_long),
          ('face_index', c_long),

          ('face_flags', c_long),
          ('style_flags', c_long),

          ('num_glyphs', c_long),
          ('family_name', c_char_p),
          ('style_name', c_char_p),

          ('num_fixed_sizes', c_int),
          ('available_sizes', POINTER(FT_Bitmap_Size)),

          ('num_charmaps', c_int),
          ('charmaps', c_void_p),

          ('generic', FT_Generic),

          ('bbox', FT_BBox),

          ('units_per_EM', c_ushort),
          ('ascender', c_short),
          ('descender', c_short),
          ('height', c_short),

          ('max_advance_width', c_short),
          ('max_advance_height', c_short),

          ('underline_position', c_short),
          ('underline_thickness', c_short),

          ('glyph', FT_GlyphSlot),
          ('size', FT_Size),
          ('charmap', c_void_p),

          ('driver', c_void_p),
          ('memory', c_void_p),
          ('stream', c_void_p),

          ('sizes_list_head', c_void_p),
          ('sizes_list_tail', c_void_p),

          ('autohint', FT_Generic),
          ('extensions', c_void_p),
          ('internal', c_void_p),
    ]

    def dump(self):
        for (name, type) in self._fields_:
            print 'FT_FaceRec', name, `getattr(self, name)`

    def has_kerning(self):
        return self.face_flags & FT_FACE_FLAG_KERNING

FT_Face = POINTER(FT_FaceRec)

class Error(Exception):
    def __init__(self, message, errcode):
        self.message = message
        self.errcode = errcode

    def __str__(self):
        return '%s: %s (%s)'%(self.__class__.__name__, self.message,
            self._ft_errors.get(self.errcode, 'unknown error'))
    _ft_errors = {
        0x00: "no error" ,
        0x01: "cannot open resource" ,
        0x02: "unknown file format" ,
        0x03: "broken file" ,
        0x04: "invalid FreeType version" ,
        0x05: "module version is too low" ,
        0x06: "invalid argument" ,
        0x07: "unimplemented feature" ,
        0x08: "broken table" ,
        0x09: "broken offset within table" ,
        0x10: "invalid glyph index" ,
        0x11: "invalid character code" ,
        0x12: "unsupported glyph image format" ,
        0x13: "cannot render this glyph format" ,
        0x14: "invalid outline" ,
        0x15: "invalid composite glyph" ,
        0x16: "too many hints" ,
        0x17: "invalid pixel size" ,
        0x20: "invalid object handle" ,
        0x21: "invalid library handle" ,
        0x22: "invalid module handle" ,
        0x23: "invalid face handle" ,
        0x24: "invalid size handle" ,
        0x25: "invalid glyph slot handle" ,
        0x26: "invalid charmap handle" ,
        0x27: "invalid cache manager handle" ,
        0x28: "invalid stream handle" ,
        0x30: "too many modules" ,
        0x31: "too many extensions" ,
        0x40: "out of memory" ,
        0x41: "unlisted object" ,
        0x51: "cannot open stream" ,
        0x52: "invalid stream seek" ,
        0x53: "invalid stream skip" ,
        0x54: "invalid stream read" ,
        0x55: "invalid stream operation" ,
        0x56: "invalid frame operation" ,
        0x57: "nested frame access" ,
        0x58: "invalid frame read" ,
        0x60: "raster uninitialized" ,
        0x61: "raster corrupted" ,
        0x62: "raster overflow" ,
        0x63: "negative height while rastering" ,
        0x70: "too many registered caches" ,
        0x80: "invalid opcode" ,
        0x81: "too few arguments" ,
        0x82: "stack overflow" ,
        0x83: "code overflow" ,
        0x84: "bad argument" ,
        0x85: "division by zero" ,
        0x86: "invalid reference" ,
        0x87: "found debug opcode" ,
        0x88: "found ENDF opcode in execution stream" ,
        0x89: "nested DEFS" ,
        0x8A: "invalid code range" ,
        0x8B: "execution context too long" ,
        0x8C: "too many function definitions" ,
        0x8D: "too many instruction definitions" ,
        0x8E: "SFNT font table missing" ,
        0x8F: "horizontal header (hhea, table missing" ,
        0x90: "locations (loca, table missing" ,
        0x91: "name table missing" ,
        0x92: "character map (cmap, table missing" ,
        0x93: "horizontal metrics (hmtx, table missing" ,
        0x94: "PostScript (post, table missing" ,
        0x95: "invalid horizontal metrics" ,
        0x96: "invalid character map (cmap, format" ,
        0x97: "invalid ppem value" ,
        0x98: "invalid vertical metrics" ,
        0x99: "could not find context" ,
        0x9A: "invalid PostScript (post, table format" ,
        0x9B: "invalid PostScript (post, table" ,
        0xA0: "opcode syntax error" ,
        0xA1: "argument stack underflow" ,
        0xA2: "ignore" ,
        0xB0: "`STARTFONT' field missing" ,
        0xB1: "`FONT' field missing" ,
        0xB2: "`SIZE' field missing" ,
        0xB3: "`CHARS' field missing" ,
        0xB4: "`STARTCHAR' field missing" ,
        0xB5: "`ENCODING' field missing" ,
        0xB6: "`BBX' field missing" ,
        0xB7: "`BBX' too big" ,
    }

FT_LOAD_RENDER = 0x4

FT_F26Dot6 = c_long

FT_Init_FreeType = _get_function('FT_Init_FreeType',
    [POINTER(FT_Library)], c_int)
FT_New_Memory_Face = _get_function('FT_New_Memory_Face',
    [FT_Library, POINTER(c_byte), c_long, c_long, POINTER(FT_Face)], c_int)
FT_New_Face = _get_function('FT_New_Face',
    [FT_Library, c_char_p, c_long, POINTER(FT_Face)], c_int)
FT_Set_Pixel_Sizes = _get_function('FT_Set_Pixel_Sizes',
    [FT_Face, c_uint, c_uint], c_int)
FT_Set_Char_Size = _get_function('FT_Set_Char_Size',
    [FT_Face, FT_F26Dot6, FT_F26Dot6, c_uint, c_uint], c_int)
FT_Load_Glyph = _get_function('FT_Load_Glyph',
    [FT_Face, c_uint, c_int32], c_int)
FT_Get_Char_Index = _get_function('FT_Get_Char_Index',
    [FT_Face, c_ulong], c_uint)
FT_Load_Char = _get_function('FT_Load_Char',
    [FT_Face, c_ulong, c_int], c_int)
FT_Get_Kerning = _get_function('FT_Get_Kerning',
    [FT_Face, c_uint, c_uint, c_uint, POINTER(FT_Vector)], c_int)


_library = None
def ft_get_library():
    global _library
    if not _library:
        _library = FT_Library()
        error = FT_Init_FreeType(byref(_library))
        if error:
            raise FontException(
                'an error occurred during library initialization', error)
    return _library
