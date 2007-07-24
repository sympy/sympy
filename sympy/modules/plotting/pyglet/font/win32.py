#!/usr/bin/env python
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

# TODO Windows Vista: need to call SetProcessDPIAware?  May affect GDI+ calls
# as well as font.

from ctypes import *

from pyglet.font import base
import pyglet.image
from pyglet.window.win32.constants import *
from pyglet.window.win32.types import *
from pyglet.window.win32 import _gdi32 as gdi32, _user32 as user32
from pyglet.window.win32 import _kernel32 as kernel32
from pyglet.window.win32 import _check

HFONT = HANDLE
HBITMAP = HANDLE
HDC = HANDLE
HGDIOBJ = HANDLE
gdi32.CreateFontIndirectA.restype = HFONT
gdi32.CreateCompatibleBitmap.restype = HBITMAP
gdi32.CreateCompatibleDC.restype = HDC
user32.GetDC.restype = HDC
gdi32.GetStockObject.restype = HGDIOBJ
gdi32.CreateDIBSection.restype = HBITMAP

class LOGFONT(Structure):
    _fields_ = [
        ('lfHeight', c_long),
        ('lfWidth', c_long),
        ('lfEscapement', c_long),
        ('lfOrientation', c_long),
        ('lfWeight', c_long),
        ('lfItalic', c_byte),
        ('lfUnderline', c_byte),
        ('lfStrikeOut', c_byte),
        ('lfCharSet', c_byte),
        ('lfOutPrecision', c_byte),
        ('lfClipPrecision', c_byte),
        ('lfQuality', c_byte),
        ('lfPitchAndFamily', c_byte),
        ('lfFaceName', (c_char * LF_FACESIZE))  # Use ASCII
    ]
    __slots__ = [f[0] for f in _fields_]

class TEXTMETRIC(Structure):
    _fields_ = [
        ('tmHeight', c_long),
        ('tmAscent', c_long),
        ('tmDescent', c_long),
        ('tmInternalLeading', c_long),
        ('tmExternalLeading', c_long),
        ('tmAveCharWidth', c_long),
        ('tmMaxCharWidth', c_long),
        ('tmWeight', c_long),
        ('tmOverhang', c_long),
        ('tmDigitizedAspectX', c_long),
        ('tmDigitizedAspectY', c_long),
        ('tmFirstChar', c_char),  # Use ASCII 
        ('tmLastChar', c_char),
        ('tmDefaultChar', c_char),
        ('tmBreakChar', c_char),
        ('tmItalic', c_byte),
        ('tmUnderlined', c_byte),
        ('tmStruckOut', c_byte),
        ('tmPitchAndFamily', c_byte),
        ('tmCharSet', c_byte)
    ]
    __slots__ = [f[0] for f in _fields_]

class ABC(Structure):
    _fields_ = [
        ('abcA', c_int),
        ('abcB', c_uint),
        ('abcC', c_int)
    ]
    __slots__ = [f[0] for f in _fields_]

class BITMAPINFOHEADER(Structure):
    _fields_ = [
        ('biSize', c_uint32),
        ('biWidth', c_int),
        ('biHeight', c_int),
        ('biPlanes', c_short),
        ('biBitCount', c_short),
        ('biCompression', c_uint32),
        ('biSizeImage', c_uint32),
        ('biXPelsPerMeter', c_long),
        ('biYPelsPerMeter', c_long),
        ('biClrUsed', c_uint32),
        ('biClrImportant', c_uint32)
    ]
    __slots__ = [f[0] for f in _fields_]

class RGBQUAD(Structure):
    _fields_ = [
        ('rgbBlue', c_byte),
        ('rgbGreen', c_byte),
        ('rgbRed', c_byte),
        ('rgbReserved', c_byte)
    ]

    def __init__(self, r, g, b):
        self.rgbRed = r
        self.rgbGreen = g
        self.rgbBlue = b


class BITMAPINFO(Structure):
    _fields_ = [
        ('bmiHeader', BITMAPINFOHEADER),
        ('bmiColors', c_ulong * 3)
    ]

def str_ucs2(text):
    if byteorder == 'big':
        text = text.encode('utf_16_be')
    else:
        text = text.encode('utf_16_le')   # explicit endian avoids BOM
    return create_string_buffer(text + '\0')

class Win32GlyphRenderer(base.GlyphRenderer):
    _bitmap = None
    _bitmap_dc = None
    _bitmap_rect = None

    def __init__(self, font):
        super(Win32GlyphRenderer, self).__init__(font)
        self.font = font
        self._create_bitmap_dc(font.max_glyph_width, font.ascent - font.descent) 
        self._black = gdi32.GetStockObject(BLACK_BRUSH)
        self._white = gdi32.GetStockObject(WHITE_BRUSH)

    def __del__(self):
        if self._bitmap_dc:
            gdi32.DeleteDC(self._bitmap_dc)
        if self._bitmap:
            gdi32.DeleteObject(self._bitmap)

    def render(self, text):
        # There's no such thing as a greyscale bitmap format in GDI.  We can
        # create an 8-bit palette bitmap with 256 shades of grey, but
        # unfortunately antialiasing will not work on such a bitmap.  So, we
        # use a 32-bit bitmap with a custom bit-mask (ABGR) to push the red
        # channel into where GL expects the alpha channel to be (using RGBA).  
        # GL ignores the remaining color components.
    
        gdi32.SelectObject(self._bitmap_dc, self._bitmap)
        gdi32.SelectObject(self._bitmap_dc, self.font.hfont)
        gdi32.SetBkColor(self._bitmap_dc, 0)
        gdi32.SetTextColor(self._bitmap_dc, 0x00ffffff)

        # Attempt to get ABC widths (only for TrueType)
        abc = ABC()
        if gdi32.GetCharABCWidthsW(self._bitmap_dc, 
            ord(text), ord(text), byref(abc)):
            width = abc.abcB 
            lsb = abc.abcA
            advance = abc.abcA + abc.abcB + abc.abcC
        else:
            width_buf = c_int()
            gdi32.GetCharWidth32(self._bitmap_dc, 
                ord(text), ord(text), byref(width_buf))
            width = width_buf.value
            lsb = 0
            advance = width

        # Can't get glyph-specific dimensions, use whole line-height.
        height = self._bitmap_rect.bottom

        # Draw to DC
        user32.FillRect(self._bitmap_dc, byref(self._bitmap_rect), self._black)
        gdi32.ExtTextOutA(self._bitmap_dc, -lsb, 0, 0, c_void_p(), text,
            len(text), c_void_p())
        gdi32.GdiFlush()

        # Create glyph object and copy bitmap data to texture
        image = pyglet.image.ImageData(width, height, 
            'RGBA', self._bitmap_data, self._bitmap_rect.right * 4)
        glyph = self.font.create_glyph(image)
        glyph.set_bearings(-self.font.descent, lsb, advance)

        return glyph
        
    def _create_bitmap_dc(self, width, height):
        if self._bitmap_dc:
            gdi32.ReleaseDC(self._bitmap_dc)
        if self._bitmap:
            gdi32.DeleteObject(self._bitmap)

        pitch = width * 4
        data = POINTER(c_byte * (height * pitch))()
        info = BITMAPINFO()
        info.bmiHeader.biSize = sizeof(info.bmiHeader)
        info.bmiHeader.biWidth = width
        info.bmiHeader.biHeight = height
        info.bmiHeader.biPlanes = 1
        info.bmiHeader.biBitCount = 32 
        info.bmiHeader.biCompression = BI_BITFIELDS
        info.bmiColors[:] = [0xff000000, 0x00ff0000, 0x0000ff00]

        self._bitmap_dc = gdi32.CreateCompatibleDC(c_void_p())
        self._bitmap = gdi32.CreateDIBSection(c_void_p(),
            byref(info), DIB_RGB_COLORS, byref(data), c_void_p(),
            0)
        # Spookiness: the above line causes a "not enough storage" error,
        # even though that error cannot be generated according to docs,
        # and everything works fine anyway.  Call SetLastError to clear it.
        kernel32.SetLastError(0)

        self._bitmap_data = data.contents
        self._bitmap_rect = RECT()
        self._bitmap_rect.left = 0
        self._bitmap_rect.right = width
        self._bitmap_rect.top = 0
        self._bitmap_rect.bottom = height

class Win32Font(base.Font):
    glyph_renderer_class = Win32GlyphRenderer

    def __init__(self, name, size, bold=False, italic=False, dpi=None):
        super(Win32Font, self).__init__()

        self.hfont = self.get_hfont(name, size, bold, italic, dpi)

        # Create a dummy DC for coordinate mapping
        dc = user32.GetDC(0)
        metrics = TEXTMETRIC()
        gdi32.SelectObject(dc, self.hfont)
        gdi32.GetTextMetricsA(dc, byref(metrics))
        self.ascent = metrics.tmAscent
        self.descent = -metrics.tmDescent
        self.max_glyph_width = metrics.tmMaxCharWidth

    @staticmethod
    def get_hfont(name, size, bold, italic, dpi):
        # Create a dummy DC for coordinate mapping
        dc = user32.GetDC(0)
        if dpi is None:
            logpixelsy = gdi32.GetDeviceCaps(dc, LOGPIXELSY)
        else:
            logpixelsy = dpi

        logfont = LOGFONT()
        # Conversion of point size to device pixels
        logfont.lfHeight = int(-size * logpixelsy / 72)
        if bold:
            logfont.lfWeight = FW_BOLD
        else:
            logfont.lfWeight = FW_NORMAL
        logfont.lfItalic = italic
        logfont.lfFaceName = name
        logfont.lfQuality = ANTIALIASED_QUALITY
        hfont = gdi32.CreateFontIndirectA(byref(logfont))
        return hfont

    @classmethod
    def have_font(cls, name):
        # CreateFontIndirect always returns a font... have to work out
        # something with EnumFontFamily... TODO
        return True

    @classmethod
    def add_font_data(cls, data):
        numfonts = c_uint32()
        gdi32.AddFontMemResourceEx(data, len(data), 0, byref(numfonts))
