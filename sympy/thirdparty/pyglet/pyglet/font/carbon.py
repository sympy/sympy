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
__version__ = '$Id: $'

# TODO Tiger and later: need to set kWindowApplicationScaledAttribute for DPI
# independence?

from ctypes import *
import string
import math
from sys import byteorder

from pyglet.font import base
import pyglet.image
from pyglet.window.carbon import carbon, _oscheck
from pyglet.window.carbon import _create_cfstring
from pyglet.window.carbon.types import *


class FixedPoint(Structure):
    _fields_ = [
        ('x', Fixed),
        ('y', Fixed)
    ]

class ATSTrapezoid(Structure):
    _fields_ = [
        ('upperLeft', FixedPoint),
        ('upperRight', FixedPoint),
        ('lowerRight', FixedPoint),
        ('lowerLeft', FixedPoint)
    ]

# TODO: most of the ATS and CG here not used any more.

CGGlyph = c_ushort
ATSUFontID = c_uint32
RGBColor = c_short * 3
ATSURGBAlphaColor = c_float * 4

kCGImageAlphaNone = 0
kCGImageAlphaPremultipliedLast = 1
kCGTextFill = 0

kATSUInvalidFontErr = -8796

kATSFontContextUnspecified = 0
kATSFontContextGlobal = 1
kATSFontContextLocal = 2

kATSFontFilterSelectorUnspecified = 0
kATSFontFilterSelectorGeneration = 3
kATSFontFilterSelectorFontFamily = 7
kATSFontFilterSelectorFontFamilyApplierFunction = 8
kATSFontFilterSelectorFontApplierFunction = 9

kATSOptionFlagsDoNotNotify = 0x00000001 << 8
kATSOptionFlagsIterationScopeMask = 0x00000007 << 12
kATSOptionFlagsDefaultScope = 0x00000000 << 12
kATSOptionFlagsUnRestrictedScope = 0x00000001 << 12
kATSOptionFlagsRestrictedScope = 0x00000002 << 12
kATSOptionFlagsProcessSubdirectories = 0x00000001 << 6

kATSUFromTextBeginning = c_ulong(0xFFFFFFFF)
kATSUToTextEnd = c_ulong(0xFFFFFFFF)

kATSULineAscentTag = 8
kATSULineDescentTag = 9
ATSUTextMeasurement = Fixed

kATSUQDBoldfaceTag            = 256
kATSUQDItalicTag              = 257
kATSUFontTag                  = 261
kATSUSizeTag                  = 262
kATSUCGContextTag             = 32767
kATSUColorTag                 = 263
kATSURGBAlphaColorTag         = 288

kATSULineWidthTag             = 1

kFontFullName                 = 4
kFontNoPlatformCode           = c_ulong(-1)
kFontNoScriptCode             = c_ulong(-1)
kFontNoLanguageCode           = c_ulong(-1)

kATSUseDeviceOrigins          = 1

kATSFontFormatUnspecified     = 0
kATSFontContextLocal          = 2

carbon.CGColorSpaceCreateWithName.restype = c_void_p
carbon.CGBitmapContextCreate.restype = POINTER(c_void_p)

UniCharArrayOffset  = c_uint32
UniCharCount = c_uint32

def fixed(value):
    # This is a guess... could easily be wrong
    #return c_int32(int(value) * (1 << 16))

    return c_int32(carbon.Long2Fix(c_long(int(value))))

carbon.Fix2X.restype = c_double
def fix2float(value):
    return carbon.Fix2X(value)

def create_atsu_style(attributes):
    # attributes is a dict of ATSUAttributeTag => ctypes value
    tags, values = zip(*attributes.items())
    tags = (c_int * len(tags))(*tags)
    sizes = (c_uint * len(values))(*[sizeof(v) for v in values])
    values = (c_void_p * len(values))(*[cast(pointer(v), c_void_p) \
                                        for v in values])

    style = c_void_p()
    carbon.ATSUCreateStyle(byref(style))
    carbon.ATSUSetAttributes(style, len(tags), tags, sizes, values)
    return style

def set_layout_attributes(layout, attributes):
    if attributes:
        # attributes is a dict of ATSUAttributeTag => ctypes value
        tags, values = zip(*attributes.items())
        tags = (c_int * len(tags))(*tags)
        sizes = (c_uint * len(values))(*[sizeof(v) for v in values])
        values = (c_void_p * len(values))(*[cast(pointer(v), c_void_p) \
                                            for v in values])

        carbon.ATSUSetLayoutControls(layout, len(tags), tags, sizes, values)

def str_ucs2(text):
    if byteorder == 'big':
        text = text.encode('utf_16_be')
    else:
        text = text.encode('utf_16_le')   # explicit endian avoids BOM
    return create_string_buffer(text + '\0')

class CarbonGlyphRenderer(base.GlyphRenderer):
    _bitmap = None
    _bitmap_context = None
    _bitmap_rect = None

    def __init__(self, font):
        super(CarbonGlyphRenderer, self).__init__(font)
        self._create_bitmap_context(256, 256)
        self.font = font

    def __del__(self):
        if self._bitmap_context:
            carbon.CGContextRelease(self._bitmap_context)

    def render(self, text):
        # Convert text to UCS2
        text_len = len(text)
        text = str_ucs2(text)

        # Create ATSU text layout for this text and font
        layout = c_void_p()
        carbon.ATSUCreateTextLayout(byref(layout))
        set_layout_attributes(layout, {
            kATSUCGContextTag: self._bitmap_context})
        carbon.ATSUSetTextPointerLocation(layout,
            text,
            kATSUFromTextBeginning,
            kATSUToTextEnd,
            text_len)
        carbon.ATSUSetRunStyle(layout, self.font.atsu_style, 
            kATSUFromTextBeginning, kATSUToTextEnd)

        # Turning on transient font matching screws up font layout
        # predictability when strange fonts are installed
        carbon.ATSUSetTransientFontMatching(layout, False)

        # Get bitmap dimensions required
        rect = Rect()
        carbon.ATSUMeasureTextImage(layout, 
            kATSUFromTextBeginning,
            kATSUToTextEnd,
            0, 0,
            byref(rect))
        image_width = rect.right - rect.left + 2
        image_height = rect.bottom - rect.top + 2
        baseline = rect.bottom + 1
        lsb = rect.left
        
        # Resize Quartz context if necessary
        if (image_width > self._bitmap_rect.size.width or
            image_height > self._bitmap_rect.size.height):
            self._create_bitmap_context(
                int(max(image_width, self._bitmap_rect.size.width)),
                int(max(image_height, self._bitmap_rect.size.height)))
            
            set_layout_attributes(layout, {
                kATSUCGContextTag: self._bitmap_context})

        # Get typographic box, which gives advance.
        bounds_actual = c_uint32()
        bounds = ATSTrapezoid()
        carbon.ATSUGetGlyphBounds(
            layout,
            0, 0,
            kATSUFromTextBeginning,
            kATSUToTextEnd,
            kATSUseDeviceOrigins,
            1,
            byref(bounds),
            byref(bounds_actual))
        advance = fix2float(bounds.lowerRight.x) - fix2float(bounds.lowerLeft.x)

        # Draw to the bitmap
        carbon.CGContextClearRect(self._bitmap_context, self._bitmap_rect)
        carbon.ATSUDrawText(layout,
            0,
            kATSUToTextEnd,
            fixed(-lsb + 1), fixed(baseline)) 
        
        # A negative pitch is required, but it is much faster to load the
        # glyph upside-down and flip the tex_coords.  Note region used
        # to start at top of glyph image.
        pitch = int(4 * self._bitmap_rect.size.width)
        image = pyglet.image.ImageData(image_width, 
            self._bitmap_rect.size.height, 'RGBA', self._bitmap, pitch)
        skip_rows = int(self._bitmap_rect.size.height - image_height)
        image = image.get_region(0, skip_rows, image.width, image_height)
        glyph = self.font.create_glyph(image)
        glyph.set_bearings(baseline, lsb - 1, int(advance))
        t = list(glyph.tex_coords)
        glyph.tex_coords = t[9:12] + t[6:9] + t[3:6] + t[:3]
        
        return glyph

    def _create_bitmap_context(self, width, height):
        '''Create or recreate bitmap and Quartz context.'''
        if self._bitmap_context:
            carbon.CGContextRelease(self._bitmap_context)
        components = 4
        pitch = width * components
        self._bitmap = (c_ubyte * (pitch * height))()
        color_space = carbon.CGColorSpaceCreateDeviceRGB()
        context = carbon.CGBitmapContextCreate(self._bitmap, 
            width, height, 8, pitch, 
            color_space, kCGImageAlphaPremultipliedLast)
        carbon.CGColorSpaceRelease(color_space)

        # Disable RGB decimated antialiasing, use standard
        # antialiasing which won't break alpha.
        carbon.CGContextSetShouldSmoothFonts(context, False)
        carbon.CGContextSetShouldAntialias(context, True)

        self._bitmap_context = context 
        self._bitmap_rect = CGRect()
        self._bitmap_rect.origin.x = 0
        self._bitmap_rect.origin.y = 0
        self._bitmap_rect.size.width = width
        self._bitmap_rect.size.height = height
        

class CarbonFont(base.Font):
    glyph_renderer_class = CarbonGlyphRenderer

    def __init__(self, name, size, bold=False, italic=False, dpi=None):
        super(CarbonFont, self).__init__()

        if not name:
            name = 'Helvetica'

        if dpi is not None:
            # If application is not DPI-aware, DPI is fixed at 72.  Scale
            # font size to emulate other DPI if necessary.  This will need
            # to be fixed if issue #87 is implemented.
            size = size * dpi / 72.

        font_id = ATSUFontID()
        carbon.ATSUFindFontFromName(
            name,
            len(name),
            kFontFullName,
            kFontNoPlatformCode,
            kFontNoScriptCode,
            kFontNoLanguageCode,
            byref(font_id))

        attributes = {
            kATSUSizeTag: fixed(size),
            kATSUFontTag: font_id,
            kATSURGBAlphaColorTag: ATSURGBAlphaColor(1, 1, 1, 1),
            kATSUQDBoldfaceTag: c_byte(bold),
            kATSUQDItalicTag: c_byte(italic)
        }
        self.atsu_style = create_atsu_style(attributes)

        self.calculate_metrics()

    @classmethod
    def have_font(cls, name):
        font_id = ATSUFontID()
        r = carbon.ATSUFindFontFromName(
            name,
            len(name),
            kFontFullName,
            kFontNoPlatformCode,
            kFontNoScriptCode,
            kFontNoLanguageCode,
            byref(font_id)) 
        return r != kATSUInvalidFontErr

    def calculate_metrics(self):
        # It seems the only way to get the font's ascent and descent is to lay
        # out some glyphs and measure them.

        # fake ucs2 string
        text = '\0a'

        layout = c_void_p()
        carbon.ATSUCreateTextLayout(byref(layout))
        carbon.ATSUSetTextPointerLocation(layout, text,
            kATSUFromTextBeginning, kATSUToTextEnd, 1)
        carbon.ATSUSetRunStyle(layout, self.atsu_style, 
            kATSUFromTextBeginning, kATSUToTextEnd)

        # determine the metrics for this font only
        carbon.ATSUSetTransientFontMatching(layout, False)

        value = ATSUTextMeasurement()
        carbon.ATSUGetLineControl(layout, 0, kATSULineAscentTag, 
            sizeof(value), byref(value), None)
        self.ascent = int(math.ceil(fix2float(value)))
        carbon.ATSUGetLineControl(layout, 0, kATSULineDescentTag,
            sizeof(value), byref(value), None)
        self.descent = -int(math.ceil(fix2float(value)))

    @classmethod
    def add_font_data(cls, data):
        container = c_void_p()
        r = carbon.ATSFontActivateFromMemory(data, len(data),
            kATSFontContextLocal, kATSFontFormatUnspecified, None, 0,
            byref(container))
        _oscheck(r)
