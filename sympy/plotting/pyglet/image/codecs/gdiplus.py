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
__version__ = '$Id: pil.py 163 2006-11-13 04:15:46Z Alex.Holkner $'

from ctypes import *

from pyglet.gl import *
from pyglet.image import *
from pyglet.image.codecs import *
from pyglet.window.win32.constants import *
from pyglet.window.win32.types import *

ole32 = windll.ole32
kernel32 = windll.kernel32
gdiplus = windll.gdiplus
mapi32 = windll.mapi32

LPSTREAM = c_void_p
REAL = c_float

PixelFormat1bppIndexed    = 196865
PixelFormat4bppIndexed    = 197634
PixelFormat8bppIndexed    = 198659
PixelFormat16bppGrayScale = 1052676
PixelFormat16bppRGB555    = 135173
PixelFormat16bppRGB565    = 135174
PixelFormat16bppARGB1555  = 397319
PixelFormat24bppRGB       = 137224
PixelFormat32bppRGB       = 139273
PixelFormat32bppARGB      = 2498570
PixelFormat32bppPARGB     = 925707
PixelFormat48bppRGB       = 1060876
PixelFormat64bppARGB      = 3424269
PixelFormat64bppPARGB     = 29622286
PixelFormatMax            = 15

ImageLockModeRead = 1
ImageLockModeWrite = 2
ImageLockModeUserInputBuf = 4

class GdiplusStartupInput(Structure):
    _fields_ = [
        ('GdiplusVersion', c_uint32),
        ('DebugEventCallback', c_void_p),
        ('SuppressBackgroundThread', BOOL),
        ('SuppressExternalCodecs', BOOL)
    ]

class GdiplusStartupOutput(Structure):
    _fields = [
        ('NotificationHookProc', c_void_p),
        ('NotificationUnhookProc', c_void_p)
    ]

class BitmapData(Structure):
    _fields_ = [
        ('Width', c_uint),
        ('Height', c_uint),
        ('Stride', c_int),
        ('PixelFormat', c_int),
        ('Scan0', POINTER(c_byte)),
        ('Reserved', POINTER(c_uint))
    ]

class Rect(Structure):
    _fields_ = [
        ('X', c_int),
        ('Y', c_int),
        ('Width', c_int),
        ('Height', c_int)
    ]

kernel32.GlobalAlloc.restype = HGLOBAL
kernel32.GlobalLock.restype = c_void_p

class GDIPlusDecoder(ImageDecoder):
    def get_file_extensions(self):
        return ['.bmp', '.gif', '.jpg', '.jpeg', '.exif', '.png', '.tif', 
                '.tiff']

    def decode(self, file, filename):
        data = file.read()

        # Create a HGLOBAL with image data
        hglob = kernel32.GlobalAlloc(GMEM_MOVEABLE, len(data))
        ptr = kernel32.GlobalLock(hglob)
        memmove(ptr, data, len(data))
        kernel32.GlobalUnlock(hglob)

        # Create IStream for the HGLOBAL
        stream = LPSTREAM()
        ole32.CreateStreamOnHGlobal(hglob, True, byref(stream))

        # Load image from stream
        bitmap = c_void_p()
        status = gdiplus.GdipCreateBitmapFromStream(stream, byref(bitmap))
        if status != 0:
            # TODO release stream
            raise ImageDecodeException(
                'GDI+ cannot load %r' % (filename or file))

        # Get size of image (Bitmap subclasses Image)
        width = REAL()
        height = REAL()
        gdiplus.GdipGetImageDimension(bitmap, byref(width), byref(height))
        width = int(width.value)
        height = int(height.value)

        # Get image pixel format
        pf = c_int()
        gdiplus.GdipGetImagePixelFormat(bitmap, byref(pf))
        pf = pf.value

        # Reverse from what's documented because of Intel little-endianness.
        format = 'BGRA'
        if pf == PixelFormat24bppRGB:
            format = 'BGR'
        elif pf == PixelFormat32bppRGB:
            pass
        elif pf == PixelFormat32bppARGB:
            pass
        elif pf in (PixelFormat16bppARGB1555, PixelFormat32bppPARGB,
                    PixelFormat64bppARGB, PixelFormat64bppPARGB):
            pf = PixelFormat32bppARGB
        else:
            format = 'BGR'
            pf = PixelFormat24bppRGB

        # Lock pixel data in best format
        rect = Rect()
        rect.X = 0
        rect.Y = 0
        rect.Width = width
        rect.Height = height
        bitmap_data = BitmapData()
        gdiplus.GdipBitmapLockBits(bitmap, 
            byref(rect), ImageLockModeRead, pf, byref(bitmap_data))
        
        # Create buffer for RawImage
        buffer = create_string_buffer(bitmap_data.Stride * height)
        memmove(buffer, bitmap_data.Scan0, len(buffer))
        
        # Unlock data
        gdiplus.GdipBitmapUnlockBits(bitmap, byref(bitmap_data))

        # Release image and stream
        gdiplus.GdipDisposeImage(bitmap)
        # TODO: How to call IUnknown::Release on stream?

        return ImageData(width, height, format, buffer, -bitmap_data.Stride)

def get_decoders():
    return [GDIPlusDecoder()]

def get_encoders():
    return []

def init():
    token = c_ulong()
    startup_in = GdiplusStartupInput()
    startup_in.GdiplusVersion = 1
    startup_out = GdiplusStartupOutput()
    gdiplus.GdiplusStartup(byref(token), byref(startup_in), byref(startup_out))

    # Shutdown later?
    # gdiplus.GdiplusShutdown(token)

init()
