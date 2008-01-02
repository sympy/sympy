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

import sys

from ctypes import *

from pyglet.gl import *
from pyglet.image import *
from pyglet.image.codecs import *

from pyglet.window.carbon import carbon, quicktime
from pyglet.window.carbon.constants import _name
from pyglet.window.carbon.types import *

Handle = POINTER(POINTER(c_byte))
GWorldPtr = c_void_p
carbon.NewHandle.restype = Handle
HandleDataHandlerSubType = _name('hndl')
ComponentInstance = c_void_p

k1MonochromePixelFormat       = 0x00000001
k2IndexedPixelFormat          = 0x00000002
k4IndexedPixelFormat          = 0x00000004
k8IndexedPixelFormat          = 0x00000008
k16BE555PixelFormat           = 0x00000010
k24RGBPixelFormat             = 0x00000018
k32ARGBPixelFormat            = 0x00000020
k32BGRAPixelFormat            = _name('BGRA')
k1IndexedGrayPixelFormat      = 0x00000021
k2IndexedGrayPixelFormat      = 0x00000022
k4IndexedGrayPixelFormat      = 0x00000024
k8IndexedGrayPixelFormat      = 0x00000028
kNativeEndianPixMap           = 1 << 8

class QuickTimeImageDecoder(ImageDecoder):
    def get_file_extensions(self):
        # Only most common ones shown here
        return ['.bmp', '.cur', '.gif', '.ico', '.jpg', '.jpeg', '.pcx', '.png',
                '.tga', '.tif', '.tiff', '.xbm', '.xpm']

    def decode(self, file, filename):
        data = file.read()
        handle = Handle()
        dataref = Handle()
        carbon.PtrToHand(data, byref(handle), len(data))
        carbon.PtrToHand(byref(handle), byref(dataref), sizeof(Handle))
        importer = ComponentInstance()
        quicktime.GetGraphicsImporterForDataRef(dataref, 
            HandleDataHandlerSubType, byref(importer))

        rect = Rect()
        quicktime.GraphicsImportGetNaturalBounds(importer, byref(rect))
        width = rect.right
        height = rect.bottom

        # TODO choose 24 bit where appropriate.
        if sys.byteorder == 'big':
            format = 'ARGB'
            qtformat = k32ARGBPixelFormat
        else:
            format = 'BGRA'
            qtformat = k32BGRAPixelFormat

        buffer = (c_byte * (width * height * len(format)))()
        world = GWorldPtr()
        quicktime.QTNewGWorldFromPtr(byref(world), qtformat,
            byref(rect), c_void_p(), c_void_p(), 0, buffer,
            len(format) * width)

        quicktime.GraphicsImportSetGWorld(importer, world, c_void_p())
        quicktime.GraphicsImportDraw(importer)
        quicktime.DisposeGWorld(world)

        pitch = len(format) * width

        return ImageData(width, height, format, buffer, -pitch)
        
def get_decoders():
    return [QuickTimeImageDecoder()]

def get_encoders():
    return []
