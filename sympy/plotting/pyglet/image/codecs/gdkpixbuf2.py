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
__version__ = '$Id: gdkpixbuf2.py 1322 2007-10-23 12:58:03Z Alex.Holkner $'

from ctypes import *

from pyglet.gl import *
from pyglet.image import *
from pyglet.image.codecs import *

import pyglet.lib
import pyglet.window

gdk = pyglet.lib.load_library('gdk-x11-2.0')
gdkpixbuf = pyglet.lib.load_library('gdk_pixbuf-2.0')

GdkPixbufLoader = c_void_p
GdkPixbuf = c_void_p
gdkpixbuf.gdk_pixbuf_loader_new.restype = GdkPixbufLoader
gdkpixbuf.gdk_pixbuf_loader_get_pixbuf.restype = GdkPixbuf
gdkpixbuf.gdk_pixbuf_get_pixels.restype = c_void_p

class GdkPixbuf2ImageDecoder(ImageDecoder):
    def get_file_extensions(self):
        return ['.png', '.xpm', '.jpg', '.jpeg', '.tif', '.tiff', '.pnm',
                '.ras', '.bmp', '.gif']

    def decode(self, file, filename):
        data = file.read()

        # Load into pixbuf
        err = c_int()
        loader = gdkpixbuf.gdk_pixbuf_loader_new()
        gdkpixbuf.gdk_pixbuf_loader_write(loader, data, len(data), byref(err))
        pixbuf = gdkpixbuf.gdk_pixbuf_loader_get_pixbuf(loader)
        if not gdkpixbuf.gdk_pixbuf_loader_close(loader, byref(err)):
            raise ImageDecodeException(filename)
        if not pixbuf:
            raise ImageDecodeException('Unable to load pixbuf: %s' % filename)
        
        # Get format and dimensions
        width = gdkpixbuf.gdk_pixbuf_get_width(pixbuf)
        height = gdkpixbuf.gdk_pixbuf_get_height(pixbuf)
        channels = gdkpixbuf.gdk_pixbuf_get_n_channels(pixbuf)
        rowstride = gdkpixbuf.gdk_pixbuf_get_rowstride(pixbuf)
        has_alpha = gdkpixbuf.gdk_pixbuf_get_has_alpha(pixbuf)
        pixels = gdkpixbuf.gdk_pixbuf_get_pixels(pixbuf)

        # Copy pixel data.
        buffer = (c_ubyte * (rowstride * height))()
        memmove(buffer, pixels, rowstride * (height - 1) + width * channels)

        # Release pixbuf
        gdk.g_object_unref(pixbuf)

        # Determine appropriate GL type
        if channels == 3:
            format = 'RGB'
        else:
            format = 'RGBA'

        return ImageData(width, height, format, buffer, -rowstride)

def get_decoders():
    return [GdkPixbuf2ImageDecoder()]

def get_encoders():
    return []

def init():
    gdk.g_type_init()

init()

