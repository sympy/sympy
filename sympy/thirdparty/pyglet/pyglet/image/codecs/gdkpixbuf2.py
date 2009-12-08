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
__version__ = '$Id: gdkpixbuf2.py 1939 2008-03-21 11:48:12Z Alex.Holkner $'

from ctypes import *

from pyglet.gl import *
from pyglet.image import *
from pyglet.image.codecs import *
from pyglet.image.codecs import gif

import pyglet.lib
import pyglet.window

gdk = pyglet.lib.load_library('gdk-x11-2.0')
gdkpixbuf = pyglet.lib.load_library('gdk_pixbuf-2.0')

GdkPixbufLoader = c_void_p
GdkPixbuf = c_void_p
gdkpixbuf.gdk_pixbuf_loader_new.restype = GdkPixbufLoader
gdkpixbuf.gdk_pixbuf_loader_get_pixbuf.restype = GdkPixbuf
gdkpixbuf.gdk_pixbuf_get_pixels.restype = c_void_p
gdkpixbuf.gdk_pixbuf_loader_get_animation.restype = c_void_p
gdkpixbuf.gdk_pixbuf_animation_get_iter.restype = c_void_p
gdkpixbuf.gdk_pixbuf_animation_iter_get_pixbuf.restype = GdkPixbuf

class GTimeVal(Structure):
    _fields_ = [
        ('tv_sec', c_long),
        ('tv_usec', c_long)
    ]

class GdkPixbuf2ImageDecoder(ImageDecoder):
    def get_file_extensions(self):
        return ['.png', '.xpm', '.jpg', '.jpeg', '.tif', '.tiff', '.pnm',
                '.ras', '.bmp', '.gif']

    def get_animation_file_extensions(self):
        return ['.gif', '.ani']

    def _load(self, file, filename, load_func):
        data = file.read()
        err = c_int()
        loader = gdkpixbuf.gdk_pixbuf_loader_new()
        gdkpixbuf.gdk_pixbuf_loader_write(loader, data, len(data), byref(err))
        result = load_func(loader)
        if not gdkpixbuf.gdk_pixbuf_loader_close(loader, byref(err)):
            raise ImageDecodeException(filename)
        if not result:
            raise ImageDecodeException('Unable to load: %s' % filename)
        return result

    def _pixbuf_to_image(self, pixbuf):
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

    def decode(self, file, filename):
        pixbuf = self._load(file, filename, 
                            gdkpixbuf.gdk_pixbuf_loader_get_pixbuf)
       
        return self._pixbuf_to_image(pixbuf)

    def decode_animation(self, file, filename):
        # Extract GIF control data.  If it's not a GIF, this method will
        # raise.
        gif_stream = gif.read(file)
        delays = [image.delay for image in gif_stream.images]

        # Get GDK animation iterator
        file.seek(0)
        anim = self._load(file, filename, 
                          gdkpixbuf.gdk_pixbuf_loader_get_animation)
        time = GTimeVal(0, 0)
        iter = gdkpixbuf.gdk_pixbuf_animation_get_iter(anim, byref(time))

        frames = []

        # Extract each image   
        for control_delay in delays:
            pixbuf = gdkpixbuf.gdk_pixbuf_animation_iter_get_pixbuf(iter)
            image = self._pixbuf_to_image(pixbuf)
            frames.append(AnimationFrame(image, control_delay))

            gdk_delay = gdkpixbuf.gdk_pixbuf_animation_iter_get_delay_time(iter)
            gdk_delay *= 1000 # milliseconds to microseconds
            # Compare gdk_delay to control_delay for interest only.
            #print control_delay, gdk_delay / 1000000.

            if gdk_delay == -1:
                break

            us = time.tv_usec + gdk_delay
            time.tv_sec += us // 1000000
            time.tv_usec = us % 1000000
            gdkpixbuf.gdk_pixbuf_animation_iter_advance(iter, byref(time))

        return Animation(frames)

def get_decoders():
    return [GdkPixbuf2ImageDecoder()]

def get_encoders():
    return []

def init():
    gdk.g_type_init()

init()

