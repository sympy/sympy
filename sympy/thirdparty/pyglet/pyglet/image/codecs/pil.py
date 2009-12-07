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
__version__ = '$Id: pil.py 1768 2008-02-17 12:45:54Z Alex.Holkner $'

import os.path

from pyglet.gl import *
from pyglet.image import *
from pyglet.image.codecs import *

import Image

class PILImageDecoder(ImageDecoder):
    def get_file_extensions(self):
        # Only most common ones shown here
        return ['.bmp', '.cur', '.gif', '.ico', '.jpg', '.jpeg', '.pcx', '.png',
                '.tga', '.tif', '.tiff', '.xbm', '.xpm']

    def decode(self, file, filename):
        try:
            image = Image.open(file)
        except Exception, e:
            raise ImageDecodeException(
                'PIL cannot read %r: %s' % (filename or file, e))

        image = image.transpose(Image.FLIP_TOP_BOTTOM)

        # Convert bitmap and palette images to component
        if image.mode in ('1', 'P'):
            image = image.convert()

        if image.mode not in ('L', 'LA', 'RGB', 'RGBA'):
            raise ImageDecodeException('Unsupported mode "%s"' % image.mode)
        type = GL_UNSIGNED_BYTE
        width, height = image.size

        return ImageData(width, height, image.mode, image.tostring())

class PILImageEncoder(ImageEncoder):
    def get_file_extensions(self):
        # Most common only
        return ['.bmp', '.eps', '.gif', '.jpg', '.jpeg',
                '.pcx', '.png', '.ppm', '.tiff', '.xbm']

    def encode(self, image, file, filename):
        # File format is guessed from filename extension, otherwise defaults
        # to PNG.
        pil_format = (filename and os.path.splitext(filename)[1][1:]) or 'png'

        if pil_format.lower() == 'jpg':
            pil_format = 'JPEG'

        image = image.get_image_data()
        format = image.format
        if format != 'RGB':
            # Only save in RGB or RGBA formats.
            format = 'RGBA'
        pitch = -(image.width * len(format))

        # Note: Don't try and use frombuffer(..); different versions of
        # PIL will orient the image differently.
        pil_image = Image.fromstring(
            format, (image.width, image.height), image.get_data(format, pitch))

        try:
            pil_image.save(file, pil_format)
        except Exception, e:
            raise ImageEncodeException(e)

def get_decoders():
    return [PILImageDecoder()]

def get_encoders():
    return [PILImageEncoder()]
