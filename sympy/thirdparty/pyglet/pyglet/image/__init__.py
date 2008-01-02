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

'''Image load, capture and high-level texture functions.

Only basic functionality is described here; for full reference see the
accompanying documentation.

To load an image::

    from pyglet import image
    pic = image.load('picture.png')

The supported image file types include PNG, BMP, GIF, JPG, and many more,
somewhat depending on the operating system.  To load an image from a file-like
object instead of a filename::

    pic = image.load('hint.jpg', file=fileobj)

The hint helps the module locate an appropriate decoder to use based on the
file extension.  It is optional.

Once loaded, images can be used directly by most other modules of pyglet (for
example, pyglet.ext.scene2d).  All images have a width and height you can
access::

    width, height = pic.width, pic.height

You can extract a region of an image (this keeps the original image intact;
the memory is shared efficiently)::

    subimage = pic.get_region(x, y, width, height)

Remember that y-coordinates are always increasing upwards.

Drawing images
--------------

To draw an image at some point on the screen::

    pic.blit(x, y, z)

This assumes an appropriate view transform and projection have been applied.

Texture access
--------------

If you are using OpenGL directly, you can access the image as a texture::

    texture = pic.texture

(This is the most efficient way to obtain a texture; some images are
immediately loaded as textures, whereas others go through an intermediate
form).  To use a texture with pyglet.gl::

    from pyglet.gl import *
    glEnable(texture.target)        # typically target is GL_TEXTURE_2D
    glBindTexture(texture.target, texture.id)
    # ... draw with the texture

Pixel access
------------

To access raw pixel data of an image::

    rawimage = pic.image_data

(If the image has just been loaded this will be a very quick operation;
however if the image is a texture a relatively expensive readback operation
will occur).  The pixels can be accessed as a string::

    pixels = rawimage.data

You determine the format of these pixels by examining rawimage.pitch,
rawimage.format.  The "pitch" of an image is the number of bytes in a row
(this may validly be more than the number required to make up the width of the
image, it is common to see this for word alignment).  If "pitch" is negative
the rows of the image are ordered from top to bottom, otherwise they are
ordedred from bottom to top.

"format" strings consist of characters that give the byte
order of each color component.  For example, if rawimage.format is 'RGBA',
there are four color components: red, green, blue and alpha, in that order.
Other common format strings are 'RGB', 'LA' (luminance, alpha) and 'I'
(intensity).

You can also set the format and pitch of an image.  This affects how "data"
will be presented the next time it is read or written to::

    # Read the alpha channel only, with rows tightly packed from bottom-to-top.
    rawimage.format = 'A'
    rawimage.pitch = rawimage.width
    data = rawimage.data

'''

__docformat__ = 'restructuredtext'
__version__ = '$Id: __init__.py 1502 2007-12-09 00:03:55Z Alex.Holkner $'

import sys
import re
import warnings
import weakref

from ctypes import *
from math import ceil
from StringIO import StringIO

from pyglet.gl import *
from pyglet.gl import gl_info
from pyglet.window import *

class ImageException(Exception):
    pass

def load(filename, file=None, decoder=None):
    '''Load an image from a file.

    :note: You can make no assumptions about the return type; usually it will
        be ImageData or CompressedImageData, but decoders are free to return
        any subclass of AbstractImage.

    :Parameters:
        `filename` : str
            Used to guess the image format, and to load the file if `file` is
            unspecified.
        `file` : file-like object or None
            Source of image data in any supported format.        
        `decoder` : ImageDecoder or None
            If unspecified, all decoders that are registered for the filename
            extension are tried.  If none succeed, the exception from the
            first decoder is raised.

    :rtype: AbstractImage
    '''

    if not file:
        file = open(filename, 'rb')
    if not hasattr(file, 'seek'):
        file = StringIO(file.read())

    if decoder:
        return decoder.decode(file, filename)
    else:
        first_exception = None
        for decoder in codecs.get_decoders(filename):
            try:
                image = decoder.decode(file, filename)
                return image
            except codecs.ImageDecodeException, e:
                first_exception = first_exception or e
                file.seek(0)

        if not first_exception:
            raise codecs.ImageDecodeException('No image decoders are available')
        raise first_exception 

def create(width, height, pattern=None):
    '''Create an image optionally filled with the given pattern.

    :note: You can make no assumptions about the return type; usually it will
        be ImageData or CompressedImageData, but patterns are free to return
        any subclass of AbstractImage.

    :Parameters:
        `width` : int
            Width of image to create
        `height` : int
            Height of image to create
        `pattern` : ImagePattern or None
            Pattern to fill image with.  If unspecified, the image will
            intially be transparent.

    :rtype: AbstractImage
    '''
    if not pattern:
        pattern = SolidColorImagePattern()
    return pattern.create_image(width, height)

class ImagePattern(object):
    '''Abstract image creation class.'''
    def create_image(width, height):
        '''Create an image of the given size.

        :Parameters:
            `width` : int
                Width of image to create
            `height` : int
                Height of image to create
        
        :rtype: AbstractImage
        '''
        raise NotImplementedError('abstract')

class SolidColorImagePattern(ImagePattern):
    '''Creates an image filled with a solid color.'''

    def __init__(self, color=(0, 0, 0, 0)):
        '''Create a solid image pattern with the given color.

        :Parameters:
            `color` : (int, int, int, int)
                4-tuple of ints in range [0,255] giving RGBA components of
                color to fill with.

        '''
        self.color = '%c%c%c%c' % color

    def create_image(self, width, height):
        data = self.color * width * height
        return ImageData(width, height, 'RGBA', data)

class CheckerImagePattern(ImagePattern):
    '''Create an image with a tileable checker image.
    ''' 

    def __init__(self, color1=(150,150,150,255), color2=(200,200,200,255)):
        '''Initialise with the given colors.

        :Parameters:
            `color1` : (int, int, int, int)
                4-tuple of ints in range [0,255] giving RGBA components of
                color to fill with.  This color appears in the top-left and
                bottom-right corners of the image.
            `color2` : (int, int, int, int)
                4-tuple of ints in range [0,255] giving RGBA components of
                color to fill with.  This color appears in the top-right and
                bottom-left corners of the image.

        '''
        self.color1 = '%c%c%c%c' % color1
        self.color2 = '%c%c%c%c' % color2

    def create_image(self, width, height):
        hw = width/2
        hh = height/2
        row1 = self.color1 * hw + self.color2 * hw
        row2 = self.color2 * hw + self.color1 * hw
        data = row1 * hh + row2 * hh
        return ImageData(width, height, 'RGBA', data)

class AbstractImage(object):
    '''Abstract class representing an image.

    :Ivariables:
        `width` : int
            Width of image
        `height` : int
            Height of image
    '''

    def __init__(self, width, height):
        self.width = width
        self.height = height

    def __repr__(self):
        return '<%s %dx%d>' % (self.__class__.__name__, self.width, self.height)

    def _get_image_data(self):
        '''Retrieve an ImageData instance for this image.'''
        raise ImageException('Cannot retrieve image data for %r' % self)

    image_data = property(_get_image_data,
        doc='''An ImageData view of this image.  
        
        Changes to the returned instance may or may not be reflected in this
        image.  Read-only.
        
        :type: `ImageData`
        ''')

    def _get_texture(self):
        '''Retrieve a Texture instance for this image.'''
        raise ImageException('Cannot retrieve texture for %r' % self)

    texture = property(_get_texture,
        doc=''' A Texture view of this image.  
        
        Changes to the returned instance may or may not be reflected in this
        image.  Read-only.

        :type: `Texture`
        ''')

    def _get_mipmapped_texture(self):
        '''Retrieve a Texture instance with all mipmap levels filled in.'''
        raise ImageException('Cannot retrieve mipmapped texture for %r' % self)

    mipmapped_texture = property(_get_mipmapped_texture,
        doc='''A Texture view of this image.  
        
        The returned Texture will have mipmaps filled in for all levels.
        Requires that image dimensions be powers of 2.  Read-only.

        :type: `Texture`
        ''')

    def get_region(self, x, y, width, height):
        '''Retrieve a rectangular region of this image.

        :Parameters:
            `x` : int
                Left edge of region.
            `y` : int
                Bottom edge of region.
            `width` : int
                Width of region.
            `height` : int
                Height of region.

        :rtype: AbstractImage
        '''
        raise ImageException('Cannot get region for %r' % self)

    def save(self, filename=None, file=None, encoder=None):
        '''Save this image to a file.

        :Parameters:
            `filename` : str
                Used to set the image file format, and to open the output file
                if `file` is unspecified.
            `file` : file-like object or None
                File to write image data to.
            `encoder` : ImageEncoder or None
                If unspecified, all encoders matching the filename extension
                are tried.  If all fail, the exception from the first one
                attempted is raised.

        '''
        if not file:
            file = open(filename, 'wb')

        if encoder:
            encoder.encode(self, file, filename)
        else:
            first_exception = None
            for encoder in codecs.get_encoders(filename):
                try:
                    encoder.encode(self, file, filename)
                    return
                except codecs.ImageDecodeException, e:
                    first_exception = first_exception or e
                    file.seek(0)

            if not first_exception:
                raise codecs.ImageEncodeException(
                        'No image encoders are available')
            raise first_exception

    def blit(self, x, y, z=0):
        '''Draw this image to the active framebuffers.'''
        raise ImageException('Cannot blit %r.' % self)

    def blit_into(self, source, x, y, z):
        '''Draw `source` on this image.
        
        Note that if `source` is larger than this image (or the positioning
        would cause the copy to go out of bounds) then you must pass a
        region of `source` to this method, typically using get_region().
        '''
        raise ImageException('Cannot blit images onto %r.' % self)

    def blit_to_texture(self, target, level, x, y, z=0):
        '''Draw this image on the currently bound texture at `target`.'''
        raise ImageException('Cannot blit %r to a texture.' % self)

class AbstractImageSequence(object):
    '''Abstract sequence of images.

    The sequence is useful for storing image animations or slices of a volume.
    For efficient access, use the `texture_sequence` member.  The class
    also implements the sequence interface (`__len__`, `__getitem__`,
    `__setitem__`).
    
    :Ivariables:
        `texture_sequence` : `TextureSequence`
            Access this image sequence as a texture sequence.
    '''
    def _get_texture_sequence(self):
        raise NotImplementedError('abstract')
    texture_sequence = property(_get_texture_sequence)

    def __getitem__(self, slice):
        '''Retrieve a (list of) image.
        
        :rtype: AbstractImage
        '''
        raise NotImplementedError('abstract')

    def __setitem__(self, slice, image):
        '''Replace one or more images in the sequence.
        
        :Parameters:
            `image` : `AbstractImage`
                The replacement image.  The actual instance may not be used,
                depending on this implementation.

        '''
        raise NotImplementedError('abstract')

    def __len__(self):
        raise NotImplementedError('abstract')


class TextureSequence(AbstractImageSequence):
    '''Interface for a sequence of textures.

    Typical implementations store multiple `TextureRegion` s within one
    `Texture` so as to minimise state changes.
    '''
    texture_sequence = property(lambda self: self)

class UniformTextureSequence(TextureSequence):
    '''Interface for a sequence of textures, each with the same dimensions.

    :Ivariables:
        `item_width` : int
            Width of each texture in the sequence.
        `item_height` : int
            Height of each texture in the sequence.
    
    '''
    def _get_item_width(self):
        raise NotImplementedError('abstract')
    item_width = property(_get_item_width)

    def _get_item_height(self):
        raise NotImplementedError('abstract')
    item_height = property(_get_item_height)

class ImageData(AbstractImage):
    '''An image represented as a string of unsigned bytes.

    :Ivariables:
        `data` : str
            Pixel data, encoded according to `format` and `pitch`.
        `format` : str
            The format string to use when reading or writing `data`.
        `pitch` : int
            Number of bytes per row.  Negative values indicate a top-to-bottom
            arrangement.

    '''

    _swap1_pattern = re.compile('(.)', re.DOTALL)
    _swap2_pattern = re.compile('(.)(.)', re.DOTALL)
    _swap3_pattern = re.compile('(.)(.)(.)', re.DOTALL)
    _swap4_pattern = re.compile('(.)(.)(.)(.)', re.DOTALL)

    _current_texture = None
    _current_mipmap_texture = None

    def __init__(self, width, height, format, data, pitch=None):
        '''Initialise image data.

        :Parameters:
            `width` : int
                Width of image data
            `height` : int
                Height of image data
            `format` : str
                A valid format string, such as 'RGB', 'RGBA', 'ARGB', etc.
            `data` : sequence
                String or array/list of bytes giving the decoded data.
            `pitch` : int or None
                If specified, the number of bytes per row.  Negative values
                indicate a top-to-bottom arrangement.  Defaults to 
                ``width * len(format)``.

        '''
        super(ImageData, self).__init__(width, height)

        self._current_format = self._desired_format = format.upper()
        self._current_data = data
        if not pitch:
            pitch = width * len(format)
        self._current_pitch = self.pitch = pitch
        self.mipmap_images = []

    def __getstate__(self):
        return {
            'width': self.width, 
            'height': self.height, 
            '_current_data': self.data, 
            '_current_format': self.format,
            '_desired_format': self.format,
            '_current_pitch': self.pitch,
            'pitch': self.pitch,
            'mipmap_images': self.mipmap_images
        }

    image_data = property(lambda self: self)

    def _set_format(self, format):
        self._desired_format = format.upper()
        self._current_texture = None

    format = property(lambda self: self._desired_format, _set_format,
        doc='''Format string of the data.  Read-write.
        
        :type: str
        ''')

    def _get_data(self):
        if self._current_pitch != self.pitch or \
           self._current_format != self.format:
            self._current_data = self._convert(self.format, self.pitch)
            self._current_format = self.format
            self._current_pitch = self.pitch

        self._ensure_string_data()
        return self._current_data

    def _set_data(self, data):
        self._current_data = data
        self._current_format = self.format
        self._current_pitch = self.pitch
        self._current_texture = None
        self._current_mipmapped_texture = None

    data = property(_get_data, _set_data, 
        doc='''The byte data of the image.  Read-write.
        
        :type: sequence of bytes, or str
        ''')

    def set_mipmap_image(self, level, image):
        '''Set a mipmap image for a particular level.

        The mipmap image will be applied to textures obtained via the
        `mipmapped_texture` property

        :Parameters:
            `level` : int
                Mipmap level to set image at, must be >= 1.
            `image` : AbstractImage
                Image to set.  Must have correct dimensions for that mipmap
                level (i.e., width >> level, height >> level)
        '''

        if level == 0:
            raise ImageException(
                'Cannot set mipmap image at level 0 (it is this image)')

        if not _is_pow2(self.width) or not _is_pow2(self.height):
            raise ImageException(
                'Image dimensions must be powers of 2 to use mipmaps.')

        # Check dimensions of mipmap
        width, height = self.width, self.height
        for i in range(level):
            width >>= 1
            height >>= 1
        if width != image.width or height != image.height:
            raise ImageException(
                'Mipmap image has wrong dimensions for level %d' % level)

        # Extend mipmap_images list to required level
        self.mipmap_images += [None] * (level - len(self.mipmap_images))
        self.mipmap_images[level - 1] = data

    def create_texture(self, cls):
        '''Create a texture containing this image.

        If the image's dimensions are not powers of 2, a TextureRegion of
        a larger Texture will be returned that matches the dimensions of this
        image.

        :Parameters:
            `cls` : class (subclass of Texture)
                Class to construct.

        :rtype: cls or cls.region_class
        '''

        texture = cls.create_for_size(
            GL_TEXTURE_2D, self.width, self.height)
        subimage = False
        if texture.width != self.width or texture.height != self.height:
            texture = texture.get_region(0, 0, self.width, self.height)
            subimage = True

        internalformat = self._get_internalformat(self.format)

        glBindTexture(texture.target, texture.id)
        glTexParameteri(texture.target, GL_TEXTURE_MIN_FILTER, GL_LINEAR)

        if subimage:
            width = texture.owner.width
            height = texture.owner.height
            blank = (c_ubyte * (width * height * 4))()
            glTexImage2D(texture.target, texture.level,
                         internalformat,
                         width, height,
                         0,
                         GL_RGBA, GL_UNSIGNED_BYTE,
                         blank) 
            internalformat = None

        self.blit_to_texture(texture.target, texture.level, 
            0, 0, 0, internalformat)
        
        return texture 

    def _get_texture(self):
        if not self._current_texture:
            self._current_texture = self.create_texture(Texture)
        return self._current_texture

    texture = property(_get_texture)

    def _get_mipmapped_texture(self):
        '''Return a Texture with mipmaps.  
        
        If `set_mipmap_image` has been called with at least one image, the set
        of images defined will be used.  Otherwise, mipmaps will be
        automatically generated.

        The texture dimensions must be powers of 2 to use mipmaps.
        '''
        if self._current_mipmap_texture:
            return self._current_mipmap_texture

        if not _is_pow2(self.width) or not _is_pow2(self.height):
            raise ImageException(
                'Image dimensions must be powers of 2 to use mipmaps.')
        
        texture = Texture.create_for_size(
            GL_TEXTURE_2D, self.width, self.height)
        internalformat = self._get_internalformat(self.format)

        glBindTexture(texture.target, texture.id)
        glTexParameteri(texture.target, GL_TEXTURE_MIN_FILTER,
                        GL_LINEAR_MIPMAP_LINEAR)

        if self.mipmap_images:
            self.blit_to_texture(texture.target, texture.level, 
                0, 0, 0, internalformat)
            level = 0
            for image in self.mipmap_images:
                level += 1
                if image:
                    image.blit_to_texture(texture.target, level, 
                        0, 0, 0, internalformat)
            # TODO: should set base and max mipmap level if some mipmaps
            # are missing.
        elif gl_info.have_version(1, 4):
            glTexParameteri(texture.target, GL_GENERATE_MIPMAP, GL_TRUE)
            self.blit_to_texture(texture.target, texture.level, 
                0, 0, 0, internalformat)
        else:
            raise NotImplementedError('TODO: gluBuild2DMipmaps')

        self._current_mipmap_texture = texture
        return texture

    mipmapped_texture = property(_get_mipmapped_texture)

    def get_region(self, x, y, width, height):
        '''Retrieve a rectangular region of this image data.

        :Parameters:
            `x` : int
                Left edge of region.
            `y` : int
                Bottom edge of region.
            `width` : int
                Width of region.
            `height` : int
                Height of region.

        :rtype: ImageDataRegion
        '''
        return ImageDataRegion(x, y, width, height, self)

    def blit(self, x, y, z=0, width=None, height=None):
        self.texture.blit(x, y, z, width, height)

    def blit_to_texture(self, target, level, x, y, z, internalformat=None):
        '''Draw this image to to the currently bound texture at `target`.

        If `internalformat` is specified, glTexImage is used to initialise
        the texture; otherwise, glTexSubImage is used to update a region.
        '''

        data_format = self.format
        data_pitch = abs(self._current_pitch)

        # Determine pixel format from format string
        matrix = None
        format, type = self._get_gl_format_and_type(data_format)
        if format is None:
            if (len(data_format) in (3, 4) and 
                gl_info.have_extension('GL_ARB_imaging')):
                # Construct a color matrix to convert to GL_RGBA
                def component_column(component):
                    try:
                        pos = 'RGBA'.index(component)
                        return [0] * pos + [1] + [0] * (3 - pos)
                    except ValueError:
                        return [0, 0, 0, 0]
                # pad to avoid index exceptions
                lookup_format = data_format + 'XXX'
                matrix = (component_column(lookup_format[0]) +
                          component_column(lookup_format[1]) +
                          component_column(lookup_format[2]) + 
                          component_column(lookup_format[3]))
                format = {
                    3: GL_RGB,
                    4: GL_RGBA}.get(len(data_format))
                type = GL_UNSIGNED_BYTE

                glMatrixMode(GL_COLOR)
                glPushMatrix()
                glLoadMatrixf((GLfloat * 16)(*matrix))
            else:
                # Need to convert data to a standard form
                data_format = {
                    1: 'L',
                    2: 'LA',
                    3: 'RGB',
                    4: 'RGBA'}.get(len(data_format))
                format, type = self._get_gl_format_and_type(data_format)

        # Workaround: don't use GL_UNPACK_ROW_LENGTH
        if gl._current_context._workaround_unpack_row_length:
            data_pitch = self.width * len(data_format)

        # Get data in required format (hopefully will be the same format it's
        # already in, unless that's an obscure format, upside-down or the
        # driver is old).
        data = self._convert(data_format, data_pitch)

        if data_pitch & 0x1:
            alignment = 1
        elif data_pitch & 0x2:
            alignment = 2
        else:
            alignment = 4
        row_length = data_pitch / len(data_format)
        glPushClientAttrib(GL_CLIENT_PIXEL_STORE_BIT)
        glPixelStorei(GL_UNPACK_ALIGNMENT, alignment)
        glPixelStorei(GL_UNPACK_ROW_LENGTH, row_length)
        self._apply_region_unpack()

        if target == GL_TEXTURE_3D:
            assert not internalformat
            glTexSubImage3D(target, level,
                            x, y, z,
                            self.width, self.height, 1,
                            format, type,
                            data)
        elif internalformat:
            glTexImage2D(target, level,
                         internalformat,
                         self.width, self.height,
                         0,
                         format, type,
                         data)
        else:
            glTexSubImage2D(target, level,
                            x, y,
                            self.width, self.height,
                            format, type,
                            data)
        glPopClientAttrib()

        if matrix:
            glPopMatrix()
            glMatrixMode(GL_MODELVIEW)

    def _apply_region_unpack(self):
        pass
   
    def _convert(self, format, pitch):
        '''Return data in the desired format; does not alter this instance's
        current format or pitch.
        '''
        if format == self._current_format and pitch == self._current_pitch:
            return self._current_data

        self._ensure_string_data()
        data = self._current_data
        current_pitch = self._current_pitch
        current_format = self._current_format
        sign_pitch = current_pitch / abs(current_pitch)
        if format != self._current_format:
            # Create replacement string, e.g. r'\4\1\2\3' to convert RGBA to
            # ARGB
            repl = ''
            for c in format:
                try:
                    idx = current_format.index(c) + 1
                except ValueError:
                    idx = 1
                repl += r'\%d' % idx

            if len(current_format) == 1:
                swap_pattern = self._swap1_pattern
            elif len(current_format) == 2:
                swap_pattern = self._swap2_pattern
            elif len(current_format) == 3:
                swap_pattern = self._swap3_pattern
            elif len(current_format) == 4:
                swap_pattern = self._swap4_pattern
            else:
                raise ImageException(
                    'Current image format is wider than 32 bits.')

            packed_pitch = self.width * len(current_format)
            if abs(self._current_pitch) != packed_pitch:
                # Pitch is wider than pixel data, need to go row-by-row.
                rows = re.findall(
                    '.' * abs(self._current_pitch), data, re.DOTALL)
                rows = [swap_pattern.sub(repl, r[:packed_pitch]) for r in rows]
                data = ''.join(rows)
            else:
                # Rows are tightly packed, apply regex over whole image.
                data = swap_pattern.sub(repl, data)

            # After conversion, rows will always be tightly packed
            current_pitch = sign_pitch * (len(format) * self.width)

        if pitch != current_pitch:
            diff = abs(current_pitch) - abs(pitch)
            if diff > 0:
                # New pitch is shorter than old pitch, chop bytes off each row
                pattern = re.compile(
                    '(%s)%s' % ('.' * abs(pitch), '.' * diff), re.DOTALL)
                data = pattern.sub(r'\1', data)    
            elif diff < 0:
                # New pitch is longer than old pitch, add '0' bytes to each row
                pattern = re.compile(
                    '(%s)' % ('.' * abs(current_pitch)), re.DOTALL)
                pad = '.' * -diff
                data = pattern.sub(r'\1%s' % pad, data)

            if current_pitch * pitch < 0:
                # Pitch differs in sign, swap row order
                rows = re.findall('.' * abs(pitch), data, re.DOTALL)
                rows.reverse()
                data = ''.join(rows)

        return data

    def _ensure_string_data(self):
        if type(self._current_data) is not str:
            buf = create_string_buffer(len(self._current_data))
            memmove(buf, self._current_data, len(self._current_data))
            self._current_data = buf.raw

    def _get_gl_format_and_type(self, format):
        if format == 'I':
            return GL_LUMINANCE, GL_UNSIGNED_BYTE
        elif format == 'L':
            return GL_LUMINANCE, GL_UNSIGNED_BYTE
        elif format == 'LA':
            return GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE
        elif format == 'R':
            return GL_RED, GL_UNSIGNED_BYTE
        elif format == 'G':
            return GL_GREEN, GL_UNSIGNED_BYTE
        elif format == 'B':
            return GL_BLUE, GL_UNSIGNED_BYTE
        elif format == 'A':
            return GL_ALPHA, GL_UNSIGNED_BYTE
        elif format == 'RGB':
            return GL_RGB, GL_UNSIGNED_BYTE
        elif format == 'RGBA':
            return GL_RGBA, GL_UNSIGNED_BYTE
        elif (format == 'ARGB' and
              gl_info.have_extension('GL_EXT_bgra') and
              gl_info.have_extension('GL_APPLE_packed_pixels')):
            return GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV
        elif (format == 'ABGR' and
              gl_info.have_extension('GL_EXT_abgr')):
            return GL_ABGR_EXT, GL_UNSIGNED_BYTE
        elif (format == 'BGR' and
              gl_info.have_extension('GL_EXT_bgra')):
            return GL_BGR, GL_UNSIGNED_BYTE
        elif (format == 'BGRA' and
              gl_info.have_extension('GL_EXT_bgra')):
            return GL_BGRA, GL_UNSIGNED_BYTE

        return None, None

    def _get_internalformat(self, format):
        if len(format) == 4:
            return GL_RGBA
        elif len(format) == 3:
            return GL_RGB
        elif len(format) == 2:
            return GL_LUMINANCE_ALPHA
        elif format == 'A':
            return GL_ALPHA
        elif format == 'L':
            return GL_LUMINANCE
        elif format == 'I':
            return GL_INTENSITY
        return GL_RGBA

class ImageDataRegion(ImageData):
    def __init__(self, x, y, width, height, image_data):
        super(ImageDataRegion, self).__init__(width, height,
            image_data._current_format, image_data._current_data, 
            image_data._current_pitch)
        self.x = x
        self.y = y

    def __getstate__(self):
        return {
            'width': self.width, 
            'height': self.height, 
            '_current_data': self.data, 
            '_current_format': self.format,
            '_desired_format': self.format,
            '_current_pitch': self.pitch,
            'pitch': self.pitch,
            'mipmap_images': self.mipmap_images,
            'x': self.x,
            'y': self.y
        }

    def _get_data(self):
        # Crop the data first
        x1 = len(self._current_format) * self.x
        x2 = len(self._current_format) * (self.x + self.width)

        self._ensure_string_data()
        data = self._convert(self._current_format, abs(self._current_pitch))
        rows = re.findall('.' * abs(self._current_pitch), data, re.DOTALL)
        rows = [row[x1:x2] for row in rows[self.y:self.y+self.height]]
        self._current_data = ''.join(rows)
        self._current_pitch = self.width * len(self._current_format)
        self._current_texture = None
        self.x = 0
        self.y = 0

        return super(ImageDataRegion, self)._get_data()

    def _set_data(self, data):
        self.x = 0
        self.y = 0
        super(ImageDataRegion, self)._set_data(data)
 
    data = property(_get_data, _set_data)
    
    def _apply_region_unpack(self):
        glPixelStorei(GL_UNPACK_SKIP_PIXELS, self.x)
        glPixelStorei(GL_UNPACK_SKIP_ROWS, self.y)

    def _ensure_string_data(self):
        super(ImageDataRegion, self)._ensure_string_data()

    def get_region(self, x, y, width, height):
        x += self.x
        y += self.y
        return super(ImageDataRegion, self).get_region(x, y, width, height)

class CompressedImageData(AbstractImage):
    '''Image representing some compressed data suitable for direct uploading
    to driver.
    '''

    _current_texture = None
    _current_mipmapped_texture = None

    def __init__(self, width, height, gl_format, data, 
                 extension=None, decoder=None):
        '''Construct a CompressedImageData with the given compressed data.

        :Parameters:
            `width` : int
                Width of image
            `height` : int
                Height of image
            `gl_format` : int
                GL constant giving format of compressed data; for example,
                ``GL_COMPRESSED_RGBA_S3TC_DXT5_EXT``.
            `data` : sequence
                String or array/list of bytes giving compressed image data.
            `extension` : str or None
                If specified, gives the name of a GL extension to check for
                before creating a texture.
            `decoder` : function(data, width, height) -> AbstractImage
                A function to decode the compressed data, to be used if the
                required extension is not present.
                
        '''
        if not _is_pow2(width) or not _is_pow2(height):
            raise ImageException('Dimensions of %r must be powers of 2' % self)

        super(CompressedImageData, self).__init__(width, height)
        self.data = data
        self.gl_format = gl_format
        self.extension = extension
        self.decoder = decoder
        self.mipmap_data = []

    def set_mipmap_data(self, level, data):
        '''Set data for a mipmap level.

        Supplied data gives a compressed image for the given mipmap level.
        The image must be of the correct dimensions for the level 
        (i.e., width >> level, height >> level); but this is not checked.  If
        any mipmap levels are specified, they are used; otherwise, mipmaps for
        `mipmapped_texture` are generated automatically.

        :Parameters:
            `level` : int
                Level of mipmap image to set.
            `data` : sequence
                String or array/list of bytes giving compressed image data.
                Data must be in same format as specified in constructor.

        '''
        # Extend mipmap_data list to required level
        self.mipmap_data += [None] * (level - len(self.mipmap_data))
        self.mipmap_data[level - 1] = data

    def _have_extension(self):
        return self.extension is None or gl_info.have_extension(self.extension)

    def _verify_driver_supported(self):
        '''Assert that the extension required for this image data is
        supported.

        Raises `ImageException` if not.
        '''

        if not self._have_extension():
            raise ImageException('%s is required to decode %r' % \
                (self.extension, self))

    def _get_texture(self):
        if self._current_texture:
            return self._current_texture

        texture = Texture.create_for_size(
            GL_TEXTURE_2D, self.width, self.height)
        glBindTexture(texture.target, texture.id)
        glTexParameteri(texture.target, GL_TEXTURE_MIN_FILTER, GL_LINEAR)

        if self._have_extension():
            glCompressedTexImage2DARB(texture.target, texture.level,
                self.gl_format,
                self.width, self.height, 0,
                len(self.data), self.data)
        else:
            image = self.decoder(self.data, self.width, self.height)
            texture = image.texture
            assert texture.width == self.width
            assert texture.height == self.height
                
        self._current_texture = texture
        return texture

    texture = property(_get_texture)

    def _get_mipmapped_texture(self):
        if self._current_mipmap_texture:
            return self._current_mipmap_texture

        if not self._have_extension():
            # TODO mip-mapped software decoded compressed textures.  For now,
            # just return a non-mipmapped texture.
            return self.texture

        texture = Texture.create_for_size(
            GL_TEXTURE_2D, self.width, self.height)
        glBindTexture(texture.target, texture.id)

        glTexParameteri(texture.target, GL_TEXTURE_MIN_FILTER,
                        GL_LINEAR_MIPMAP_LINEAR)

        if not self.mipmap_data:
            if not gl_info.have_version(1, 4):
                raise ImageException(
                  'Require GL 1.4 to generate mipmaps for compressed textures')
            glTexParameteri(texture.target, GL_GENERATE_MIPMAP, GL_TRUE)

        glCompressedTexImage2DARB(texture.target, texture.level,
            self.gl_format,
            self.width, self.height, 0,
            len(self.data), self.data) 

        width, height = self.width, self.height
        level = 0
        for data in self.mipmap_data:
            width >>= 1
            height >>= 1
            level += 1
            glCompressedTexImage2DARB(texture.target, level,
                self.gl_format,
                width, height, 0,
                len(data), data)

        self._current_mipmap_texture = texture
        return texture

    mipmapped_texture = property(_get_mipmapped_texture)

    def blit_to_texture(self, target, level, x, y, z):
        self._verify_driver_supported()

        if target == GL_TEXTURE_3D:
            glCompressedTexSubImage3DARB(target, level,
                x, y, z,
                self.width, self.height, 1,
                self.gl_format,
                len(self.data), self.data)
        else:
            glCompressedTexSubImage2DARB(target, level, 
                x, y,
                self.width, self.height,
                self.gl_format,
                len(self.data), self.data)
        
def _nearest_pow2(v):
    # From http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
    # Credit: Sean Anderson
    v -= 1
    v |= v >> 1
    v |= v >> 2
    v |= v >> 4
    v |= v >> 8
    v |= v >> 16
    return v + 1

def _is_pow2(v):
    # http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
    return (v & (v - 1)) == 0

class Texture(AbstractImage):
    '''An image loaded into video memory that can be efficiently drawn
    to the framebuffer.

    Typically you will get an instance of Texture by accessing the `texture`
    member of any other AbstractImage.

    :Ivariables:
        `region_class` : class (subclass of TextureRegion)
            Class to use when constructing regions of this texture.
        `tex_coords` : tuple
            12-tuple of float, named (u1, v1, r1, u2, v2, r2, ...).  u, v, r
            give the 3D texture coordinates for vertices 1-4.  The vertices
            are specified in the order bottom-left, bottom-right, top-right
            and top-left.
        `target` : int
            The GL texture target (e.g., ``GL_TEXTURE_2D``).
        `level` : int
            The mipmap level of this texture.

    '''

    region_class = None # Set to TextureRegion after it's defined
    tex_coords = (0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0.)
    level = 0
    images = 1
    x = y = z = 0

    def __init__(self, width, height, target, id):
        super(Texture, self).__init__(width, height)
        self.target = target
        self.id = id
        self._context = get_current_context()

    def delete(self):
        warnings.warn(
            'Texture.delete() is deprecated; textures are '
            'released through GC now')
        self._context.delete_texture(self.id)
        self.id = 0

    def __del__(self):
        try:
            self._context.delete_texture(self.id)
        except AttributeError:
            pass

    @classmethod
    def create_for_size(cls, target, min_width, min_height,
                        internalformat=None):
        '''Create a Texture with dimensions at least min_width, min_height.
        On return, the texture will be bound.

        :Parameters:
            `target` : int
                GL constant giving texture target to use, typically
                ``GL_TEXTURE_2D``.
            `min_width` : int
                Minimum width of texture (may be increased to create a power
                of 2).
            `min_height` : int
                Minimum height of texture (may be increased to create a power
                of 2).
            `internalformat` : int
                GL constant giving internal format of texture; for example,
                ``GL_RGBA``.  If unspecified, the texture will not be
                initialised (only the texture name will be created on the
                instance).   If specified, the image will be initialised
                to this format with zero'd data.
        '''
        if target not in (GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_RECTANGLE_ARB):
            width = _nearest_pow2(min_width)
            height = _nearest_pow2(min_height)
            tex_coords = cls.tex_coords
        else:
            width = min_width
            height = min_height
            tex_coords = (0., 0., 0., 
                          width, 0., 0., 
                          width, height, 0., 
                          0., height, 0.)
        id = GLuint()
        glGenTextures(1, byref(id))
        glBindTexture(target, id.value)
        glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR)

        if internalformat is not None:
            blank = (GLubyte * (width * height * 4))()
            glTexImage2D(target, 0,
                         internalformat,
                         width, height,
                         0,
                         GL_RGBA, GL_UNSIGNED_BYTE,
                         blank)
                         
        texture = cls(width, height, target, id.value)
        texture.tex_coords = tex_coords
        return texture

    def get_image_data(self, z=0):
        '''Get the image data of this texture.

        :Parameters:
            `z` : int
                For 3D textures, the image slice to retrieve.

        :rtype: `ImageData`
        '''
        glBindTexture(self.target, self.id)

        # Always extract complete RGBA data.  Could check internalformat
        # to only extract used channels. XXX
        format = 'RGBA'
        gl_format = GL_RGBA

        glPushClientAttrib(GL_CLIENT_PIXEL_STORE_BIT)
        glPixelStorei(GL_PACK_ALIGNMENT, 1)
        buffer = \
            (GLubyte * (self.width * self.height * self.images * len(format)))()
        glGetTexImage(self.target, self.level, 
                      gl_format, GL_UNSIGNED_BYTE, buffer)
        glPopClientAttrib()

        data = ImageData(self.width, self.height, format, buffer)
        if self.images > 1:
            data = data.get_region(0, z * self.height, self.width, self.height)
        return data

    image_data = property(get_image_data,
        doc='''An ImageData view of this texture.  
        
        Changes to the returned instance will not be reflected in this
        texture.  If the texture is a 3D texture, the first image will be 
        returned.  See also `get_image_data`.  Read-only.
        
        :type: `ImageData`
        ''')


    texture = property(lambda self: self)

    # no implementation of blit_to_texture yet (could use aux buffer)

    def blit(self, x, y, z=0, width=None, height=None):
        # Create interleaved array in T4F_V4F format
        t = self.tex_coords
        w = width is None and self.width or width
        h = height is None and self.height or height
        array = (GLfloat * 32)(
             t[0],  t[1],  t[2],  1.,
             x,     y,     z,     1.,
             t[3],  t[4],  t[5],  1., 
             x + w, y,     z,     1.,
             t[6],  t[7],  t[8],  1., 
             x + w, y + h, z,     1.,
             t[9],  t[10], t[11], 1., 
             x,     y + h, z,     1.)

        glPushAttrib(GL_ENABLE_BIT)
        glEnable(self.target)
        glBindTexture(self.target, self.id)
        glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT)
        glInterleavedArrays(GL_T4F_V4F, 0, array)
        glDrawArrays(GL_QUADS, 0, 4)
        glPopClientAttrib()
        glPopAttrib()

    def blit_into(self, source, x, y, z):
        glBindTexture(self.target, self.id)
        source.blit_to_texture(self.target, self.level, x, y, z)

    def get_region(self, x, y, width, height):
        return self.region_class(x, y, 0, width, height, self)


class TextureRegion(Texture):
    '''A rectangular region of a texture, presented as if it were
    a separate texture.
    '''

    def __init__(self, x, y, z, width, height, owner):
        super(TextureRegion, self).__init__(
            width, height, owner.target, owner.id)
        
        self.x = x
        self.y = y
        self.z = z
        self.owner = owner
        u1 = x / float(owner.width)
        v1 = y / float(owner.height)
        u2 = (x + width) / float(owner.width)
        v2 = (y + height) / float(owner.height)
        r = z / float(owner.images)
        self.tex_coords = (u1, v1, r, u2, v1, r, u2, v2, r, u1, v2, r)


    def _get_image_data(self):
        image_data = self.owner.get_image_data(self.z)
        return image_data.get_region(self.x, self.y, self.width, self.height)

    image_data = property(_get_image_data)

    def get_region(self, x, y, width, height):
        x += self.x
        y += self.y
        return self.region_class(x, y, self.z, width, height, self.owner)

    def blit_into(self, source, x, y, z):
        self.owner.blit_into(source, x + self.x, y + self.y, z + self.z)

    def __del__(self):
        # only the owner Texture should handle deletion
        pass

Texture.region_class = TextureRegion

class Texture3D(Texture, UniformTextureSequence):
    '''A texture with more than one image slice.

    Use `create_for_images` or `create_for_image_grid` classmethod to
    construct.
    '''
    item_width = 0
    item_height = 0
    items = ()

    @classmethod
    def create_for_images(cls, images, internalformat=GL_RGBA):
        item_width = images[0].width
        item_height = images[0].height
        for image in images:
            if image.width != item_width or image.height != item_height:
                raise ImageException('Images do not have same dimensions.')

        depth = len(images)
        if not gl_info.have_version(2,0):
            depth = _nearest_pow2(depth)

        texture = cls.create_for_size(GL_TEXTURE_3D, item_width, item_height)
        texture.images = depth
        
        blank = (GLubyte * (texture.width * texture.height * texture.images))()
        glBindTexture(texture.target, texture.id)
        glTexImage3D(texture.target, texture.level,
                     internalformat,
                     texture.width, texture.height, texture.images, 0,
                     GL_ALPHA, GL_UNSIGNED_BYTE,
                     blank)

        items = []
        for i, image in enumerate(images):
            item = cls.region_class(0, 0, i, item_width, item_height, texture)
            items.append(item)
            image.blit_to_texture(texture.target, texture.level, 0, 0, i)

        texture.items = items
        texture.item_width = item_width
        texture.item_height = item_height
        return texture

    @classmethod
    def create_for_image_grid(cls, grid, internalformat=GL_RGBA):
        return cls.create_for_images(grid[:], internalformat)

    def __len__(self):
        return len(self.items)

    def __getitem__(self, index):
        return self.items[index]

    def __setitem__(self, index, value):
        if type(index) is slice:
            for item, image in zip(self[index], value):
                image.blit_to_texture(self.target, self.level, 0, 0, item.z)
        else:
            value.blit_to_texture(self.target, self.level, 0, 0, self[index].z)

class TileableTexture(Texture):
    '''A texture that can be tiled efficiently.

    Use `create_for_image` classmethod to construct.
    '''
    def __init__(self, width, height, target, id):
        if not _is_pow2(width) or not _is_pow2(height):
            raise ImageException(
                'TileableTexture requires dimensions that are powers of 2')
        super(TileableTexture, self).__init__(width, height, target, id)
        
    def get_region(self, x, y, width, height):
        raise ImageException('Cannot get region of %r' % self)

    def blit_tiled(self, x, y, z, width, height):
        '''Blit this texture tiled over the given area.'''
        u1 = v1 = 0
        u2 = width / float(self.width)
        v2 = height / float(self.height)
        w, h = width, height
        t = self.tex_coords
        array = (GLfloat * 32)(
             u1,      v1,      t[2],  1.,
             x,       y,       z,     1.,
             u2,      v1,      t[5],  1., 
             x + w,   y,       z,     1.,
             u2,      v2,      t[8],  1., 
             x + w,   y + h,   z,     1.,
             u1,      v2,      t[11], 1., 
             x,       y + h,   z,     1.)

        glPushAttrib(GL_ENABLE_BIT)
        glEnable(self.target)
        glBindTexture(self.target, self.id)
        glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT)
        glInterleavedArrays(GL_T4F_V4F, 0, array)
        glDrawArrays(GL_QUADS, 0, 4)
        glPopClientAttrib()
        glPopAttrib()
        

    @classmethod
    def create_for_image(cls, image):
        if not _is_pow2(image.width) or not _is_pow2(image.height):
            # Potentially unnecessary conversion if a GL format exists.
            image = image.image_data
            image.format = 'RGBA'
            image.pitch = image.width * len(image.format)
            texture_width = _nearest_pow2(image.width)
            texture_height = _nearest_pow2(image.height)
            newdata = c_buffer(texture_width * texture_height *
                               len(image.format))
            gluScaleImage(GL_RGBA,
                          image.width, image.height,
                          GL_UNSIGNED_BYTE,
                          image.data,
                          texture_width,
                          texture_height,
                          GL_UNSIGNED_BYTE,
                          newdata)
            image = ImageData(texture_width, texture_height, image.format,
                              newdata)

        image = image.image_data
        return image.create_texture(cls)

class DepthTexture(Texture):
    '''A texture with depth samples (typically 24-bit).'''
    def blit_into(self, source, x, y, z):
        glBindTexture(self.target, self.id)
        source.blit_to_texture(self.level, x, y, z)

class BufferManager(object):
    '''Manages the set of framebuffers for a context.

    Use `get_buffer_manager` to obtain the singleton instance of this class.
    '''
    def __init__(self):
        self.color_buffer = None
        self.depth_buffer = None

        aux_buffers = GLint()
        glGetIntegerv(GL_AUX_BUFFERS, byref(aux_buffers))
        self.free_aux_buffers = [GL_AUX0, 
                                 GL_AUX1, 
                                 GL_AUX2,
                                 GL_AUX3][:aux_buffers.value]

        stencil_bits = GLint()
        glGetIntegerv(GL_STENCIL_BITS, byref(stencil_bits))
        self.free_stencil_bits = range(stencil_bits.value)

        self.refs = []

    def get_viewport(self):
        '''Get the current OpenGL viewport dimensions.

        :rtype: 4-tuple of float.
        :return: Left, top, right and bottom dimensions.
        '''
        viewport = (GLint * 4)()
        glGetIntegerv(GL_VIEWPORT, viewport)
        return viewport
    
    def get_color_buffer(self):
        '''Get the color buffer.

        :rtype: `ColorBufferImage`
        '''
        if not self.color_buffer:
            viewport = self.get_viewport()
            self.color_buffer = ColorBufferImage(*viewport)
        return self.color_buffer

    def get_aux_buffer(self):
        '''Get a free auxilliary buffer.

        If not aux buffers are available, `ImageException` is raised.  Buffers
        are released when they are garbage collected.
        
        :rtype: `ColorBufferImage`
        '''
        if not self.free_aux_buffers:
            raise ImageException('No free aux buffer is available.')

        gl_buffer = self.free_aux_buffers.pop(0)
        viewport = self.get_viewport()
        buffer = ColorBufferImage(*viewport)
        buffer.gl_buffer = gl_buffer

        def release_buffer(ref, self=self):
            self.free_aux_buffers.insert(0, gl_buffer)
        self.refs.append(weakref.ref(buffer, release_buffer))
            
        return buffer

    def get_depth_buffer(self):
        '''Get the depth buffer.

        :rtype: `DepthBufferImage`
        '''
        if not self.depth_buffer:
            viewport = self.get_viewport()
            self.depth_buffer = DepthBufferImage(*viewport)
        return self.depth_buffer

    def get_buffer_mask(self):
        '''Get a free bitmask buffer.

        A bitmask buffer is a buffer referencing a single bit in the stencil
        buffer.  If no bits are free, `ImageException` is raised.  Bits are
        released when the bitmask buffer is garbage collected.

        :rtype: `BufferImageMask`
        '''
        if not self.free_stencil_bits:
            raise ImageException('No free stencil bits are available.')

        stencil_bit = self.free_stencil_bits.pop(0)
        viewport = self.get_viewport()
        buffer = BufferImageMask(x, y, width, height)
        buffer.stencil_bit = stencil_bit

        def release_buffer(ref, self=self):
            self.free_stencil_bits.insert(0, stencil_bit)
        self.refs.append(weakref.ref(buffer, release_buffer))

        return buffer

def get_buffer_manager():
    '''Get the singleton buffer manager.
    
    :rtype: `BufferManager`
    '''
    context = get_current_context()
    if not hasattr(context, 'image_buffer_manager'):
        context.image_buffer_manager = BufferManager()
    return context.image_buffer_manager

# XXX BufferImage could be generalised to support EXT_framebuffer_object's
# renderbuffer.
class BufferImage(AbstractImage):
    '''An abstract framebuffer.
    '''
    #: The OpenGL read and write target for this buffer.
    gl_buffer = GL_BACK

    #: The OpenGL format constant for image data.
    gl_format = 0

    #: The format string used for image data.
    format = ''

    owner = None

    # TODO: enable methods

    def __init__(self, x, y, width, height):
        self.x = x
        self.y = y
        self.width = width
        self.height = height

    def _get_image_data(self):
        buffer = (GLubyte * (len(self.format) * self.width * self.height))()

        x = self.x
        y = self.y
        if self.owner:
            x += self.owner.x
            y += self.owner.y

        glReadBuffer(self.gl_buffer)
        glPushClientAttrib(GL_CLIENT_PIXEL_STORE_BIT)
        glPixelStorei(GL_PACK_ALIGNMENT, 1)
        glReadPixels(x, y, self.width, self.height, 
                     self.gl_format, GL_UNSIGNED_BYTE, buffer)
        glPopClientAttrib()

        return ImageData(self.width, self.height, self.format, buffer)

    image_data = property(_get_image_data)

    def get_region(self, x, y, width, height):
        if self.owner:
            return self.owner.get_region(x + self.x, y + self.y, width, height)

        region = self.__class__(x + self.x, y + self.y, width, height)
        region.owner = self
        return region

class ColorBufferImage(BufferImage):
    '''A color framebuffer.

    This class is used to wrap both the primary color buffer (i.e., the back
    buffer) or any one of the auxilliary buffers.
    '''
    gl_format = GL_RGBA
    format = 'RGBA'

    def _get_texture(self):
        texture = Texture.create_for_size(GL_TEXTURE_2D, 
            self.width, self.height)

        if texture.width != self.width or texture.height != self.height:
            texture = texture.get_region(0, 0, self.width, self.height)
            width = texture.owner.width
            height = texture.owner.height
            blank = (c_ubyte * (width * height * 4))()
            glTexImage2D(texture.target, texture.level,
                         GL_RGBA,
                         width, height,
                         0,
                         GL_RGBA, GL_UNSIGNED_BYTE,
                         blank)
            self.blit_to_texture(texture.target, texture.level, 0, 0, 0)
        else:
            glReadBuffer(self.gl_buffer)
            glCopyTexImage2D(texture.target, texture.level,
                             GL_RGBA,
                             self.x, self.y, self.width, self.height,
                             0)
        return texture

    texture = property(_get_texture)

    def blit_to_texture(self, target, level, x, y, z):
        glReadBuffer(self.gl_buffer)
        glCopyTexSubImage2D(target, level, 
                            x, y,
                            self.x, self.y, self.width, self.height) 

class DepthBufferImage(BufferImage):
    '''The depth buffer.
    '''
    gl_format = GL_DEPTH_COMPONENT
    format = 'L'

    def _get_texture(self):
        if not _is_pow2(self.width) or not _is_pow2(self.height):
            raise ImageException(
                'Depth texture requires that buffer dimensions be powers of 2')
        
        texture = DepthTexture.create_for_size(GL_TEXTURE_2D,
            self.width, self.height)
        glReadBuffer(self.gl_buffer)
        glCopyTexImage2D(texture.target, 0,
                         GL_DEPTH_COMPONENT,
                         self.x, self.y, self.width, self.height,
                         0)
        return texture

    texture = property(_get_texture)

    def blit_to_texture(self, target, level, x, y, z):
        glReadBuffer(self.gl_buffer)
        glCopyTexSubImage2D(target, level,
                            x, y,
                            self.x, self.y, self.width, self.height)


class BufferImageMask(BufferImage):
    '''A single bit of the stencil buffer.
    '''
    gl_format = GL_STENCIL_INDEX
    format = 'L'

    # TODO mask methods

class ImageGrid(AbstractImage, AbstractImageSequence):
    '''An imaginary grid placed over an image allowing easy access to
    regular regions of that image.

    The grid can be accessed either as a complete image, or as a sequence
    of images.  The most useful applications are to access the grid
    as a `TextureGrid`::

        image_grid = ImageGrid(...)
        texture_grid = image_grid.texture_sequence

    or as a `Texture3D`::

        image_grid = ImageGrid(...)
        texture_3d = Texture3D.create_for_image_grid(image_grid)

    '''
    _items = ()
    _texture_grid = None

    def __init__(self, image, rows, columns, 
                 item_width=None, item_height=None,
                 row_padding=0, column_padding=0):
        '''Construct a grid for the given image.

        You can specify parameters for the grid, for example setting
        the padding between cells.  Grids are always aligned to the 
        bottom-left corner of the image.

        :Parameters:
            `image` : AbstractImage
                Image over which to construct the grid.
            `rows` : int
                Number of rows in the grid.
            `columns` : int
                Number of columns in the grid.
            `item_width` : int
                Width of each column.  If unspecified, is calculated such
                that the entire image width is used.
            `item_height` : int
                Height of each row.  If unspecified, is calculated such that
                the entire image height is used.
            `row_padding` : int
                Pixels separating adjacent rows.  The padding is only
                inserted between rows, not at the edges of the grid.
            `column_padding` : int
                Pixels separating adjacent columns.  The padding is only 
                inserted between columns, not at the edges of the grid.
        '''
        super(ImageGrid, self).__init__(image.width, image.height)

        if item_width is None:
            item_width = \
                int((image.width - column_padding * (columns - 1)) / columns)
        if item_height is None:
            item_height = \
                int((image.height - row_padding * (rows - 1)) / rows) 
        self.image = image
        self.rows = rows
        self.columns = columns
        self.item_width = item_width
        self.item_height = item_height
        self.row_padding = row_padding
        self.column_padding = column_padding

    def _get_texture(self):
        return self.image.texture

    texture = property(_get_texture)

    def _get_image_data(self):
        return self.image.image_data

    image_data = property(_get_image_data)

    def _get_texture_sequence(self):
        if not self._texture_grid:
            self._texture_grid = TextureGrid(self)
        return self._texture_grid

    texture_sequence = property(_get_texture_sequence)

    def __len__(self):
        return self.rows * self.columns

    def __getitem__(self, index):
        if not self._items:
            self._items = []
            y = 0
            for row in range(self.rows):
                x = 0
                for col in range(self.columns):
                    self._items.append(self.image.get_region(
                        x, y, self.item_width, self.item_height))
                    x += self.item_width + self.column_padding
                y += self.item_height + self.row_padding

        # TODO tuples
        return self._items[index]

class TextureGrid(TextureRegion, UniformTextureSequence):
    '''A texture containing a regular grid of texture regions.

    To construct, create an `ImageGrid` first::

        image_grid = ImageGrid(...)
        texture_grid = TextureGrid(image_grid)

    The texture grid can be accessed as a single texture, or as a sequence
    of `TextureRegion`.  When accessing as a sequence, you can specify
    integer indexes, in which the images are arranged in rows from the
    bottom-left to the top-right::

        # assume the texture_grid is 3x3:
        current_texture = texture_grid[3] # get the middle-left image

    You can also specify tuples in the sequence methods, which are addressed
    as ``row, column``::

        # equivalent to the previous example:
        current_texture = texture_grid[1, 0]

    When using tuples in a slice, the returned sequence is over the
    rectangular region defined by the slice::

        # returns center, center-right, center-top, top-right images in that
        # order:
        images = texture_grid[(1,1):]
        # equivalent to
        images = texture_grid[(1,1):(3,3)]

    '''
    items = ()
    rows = 1
    columns = 1
    item_width = 0
    item_height = 0

    def __init__(self, grid):
        image = grid.texture
        if isinstance(image, TextureRegion):
            owner = image.owner
        else:
            owner = image

        super(TextureGrid, self).__init__(
            image.x, image.y, image.z, image.width, image.height, owner)
        
        items = []
        y = 0
        for row in range(grid.rows):
            x = 0
            for col in range(grid.columns):
                items.append(
                    self.get_region(x, y, grid.item_width, grid.item_height))
                x += grid.item_width + grid.column_padding
            y += grid.item_height + grid.row_padding

        self.items = items
        self.rows = grid.rows
        self.columns = grid.columns
        self.item_width = grid.item_width
        self.item_height = grid.item_height
        
    def get(self, row, column):
        return self[(row, column)]

    def __getitem__(self, index):
        if type(index) is slice:
            if type(index.start) is not tuple and \
               type(index.stop) is not tuple:
                return self.items[index]
            else:
                row1 = 0
                col1 = 0
                row2 = self.rows
                col2 = self.columns
                if type(index.start) is tuple:
                    row1, col1 = index.start
                elif type(index.start) is int:
                    row1 = index.start / self.columns
                    col1 = index.start % self.columns
                assert row1 >= 0 and col1 >= 0 and \
                       row1 < self.rows and col1 < self.columns

                if type(index.stop) is tuple:
                    row2, col2 = index.stop
                elif type(index.stop) is int:
                    row2 = index.stop / self.columns
                    col2 = index.stop % self.columns
                assert row2 >= 0 and col2 >= 0 and \
                       row2 < self.rows and col2 < self.columns

                result = []
                i = row1 * self.columns
                for row in range(row1, row2):
                    result += self.items[i+col1:i+col2]
                    i += self.columns
                return result
        else:
            if type(index) is tuple:
                row, column = index
                assert row >= 0 and column >= 0 and \
                       row < self.rows and column < self.columns
                return self.items[row * self.columns + column]
            elif type(index) is int:
                return self.items[index]

    def __setitem__(self, index, value):
        if type(index) is slice:
            for region, image in zip(self[index], value):
                if image.width != self.item_width or \
                   image.height != self.item_height:
                    raise ImageException('Image has incorrect dimensions')
                image.blit_into(region, 0, 0, 0)
        else:
            image = value
            if image.width != self.item_width or \
               image.height != self.item_height:
                raise ImageException('Image has incorrect dimensions')
            image.blit_into(self[index], 0, 0, 0)

    def __len__(self):
        return len(self.items)

# Initialise default codecs
from pyglet.image import codecs as _codecs
_codecs.add_default_image_codecs()
