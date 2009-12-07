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
# png.py - PNG encoder in pure Python
# Copyright (C) 2006 Johann C. Rocholl <johann@browsershots.org>
# <ah> Modifications for pyglet by Alex Holkner <alex.holkner@gmail.com> 
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Contributors (alphabetical):
# Nicko van Someren <nicko@nicko.org>
#
# Changelog (recent first):
# 2006-06-17 Nicko: Reworked into a class, faster interlacing.
# 2006-06-17 Johann: Very simple prototype PNG decoder.
# 2006-06-17 Nicko: Test suite with various image generators.
# 2006-06-17 Nicko: Alpha-channel, grey-scale, 16-bit/plane support.
# 2006-06-15 Johann: Scanline iterator interface for large input files.
# 2006-06-09 Johann: Very simple prototype PNG encoder.


"""
Pure Python PNG Reader/Writer

This is an implementation of a subset of the PNG specification at
http://www.w3.org/TR/2003/REC-PNG-20031110 in pure Python. It reads
and writes PNG files with 8/16/24/32/48/64 bits per pixel (greyscale,
RGB, RGBA, with 8 or 16 bits per layer), with a number of options. For
help, type "import png; help(png)" in your python interpreter.

This file can also be used as a command-line utility to convert PNM
files to PNG. The interface is similar to that of the pnmtopng program
from the netpbm package. Type "python png.py --help" at the shell
prompt for usage and a list of options.
"""


__revision__ = '$Rev$'
__date__ = '$Date$'
__author__ = '$Author$'


import sys
import zlib
import struct
import math
from array import array


_adam7 = ((0, 0, 8, 8),
          (4, 0, 8, 8),
          (0, 4, 4, 8),
          (2, 0, 4, 4),
          (0, 2, 2, 4),
          (1, 0, 2, 2),
          (0, 1, 1, 2))


def interleave_planes(ipixels, apixels, ipsize, apsize):
    """
    Interleave color planes, e.g. RGB + A = RGBA.

    Return an array of pixels consisting of the ipsize bytes of data
    from each pixel in ipixels followed by the apsize bytes of data
    from each pixel in apixels, for an image of size width x height.
    """
    itotal = len(ipixels)
    atotal = len(apixels)
    newtotal = itotal + atotal
    newpsize = ipsize + apsize
    # Set up the output buffer
    out = array('B')
    # It's annoying that there is no cheap way to set the array size :-(
    out.extend(ipixels)
    out.extend(apixels)
    # Interleave in the pixel data
    for i in range(ipsize):
        out[i:newtotal:newpsize] = ipixels[i:itotal:ipsize]
    for i in range(apsize):
        out[i+ipsize:newtotal:newpsize] = apixels[i:atotal:apsize]
    return out

class Error(Exception):
    pass

class Writer:
    """
    PNG encoder in pure Python.
    """

    def __init__(self, width, height,
                 transparent=None,
                 background=None,
                 gamma=None,
                 greyscale=False,
                 has_alpha=False,
                 bytes_per_sample=1,
                 compression=None,
                 interlaced=False,
                 chunk_limit=2**20):
        """
        Create a PNG encoder object.

        Arguments:
        width, height - size of the image in pixels
        transparent - create a tRNS chunk
        background - create a bKGD chunk
        gamma - create a gAMA chunk
        greyscale - input data is greyscale, not RGB
        has_alpha - input data has alpha channel (RGBA)
        bytes_per_sample - 8-bit or 16-bit input data
        compression - zlib compression level (1-9)
        chunk_limit - write multiple IDAT chunks to save memory

        If specified, the transparent and background parameters must
        be a tuple with three integer values for red, green, blue, or
        a simple integer (or singleton tuple) for a greyscale image.

        If specified, the gamma parameter must be a float value.

        """
        if width <= 0 or height <= 0:
            raise ValueError("width and height must be greater than zero")

        if has_alpha and transparent is not None:
            raise ValueError(
                "transparent color not allowed with alpha channel")

        if bytes_per_sample < 1 or bytes_per_sample > 2:
            raise ValueError("bytes per sample must be 1 or 2")

        if transparent is not None:
            if greyscale:
                if type(transparent) is not int:
                    raise ValueError(
                        "transparent color for greyscale must be integer")
            else:
                if not (len(transparent) == 3 and
                        type(transparent[0]) is int and
                        type(transparent[1]) is int and
                        type(transparent[2]) is int):
                    raise ValueError(
                        "transparent color must be a triple of integers")

        if background is not None:
            if greyscale:
                if type(background) is not int:
                    raise ValueError(
                        "background color for greyscale must be integer")
            else:
                if not (len(background) == 3 and
                        type(background[0]) is int and
                        type(background[1]) is int and
                        type(background[2]) is int):
                    raise ValueError(
                        "background color must be a triple of integers")

        self.width = width
        self.height = height
        self.transparent = transparent
        self.background = background
        self.gamma = gamma
        self.greyscale = greyscale
        self.has_alpha = has_alpha
        self.bytes_per_sample = bytes_per_sample
        self.compression = compression
        self.chunk_limit = chunk_limit
        self.interlaced = interlaced

        if self.greyscale:
            self.color_depth = 1
            if self.has_alpha:
                self.color_type = 4
                self.psize = self.bytes_per_sample * 2
            else:
                self.color_type = 0
                self.psize = self.bytes_per_sample
        else:
            self.color_depth = 3
            if self.has_alpha:
                self.color_type = 6
                self.psize = self.bytes_per_sample * 4
            else:
                self.color_type = 2
                self.psize = self.bytes_per_sample * 3

    def write_chunk(self, outfile, tag, data):
        """
        Write a PNG chunk to the output file, including length and checksum.
        """
        # http://www.w3.org/TR/PNG/#5Chunk-layout
        outfile.write(struct.pack("!I", len(data)))
        outfile.write(tag)
        outfile.write(data)
        checksum = zlib.crc32(tag)
        checksum = zlib.crc32(data, checksum)
        # <ah> Avoid DeprecationWarning: struct integer overflow masking
        #      with Python2.5/Windows.
        checksum = checksum & 0xffffffff
        outfile.write(struct.pack("!I", checksum))

    def write(self, outfile, scanlines):
        """
        Write a PNG image to the output file.
        """
        # http://www.w3.org/TR/PNG/#5PNG-file-signature
        outfile.write(struct.pack("8B", 137, 80, 78, 71, 13, 10, 26, 10))

        # http://www.w3.org/TR/PNG/#11IHDR
        if self.interlaced:
            interlaced = 1
        else:
            interlaced = 0
        self.write_chunk(outfile, 'IHDR',
                         struct.pack("!2I5B", self.width, self.height,
                                     self.bytes_per_sample * 8,
                                     self.color_type, 0, 0, interlaced))

        # http://www.w3.org/TR/PNG/#11tRNS
        if self.transparent is not None:
            if self.greyscale:
                self.write_chunk(outfile, 'tRNS',
                                 struct.pack("!1H", *self.transparent))
            else:
                self.write_chunk(outfile, 'tRNS',
                                 struct.pack("!3H", *self.transparent))

        # http://www.w3.org/TR/PNG/#11bKGD
        if self.background is not None:
            if self.greyscale:
                self.write_chunk(outfile, 'bKGD',
                                 struct.pack("!1H", *self.background))
            else:
                self.write_chunk(outfile, 'bKGD',
                                 struct.pack("!3H", *self.background))

        # http://www.w3.org/TR/PNG/#11gAMA
        if self.gamma is not None:
            self.write_chunk(outfile, 'gAMA',
                             struct.pack("!L", int(self.gamma * 100000)))

        # http://www.w3.org/TR/PNG/#11IDAT
        if self.compression is not None:
            compressor = zlib.compressobj(self.compression)
        else:
            compressor = zlib.compressobj()

        data = array('B')
        for scanline in scanlines:
            data.append(0)
            data.extend(scanline)
            if len(data) > self.chunk_limit:
                compressed = compressor.compress(data.tostring())
                if len(compressed):
                    # print >> sys.stderr, len(data), len(compressed)
                    self.write_chunk(outfile, 'IDAT', compressed)
                data = array('B')
        if len(data):
            compressed = compressor.compress(data.tostring())
        else:
            compressed = ''
        flushed = compressor.flush()
        if len(compressed) or len(flushed):
            # print >> sys.stderr, len(data), len(compressed), len(flushed)
            self.write_chunk(outfile, 'IDAT', compressed + flushed)

        # http://www.w3.org/TR/PNG/#11IEND
        self.write_chunk(outfile, 'IEND', '')

    def write_array(self, outfile, pixels):
        """
        Encode a pixel array to PNG and write output file.
        """
        if self.interlaced:
            self.write(outfile, self.array_scanlines_interlace(pixels))
        else:
            self.write(outfile, self.array_scanlines(pixels))

    def convert_ppm(self, ppmfile, outfile):
        """
        Convert a PPM file containing raw pixel data into a PNG file
        with the parameters set in the writer object.
        """
        if self.interlaced:
            pixels = array('B')
            pixels.fromfile(ppmfile,
                            self.bytes_per_sample * self.color_depth *
                            self.width * self.height)
            self.write(outfile, self.array_scanlines_interlace(pixels))
        else:
            self.write(outfile, self.file_scanlines(ppmfile))

    def convert_ppm_and_pgm(self, ppmfile, pgmfile, outfile):
        """
        Convert a PPM and PGM file containing raw pixel data into a
        PNG outfile with the parameters set in the writer object.
        """
        pixels = array('B')
        pixels.fromfile(ppmfile,
                        self.bytes_per_sample * self.color_depth *
                        self.width * self.height)
        apixels = array('B')
        apixels.fromfile(pgmfile,
                         self.bytes_per_sample *
                         self.width * self.height)
        pixels = interleave_planes(pixels, apixels,
                                   self.bytes_per_sample * self.color_depth,
                                   self.bytes_per_sample)
        if self.interlaced:
            self.write(outfile, self.array_scanlines_interlace(pixels))
        else:
            self.write(outfile, self.array_scanlines(pixels))

    def file_scanlines(self, infile):
        """
        Generator for scanlines from an input file.
        """
        row_bytes = self.psize * self.width
        for y in range(self.height):
            scanline = array('B')
            scanline.fromfile(infile, row_bytes)
            yield scanline

    def array_scanlines(self, pixels):
        """
        Generator for scanlines from an array.
        """
        row_bytes = self.width * self.psize
        stop = 0
        for y in range(self.height):
            start = stop
            stop = start + row_bytes
            yield pixels[start:stop]

    def old_array_scanlines_interlace(self, pixels):
        """
        Generator for interlaced scanlines from an array.
        http://www.w3.org/TR/PNG/#8InterlaceMethods
        """
        row_bytes = self.psize * self.width
        for xstart, ystart, xstep, ystep in _adam7:
            for y in range(ystart, self.height, ystep):
                if xstart < self.width:
                    if xstep == 1:
                        offset = y*row_bytes
                        yield pixels[offset:offset+row_bytes]
                    else:
                        row = array('B')
                        offset = y*row_bytes + xstart* self.psize
                        skip = self.psize * xstep
                        for x in range(xstart, self.width, xstep):
                            row.extend(pixels[offset:offset + self.psize])
                            offset += skip
                        yield row

    def array_scanlines_interlace(self, pixels):
        """
        Generator for interlaced scanlines from an array.
        http://www.w3.org/TR/PNG/#8InterlaceMethods
        """
        row_bytes = self.psize * self.width
        for xstart, ystart, xstep, ystep in _adam7:
            for y in range(ystart, self.height, ystep):
                if xstart >= self.width:
                    continue
                if xstep == 1:
                    offset = y * row_bytes
                    yield pixels[offset:offset+row_bytes]
                else:
                    row = array('B')
                    # Note we want the ceiling of (self.width - xstart) / xtep
                    row_len = self.psize * (
                        (self.width - xstart + xstep - 1) / xstep)
                    # There's no easier way to set the length of an array
                    row.extend(pixels[0:row_len])
                    offset = y * row_bytes + xstart * self.psize
                    end_offset = (y+1) * row_bytes
                    skip = self.psize * xstep
                    for i in range(self.psize):
                        row[i:row_len:self.psize] = \
                            pixels[offset+i:end_offset:skip]
                    yield row

class _readable:
    """
    A simple file-like interface for strings and arrays.
    """

    def __init__(self, buf):
        self.buf = buf
        self.offset = 0

    def read(self, n):
        r = buf[offset:offset+n]
        if isinstance(r, array):
            r = r.tostring()
        offset += n
        return r

class Reader:
    """
    PNG decoder in pure Python.
    """

    def __init__(self, _guess=None, **kw):
        """
        Create a PNG decoder object.

        The constructor expects exactly one keyword argument. If you
        supply a positional argument instead, it will guess the input
        type. You can choose among the following arguments:
        filename - name of PNG input file
        file - object with a read() method
        pixels - array or string with PNG data

        """
        if ((_guess is not None and len(kw) != 0) or
            (_guess is None and len(kw) != 1)):
            raise TypeError("Reader() takes exactly 1 argument")

        if _guess is not None:
            if isinstance(_guess, array):
                kw["pixels"] = _guess
            elif isinstance(_guess, str):
                kw["filename"] = _guess
            elif isinstance(_guess, file):
                kw["file"] = _guess

        if "filename" in kw:
            self.file = file(kw["filename"])
        elif "file" in kw:
            self.file = kw["file"]
        elif "pixels" in kw:
            self.file = _readable(kw["pixels"])
        else:
            raise TypeError("expecting filename, file or pixels array")

    def read_chunk(self):
        """
        Read a PNG chunk from the input file, return tag name and data.
        """
        # http://www.w3.org/TR/PNG/#5Chunk-layout
        try:
            data_bytes, tag = struct.unpack('!I4s', self.file.read(8))
        except struct.error:
            raise ValueError('Chunk too short for header')
        data = self.file.read(data_bytes)
        if len(data) != data_bytes:
            raise ValueError('Chunk %s too short for required %i data octets'
                             % (tag, data_bytes))
        checksum = self.file.read(4)
        if len(checksum) != 4:
            raise ValueError('Chunk %s too short for checksum', tag)
        verify = zlib.crc32(tag)
        verify = zlib.crc32(data, verify)
        verify = struct.pack('!i', verify)
        if checksum != verify:
            # print repr(checksum)
            (a,) = struct.unpack('!I', checksum)
            (b,) = struct.unpack('!I', verify)
            raise ValueError("Checksum error in %s chunk: 0x%X != 0x%X"
                             % (tag, a, b))
        return tag, data

    def _reconstruct_sub(self, offset, xstep, ystep):
        """
        Reverse sub filter.
        """
        pixels = self.pixels
        a_offset = offset
        offset += self.psize * xstep
        if xstep == 1:
            for index in range(self.psize, self.row_bytes):
                x = pixels[offset]
                a = pixels[a_offset]
                pixels[offset] = (x + a) & 0xff
                offset += 1
                a_offset += 1
        else:
            byte_step = self.psize * xstep
            for index in range(byte_step, self.row_bytes, byte_step):
                for i in range(self.psize):
                    x = pixels[offset + i]
                    a = pixels[a_offset + i]
                    pixels[offset + i] = (x + a) & 0xff
                offset += self.psize * xstep
                a_offset += self.psize * xstep

    def _reconstruct_up(self, offset, xstep, ystep):
        """
        Reverse up filter.
        """
        pixels = self.pixels
        b_offset = offset - (self.row_bytes * ystep)
        if xstep == 1:
            for index in range(self.row_bytes):
                x = pixels[offset]
                b = pixels[b_offset]
                pixels[offset] = (x + b) & 0xff
                offset += 1
                b_offset += 1
        else:
            for index in range(0, self.row_bytes, xstep * self.psize):
                for i in range(self.psize):
                    x = pixels[offset + i]
                    b = pixels[b_offset + i]
                    pixels[offset + i] = (x + b) & 0xff
                offset += self.psize * xstep
                b_offset += self.psize * xstep

    def _reconstruct_average(self, offset, xstep, ystep):
        """
        Reverse average filter.
        """
        pixels = self.pixels
        a_offset = offset - (self.psize * xstep)
        b_offset = offset - (self.row_bytes * ystep)
        if xstep == 1:
            for index in range(self.row_bytes):
                x = pixels[offset]
                if index < self.psize:
                    a = 0
                else:
                    a = pixels[a_offset]
                if b_offset < 0:
                    b = 0
                else:
                    b = pixels[b_offset]
                pixels[offset] = (x + ((a + b) >> 1)) & 0xff
                offset += 1
                a_offset += 1
                b_offset += 1
        else:
            for index in range(0, self.row_bytes, self.psize * xstep):
                for i in range(self.psize):
                    x = pixels[offset+i]
                    if index < self.psize:
                        a = 0
                    else:
                        a = pixels[a_offset + i]
                    if b_offset < 0:
                        b = 0
                    else:
                        b = pixels[b_offset + i]
                    pixels[offset + i] = (x + ((a + b) >> 1)) & 0xff
                offset += self.psize * xstep
                a_offset += self.psize * xstep
                b_offset += self.psize * xstep

    def _reconstruct_paeth(self, offset, xstep, ystep):
        """
        Reverse Paeth filter.
        """
        pixels = self.pixels
        a_offset = offset - (self.psize * xstep)
        b_offset = offset - (self.row_bytes * ystep)
        c_offset = b_offset - (self.psize * xstep)
        # There's enough inside this loop that it's probably not worth
        # optimising for xstep == 1
        for index in range(0, self.row_bytes, self.psize * xstep):
            for i in range(self.psize):
                x = pixels[offset+i]
                if index < self.psize:
                    a = c = 0
                    b = pixels[b_offset+i]
                else:
                    a = pixels[a_offset+i]
                    b = pixels[b_offset+i]
                    c = pixels[c_offset+i]
                p = a + b - c
                pa = abs(p - a)
                pb = abs(p - b)
                pc = abs(p - c)
                if pa <= pb and pa <= pc:
                    pr = a
                elif pb <= pc:
                    pr = b
                else:
                    pr = c
                pixels[offset+i] = (x + pr) & 0xff
            offset += self.psize * xstep
            a_offset += self.psize * xstep
            b_offset += self.psize * xstep
            c_offset += self.psize * xstep

    # N.B. PNG files with 'up', 'average' or 'paeth' filters on the
    # first line of a pass are legal. The code above for 'average'
    # deals with this case explicitly. For up we map to the null
    # filter and for paeth we map to the sub filter.

    def reconstruct_line(self, filter_type, first_line, offset, xstep, ystep):
        # print >> sys.stderr, "Filter type %s, first_line=%s" % (
        #                      filter_type, first_line)
        filter_type += (first_line << 8)
        if filter_type == 1 or filter_type == 0x101 or filter_type == 0x104:
            self._reconstruct_sub(offset, xstep, ystep)
        elif filter_type == 2:
            self._reconstruct_up(offset, xstep, ystep)
        elif filter_type == 3 or filter_type == 0x103:
            self._reconstruct_average(offset, xstep, ystep)
        elif filter_type == 4:
            self._reconstruct_paeth(offset, xstep, ystep)
        return

    def deinterlace(self, scanlines):
        # print >> sys.stderr, ("Reading interlaced, w=%s, r=%s, planes=%s," +
        #     " bpp=%s") % (self.width, self.height, self.planes, self.bps)
        a = array('B')
        self.pixels = a
        # Make the array big enough
        temp = scanlines[0:self.width*self.height*self.psize]
        a.extend(temp)
        source_offset = 0
        for xstart, ystart, xstep, ystep in _adam7:
            # print >> sys.stderr, "Adam7: start=%s,%s step=%s,%s" % (
            #     xstart, ystart, xstep, ystep)
            filter_first_line = 1
            for y in range(ystart, self.height, ystep):
                if xstart >= self.width:
                    continue
                filter_type = scanlines[source_offset]
                source_offset += 1
                if xstep == 1:
                    offset = y * self.row_bytes
                    a[offset:offset+self.row_bytes] = \
                        scanlines[source_offset:source_offset + self.row_bytes]
                    source_offset += self.row_bytes
                else:
                    # Note we want the ceiling of (width - xstart) / xtep
                    row_len = self.psize * (
                        (self.width - xstart + xstep - 1) / xstep)
                    offset = y * self.row_bytes + xstart * self.psize
                    end_offset = (y+1) * self.row_bytes
                    skip = self.psize * xstep
                    for i in range(self.psize):
                        a[offset+i:end_offset:skip] = \
                            scanlines[source_offset + i:
                                      source_offset + row_len:
                                      self.psize]
                    source_offset += row_len
                if filter_type:
                    self.reconstruct_line(filter_type, filter_first_line,
                                          offset, xstep, ystep)
                filter_first_line = 0
        return a

    def read_flat(self, scanlines):
        a = array('B')
        self.pixels = a
        offset = 0
        source_offset = 0
        filter_first_line = 1
        for y in range(self.height):
            filter_type = scanlines[source_offset]
            source_offset += 1
            a.extend(scanlines[source_offset: source_offset + self.row_bytes])
            if filter_type:
                self.reconstruct_line(filter_type, filter_first_line,
                                      offset, 1, 1)
            filter_first_line = 0
            offset += self.row_bytes
            source_offset += self.row_bytes
        return a

    def read(self):
        """
        Read a simple PNG file, return width, height, pixels and image metadata

        This function is a very early prototype with limited flexibility
        and excessive use of memory.
        """
        signature = self.file.read(8)
        if (signature != struct.pack("8B", 137, 80, 78, 71, 13, 10, 26, 10)):
            raise Error("PNG file has invalid header")
        compressed = []
        image_metadata = {}
        while True:
            try:
                tag, data = self.read_chunk()
            except ValueError, e:
                raise Error('Chunk error: ' + e.args[0])

            # print >> sys.stderr, tag, len(data)
            if tag == 'IHDR': # http://www.w3.org/TR/PNG/#11IHDR
                (width, height, bits_per_sample, color_type,
                 compression_method, filter_method,
                 interlaced) = struct.unpack("!2I5B", data)
                bps = bits_per_sample / 8
                if bps == 0:
                    raise Error("unsupported pixel depth")
                if bps > 2 or bits_per_sample != (bps * 8):
                    raise Error("invalid pixel depth")
                if color_type == 0:
                    greyscale = True
                    has_alpha = False
                    planes = 1
                elif color_type == 2:
                    greyscale = False
                    has_alpha = False
                    planes = 3
                elif color_type == 4:
                    greyscale = True
                    has_alpha = True
                    planes = 2
                elif color_type == 6:
                    greyscale = False
                    has_alpha = True
                    planes = 4
                else:
                    raise Error("unknown PNG colour type %s" % color_type)
                if compression_method != 0:
                    raise Error("unknown compression method")
                if filter_method != 0:
                    raise Error("unknown filter method")
                self.bps = bps
                self.planes = planes
                self.psize = bps * planes
                self.width = width
                self.height = height
                self.row_bytes = width * self.psize
            elif tag == 'IDAT': # http://www.w3.org/TR/PNG/#11IDAT
                compressed.append(data)
            elif tag == 'bKGD':
                if greyscale:
                    image_metadata["background"] = struct.unpack("!1H", data)
                else:
                    image_metadata["background"] = struct.unpack("!3H", data)
            elif tag == 'tRNS':
                if greyscale:
                    image_metadata["transparent"] = struct.unpack("!1H", data)
                else:
                    image_metadata["transparent"] = struct.unpack("!3H", data)
            elif tag == 'gAMA':
                image_metadata["gamma"] = (
                    struct.unpack("!L", data)[0]) / 100000.0
            elif tag == 'IEND': # http://www.w3.org/TR/PNG/#11IEND
                break
        scanlines = array('B', zlib.decompress(''.join(compressed)))
        if interlaced:
            pixels = self.deinterlace(scanlines)
        else:
            pixels = self.read_flat(scanlines)
        image_metadata["greyscale"] = greyscale
        image_metadata["has_alpha"] = has_alpha
        image_metadata["bytes_per_sample"] = bps
        image_metadata["interlaced"] = interlaced
        return width, height, pixels, image_metadata


def test_suite(options):
    """
    Run regression test and write PNG file to stdout.
    """

    # Below is a big stack of test image generators

    def test_gradient_horizontal_lr(x, y):
        return x

    def test_gradient_horizontal_rl(x, y):
        return 1-x

    def test_gradient_vertical_tb(x, y):
        return y

    def test_gradient_vertical_bt(x, y):
        return 1-y

    def test_radial_tl(x, y):
        return max(1-math.sqrt(x*x+y*y), 0.0)

    def test_radial_center(x, y):
        return test_radial_tl(x-0.5, y-0.5)

    def test_radial_tr(x, y):
        return test_radial_tl(1-x, y)

    def test_radial_bl(x, y):
        return test_radial_tl(x, 1-y)

    def test_radial_br(x, y):
        return test_radial_tl(1-x, 1-y)

    def test_stripe(x, n):
        return 1.0*(int(x*n) & 1)

    def test_stripe_h_2(x, y):
        return test_stripe(x, 2)

    def test_stripe_h_4(x, y):
        return test_stripe(x, 4)

    def test_stripe_h_10(x, y):
        return test_stripe(x, 10)

    def test_stripe_v_2(x, y):
        return test_stripe(y, 2)

    def test_stripe_v_4(x, y):
        return test_stripe(y, 4)

    def test_stripe_v_10(x, y):
        return test_stripe(y, 10)

    def test_stripe_lr_10(x, y):
        return test_stripe(x+y, 10)

    def test_stripe_rl_10(x, y):
        return test_stripe(x-y, 10)

    def test_checker(x, y, n):
        return 1.0*((int(x*n) & 1) ^ (int(y*n) & 1))

    def test_checker_8(x, y):
        return test_checker(x, y, 8)

    def test_checker_15(x, y):
        return test_checker(x, y, 15)

    def test_zero(x, y):
        return 0

    def test_one(x, y):
        return 1

    test_patterns = {
        "GLR": test_gradient_horizontal_lr,
        "GRL": test_gradient_horizontal_rl,
        "GTB": test_gradient_vertical_tb,
        "GBT": test_gradient_vertical_bt,
        "RTL": test_radial_tl,
        "RTR": test_radial_tr,
        "RBL": test_radial_bl,
        "RBR": test_radial_br,
        "RCTR": test_radial_center,
        "HS2": test_stripe_h_2,
        "HS4": test_stripe_h_4,
        "HS10": test_stripe_h_10,
        "VS2": test_stripe_v_2,
        "VS4": test_stripe_v_4,
        "VS10": test_stripe_v_10,
        "LRS": test_stripe_lr_10,
        "RLS": test_stripe_rl_10,
        "CK8": test_checker_8,
        "CK15": test_checker_15,
        "ZERO": test_zero,
        "ONE": test_one,
        }

    def test_pattern(width, height, depth, pattern):
        a = array('B')
        fw = float(width)
        fh = float(height)
        pfun = test_patterns[pattern]
        if depth == 1:
            for y in range(height):
                for x in range(width):
                    a.append(int(pfun(float(x)/fw, float(y)/fh) * 255))
        elif depth == 2:
            for y in range(height):
                for x in range(width):
                    v = int(pfun(float(x)/fw, float(y)/fh) * 65535)
                    a.append(v >> 8)
                    a.append(v & 0xff)
        return a

    def test_rgba(size=256, depth=1,
                    red="GTB", green="GLR", blue="RTL", alpha=None):
        r = test_pattern(size, size, depth, red)
        g = test_pattern(size, size, depth, green)
        b = test_pattern(size, size, depth, blue)
        if alpha:
            a = test_pattern(size, size, depth, alpha)
        i = interleave_planes(r, g, depth, depth)
        i = interleave_planes(i, b, 2 * depth, depth)
        if alpha:
            i = interleave_planes(i, a, 3 * depth, depth)
        return i

    # The body of test_suite()
    size = 256
    if options.test_size:
        size = options.test_size
    depth = 1
    if options.test_deep:
        depth = 2

    kwargs = {}
    if options.test_red:
        kwargs["red"] = options.test_red
    if options.test_green:
        kwargs["green"] = options.test_green
    if options.test_blue:
        kwargs["blue"] = options.test_blue
    if options.test_alpha:
        kwargs["alpha"] = options.test_alpha
    pixels = test_rgba(size, depth, **kwargs)

    writer = Writer(size, size,
                    bytes_per_sample=depth,
                    transparent=options.transparent,
                    background=options.background,
                    gamma=options.gamma,
                    has_alpha=options.test_alpha,
                    compression=options.compression,
                    interlaced=options.interlace)
    writer.write_array(sys.stdout, pixels)


def read_pnm_header(infile, supported='P6'):
    """
    Read a PNM header, return width and height of the image in pixels.
    """
    header = []
    while len(header) < 4:
        line = infile.readline()
        sharp = line.find('#')
        if sharp > -1:
            line = line[:sharp]
        header.extend(line.split())
        if len(header) == 3 and header[0] == 'P4':
            break # PBM doesn't have maxval
    if header[0] not in supported:
        raise NotImplementedError('file format %s not supported' % header[0])
    if header[0] != 'P4' and header[3] != '255':
        raise NotImplementedError('maxval %s not supported' % header[3])
    return int(header[1]), int(header[2])


def color_triple(color):
    """
    Convert a command line color value to a RGB triple of integers.
    FIXME: Somewhere we need support for greyscale backgrounds etc.
    """
    if color.startswith('#') and len(color) == 4:
        return (int(color[1], 16),
                int(color[2], 16),
                int(color[3], 16))
    if color.startswith('#') and len(color) == 7:
        return (int(color[1:3], 16),
                int(color[3:5], 16),
                int(color[5:7], 16))
    elif color.startswith('#') and len(color) == 13:
        return (int(color[1:5], 16),
                int(color[5:9], 16),
                int(color[9:13], 16))


def _main():
    """
    Run the PNG encoder with options from the command line.
    """
    # Parse command line arguments
    from optparse import OptionParser
    version = '%prog ' + __revision__.strip('$').replace('Rev: ', 'r')
    parser = OptionParser(version=version)
    parser.set_usage("%prog [options] [pnmfile]")
    parser.add_option("-i", "--interlace",
                      default=False, action="store_true",
                      help="create an interlaced PNG file (Adam7)")
    parser.add_option("-t", "--transparent",
                      action="store", type="string", metavar="color",
                      help="mark the specified color as transparent")
    parser.add_option("-b", "--background",
                      action="store", type="string", metavar="color",
                      help="save the specified background color")
    parser.add_option("-a", "--alpha",
                      action="store", type="string", metavar="pgmfile",
                      help="alpha channel transparency (RGBA)")
    parser.add_option("-g", "--gamma",
                      action="store", type="float", metavar="value",
                      help="save the specified gamma value")
    parser.add_option("-c", "--compression",
                      action="store", type="int", metavar="level",
                      help="zlib compression level (0-9)")
    parser.add_option("-T", "--test",
                      default=False, action="store_true",
                      help="create a test image")
    parser.add_option("-R", "--test-red",
                      action="store", type="string", metavar="pattern",
                      help="test pattern for the red image layer")
    parser.add_option("-G", "--test-green",
                      action="store", type="string", metavar="pattern",
                      help="test pattern for the green image layer")
    parser.add_option("-B", "--test-blue",
                      action="store", type="string", metavar="pattern",
                      help="test pattern for the blue image layer")
    parser.add_option("-A", "--test-alpha",
                      action="store", type="string", metavar="pattern",
                      help="test pattern for the alpha image layer")
    parser.add_option("-D", "--test-deep",
                      default=False, action="store_true",
                      help="use test patterns with 16 bits per layer")
    parser.add_option("-S", "--test-size",
                      action="store", type="int", metavar="size",
                      help="width and height of the test image")
    (options, args) = parser.parse_args()

    # Convert options
    if options.transparent is not None:
        options.transparent = color_triple(options.transparent)
    if options.background is not None:
        options.background = color_triple(options.background)

    # Run regression tests
    if options.test:
        return test_suite(options)

    # Prepare input and output files
    if len(args) == 0:
        ppmfilename = '-'
        ppmfile = sys.stdin
    elif len(args) == 1:
        ppmfilename = args[0]
        ppmfile = open(ppmfilename, 'rb')
    else:
        parser.error("more than one input file")
    outfile = sys.stdout

    # Encode PNM to PNG
    width, height = read_pnm_header(ppmfile)
    writer = Writer(width, height,
                    transparent=options.transparent,
                    background=options.background,
                    has_alpha=options.alpha is not None,
                    gamma=options.gamma,
                    compression=options.compression)
    if options.alpha is not None:
        pgmfile = open(options.alpha, 'rb')
        awidth, aheight = read_pnm_header(pgmfile, 'P5')
        if (awidth, aheight) != (width, height):
            raise ValueError("alpha channel image size mismatch" +
                             " (%s has %sx%s but %s has %sx%s)"
                             % (ppmfilename, width, height,
                                options.alpha, awidth, aheight))
        writer.convert_ppm_and_pgm(ppmfile, pgmfile, outfile,
                           interlace=options.interlace)
    else:
        writer.convert_ppm(ppmfile, outfile,
                           interlace=options.interlace)


if __name__ == '__main__':
    _main()
