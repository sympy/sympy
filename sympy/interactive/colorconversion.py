"""
Functions to convert between different LaTeX color spaces.

Main purpose is to be able to find one valid color space among those provided.

Will always return the original string if not possible to convert and not say
anything about it.
"""

# -*- coding: utf-8 -*-

from enum import Enum
from collections import namedtuple

ColorType = namedtuple('ColorType', ['value', 'int_type', 'to_str'])

class ColorTypes(Enum):
    @property
    def int_type(self):
        return self.value.int_type

    def to_str(self, s):
        return self.value.to_str(s)

    STRING = ColorType(1, None, lambda x: 'unknown')
    rgb = ColorType(2, False, lambda x: 'rgb {}'.format(" ".join(['{:.5g}'.format(f) for f in x])))
    RGB = ColorType(3, True, lambda x: 'RGB {}'.format(" ".join(['{}'.format(f) for f in x])))
    HTML = ColorType(4, True, lambda x: 'HTML {}'.format("".join(['{:02X}'.format(f) for f in x])))
    cmyk = ColorType(5, False, lambda x: 'rgb {}'.format(" ".join(['{:.5g}'.format(f) for f in rgb_to_cmyk(x)])))
    rgbhex = ColorType(6, True, lambda x: '#{}'.format("".join(['{:02X}'.format(f) for f in x])))
    gray = ColorType(7, False, lambda x: 'gray {:.5g}'.format(rgb_to_cmyk(x)[3]))


def get_color_type(s):
    s = s.strip()
    if s.startswith("rgb"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 4:
            return ColorTypes.rgb
    if s.startswith("RGB"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 4:
            return ColorTypes.RGB
    if s.startswith("cmyk"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 5:
            return ColorTypes.cmyk
    if s.startswith("HTML"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 2 and len(p[1]) == 6:
            return ColorTypes.HTML
    if s.startswith("#"):
        if len(s) == 7:
            return ColorTypes.rgbhex
    if s.startswith("gray"):
        p = [x for x in s.split(" ") if x]
        if len(p) == 2:
            return ColorTypes.gray
    return ColorTypes.STRING


def convert_color_to(s, valid_types):
    t = get_color_type(s)
    if t == ColorTypes.STRING or t in valid_types:
        return s

    val, int_type = _get_color_value(s, t)
    int_types = [x for x in valid_types if x.int_type == True]
    float_types = [x for x in valid_types if x.int_type == False]
    if int_type == True:
        if int_types:
            return int_types[0].to_str(val)
        elif float_types:
            return float_types[0].to_str([x/255 for x in val])
    elif int_type == False:
        if float_types:
            return float_types[0].to_str(val)
        elif int_types:
            return int_types[0].to_str([round(x*255) for x in val])
    return s


def _get_color_value(s, t=None):
    if t is None:
        t = get_color_type(s)
        if t == ColorTypes.STRING:
            return (None, None)
    if t == ColorTypes.gray:
        p = [x for x in s.split(" ") if x]
        if len(p) == 2:
            try:
                f = float(p[1])
                return ([f, f, f], t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.rgb:
        p = [x for x in s.split(" ") if x]
        if len(p) == 4:
            try:
                r = float(p[1])
                g = float(p[2])
                b = float(p[3])
                return ([r, g, b], t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.RGB:
        p = [x for x in s.split(" ") if x]
        if len(p) == 4:
            try:
                R = int(p[1])
                G = int(p[2])
                B = int(p[3])
                return ([R, G, B], t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.cmyk:
        p = [x for x in s.split(" ") if x]
        if len(p) == 5:
            try:
                c = float(p[1])
                m = float(p[2])
                y = float(p[3])
                k = float(p[4])
                return (cmyk_to_rgb([c, m, y, k]), t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.HTML:
        p = [x for x in s.split(" ") if x]
        if len(p) == 2 and len(p[1]) == 6:
            try:
                R = int(p[1][0:2], 16)
                G = int(p[1][2:4], 16)
                B = int(p[1][4:], 16)
                return ([R, G, B], t.int_type)
            except ValueError:
                pass
    elif t == ColorTypes.rgbhex:
        p = s.strip()
        if len(p) == 7:
            try:
                R = int(p[1:3], 16)
                G = int(p[3:5], 16)
                B = int(p[5:], 16)
                return ([R, G, B], t.int_type)
            except ValueError:
                pass
    return (None, None)


def rgb_to_cmyk(rgb):
    """
    Convert from rgb to cmyk.

    Assume that rgb values are floating-point between 0 and 1.
    """
    r = rgb[0]
    g = rgb[1]
    b = rgb[2]
    k = 1 - max(r, g, b)
    if k != 1:
        c = (1 - r - k)/(1 - k)
        m = (1 - g - k)/(1 - k)
        y = (1 - b - k)/(1 - k)
    else:
        c = 0
        m = 0
        y = 0
    return [c, m, y, k]


def cmyk_to_rgb(cmyk):
    """
    Convert from cmyk to rgb.

    Assume that cmyk values are floating-point between 0 and 1.
    """
    c = cmyk[0]
    m = cmyk[1]
    y = cmyk[2]
    k = cmyk[3]
    r = (1 - c)*(1 - k)
    g = (1 - m)*(1 - k)
    b = (1 - y)*(1 - k)
    return [r, g, b]
