"""Tests the color conversion functions."""

from sympy.interactive.colorconversion import ColorTypes, convert_color_to

def test_conversion():
    assert convert_color_to('Transparent', [ColorTypes.RGB]) == 'Transparent'
    assert convert_color_to('cmyk 0.1 0.4 0.5 0.5', [ColorTypes.RGB]) == 'RGB 115 76 64'
    assert convert_color_to('#aa0011', [ColorTypes.RGB, ColorTypes.rgb, ColorTypes.HTML]) == 'RGB 170 0 17'

def test_conversion_for_svg():
    svg_colors = [ColorTypes.rgb, ColorTypes.cmyk, ColorTypes.gray]
    assert convert_color_to('Transparent', svg_colors) == 'Transparent'
    assert convert_color_to('Black', svg_colors) == 'Black'
    assert convert_color_to('rgb 0.1 0 1', svg_colors) == 'rgb 0.1 0 1'
    assert convert_color_to('RGB 255 100 0', svg_colors) ==  'rgb 1 0.39216 0'
    assert convert_color_to('cmyk 0.1 0 1 0.5', svg_colors) == 'cmyk 0.1 0 1 0.5'
    assert convert_color_to('#00FF40', svg_colors) == 'rgb 0 1 0.25098'
    assert convert_color_to('HTML 00FF40', svg_colors) == 'rgb 0 1 0.25098'
    assert convert_color_to('gray 0.4', svg_colors) == 'gray 0.4'


def test_conversion_for_png():
    png_colors = [ColorTypes.rgb, ColorTypes.cmyk, ColorTypes.RGB, ColorTypes.HTML, ColorTypes.gray]
    assert convert_color_to('Transparent', png_colors) == 'Transparent'
    assert convert_color_to('Black', png_colors) == 'Black'
    assert convert_color_to('rgb 0.1 0 1', png_colors) == 'rgb 0.1 0 1'
    assert convert_color_to('RGB 255 100 0', png_colors) ==  'RGB 255 100 0'
    assert convert_color_to('cmyk 0.1 0 1 0.5', png_colors) == 'cmyk 0.1 0 1 0.5'
    assert convert_color_to('#00FF40', png_colors) == 'RGB 0 255 64'
    assert convert_color_to('HTML 00FF40', png_colors) == 'HTML 00FF40'
    assert convert_color_to('gray 0.4', png_colors) == 'gray 0.4'


def test_converstion_for_mpl():
    mpl_colors = [ColorTypes.rgbhex]
    assert convert_color_to('Transparent', mpl_colors) == 'Transparent'
    assert convert_color_to('Black', mpl_colors) == 'Black'
    assert convert_color_to('rgb 0.1 0 1', mpl_colors) == '#1A00FF'
    assert convert_color_to('RGB 255 100 0', mpl_colors) ==  '#FF6400'
    assert convert_color_to('cmyk 0.1 0 1 0.5', mpl_colors) == '#738000'
    assert convert_color_to('#00FF40', mpl_colors) == '#00FF40'
    assert convert_color_to('HTML 00FF40', mpl_colors) == '#00FF40'
    assert convert_color_to('gray 0.4', mpl_colors) == '#666666'
