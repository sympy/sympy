#!/usr/bin/env python

"""
This script creates logos of different formats from the source "sympy.svg"

Requirements:
    rsvg-convert    - for converting to *.png format (librsvg2-bin deb package)
    imagemagick     - for converting to *.ico favicon format
"""

from optparse import OptionParser
import xml.dom.minidom
import os.path
import logging
import subprocess
import sys


default_source_dir = os.path.join(os.path.dirname(__file__), "src/logo")
default_source_svg = "sympy.svg"
default_output_dir = os.path.join(os.path.dirname(__file__), "_build/logo")

# those are the options for resizing versions without tail or text
svg_sizes = {}
svg_sizes['notail'] = {"prefix":"notail", "dx":-70, "dy":-20, "size":690, "title":"SymPy Logo, with no tail"}
svg_sizes['notail-notext'] = {"prefix":"notailtext", "dx":-70, "dy":60, "size":690, "title":"SymPy Logo, with no tail, no text"}
svg_sizes['notext'] = {"prefix":"notext", "dx":-7, "dy":90, "size":750, "title":"SymPy Logo, with no text"}

# The list of identifiers of various versions
versions = ['notail', 'notail-notext', 'notext']

parser = OptionParser()
parser.set_usage("%s [options ...]")

parser.add_option("--source-dir", type="string", dest="source_dir",
    help="Directory of the source *.svg file [default: %default]",
    default=default_source_dir)

parser.add_option("--source-svg", type="string", dest="source_svg",
    help="File name of the source *.svg file [default: %default]",
    default=default_source_svg)

parser.add_option("--svg", action="store_true", dest="generate_svg",
    help="Generate *.svg versions without tails and without text 'SymPy' [default: %default]",
    default=False)

parser.add_option("--png", action="store_true", dest="generate_png",
    help="Generate *.png versions [default: %default]",
    default=False)

parser.add_option("--ico", action="store_true", dest="generate_ico",
    help="Generate *.ico versions [default: %default]",
    default=False)

parser.add_option("--clear", action="store_true", dest="clear",
    help="Remove temporary files [default: %default]",
    default=False)

parser.add_option("-a", "--all", action="store_true", dest="generate_all",
    help="Shorthand for '--svg --png --ico --clear' options [default: %default]",
    default=True)

parser.add_option("-s", "--sizes", type="string", dest="sizes",
    help="Sizes of png pictures [default: %default]",
    default="160,500")

parser.add_option("--icon-sizes", type="string", dest="icon_sizes",
    help="Sizes of icons embedded in favicon file [default: %default]",
    default="16,32,48,64")

parser.add_option("--output-dir", type="string", dest="output_dir",
    help="Outpu dir [default: %default]",
    default=default_output_dir)

parser.add_option("-d", "--debug", action="store_true", dest="debug",
    default=False)

def main():
    options, args = parser.parse_args()
    if options.debug:
        logging.basicConfig(level=logging.DEBUG)

    fn_source = os.path.join(options.source_dir, options.source_svg)

    if options.generate_svg or options.generate_all:
        generate_notail_notext_versions(fn_source, options.output_dir)

    if options.generate_png or options.generate_all:
        sizes = options.sizes.split(",")
        sizes = [int(s) for s in sizes]
        convert_to_png(fn_source, options.output_dir, sizes)

    if options.generate_ico or options.generate_all:
        sizes = options.icon_sizes.split(",")
        sizes = [int(s) for s in sizes]
        convert_to_ico(fn_source, options.output_dir, sizes)

def generate_notail_notext_versions(fn_source, output_dir):
    for ver in versions:
        properties = svg_sizes[ver]

        doc = load_svg(fn_source)

        (notail, notext) = versionkey_to_boolean_tuple(ver)

        g_tail = searchElementById(doc, "SnakeTail", "g")
        if notail:
            g_tail.setAttribute("display", "none")

        g_text = searchElementById(doc, "SymPy_text", "g")
        if notext:
            g_text.setAttribute("display", "none")

        g_logo = searchElementById(doc, "SympyLogo", "g")
        dx = properties["dx"]
        dy = properties["dy"]
        transform = "translate(%d,%d)" % (dx, dy)
        g_logo.setAttribute("transform", transform)

        svg = searchElementById(doc, "svg_SympyLogo", "svg")
        newsize = properties["size"]
        svg.setAttribute("width", "%d" % newsize)
        svg.setAttribute("height", "%d" % newsize)

        title = svg.getElementsByTagName("title")[0]
        title.firstChild.data = properties["title"]

        desc = svg.getElementsByTagName("desc")[0]
        desc.appendChild(doc.createTextNode("\n\nThis file is generated from %s !" % fn_source))

        fn_out = get_svg_filename_from_versionkey(fn_source, ver)
        fn_out = os.path.join(output_dir, fn_out)
        save_svg(fn_out, doc)

def convert_to_png(fn_source, output_dir, sizes):
    svgs = list(versions)
    svgs.insert(0, '')

    cmd = "rsvg-convert"
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    p.communicate()
    if p.returncode == 127:
        logging.error("%s: command not found" % cmd)
        sys.exit(p.returncode)

    for ver in svgs:
        if ver == '':
            fn_svg = fn_source
        else:
            fn_svg = get_svg_filename_from_versionkey(fn_source, ver)
            fn_svg = os.path.join(output_dir, fn_svg)

        basename = os.path.basename(fn_svg)
        name, ext = os.path.splitext(basename)
        for size in sizes:
            fn_out = "%s-%dpx.png" % (name, size)
            fn_out = os.path.join(output_dir, fn_out)

            cmd = "rsvg-convert %s -f png -o %s -h %d -w %d" % (fn_svg, fn_out,
                size, size)

            p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
            p.communicate()
            if p.returncode != 0:
                logging.error("Return code is not 0: Command: %s" % cmd)
                logging.error("return code: %s" % p.returncode)
                sys.exit(p.returncode)
            else:
                logging.debug("command: %s" % cmd)
                logging.debug("return code: %s" % p.returncode)

def convert_to_ico(fn_source, output_dir, sizes):
    # firstly prepare *.png files, which will be embeded
    # into the *.ico files.
    convert_to_png(fn_source, output_dir, sizes)

    svgs = list(versions)
    svgs.insert(0, '')

    cmd = "convert"
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    p.communicate()
    if p.returncode == 127:
        logging.error("%s: command not found" % cmd)
        sys.exit(p.returncode)

    for ver in svgs:
        if ver == '':
            fn_svg = fn_source
        else:
            fn_svg = get_svg_filename_from_versionkey(fn_source, ver)
            fn_svg = os.path.join(output_dir, fn_svg)
        basename = os.path.basename(fn_svg)
        name, ext = os.path.splitext(basename)

        # calculate the list of *.png files
        pngs = []
        for size in sizes:
            fn_png= "%s-%dpx.png" % (name, size)
            fn_png = os.path.join(output_dir, fn_png)
            pngs.append(fn_png)

        # convert them to *.ico
        fn_out = "%s-favicon.ico" % name
        fn_out = os.path.join(output_dir, fn_out)

        cmd = "convert %s %s" % (" ".join(pngs), fn_out)

        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        p.communicate()
        if p.returncode != 0:
            logging.error("Return code is not 0: Command: %s" % cmd)
            logging.error("return code: %s" % p.returncode)
            sys.exit(p.returncode)
        else:
            logging.debug("command: %s" % cmd)
            logging.debug("return code: %s" % p.returncode)

def versionkey_to_boolean_tuple(ver):
    notail = False
    notext = False
    vers = ver.split("-")
    notail = 'notail' in vers
    notext = 'notext' in vers
    return (notail, notext)

def get_svg_filename_from_versionkey(fn_source, ver):
    basename = os.path.basename(fn_source)
    if ver == '':
        return basename
    name, ext = os.path.splitext(basename)
    prefix = svg_sizes[ver]["prefix"]
    fn_out = "%s-%s.svg" % (name, prefix)
    return fn_out

def searchElementById(node, Id, tagname):
    """
    Search element by id in the all childs and descendant of node.

    id is lower case, not ID which is usually used for gelElementById
    """
    nodes = node.getElementsByTagName(tagname)
    for node in nodes:
        an = node.getAttributeNode('id')
        if an and an.nodeValue == Id:
            return node

def load_svg(fn):
    doc = xml.dom.minidom.parse(fn)
    return doc

def save_svg(fn, doc):
    with open(fn, "wb") as f:
        xmlstr = doc.toxml("utf-8")
        f.write(xmlstr)
        logging.info(" File saved: %s" % fn)

main()
