import re
import fnmatch


# XXX Python 2 unicode import test.
# May remove after deprecating python 2.7.
message_unicode_A = \
    "File contains a unicode character : %s, line %s. " \
    "But with no encoding header. " \
    "See https://www.python.org/dev/peps/pep-0263/ " \
    "and add '# coding=utf-8'"
message_unicode_B = \
    "File contains a unicode character : %s, line %s. " \
    "But not in the whitelist. " \
    "Add the file to the whitelist in " + __file__
message_unicode_C = \
    "File contains a unicode character : %s, line %s. " \
    "And is in the whitelist, but without the encoding header. " \
    "See https://www.python.org/dev/peps/pep-0263/ " \
    "and add '# coding=utf-8'."
message_unicode_D = \
    "File does not contain a unicode character : %s." \
    "but is in the whitelist. " \
    "Remove the file from the whitelist in " + __file__
message_unicode_E = \
    "File does not contain a unicode character : %s." \
    "but contains the header '# coding=utf-8' or equivalent." \
    "Remove the header."


encoding_header_re = re.compile(
    r'^[ \t\f]*#.*?coding[:=][ \t]*([-_.a-zA-Z0-9]+)')

# Whitelist pattern for files which can have unicode.
unicode_whitelist = [
    # Author names can include non-ASCII characters
    r'*/bin/authors_update.py',

    # These files have functions and test functions for unicode input and
    # output.
    r'*/sympy/utilities/tests/test_code_quality.py',
    r'*/sympy/physics/vector/tests/test_printing.py',
    r'*/physics/quantum/tests/test_printing.py',
    r'*/sympy/vector/tests/test_printing.py',
    r'*/sympy/parsing/tests/test_sympy_parser.py',
    r'*/sympy/printing/pretty/tests/test_pretty.py',
    r'*/sympy/printing/tests/test_preview.py',
    r'*/liealgebras/type_g.py',
    r'*/liealgebras/weyl_group.py',
    r'*/liealgebras/tests/test_type_G.py',

    # wigner.py and polarization.py have unicode doctests. These probably
    # don't need to be there but some of the examples that are there are
    # pretty ugly without use_unicode (matrices need to be wrapped across
    # multiple lines etc)
    r'*/sympy/physics/wigner.py',
    r'*/sympy/physics/optics/polarization.py',
]

unicode_strict_whitelist = [
    r'*/sympy/parsing/latex/_antlr/__init__.py',
]


def test_this_file_encoding(
    fname, test_file,
    unicode_whitelist=unicode_whitelist,
    unicode_strict_whitelist=unicode_strict_whitelist):
    """Test helper function for python 2 importability test

    This test checks whether the file has
    # coding=utf-8
    or
    # -*- coding: utf-8 -*-
    line if there is a unicode character in the code

    The test may have to operate on filewise manner, so it had moved
    to a separate process.
    May remove after deprecating python 2.7.
    """
    has_coding_utf8 = False
    has_unicode = False

    is_in_whitelist = False
    is_in_strict_whitelist = False
    for patt in unicode_whitelist:
        if fnmatch.fnmatch(fname, patt):
            is_in_whitelist = True
            break
    for patt in unicode_strict_whitelist:
        if fnmatch.fnmatch(fname, patt):
            is_in_strict_whitelist = True
            is_in_whitelist = True
            break

    if is_in_whitelist:
        for idx, line in enumerate(test_file):
            if idx in (0, 1):
                match = encoding_header_re.match(line)
                if match and match.group(1).lower() == 'utf-8':
                    has_coding_utf8 = True
            try:
                line.encode(encoding='ascii')
            except (UnicodeEncodeError, UnicodeDecodeError):
                has_unicode = True
                if has_coding_utf8 is False:
                    assert False, \
                        message_unicode_C % (fname, idx + 1)

        if not has_unicode and not is_in_strict_whitelist:
            assert False, message_unicode_D % fname

    else:
        for idx, line in enumerate(test_file):
            if idx in (0, 1):
                match = encoding_header_re.match(line)
                if match and match.group(1).lower() == 'utf-8':
                    has_coding_utf8 = True
            try:
                line.encode(encoding='ascii')
            except (UnicodeEncodeError, UnicodeDecodeError):
                has_unicode = True
                if has_coding_utf8:
                    assert False, \
                        message_unicode_B % (fname, idx + 1)
                else:
                    assert False, \
                        message_unicode_A % (fname, idx + 1)

        if not has_unicode and has_coding_utf8:
            assert False, \
                message_unicode_E % fname
