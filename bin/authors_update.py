#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A tool generate AUTHORS. We started tracking authors before moving to git, so
we have to do some manual rearrangement of the git history authors in order to
get the order in AUTHORS.
"""

from __future__ import unicode_literals
from __future__ import print_function

import codecs
import os
import sys

if sys.version_info < (3, 6):
    sys.exit("This script requires Python 3.6 or newer")

from subprocess import run, PIPE
from collections import OrderedDict


def yellow(text):
    return "\033[33m%s\033[0m" % text

def blue(text):
    return "\033[34m%s\033[0m" % text

mailmap_update_path = os.path.abspath(__file__)
mailmap_update_dir = os.path.dirname(mailmap_update_path)
sympy_top = os.path.split(mailmap_update_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')

if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

from sympy.utilities.misc import filldedent

git_command = ["git", "log", "--topo-order", "--reverse", "--format=%aN <%aE>"]

git_people = run(git_command, stdout=PIPE, encoding='utf-8').stdout.strip().split("\n")

# Remove duplicates, keeping the original order
git_people = list(OrderedDict.fromkeys(git_people))

from distutils.version import LooseVersion

git_ver = run(['git', '--version'], stdout=PIPE, encoding='utf-8').stdout[12:]
if LooseVersion(git_ver) < LooseVersion('1.8.4.2'):
    print(yellow("Please use a newer git version >= 1.8.4.2"))

def move(l, i1, i2):
    x = l.pop(i1)
    l.insert(i2, x)


# Do the few changes necessary in order to reproduce AUTHORS:

move(git_people, 2, 0) # Ondřej Čertík
move(git_people, 42, 1) # Fabian Pedregosa
move(git_people, 22, 2) # Jurjen N.E. Bos
git_people.insert(4, "*Marc-Etienne M.Leveille <protonyc@gmail.com>")
move(git_people, 10, 5) # Brian Jorgensen
git_people.insert(11, "*Ulrich Hecht <ulrich.hecht@gmail.com>")
git_people.pop(12) # Kirill Smelkov
move(git_people, 12, 32) # Sebastian Krämer
move(git_people, 227, 35) # Case Van Horsen
git_people.insert(43, "*Dan <coolg49964@gmail.com>")
move(git_people, 57, 59) # Aaron Meurer
move(git_people, 58, 57) # Andrew Docherty
move(git_people, 67, 66) # Chris Smith
move(git_people, 79, 76) # Kevin Goodsell
git_people.insert(84, "*Chu-Ching Huang <cchuang@mail.cgu.edu.tw>")
move(git_people, 93, 92) # James Pearson

git_people.pop(226) # Sergey B Kirpichev

header = filldedent("""
    All people who contributed to SymPy by sending at least a patch or
    more (in the order of the date of their first contribution), except
    those who explicitly didn't want to be mentioned. People with a * next
    to their names are not found in the metadata of the git history. This
    file is generated automatically by running `./bin/authors_update.py`.
    """).lstrip()
fmt = """\n\nThere are a total of {authors_count} authors.\n"""
header_extra = fmt.format(authors_count=len(git_people))

with codecs.open(os.path.realpath(os.path.join(
        __file__, os.path.pardir, os.path.pardir, "AUTHORS")),
        "w", "utf-8") as fd:
    fd.write(header)
    fd.write(header_extra)
    fd.write("\n")
    fd.write("\n".join(git_people))
    fd.write("\n")

print(blue(filldedent("""
    Please make sure that there are no duplicates in the new AUTHORS, then
    commit the changes. You may also want to run ./bin/mailmap_update.py
    to update .mailmap as well.
""")))
