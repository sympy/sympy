#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A tool generate AUTHORS. We started tracking authors before moving to git, so
we have to do some manual rearrangement of the git history authors in order to
get the order in AUTHORS.
"""

from __future__ import unicode_literals
from __future__ import print_function

import os
import sys

from fabric.api import local, env
from fabric.colors import yellow, blue, green, red
from fabric.utils import error

mailmap_update_path = os.path.abspath(__file__)
mailmap_update_dir = os.path.dirname(mailmap_update_path)
sympy_top = os.path.split(mailmap_update_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')

if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

from sympy.utilities.misc import filldedent

try:
    # Only works in newer versions of fabric
    env.colorize_errors = True
except AttributeError:
    pass

git_command = """git log --topo-order --reverse --format="%aN <%aE>" | awk ' !x[$0]++'"""

git_people = unicode(local(git_command, capture=True), 'utf-8').strip().split("\n")

from distutils.version import LooseVersion

git_ver = local('git --version', capture=True)[12:]
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
git_people.insert(35, "*Case Van Horsen <casevh@gmail.com>")
git_people.insert(43, "*Dan <coolg49964@gmail.com>")
move(git_people, 57, 59) # Aaron Meurer
move(git_people, 58, 57) # Andrew Docherty
move(git_people, 67, 66) # Chris Smith
move(git_people, 79, 76) # Kevin Goodsell
git_people.insert(84, "*Chu-Ching Huang <cchuang@mail.cgu.edu.tw>")
move(git_people, 93, 92) # James Pearson

git_people.pop(226) # Sergey B Kirpichev


header = """\
All people who contributed to SymPy by sending at least a patch or more (in the
order of the date of their first contribution), except those who explicitly
didn't want to be mentioned. People with a * next to their names are not found
in the metadata of the git history. This file is generated automatically by
running `./bin/authors_update.py`.
"""

fd = open(os.path.realpath(os.path.join(__file__, os.path.pardir,
    os.path.pardir, "AUTHORS")), "w")
fd.write(header)
fd.write("\n")
fd.write("\n".join(git_people).encode("utf8"))
fd.write("\n")
