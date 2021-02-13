#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A tool to generate AUTHORS. We started tracking authors before moving
to git, so we have to do some manual rearrangement of the git history
authors in order to get the order in AUTHORS. bin/mailmap_update.py
should be run before committing the results.
"""

from __future__ import unicode_literals
from __future__ import print_function

import codecs
import sys
import os


if sys.version_info < (3, 6):
    sys.exit("This script requires Python 3.6 or newer")

from subprocess import run, PIPE
from distutils.version import LooseVersion
from collections import OrderedDict

def red(text):
    return "\033[31m%s\033[0m" % text

def yellow(text):
    return "\033[33m%s\033[0m" % text

def green(text):
    return "\033[32m%s\033[0m" % text

# put sympy on the path
mailmap_update_path = os.path.abspath(__file__)
mailmap_update_dir = os.path.dirname(mailmap_update_path)
sympy_top = os.path.split(mailmap_update_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')
if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

from sympy.utilities.misc import filldedent

# check git version
minimal = '1.8.4.2'
git_ver = run(['git', '--version'], stdout=PIPE, encoding='utf-8').stdout[12:]
if LooseVersion(git_ver) < LooseVersion(minimal):
    print(yellow("Please use a git version >= %s" % minimal))

def author_name(line):
    assert line.count("<") == line.count(">") == 1
    assert line.endswith(">")
    return line.split("<", 1)[0].strip()

def move(l, i1, i2, who):
    x = l.pop(i1)
    # this will fail if the .mailmap is not right
    assert who == author_name(x), \
        '%s was not found at line %i' % (who, i1)
    l.insert(i2, x)

# find who git knows ahout
git_command = ["git", "log", "--topo-order", "--reverse", "--format=%aN <%aE>"]
git_people = run(git_command, stdout=PIPE, encoding='utf-8').stdout.strip().split("\n")

# remove duplicates, keeping the original order
git_people = list(OrderedDict.fromkeys(git_people))

# Do the few changes necessary in order to reproduce AUTHORS:

try:
    move(git_people, 2, 0, 'Ondřej Čertík')
    move(git_people, 42, 1, 'Fabian Pedregosa')
    move(git_people, 22, 2, 'Jurjen N.E. Bos')
    git_people.insert(4, "*Marc-Etienne M.Leveille <protonyc@gmail.com>")
    move(git_people, 10, 5, 'Brian Jorgensen')
    git_people.insert(11, "*Ulrich Hecht <ulrich.hecht@gmail.com>")
    # this will fail if the .mailmap is not right
    assert 'Kirill Smelkov' == author_name(git_people.pop(12)
        ), 'Kirill Smelkov was not found at line 12'
    move(git_people, 12, 32, 'Sebastian Krämer')
    move(git_people, 227, 35, 'Case Van Horsen')
    git_people.insert(43, "*Dan <coolg49964@gmail.com>")
    move(git_people, 57, 59, 'Aaron Meurer')
    move(git_people, 58, 57, 'Andrew Docherty')
    move(git_people, 67, 66, 'Chris Smith')
    move(git_people, 79, 76, 'Kevin Goodsell')
    git_people.insert(84, "*Chu-Ching Huang <cchuang@mail.cgu.edu.tw>")
    move(git_people, 93, 92, 'James Pearson')
    # this will fail if the .mailmap is not right
    assert 'Sergey B Kirpichev' == author_name(git_people.pop(226)
        ), 'Sergey B Kirpichev was not found at line 226.'

    index = git_people.index(
        "azure-pipelines[bot] " +
        "<azure-pipelines[bot]@users.noreply.github.com>")
    git_people.pop(index)
    index = git_people.index(
        "whitesource-bolt-for-github[bot] " +
        "<whitesource-bolt-for-github[bot]@users.noreply.github.com>")
    git_people.pop(index)

except AssertionError as msg:
    print(red(msg))
    sys.exit(1)

# define new lines for the file

header = filldedent("""
    All people who contributed to SymPy by sending at least a patch or
    more (in the order of the date of their first contribution), except
    those who explicitly didn't want to be mentioned. People with a * next
    to their names are not found in the metadata of the git history. This
    file is generated automatically by running `./bin/authors_update.py`.
    """).lstrip()
fmt = """There are a total of {authors_count} authors."""
header_extra = fmt.format(authors_count=len(git_people))
lines = header.splitlines()
lines.append('')
lines.append(header_extra)
lines.append('')
lines.extend(git_people)

# compare to old lines and stop if no changes were made

old_lines = codecs.open(os.path.realpath(os.path.join(
        __file__, os.path.pardir, os.path.pardir, "AUTHORS")),
        "r", "utf-8").read().splitlines()
if old_lines == lines:
    print(green('No changes made to AUTHORS.'))
    sys.exit(0)

# check for new additions
new_authors = []
for i in sorted(set(lines) - set(old_lines)):
    try:
        author_name(i)
        new_authors.append(i)
    except AssertionError:
        continue

# write the new file
with codecs.open(os.path.realpath(os.path.join(
        __file__, os.path.pardir, os.path.pardir, "AUTHORS")),
        "w", "utf-8") as fd:
    fd.write('\n'.join(lines))
    fd.write('\n')

# warn about additions
if new_authors:
    print(yellow(filldedent("""
        The following authors were added to AUTHORS.
        If mailmap_update.py has already been run and
        each author appears as desired and is not a
        duplicate of some other author, then the
        changes can be committed. Otherwise, see
        .mailmap for instructions on how to change
        an author's entry.""")))
    print()
    for i in sorted(new_authors, key=lambda x: x.lower()):
        print('\t%s' % i)
else:
    print(yellow("The AUTHORS file was updated."))

print(red("Changes were made in the authors file"))
sys.exit(1)
