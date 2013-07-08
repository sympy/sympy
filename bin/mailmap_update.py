#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A tool to help keep .mailmap and AUTHORS up-to-date.
"""

import os

from fabric.api import local, env
from fabric.colors import yellow
from fabric.utils import warn

try:
    # Only works in newer versions of fabric
    env.colorize_errors = True
except AttributeError:
    pass

git_command = 'git log --format="%aN <%aE>" | sort -u'

git_people = unicode(local(git_command, capture=True), 'utf-8').strip().split("\n")

with open(os.path.realpath(os.path.join(__file__, os.path.pardir,
    os.path.pardir, "AUTHORS"))) as fd:
    AUTHORS = unicode(fd.read(), 'utf-8')

firstauthor = u"Ondřej Čertík"

authors = AUTHORS[AUTHORS.find(firstauthor):].strip().split('\n')

# People who don't want to be listed in AUTHORS
authors_skip = ["Kirill Smelkov <kirr@landau.phys.spbu.ru>"]

predate_git = 0

print
print yellow("People who are in AUTHORS but not in git:")
print

for name in sorted(set(authors) - set(git_people)):
    if name.startswith("*"):
        # People who are in AUTHORS but predate git
        predate_git += 1
        continue
    print name

print
print yellow("People who are in git but not in AUTHORS:")
print

for name in sorted(set(git_people) - set(authors) - set(authors_skip)):
    print name

authors_count = AUTHORS[AUTHORS.find(firstauthor):].strip().count("\n")
adjusted_authors_count = (
    authors_count
    + 1 # Final line won't have newline due to above strip()
    - predate_git
    + len(authors_skip)
    )
git_count = len(git_people)

print
print yellow("There are {git_count} people in git, and {adjusted_authors_count} "
    "(adjusted) people from AUTHORS".format(git_count=git_count,
    adjusted_authors_count=adjusted_authors_count))

if git_count != adjusted_authors_count:
    warn("These two are not the same!")
