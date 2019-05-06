#!/usr/bin/env python
"""
A tool to help keep .mailmap up-to-date with the
current git authors.

See also bin/authors_update.py
"""

import codecs
import sys
import os


if sys.version_info < (3, 6):
    sys.exit("This script requires Python 3.6 or newer")

from subprocess import run, PIPE
from distutils.version import LooseVersion
from collections import defaultdict, OrderedDict

def red(text):
    return "\033[31m%s\033[0m" % text

def yellow(text):
    return "\033[33m%s\033[0m" % text

def blue(text):
    return "\033[34m%s\033[0m" % text

# put sympy on the path
mailmap_update_path = os.path.abspath(__file__)
mailmap_update_dir = os.path.dirname(mailmap_update_path)
sympy_top = os.path.split(mailmap_update_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')
if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

from sympy.utilities.misc import filldedent
from sympy.utilities.iterables import sift

# check git version
minimal = '1.8.4.2'
git_ver = run(['git', '--version'], stdout=PIPE, encoding='utf-8').stdout[12:]
if LooseVersion(git_ver) < LooseVersion(minimal):
    print(yellow("Please use a git version >= %s" % minimal))

def author_name(line):
    assert line.count("<") == line.count(">") == 1
    assert line.endswith(">")
    return line.split("<", 1)[0].strip()

def author_email(line):
    assert line.count("<") == line.count(">") == 1
    assert line.endswith(">")
    return line.split("<", 1)[1][:-1].strip()

sysexit = 0
print(blue("checking git authors..."))

# read git authors
git_command = ['git', 'log', '--format=%aN <%aE>']
git_people = sorted(set(run(git_command, stdout=PIPE, encoding='utf-8').stdout.strip().split("\n")))

# check for ambiguous emails

dups = defaultdict(list)
near_dups = defaultdict(list)
for i in git_people:
    k = i.split('<')[1]
    dups[k].append(i)
    near_dups[k.lower()].append((k, i))
multi = [k for k in dups if len(dups[k]) > 1]
if multi:
    print()
    print(red(filldedent("""
        Ambiguous email address error: each address should refer to a
        single author. Disambiguate the following in .mailmap.
        Then re-run this script.""")))
    for k in multi:
        print()
        for e in sorted(dups[k]):
            print('\t%s' % e)
    sysexit = 1

# warn for nearly ambiguous email addresses
dups = near_dups
# some may have been real dups, so disregard those
# for which all email addresses were the same
multi = [k for k in dups if len(dups[k]) > 1 and
    len(set([i for i, _ in dups[k]])) > 1]
if multi:
    # not fatal but make it red
    print()
    print(red(filldedent("""
        Ambiguous email address warning: git treats the
        following as distinct but .mailmap will treat them
        the same. If these are not all the same person then,
        when making an entry in .mailmap, be sure to include
        both commit name and address (not just the address).""")))
    for k in multi:
        print()
        for _, e in sorted(dups[k]):
            print('\t%s' % e)

# warn for ambiguous names
dups = defaultdict(list)
for i in git_people:
    dups[author_name(i)].append(i)
multi = [k for k in dups if len(dups[k]) > 1]
if multi:
    print()
    print(yellow(filldedent("""
        Ambiguous name warning: if a person uses more than
        one email address, entries should be added to .mailmap
        to merge them into a single canonical address.
        Then re-run this script.
        """)))
    for k in multi:
        print()
        for e in sorted(dups[k]):
            print('\t%s' % e)

bad_names = []
bad_emails = []
for i in git_people:
    name = author_name(i)
    email = author_email(i)
    if '@' in name:
        bad_names.append(i)
    elif '@' not in email:
        bad_emails.append(i)
if bad_names:
    print()
    print(yellow(filldedent("""
        The following people appear to have an email address
        listed for their name. Entries should be added to
        .mailmap so that names are formatted like
        "Name <email address>".
        """)))
    for i in bad_names:
        print("\t%s" % i)

# TODO: Should we check for bad emails as well? Some people have empty email
# addresses. The above check seems to catch people who get the name and email
# backwards, so let's leave this alone for now.

# if bad_emails:
#     print()
#     print(yellow(filldedent("""
#         The following names do not appear to have valid
#         emails. Entries should be added to .mailmap that
#         use a proper email address. If there is no email
#         address for a person, use "none@example.com".
#         """)))
#     for i in bad_emails:
#         print("\t%s" % i)

print()
print(blue("checking .mailmap..."))

# put entries in order -- this will help the user
# to see if there are already existing entries for an author
file = codecs.open(os.path.realpath(os.path.join(
        __file__, os.path.pardir, os.path.pardir, ".mailmap")),
        "r", "utf-8").read()
blankline = not file or file.endswith('\n')
lines = file.splitlines()
def key(line):
    # return lower case first address on line or
    # raise an error if not an entry
    if '#' in line:
        line = line.split('#')[0]
    L, R = line.count("<"), line.count(">")
    assert L == R and L in (1, 2)
    return line.split(">", 1)[0].split("<")[1].lower()

who = OrderedDict()
for i, line in enumerate(lines):
    try:
        who.setdefault(key(line), []).append(line)
    except AssertionError:
        who[i] = [line]

out = []
for k in who:
    # put long entries before short since if they match, the
    # short entries will be ignored. The ORDER MATTERS
    # so don't re-order the lines for a given address.
    # Other tidying up could be done but we won't do that here.
    def short_entry(line):
        if line.count('<') == 2:
            if line.split('>', 1)[1].split('<')[0].strip():
                return False
        return True
    if len(who[k]) == 1:
        line = who[k][0]
        if not line.strip():
            continue  # ignore blank lines
        out.append(line)
    else:
        uniq = list(OrderedDict.fromkeys(who[k]))
        short, long = sift(uniq, short_entry, binary=True)
        out.extend(long)
        out.extend(short)

if out != lines or not blankline:
    # write lines
    with codecs.open(os.path.realpath(os.path.join(
            __file__, os.path.pardir, os.path.pardir, ".mailmap")),
            "w", "utf-8") as fd:
        fd.write('\n'.join(out))
        fd.write('\n')
    print()
    if out != lines:
        print(yellow('.mailmap lines were re-ordered.'))
    else:
        print(yellow('blank line added to end of .mailmap'))

sys.exit(sysexit)
